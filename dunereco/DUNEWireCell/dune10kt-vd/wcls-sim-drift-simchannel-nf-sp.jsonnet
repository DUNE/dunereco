// High layer / entry point for the dune10kt-vd WC/LS sim-drift-simchannel-nf-sp
// job.  It is named to match the hosting FHiCL (wcls-sim-drift-simchannel-nf-sp.fcl).
//
// Responsibilities of THIS layer only:
//   1. resolve all std.extVar() fcl knobs into one `fcl_params` object,
//   2. build `params` (params-10kt) and `tools` (narrowed to the anodes to run),
//   3. dispatch on process_mode to a modal layer that builds the graph:
//        - "multi-sim"            -> -multi.jsonnet  (per-anode fan-out, per-anode
//                                     FrameSavers; process_anode_indices list)
//        - "single-sim/-sp/-sim-sp" and the "full" fallback -> -single.jsonnet
//   4. wrap the graph in the flow-engine app and emit.
//
// Shared detector glue lives in -common.jsonnet; the modal graph assembly lives
// in -single.jsonnet / -multi.jsonnet.  See those files for details.
//
// NOTE on extVar laziness: per-mode knobs (process_anode_index vs
// process_anode_indices / frame_tag_pattern) are referenced only on the code path
// that needs them, so a FHiCL need only supply the knobs for its chosen mode.

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local util = import 'pgrapher/experiment/protodunevd/funcs.jsonnet';
local params_maker = import 'pgrapher/experiment/dune10kt-vd/params-10kt.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';

local common_maker = import 'pgrapher/experiment/dune10kt-vd/wcls-sim-drift-simchannel-nf-sp-common.jsonnet';
local single_maker = import 'pgrapher/experiment/dune10kt-vd/wcls-sim-drift-simchannel-nf-sp-single.jsonnet';
local multi_maker = import 'pgrapher/experiment/dune10kt-vd/wcls-sim-drift-simchannel-nf-sp-multi.jsonnet';

// ---- 1. resolve the fcl knobs -------------------------------------------
// IMPORTANT: the extVar-derived knobs are forced to concrete literals ONCE, via
// std.parseJson(std.manifestJson(...)).  Without this, go-jsonnet does NOT
// memoize the extVar-valued numeric params (nticks, elecGain, response_plane,
// ...) as they thread through params_maker into the per-anode response / sim /
// sigproc construction: it REBUILDS the (nticks-sized) response arrays once per
// processed anode.  For the multi-sim fan that is ~GBs of config-compile memory
// PER CRU (the dominant per-CRU job-memory cost); forcing the values concrete
// restores normal sharing (measured: ~10 GB -> ~0.13 GB for 4 CRUs).
local fcl_params = std.parseJson(std.manifestJson({
  G4RefTime: std.extVar('G4RefTime') * wc.us,
  response_plane: std.extVar('response_plane') * wc.cm,
  nticks: std.extVar('nticks'),
  use_dnnroi: std.extVar('use_dnnroi'),
  process_mode: std.extVar('process_mode'),
  geo_planeid_labels: std.extVar('geo_planeid_labels'),
  // art InputTag ("label[:instance[:process]]") of the input SimEnergyDeposit
  // collection the WCLS depo source reads (and the DepoFluxWriter's sed_label).
  input_depos_tag: std.extVar('input_depos_tag'),
  use_hydra: std.extVar('use_hydra'),
  save_rawdigits: std.extVar('save_rawdigits'),
  adc_resolution: std.extVar('adc_resolution'),
  // Bottom drift (cold electronics) front-end gain [mV/fC]
  elecGain: std.extVar('elecGain') * wc.mV / wc.fC,
  // single-* : which one anode to process
  process_anode_index: std.extVar('process_anode_index'),
  // multi-*  : list of anode indices, and the per-anode RawDigit tag/instance pattern
  process_anode_indices: std.extVar('process_anode_indices'),
  frame_tag_pattern: std.extVar('frame_tag_pattern'),
}));

// ---- 2. params + tools ---------------------------------------------------
// The lar/files extVar overrides are forced concrete for the same reason.
local lar_over = std.parseJson(std.manifestJson({
  DL: std.extVar('DL') * wc.cm2 / wc.ns,
  DT: std.extVar('DT') * wc.cm2 / wc.ns,
  lifetime: std.extVar('lifetime') * wc.us,
  drift_speed: util.drift_velocity(std.extVar('efield'), std.extVar('temperature')) * wc.mm / wc.us,
}));
local files_over = std.parseJson(std.manifestJson({
  wires: std.extVar('files_wires'),
  // same field response for [bottom, top] drift volumes
  fields: [std.extVar('files_fields'), std.extVar('files_fields')],
}));
local params = params_maker(fcl_params) {
  lar: super.lar + lar_over,
  files: super.files + files_over,
};

// Where the anodes get their wire geometry (params.wire_schema):
//   "service" -- larwirecell's art WireSchemaService, built from the job's own
//                LArSoft geometry + channel map, so WCT and LArSoft wires and
//                channels agree by construction.  Needs a WireSchemaService in
//                the FHiCL services (grouping "pair-tpc" for this detector).
//   "file"    -- the WCT wires file named by params.files_wires.
local wire_schema = std.extVar('wire_schema');
local tools_base = tools_maker(params);
local tools_wires =
  if wire_schema == 'service'
  then tools_base { wires: { type: 'wclsWireSchemaSource', name: '', data: {} } }
  else if wire_schema == 'file'
  then tools_base
  else error 'wire_schema must be "service" or "file", got "%s"' % wire_schema;

// Where the anodes and the drifter get their "anode", "response" and "cathode"
// drift-axis planes:
//   "service" -- larwirecell's WireSchemaService, from the LArSoft geometry.
//                The response plane is params.det.response_plane from the
//                collection plane, which must match the field response.
//   "params"  -- the fixed per-volume faces of params.det.volumes.
local drift_planes = std.extVar('drift_planes');
local tools_all =
  if drift_planes == 'service'
  then tools_wires {
    drift_planes:: {
      type: 'wclsDriftPlanesSource',
      name: '',
      data: { response_distance: params.det.response_plane },
    },
  }
  else if drift_planes == 'params'
  then tools_wires
  else error 'drift_planes must be "service" or "params", got "%s"' % drift_planes;

local is_single =
  fcl_params.process_mode == 'single-sim'
  || fcl_params.process_mode == 'single-sp'
  || fcl_params.process_mode == 'single-sim-sp';
local is_multi = fcl_params.process_mode == 'multi-sim';

// The anodes this job will process, as positions into tools_all.anodes.
local anode_indices =
  if is_multi then fcl_params.process_anode_indices
  else if is_single then [fcl_params.process_anode_index]
  else std.range(0, std.length(tools_all.anodes) - 1);  // "full"

local tools = tools_all { anodes: [tools_all.anodes[i] for i in anode_indices] };

// ---- 3. dispatch to the modal layer -------------------------------------
local com = common_maker(params, tools, fcl_params);
local built =
  if is_multi
  then multi_maker(params, tools, fcl_params, com)
  else single_maker(params, tools, fcl_params, com);
local graph = built.graph;

// ---- 4. emit -------------------------------------------------------------
// Order-preserving O(N log N) dedup.  pgraph's g.uses()/g.edges() dedup via
// wc.unique_list() which is O(N^2) DEEP-object comparison -- fine for ProtoDUNE
// but impractical for a full FD-VD single graph (~5800 nodes, all uniquely
// named, so ~17M deep compares finding ZERO dups -> ~20 min just to compile).
// Nodes are uniquely identified by "type:name" (the WCT NamedFactory identity)
// and edges by their serialized form, so we dedup by a cheap string key while
// keeping the FIRST occurrence -- order matters: WCT configures components in
// array order and resolve_uses() already emits dependencies before dependents
// (e.g. AnodePlane before OmnibusSigProc).  [The proper fix belongs in pgraph's
// g.uses/g.edges upstream.]
local unique_by(lst, keyf) =
  local n = std.length(lst);
  if n <= 1 then lst
  else
    local pairs = std.makeArray(n, function(i) { k: keyf(lst[i]), i: i });
    local sorted = std.sort(pairs, function(p) [p.k, p.i]);
    local keep = [ sorted[j].i for j in std.range(0, n - 1)
                   if j == 0 || sorted[j].k != sorted[j - 1].k ];
    [ lst[i] for i in std.sort(keep) ];
local uses = unique_by(g.resolve_uses(graph.uses),
                       function(n) std.get(n, 'type', '') + ':' + std.get(n, 'name', ''));
local edges = unique_by(graph.edges, function(e) std.manifestJson(e));

// The graph-execution engine is fcl-selectable: 'TbbFlow' (multi-threaded TBB
// flow graph) or 'Pgrapher' (single-threaded serial executor).  Both consume
// the same edge list.
local engine = std.extVar('engine');
local app = {
  type: engine,
  data: { edges: edges },
};
// For TbbFlow only, a TbbDataFlowGraph carries the thread limit: max_threads
// feeds tbb::global_control(max_allowed_parallelism); 0 leaves TBB free to use
// all cores (the historical default).  Pgrapher ignores it (always serial).
local dfg = if engine == 'TbbFlow'
  then [{ type: 'TbbDataFlowGraph', data: { max_threads: std.extVar('max_threads') } }]
  else [];

uses + dfg + [app]
