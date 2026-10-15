// High layer / entry point for the dune10kt-vd WC/LS sim-drift-simchannel-nf-sp
// job.  It is named to match the hosting FHiCL (wcls-sim-drift-simchannel-nf-sp.fcl).
//
// Responsibilities of THIS layer only:
//   1. resolve all std.extVar() fcl knobs into one `fcl_params` object,
//   2. build `params` (params-10kt) and `tools` (narrowed to the anodes to run),
//   3. dispatch on process_mode to a modal layer that builds the graph:
//        - "multi-sim"            -> -multi.jsonnet  (per-TPC fan-out, per-TPC
//                                     FrameSavers; process_tpc_indices list)
//        - "single-sim/-sp/-sim-sp" and the "full" fallback -> -single.jsonnet
//   4. wrap the graph in the flow-engine app and emit.
//
// Shared detector glue lives in -common.jsonnet; the modal graph assembly lives
// in -single.jsonnet / -multi.jsonnet.  See those files for details.
//
// NOTE on extVar laziness: per-mode knobs (process_tpc_index vs
// process_tpc_indices / frame_tag_pattern) are referenced only on the code path
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
local fcl_params = {
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
  // single-* : which one TPC to process
  process_tpc_index: std.extVar('process_tpc_index'),
  // multi-*  : list of TPC ids, and the per-TPC RawDigit tag/instance pattern
  process_tpc_indices: std.extVar('process_tpc_indices'),
  frame_tag_pattern: std.extVar('frame_tag_pattern'),
};

// ---- 2. params + tools ---------------------------------------------------
local params = params_maker(fcl_params) {
  lar: super.lar {
    DL: std.extVar('DL') * wc.cm2 / wc.ns,
    DT: std.extVar('DT') * wc.cm2 / wc.ns,
    lifetime: std.extVar('lifetime') * wc.us,
    drift_speed: util.drift_velocity(std.extVar('efield'), std.extVar('temperature')) * wc.mm / wc.us,
  },
  files: super.files {
    wires: std.extVar('files_wires'),
    // same field response for [bottom, top] drift volumes
    fields: [std.extVar('files_fields'), std.extVar('files_fields')],
  },
};

local tools_all = tools_maker(params);

local is_single =
  fcl_params.process_mode == 'single-sim'
  || fcl_params.process_mode == 'single-sp'
  || fcl_params.process_mode == 'single-sim-sp';
local is_multi = fcl_params.process_mode == 'multi-sim';

// The anodes this job will process, as positions into tools_all.anodes.
local anode_indices =
  if is_multi then fcl_params.process_tpc_indices
  else if is_single then [fcl_params.process_tpc_index]
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
local app = {
  type: 'TbbFlow',  // Pgrapher, TbbFlow
  data: { edges: g.edges(graph) },
};

g.uses(graph) + [app]
