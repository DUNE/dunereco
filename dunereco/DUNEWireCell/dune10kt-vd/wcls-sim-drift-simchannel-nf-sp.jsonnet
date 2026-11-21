// Legacy entry point for the dune10kt-vd WC/LS sim + sigproc job.
//
// This name is what the wirecell_dune.fcl blocks dune10kt_vertdrift_sim_nfsp
// and dune10kt_vertdrift_sim_nfsp_singledigits give in their configs (and so
// what dunesw standard_detsim_dunevd10kt.fcl loads).  It reads the std.extVar()
// knobs of those blocks and calls wcls-fdvd-sim-sigproc.jsonnet with them, in
// one of its merged-product modes: "full" or "single-sim/-sp/-sim-sp".
// "multi-sim" is not available here: use wcls-fdvd-sim-sigproc.jsonnet directly.
//
// Knobs not given by those blocks take their legacy values here.  In
// particular the WCT wires come from the files_wires file and the drift planes
// from params-10kt.jsonnet: the wire_schema and drift_planes params of the
// blocks are not read, as "service" would need the larwirecell
// WireSchemaService, which jobs of this legacy configuration do not configure.

local entry = import 'pgrapher/experiment/dune10kt-vd/wcls-fdvd-sim-sigproc.jsonnet';

local process_mode = std.extVar('process_mode');

entry({
  // from the FHiCL block: params
  files_wires: std.extVar('files_wires'),
  files_fields: std.extVar('files_fields'),
  geo_planeid_labels: std.extVar('geo_planeid_labels'),
  process_mode:
    if process_mode == 'multi-sim'
    then error 'process_mode "multi-sim" needs wcls-fdvd-sim-sigproc.jsonnet'
    else process_mode,

  // from the FHiCL block: structs
  elecGain: std.extVar('elecGain'),
  nticks: std.extVar('nticks'),
  DL: std.extVar('DL'),
  DT: std.extVar('DT'),
  lifetime: std.extVar('lifetime'),
  efield: std.extVar('efield'),
  temperature: std.extVar('temperature'),
  G4RefTime: std.extVar('G4RefTime'),
  response_plane: std.extVar('response_plane'),
  process_anode_index: std.extVar('process_anode_index'),
  use_hydra: std.extVar('use_hydra'),
  save_rawdigits: std.extVar('save_rawdigits'),
  use_dnnroi: std.extVar('use_dnnroi'),
  adc_resolution: std.extVar('adc_resolution'),

  // legacy values
  wire_schema: 'file',
  drift_planes: 'params',
  input_depos_tag: 'IonAndScint',
  engine: 'TbbFlow',  // matches the blocks' apps
  max_threads: 0,  // TBB free to use all cores
  process_anode_indices: [],  // multi-sim only
  frame_tag_pattern: 'daq%d',  // multi-sim only
})
