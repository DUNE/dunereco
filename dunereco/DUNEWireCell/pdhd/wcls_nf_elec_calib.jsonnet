// PDHD WC/LS entry point for ELECTRONICS-RESPONSE CALIBRATION (validation).
//
// Copy of wcls-nf.jsonnet that:
//   * injects the fitted per-channel response JSON via params.files.chresp so
//     tools.perchanresp -> ParamsPerChannelResponse is instantiated, and
//   * uses nf_elec_calib.jsonnet (which adds the perChannelShaper to the NF)
//     instead of nf.jsonnet, passing `tools` into the maker.
// The saved RawDigit (frame tag "raw", art product "wclsdatahdfilter:raw") is
// the per-channel response-equalized waveform, to be re-fit with the ideal
// cold-electronics response by the downstream fitter.

local reality = std.extVar('reality');
local sigoutform = std.extVar('signal_output_form');  // eg "sparse" or "dense"


local wc = import 'wirecell.jsonnet';
local g = import 'pgraph.jsonnet';

local raw_input_label = std.extVar('raw_input_label');  // eg "daq"

// Fitted per-channel response file (basename resolved on WIRECELL_PATH).
local chresp_file = 'protodunehd-params-channel-response-30413-no_avg_3.json.bz2';

local data_params = import 'pgrapher/experiment/pdhd/params.jsonnet';
local simu_params = import 'pgrapher/experiment/pdhd/simparams.jsonnet';
local base = if reality == 'data' then data_params else simu_params;
local params = base {
    daq: super.daq {
      tick: 1.0/std.extVar('clock_speed') * wc.us,
    },
    files: super.files {
      chresp: chresp_file,   // <-- activate per-channel response correction
    },
};

local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);

local sp_maker = import 'pgrapher/experiment/pdhd/sp.jsonnet';

local use_resampler = (reality == 'data');

// Collect the WC/LS input converters for use below.
local wcls_input = {
  adc_digits: g.pnode({
    type: 'wclsRawFrameSource',
    name: '',
    data: {
      art_tag: raw_input_label,
      frame_tags: ['orig'],  // this is a WCT designator
      tick: 512*wc.ns, //Use 512ns here for input, we resample to 500ns later
    },
  }, nin=0, nout=1),

};

local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
  },
};
local wcls_output = {
  // The noise-filtered + response-equalized "ADC" values, saved as RawDigit.
  nf_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'nfsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['raw'],
      chanmaskmaps: ['bad'],
      nticks: 0,
    },
  }, nin=1, nout=1, uses=[mega_anode]),
};

local base_chndb = import 'chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: base_chndb(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in std.range(0, std.length(tools.anodes) - 1)];

local nf_maker = import 'pgrapher/experiment/pdhd/nf_elec_calib.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], tools, chndb[n], n, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];

local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % n,
    data: {
      channels: std.range(2560 * n, 2560 * (n + 1) - 1),
    },
  }, nin=1, nout=1)
  for n in std.range(0, std.length(tools.anodes) - 1)
];


local resamplers_config = import 'pgrapher/common/resamplers.jsonnet';
local load_resamplers = resamplers_config(g, wc, tools);
local resamplers = load_resamplers.resamplers;

local nfsp_pipes = [
  g.pipeline(
    [chsel_pipes[n]] +
    (if use_resampler then [resamplers[n]] else []) +
    [
     nf_pipes[n],
    ],
    'nfsp_pipe_%d' % n)
    for n in std.range(0, std.length(tools.anodes) - 1)
];

local f = import 'pgrapher/experiment/pdhd/funcs.jsonnet';
local fanpipe = f.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'sn_mag_nf');

local retagger = g.pnode({
  type: 'Retagger',
  data: {
    tag_rules: [{
      frame: {
        '.*': 'retagger',
      },
      merge: {
        'raw\\d': 'raw',
      },
    }],
  },
}, nin=1, nout=1);

local sink = g.pnode({ type: 'DumpFrames' }, nin=1, nout=0);

local graph = g.pipeline([wcls_input.adc_digits, fanpipe, retagger, wcls_output.nf_digits, sink]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

// Finally, the configuration sequence
g.uses(graph) + [app]
