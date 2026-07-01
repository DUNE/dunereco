// PDVD WC/LS entry point for ELECTRONICS-RESPONSE CALIBRATION (validation).
//
// Copy of wcls-nf.jsonnet that:
//   * injects the fitted per-channel response JSON via params.files.chresp so
//     tools.perchanresp -> ParamsPerChannelResponse is instantiated, and
//   * uses nf_elec_calib.jsonnet (which adds the perChannelShaper to the NF)
//     instead of nf.jsonnet, passing `tools` into the maker.
// The saved RawDigit (frame tag "raw", art product "wclsdatavdfilter:raw") is
// the per-channel response-equalized waveform, to be re-fit with the ideal
// cold-electronics response by the downstream fitter.

local reality = std.extVar('reality');
local sigoutform = std.extVar('signal_output_form');  // eg "sparse" or "dense"


local wc = import 'wirecell.jsonnet';
local g = import 'pgraph.jsonnet';

local raw_input_label = std.extVar('raw_input_label');  // eg "daq"

// Fitted per-channel response file (basename resolved on WIRECELL_PATH).
local chresp_file = 'protodunevd-params-channel-response-41603-no_avg_3.json.bz2';

local data_params = import 'params.jsonnet';
local simu_params = import 'simparams.jsonnet';
local base = if reality == 'data' then data_params else simu_params;
local params = base {
    files: super.files {
      chresp: chresp_file,   // <-- activate per-channel response correction
    },
};


local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);

local sp_maker = import 'pgrapher/experiment/protodunevd/sp.jsonnet';

local use_resampler = std.extVar("use_resampler");

// This run is BOTTOM-DRIFTER only (anodes 0..3, channels 0..6143).  Restrict
// every anode-indexed structure (NF pipes, MegaAnodePlane, fan rules) to the
// bottom anodes so anodes 4..7 are never referenced/instantiated.
local bot_anodes = [n for n in std.range(0, std.length(tools.anodes) - 1)
                    if tools.anodes[n].data.ident < 4];

// Collect the WC/LS input converters for use below.
local wcls_input = {
  adc_digits: g.pnode({
    type: 'wclsRawFrameSource',
    name: '',
    data: {
      art_tag: raw_input_label,
      frame_tags: ['orig'],  // this is a WCT designator
      tick: 512*wc.ns,       // native bottom-drifter sampling; no resampler
    },
  }, nin=0, nout=1),

};

local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(tools.anodes[n]) for n in bot_anodes],
  },
};

local wcls_output = {
  // The ORIGINAL (input) "ADC" values, saved as RawDigit for before/after
  // comparison.  Tapped off right after the raw-frame source, before NF.
  orig_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'origsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['orig'],
      chanmaskmaps: [],
    },
  }, nin=1, nout=1, uses=[mega_anode]),

  // The noise-filtered + response-equalized "ADC" values, saved as RawDigit.
  nf_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'nfsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['raw'],
      chanmaskmaps: ['bad'],
    },
  }, nin=1, nout=1, uses=[mega_anode]),
};

local base_chndb = import 'chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: base_chndb(params, tools.anodes[n], tools.field, n) { dft:wc.tn(tools.dft) },
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in std.range(0, std.length(tools.anodes) - 1)];

local nf_maker = import 'pgrapher/experiment/protodunevd/nf_elec_calib.jsonnet';
// The FrameFanout (fanout_tag_rules) renames the per-anode FRAME tag to
// 'orig<ident>', so the OmnibusNoiseFilter in each pipe must read that tag.
// Its default intraces='orig' matches neither a trace tag nor the renamed
// frame tag -> Aux::tagged_traces returns empty -> ONF emits an empty frame
// -> the FrameSaver writes 0-sample RawDigits.  Pass the matching per-anode tag.
local nf_pipes = [nf_maker(params, tools.anodes[n], tools, chndb[n], n, name='nf%d' % n,
                          intraces='orig%d' % tools.anodes[n].data.ident) for n in bot_anodes];

local util = import 'pgrapher/experiment/protodunevd/funcs.jsonnet';
local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % n,
    data: {
      channels: util.anode_channels(n),
    },
  }, nin=1, nout=1)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

// No resampler: the shaper operates on the native 512 ns sampling.
local nfsp_pipes = [
  g.pipeline(
             [ chsel_pipes[n] ]
             + [ nf_pipes[std.find(n, bot_anodes)[0]] ],

             'nfsp_pipe_%d' % n)
  for n in bot_anodes
];

local fanout_tag_rules = [
          {
            frame: {
              '.*': 'orig%d' % tools.anodes[n].data.ident,
            },
            trace: {
            },
          }
          for n in bot_anodes
        ];

local anode_ident = [tools.anodes[n].data.ident for n in bot_anodes];
local fanin_tag_rules = [
          {
            frame: {
              '.*': 'framefanin',
            },
            trace: {
              ['raw%d'%ind]:'raw%d'%ind,
              ['gauss%d'%ind]:'gauss%d'%ind,
              ['wiener%d'%ind]:'wiener%d'%ind,
              ['threshold%d'%ind]:'threshold%d'%ind,
              ['loose_lf%d'%ind]:'loose_lf%d'%ind,
            },

          }
          for ind in anode_ident
        ];
local fanpipe = util.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'nfsp', [], fanout_tag_rules, fanin_tag_rules);


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
// NOTE: the inline origsaver (wcls_output.orig_digits) does NOT produce a
// populated product in this WCLS version -- an inline FrameSaver placed before
// the fanout writes 0 RawDigits regardless of frame/trace tagging.  The
// uncorrected before-image is instead read directly from "tpcrawdecoder:daq"
// (identical raw ADC), so origsaver is dropped from the graph.
local graph = g.pipeline([wcls_input.adc_digits, fanpipe, retagger, wcls_output.nf_digits, sink]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

// Finally, the configuration sequence
g.uses(graph) + [app]
