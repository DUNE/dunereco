// This provides multiple MagnifySink for e.g. protoDUNE

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// multiple MagnifySink
// tagn (n = 0, 1, ... 5) for anode[n]
// FrameFanin tags configured in sim.jsonnet

function(tag, tools, outputfile) {

  local nanodes = std.length(tools.anodes),


  local multimagnify = [
    g.pnode({
      type: 'MagnifySink',
      name: 'mag%s%d' % [tag, n],  // traces from source have NO tag
      data: {
        output_filename: outputfile,
        root_file_mode: 'RECREATE',
        frames: ['%s%d' % [tag, n]],
        trace_has_tag: false 
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],

  return: {
    magnify_pipelines: [g.pipeline([multimagnify[n]], name='magnifypipes%d' % n) for n in std.range(0, nanodes - 1)],
  },


}.return
