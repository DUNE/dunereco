// Inspired to the PCBRP cli-sim-npz.jsonnet

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params = import 'pgrapher/experiment/dune-vd/simparams.jsonnet';

local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/dune-vd/sim.jsonnet';
local sim = sim_maker(params, tools);



// LOAD THE DETECTOR




local tracklist = [
  {
    time: 0,
    charge: -5000,
    ray: {
        tail: wc.point(20, -10, -10, wc.cm), // << TO EDIT
        head: wc.point(30, +10, +10, wc.cm), // << TO EDIT
    }
  },
];


local depos = sim.tracks(tracklist);
local drifter = sim.drifter;
local bagger = sim.make_bagger();








//local graph = g.pipeline([depos, drifter, bagger ]);



// End commands 

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot"],
        apps: ["Pgrapher"]
    }
};

[cmdline] + g.uses(graph) + [app]