
local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local util = import 'pgrapher/experiment/pdhd/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';

local params = import 'pgrapher/experiment/pdhd/simparams.jsonnet';

// local tools = tools_maker(params);
local btools = tools_maker(params);
local anode_index = std.extVar('anode_index');
assert (anode_index == -1 || (anode_index >= 0 && anode_index < std.length(btools.anodes))) ||
    error ("anode_index extVar must be -1 (all anodes) or less than number of anodes (%d)." +
           "\nThis can be changed by changing the wcls_main.structs.anode_index value in your fcl file.") % std.length(btools.anodes);
local tools = if anode_index == -1 then btools else
    btools {
        anodes: [btools.anodes[anode_index], ],
    };

local sim_maker = import 'pgrapher/experiment/pdhd/sim.jsonnet';
local sim = sim_maker(params, tools);
local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes-1);
local anode_indices = if anode_index == -1 then [n for n in std.range(0, nanodes-1)] else [anode_index];

local wcls_maker = import "pgrapher/ui/wcls/nodes.jsonnet";
local wcls = wcls_maker(params, tools);
local wcls_input = {
    depos: wcls.input.depos(name="", art_tag="IonAndScint"),
};

local drifter = sim.drifter;
local bagger = [sim.make_bagger("bagger%d"%n) for n in anode_iota];

local depo_fanout_1st = g.pnode({
    type:'DepoFanout',
    name:'depo_fanout_1st',
    data:{
        multiplicity:nanodes,
        tags: [],
    }}, nin=1, nout=nanodes);

// N sinks

local sinks = [g.pnode({
    name: "deposink%02d" % n,
    type: "DepoFileSink",
    data: {
        outname: "pdhd-depos-out-%d.tar.bz2" % anode_indices[n],
    }
}, nin=1, nout=0) for n in anode_iota];

// g4 sim as input
local graph = g.intern(
    innodes=[wcls_input.depos], centernodes=[drifter, depo_fanout_1st] + bagger, outnodes=sinks,
    edges = 
        [
            g.edge(wcls_input.depos, drifter, 0, 0),
            g.edge(drifter, depo_fanout_1st, 0, 0),
        ] + [
            g.edge(depo_fanout_1st, bagger[n],  n, 0) for n in anode_iota
        ] + [
            g.edge(bagger[n], sinks[n], 0, 0) for n in anode_iota
        ],
);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};
g.uses(graph) + [app]
