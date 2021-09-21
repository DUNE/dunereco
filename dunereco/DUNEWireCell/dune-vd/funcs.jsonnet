// We define here some fuunctions for the two-layers fans
// To do generalization to custom number of layers

local g = import "pgraph.jsonnet";

{

		multifanpipe :: function( fout, pipelines, fin, fanmult=6, name='twolayerfanpipe', outtags=[], tag_rules=[] ) {

			local idents = std.range(0, fanmult-1),
			local identsquare = std.range(0, fanmult*fanmult-1),


			// Two layer fanout

			local layer1 = g.pnode({
    			type: fout,
    			name: "layer1",
          data: {
            multiplicity: fanmult, 
            tag_rules: tag_rules,
          },
			}, nin=1, nout=fanmult),

			local layer2s = [g.pnode({
    			type: fout,
    			name: "layer2_%02d" % n,
          data: {
            multiplicity: fanmult, 
            tag_rules: tag_rules,
          },
			}, nin=1, nout=fanmult) for n in idents],

			local bigfanout = g.intern(innodes = [layer1],
            	outnodes = layer2s,
            	edges = [g.edge(layer1, layer2s[n], n, 0) for n in idents],
            name="bigfanout"),


      // Two layer fanin

      local layer3 = g.pnode({
    			type: fin,
    			name: "layer3",
          data: {
            multiplicity: fanmult, 
            tags: outtags,
          },
			}, nin=fanmult, nout=1),

			local layer4 = [g.pnode({
    			type: fin,
    			name: "layer4_%02d" % n,
          data: {
            multiplicity: fanmult, 
            tags: outtags,
          },
			}, nin=fanmult, nout=1) for n in idents],


			local bigfanin = g.intern( innodes=layer4,
                outnodes=[layer3],
                edges=[g.edge(layer4[n], layer3, 0, n) for n in idents],
                name="bigfanin"),


            // Pipeline

            ret: g.intern(  innodes=[bigfanout],
                      		centernodes=pipelines,
                      		outnodes = [bigfanin],
                      		edges = [g.edge(bigfanout, pipelines[n], n, 0) for n in identsquare] +
                          	[g.edge(pipelines[n], bigfanin, 0, n) for n in identsquare],
                      		name=name,
                    	),

		}.ret, 
}