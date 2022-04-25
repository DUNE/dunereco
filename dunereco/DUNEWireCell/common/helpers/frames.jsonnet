local pg = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(anode)
{
    local apaid = anode.data.ident,

    // A sink of frames which dumps info to log
    logframes :: function(name=apaid) pg.pnode({
        type: "DumpFrames",
        name: "name",
    }, nin=1, nout=0),
}

