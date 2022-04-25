local pg = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(anode)
{
    local apaid = anode.data.ident,

    // fill ROOT histograms with frames
    magnify :: function(name, filename="magnify-"+apaid+".root", frame_tags=["orig"+apaid]) ret: g.pnode({
        type: 'MagnifySink',
        name: name,
        data: {
            output_filename: filename,
            root_file_mode: 'UPDATE',
            frames: frame_tags,
            trace_has_tag: true,
            anode: wc.tn(anode),
        },
    }, nin=1, nout=1),

}
    
