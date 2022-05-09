// Help build imaging nodes.  This provides a function which returns
// an object of methods and (sub)objects on the given anode.

local pg = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(anode)
{
    local apaid = anode.data.ident,

    // The conventional slicing node 
    slicing :: function(name=apaid, tag="", span=4) pg.pnode({
        type: "SumSlices",
        name: name,
        data: {
            tag: tag,
            tick_span: span,
            anode: wc.tn(anode),
        },
    }, nin=1, nout=1, uses=[anode]),

    // A function that tries to figure out if anode is wrapped or if
    // single faced which face number is sensitive.
    tiling :: function(name=apaid)
        if std.length(std.prune(anode.data.faces)) == 2 then
            $.wrapped_tiling(name)
        else if std.type(anode.data.faces[0]) == "null" then
            $.single_tiling(name, 1)
        else
            $.single_tiling(name, 0),

    // A function that returns a tiling node for single faced anodes
    single_tiling :: function(name=apaid, face=0) pg.pnode({
        type: "GridTiling",
        name: "%s-%d"%[name, face],
        data: {
            anode: wc.tn(anode),
            face: face,
        }
    }, nin=1, nout=1, uses=[anode]),

    // A function providing a subgraph to do tiling for wrapped (two
    // faced) anodes
    wrapped_tiling :: function(name=apaid) {

        local slice_fanout = pg.pnode({
            type: "SliceFanout",
            name: name,
            data: { multiplicity: 2 },
        }, nin=1, nout=2),

        local tilings = [
            pg.pnode({
                type: "GridTiling",
                name: "%s-%d"%[name, face],
                data: {
                    anode: wc.tn(anode),
                    face: face,
                }
            }, nin=1, nout=1, uses=[anode]) for face in [0,1]],

        local blobsync = pg.pnode({
            type: "BlobSetSync",
            name: name,
            data: { multiplicity: 2 }
        }, nin=2, nout=1),

        ret: pg.intern(
            innodes=[slice_fanout],
            outnodes=[blobsync],
            centernodes=tilings,
            edges=
                [pg.edge(slice_fanout, tilings[n], n, 0) for n in [0,1]] +
                [pg.edge(tilings[n], blobsync, 0, n) for n in [0,1]],
            name=name),
    }.ret,

    // A node that clusters blobs between slices w/in the given number
    // of slice spans.  Makes blob-blob edges
    clustering :: function(name=apaid, spans=1.0) pg.pnode({
        type: "BlobClustering",
        name: name,
        data:  { spans : spans }
    }, nin=1, nout=1),

    // A node that groups wires+channels into measurements.  Makes
    // blob-measure edges.
    grouping :: function(name=apaid) pg.pnode({
        type: "BlobGrouping",
        name: name,
        data:  {
        }
    }, nin=1, nout=1),

    // A node that "solves" for blob charge given blob-measure edges.
    // Only blobs with charge above threshold are kept.  
    charge_solving :: function(name=apaid,
                               meas_val_thresh=10.0,
                               meas_err_thresh=1.0e9,
                               blob_val_thresh=0,
                               // fixme; there are more with
                               // reasonable defaults that could be
                               // filled in.
    )
        pg.pnode({
            type: "ChargeSolving",
            name: name,
            data:  {
                meas_value_threshold: meas_val_thresh,
                meas_error_threshold: meas_err_thresh,
                blob_value_threshold: blob_val_thresh,
            }
        }, nin=1, nout=1),

    // This uses a simpler, less effective solving algorithm.  Use
    // charge_solving for the big guns.
    blob_solving :: function(name=apaid, threshold=0.0) pg.pnode({
        type: "BlobSolving",
        name: name,
        data:  { threshold: threshold }
    }, nin=1, nout=1),

    // A nominal imaging omnibus pipeline starting with a frame ending
    // in solved cluster.  Note, there are many ways to reasonably
    // form a full chain.  This is but one which implements the graph
    // in the raygrid.pdf figure 12 but starting with IFrame input and
    // ending with a solving.
    simple :: function(name=apaid, tag="", span=4, spans=1.0, threshold=0.0) pg.pipeline([
        $.slicing(name, tag, span),
        $.tiling(name),
        $.clustering(name, spans),
        $.grouping(name),
        $.blob_solving(name, threshold)
    ]),

    

    // A function that projects blobs back to frames.  
    reframing :: function(name=apaid, tag="") pg.pnode({
        type: "BlobReframer",
        name: name,
        data: {
            frame_tag: tag,
        }
    }, nin=1, nout=1),

}
