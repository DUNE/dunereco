// Helpers to make configuration objects for file input/output
// components.

local pg = import "pgraph.jsonnet";
local wc = import "wirecell.jsonnet";

{
    // If a frame is written out "dense" it the sink needs to know
    // what the bounds are in channel IDs (rows) and ticks (columns).
    frame_bounds(nchannels, nticks, first_chid=0, first_tick=0) :: {
        chbeg: first_chid, chend: first_chid+nchannels,
        tbbeg: first_tick, tbend: first_tick+nticks
    },

    /// Sink a stream of frames to file. 
    frame_file_sink(filename, tags=[],
                    digitize=false, dense=null) :: 
        pg.pnode({
            type: "FrameFileSink",
            name: filename,
            data: {
                outname: filename,
                tags: tags,
                digitize: digitize,
                dense: dense,
            },
        }, nin=1, nout=0),
        
    /// Source a stream of frames from file.
    frame_file_source(filename, tags=[]) ::
        pg.pnode({
            type: "FrameFileSource",
            name: filename,
            data: {
                inname: filename,
                tags: tags,
            },
        }, nin=0, nout=1),


    // Like a frame_file_sink but pass input frames both to output
    // port 0 and to a sink.
    frame_file_tap(filename, tags=[], digitize=false, dense=null) ::
        pg.fan.tap('FrameFanout', 
                   $.frame_file_sink(filename, 
                                     tags=tags,
                                     digitize=digitize,
                                     dense=dense), filename),

    // A pass-through node for ICluster which saves as a side-effect.
    cluster_json_tap :: function(name, drift_speed=1.6*wc.mm/wc.us, filepat=null) pg.pnode({
        type: "JsonClusterTap",
        name: name,
        data: {
            filename: if std.type(filepat) == "null" then "clusters-"+name+"-%04d.json" else filepat,
            drift_speed: drift_speed
        },
    }, nin=1, nout=1),

    // Write a cluster to a graphviz file
    cluster_graphviz :: function(name, filepat=null) pg.pnode({
        type: "ClusterSink",
        name: name,
        data: {
            filename: if std.type(filepat) == "null" then "clusters-apa-"+name+"-%04d.dot" else filepat,
        }
    }, nin=1, nout=0),

    // Write out in a WCT cluster file format (json+tar[+compression])
    cluster_file_sink :: function(name, filename=null) pg.pnode({
        type: 'ClusterFileSink',
        name: name,
        data: {
            outname: if std.type(filename) == "null" then "clusters-"+name+".tar.bz2" else filename,
        },
    }, nin=1, nout=0),


}
