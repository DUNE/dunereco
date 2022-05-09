// Utility helper functions.

local pg = import "pgraph.jsonnet";

{
    // Return an ident as a string, eg for an anode
    idents(obj) :: std.toString(obj.data.ident),

    // Return ADC counts per voltage given params.adc object
    adcpermv(adc) :: 
        ((1 << adc.resolution)-1) / (adc.fullscale[1]-adc.fullscale[0]),


    local basic_plugins = [
        "WireCellSio", "WireCellAux",
        "WireCellGen", "WireCellSigProc", "WireCellImg", 
        "WireCellApps"],
    
    local app_plugins = {
        'TbbFlow': ["WireCellTbb"],
        'Pgrapher': ["WireCellPgraph"],
    },

    main(graph, app='Pgrapher', extra_plugins = []) :: 
        local plugins = std.set(basic_plugins + extra_plugins + app_plugins[app]);
        local appcfg = {
            type: app,
            data: {
                edges: pg.edges(graph)
            },
        };
        local cmdline = {
            type: "wire-cell",
            data: {
                plugins: plugins,
                apps: [appcfg.type]
            }
        };
        [cmdline] + pg.uses(graph) + [appcfg],

}


