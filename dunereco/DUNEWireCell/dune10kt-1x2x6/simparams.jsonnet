// Here we override params.jsonnet to provide simulation-specific params.

local base = import 'pgrapher/experiment/dune10kt-1x2x6/params.jsonnet';
local wc = import 'wirecell.jsonnet';

function(params) base(params) {
  // lar: super.lar {
  //     // Longitudinal diffusion constant
  //     DL : 4.0 * wc.cm2/wc.s,
  //     // Transverse diffusion constant
  //     DT : 8.8 * wc.cm2/wc.s,
  //     // Electron lifetime
  //     lifetime : 35*wc.ms,
  //     // Electron drift speed, assumes a certain applied E-field
  //     drift_speed : 1.565*wc.mm/wc.us,
  // },

  // redefine the detector volumes with the cryostat side included
  det : {

      // The "faces" is consumed by, at least, the Drifter and
      // AnodePlane.  The "wires" number is used to set
      // AnodePlane.ident used to lookup the anode in WireSchema.
      // It corresponds to the anode number.

      // between center lines
      // local apa_cpa = 3.637*wc.m, // DocDB 203
      local apa_cpa = 3.63075*wc.m, // LArSoft
      // local cpa_thick = 50.8*wc.mm, // DocDB 203
      local cpa_thick = 3.175*wc.mm, // 1/8", from Bo Yu (BNL) and confirmed with LArSoft

      // X positions from dune10kt gdml:
      // U: 39.5355 mm, V: 34.7755 mm, W: 30.0155 mm
      local apa_w2w = 60.031*wc.mm,
      local plane_gap = 4.76*wc.mm,
      local apa_g2g = apa_w2w + 6*plane_gap, // 88.591*wc.mm,

      // The "anode" cut off plane, here measured from APA
      // centerline, determines how close to the wires do we
      // consider any depo.  Anything closer will simply be
      // discarded, else it will either be drifted or "backed up" to
      // the response plane.  This is somewhat arbitrary choice.
      // Placing it w/in the response plane means any depos that are
      // "backed up" won't have proper field response.  But, the
      // tighter this is made, the less volume is simulated.
      // local apa_plane = 0.5*apa_g2g, // pick it to be at the grid wires
      local apa_plane = 0.5*apa_g2g - plane_gap, // pick it to be at the first induction wires

      // The "response" plane is where the field response functions
      // start.  Garfield calcualtions start somewhere relative to
      // something, here's where that is made concrete.  This MUST
      // match what field response functions also used.
      response_plane: 10*wc.cm, // relative to collection wires
      local res_plane = 0.5*apa_w2w + self.response_plane,

      // The cathode plane is like the anode cut off plane.  Any
      // depo not between the two is dropped prior to drifting.
      local cpa_plane = apa_cpa - 0.5*cpa_thick,


      // The volumes are then defined in terms of these above
      // numbers.  You can use "wirecell-util wires-info" or
      // "wirecell-util wires-volumes" or others to understand the
      // mapping of anode number to the 6 locations in X and Z.  For
      // Larsoft wires the numbering is column major starting at
      // small X and Z so the centerline is -/+/-/+/-/+.  Also
      // important is that the faces are listed "front" first.
      // Front is the one with the more positive X coordinates and
      // if we want to ignore a face it is made null.
      volumes: [
          {
              local sign = 2*(n%2)-1,
              local centerline = 0, // sign*apa_cpa,
              wires: n,       // anode number
              name: "apa%d"%n,
              faces:
              // top, front face is against cryo wall
              if sign > 0
              then [
                  {
                      anode: centerline + apa_plane,
                      response: centerline + res_plane,
                      cathode: centerline + cpa_plane, 
                  },
                  {
                      anode: centerline - apa_plane,
                      response: centerline - res_plane,
                      cathode: centerline - cpa_plane, 
                  }
              ]
              // bottom, back face is against cryo wall
              else [
                  {
                      anode: centerline + apa_plane,
                      response: centerline + res_plane,
                      cathode: centerline + cpa_plane, 
                  },
                  {
                      anode: centerline - apa_plane,
                      response: centerline - res_plane,
                      cathode: centerline - cpa_plane, 
                  }
              ],
          } for n in std.range(0,11)],

      // This describes some rough, overall bounding box.  It's not
      // directly needed but can be useful on the Jsonnet side, for
      // example when defining some simple kinematics.  It is
      // represented by a ray going from extreme corners of a
      // rectangular solid.  Again "wirecell-util wires-info" helps
      // to choose something.
      bounds : {
          tail: wc.point(-363.376, -600.019, -0.87625, wc.cm),
          head: wc.point( 363.376,  600.019, 1393.46, wc.cm),
      }
  },

  daq: super.daq {

    // Number of readout ticks.  See also sim.response.nticks.
    // In MB LArSoft simulation, they expect a different number of
    // ticks than acutal data.
    nticks: 6000,
  },

  // These parameters only make sense for running WCT simulation on
  // microboone in larsoft.  The "trigger" section is not
  // "standard".  This section is just a set up for use below in
  // "sim".  There is no trigger, per se, in the simulation but
  // rather a contract between the generators of energy depositions
  // (ie, LarG4) and the drift and induction simulation (WCT).  For
  // details of this contract see:
  // https://microboone-docdb.fnal.gov/cgi-bin/private/ShowDocument?docid=12290
  //trigger : {

  //    // A hardware trigger occurs at some "absolute" time but near
  //    // 0.0 for every "event".  It is measured in "simulation" time
  //    // which is the same clock used for expressing energy
  //    // deposition times.  The following values are from table 3 of
  //    // DocDB 12290.
  //    hardware: {
  //        times: {

  //            none: 0.0,

  //            // BNB hardware trigger time.  Note interactions
  //            // associated with BNB neutrinos should all be produced
  //            // starting in the beam gate which begins at 3125ns and is
  //            // 1600ns in duration.
  //            bnb : -31.25*wc.ns,


  //            // Like above but for NUMI.  It's gate is 9600ns long starting
  //            // at 4687.5ns.
  //            numi : -54.8675*wc.ns,

  //            ext : -414.0625*wc.ns,
  //
  //            mucs: -405.25*wc.ns,
  //        },

  //        // Select on of the trigger types
  //        type: "bnb",

  //        time: self.times[self.type],
  //    },

  //    // Measured relative to the hardware trigger time above is a
  //    // time offset to the time that the first tick of the readout
  //    // should sample.  This is apparently fixed for all hardware
  //    // trigger types (?).
  //    time_offset: -1.6*wc.ms,

  //    time: self.hardware.time + self.time_offset,
  //},

  sim: super.sim {

    // For running in LArSoft, the simulation must be in fixed time mode.
    fixed: true,
    continuous: false,
    fluctuate: true,

    //ductor : super.ductor {
    //    start_time: $.daq.start_time - $.elec.fields.drift_dt + $.trigger.time,
    //},


    // Additional e.g. 10 us time difference is due to the larger drift velocity
    // in Garfield field response where the collection plane peak
    // at around 81 us instead of response_plane (10 cm to Y plane) /drift_speed.
    // Assuming a constant drift velocity, this correction is needed.
    // Interplane timeoffset still holds and will be intrinsically taken into account
    // in the 2D decon.
    // ATTENTION: when response variation (sys_status) turned on, an offset is needed.
    // smearing function is centralized at t=0 instead of starting from t=0
    //reframer: super.reframer{
    //    tbin: if $.sys_status == true
    //            then (81*wc.us-($.sys_resp.start))/($.daq.tick)
    //            else (81*wc.us)/($.daq.tick),
    //    nticks: $.daq.nticks,
    //    toffset: if $.sys_status == true
    //                then $.elec.fields.drift_dt - 81*wc.us + $.sys_resp.start
    //                else $.elec.fields.drift_dt - 81*wc.us,
    //},

  },
  // This is a non-standard, MB-specific variable.  Each object
  // attribute holds an array of regions corresponding to a
  // particular set of field response functions.  A region is
  // defined as an array of trios: plane, min and max wire index.
  // Each trio defines a swath in the transverse plane bounded by
  // the min/max wires.  A region is finally the intersection or
  // overlap of all its trios in the transverse plane.
  //shorted_regions : {
  //    uv: [
  //        [ { plane:0, min:296, max:296 } ],
  //        [ { plane:0, min:298, max:315 } ],
  //        [ { plane:0, min:317, max:317 } ],
  //        [ { plane:0, min:319, max:327 } ],
  //        [ { plane:0, min:336, max:337 } ],
  //        [ { plane:0, min:343, max:345 } ],
  //        [ { plane:0, min:348, max:351 } ],
  //        [ { plane:0, min:376, max:400 } ],
  //        [ { plane:0, min:410, max:445 } ],
  //        [ { plane:0, min:447, max:484 } ],
  //        [ { plane:0, min:501, max:503 } ],
  //        [ { plane:0, min:505, max:520 } ],
  //        [ { plane:0, min:522, max:524 } ],
  //        [ { plane:0, min:536, max:559 } ],
  //        [ { plane:0, min:561, max:592 } ],
  //        [ { plane:0, min:595, max:598 } ],
  //        [ { plane:0, min:600, max:632 } ],
  //        [ { plane:0, min:634, max:652 } ],
  //        [ { plane:0, min:654, max:654 } ],
  //        [ { plane:0, min:656, max:671 } ],
  //    ],
  //    vy: [
  //        [ { plane:2, min:2336, max:2399 } ],
  //        [ { plane:2, min:2401, max:2414 } ],
  //        [ { plane:2, min:2416, max:2463 } ],
  //    ],
  //},

  //files: super.files{
  //    chresp: null,
  //},

  // This sets a relative gain at the input to the ADC.  Note, if
  // you are looking to fix SimDepoSource, you are in the wrong
  // place.  See the "scale" parameter of wcls.input.depos() defined
  // in pgrapher/common/ui/wcls/nodes.jsonnet.
  // elec: super.elec {
  //   postgain: 0.858,
  //   shaping: 2.2 * wc.us,
  // },

  sys_status: false,
  sys_resp: {
    // overall_short_padding should take into account this offset "start".
    start: -10 * wc.us,
    magnitude: 1.0,
    time_smear: 1.0 * wc.us,
  },

  rc_resp: {
    width: 1.1*wc.ms,
    rc_layers: 0, // 1
  }
}
