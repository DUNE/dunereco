// ProtoDUNE-SP specific parameters.  This file inerets from the
// generic set of parameters and overrides things specific to PDSP.

local base = import 'pgrapher/common/params_gain78.jsonnet';
local wc = import 'wirecell.jsonnet';

base {
  // This section will be overwritten in simparams.jsonnet
  det: {

    // define the 6 APAs.  This must use the coordinate system
    // defined by the wire geometry file.
    // A full drift is a box: xyz=[3.594*wc.m, 5.9*wc.m, 2.2944*wc.m].
    //
    // The "faces" is consumed by, at least, the Drifter and
    // AnodePlane.  The "wires" number is used to set
    // AnodePlane.ident used to lookup the anode in WireSchema.
    // It corresponds to the anode number.
    //
    // Also see:
    //   wirecell-util wire-volumes protodune-wires-larsoft-v3.json.bz2
    // to help with defining these parameters.

    // from DocDB 203 and assuming wires are symmetric across x=0

    // between center lines
    local apa_cpa = 3.637 * wc.m,
    local cpa_thick = 50.8 * wc.mm,
    local apa_w2w = 85.725 * wc.mm,  // DocDB 203 calls "W" as "X"
    local plane_gap = 4.76 * wc.mm,
    local apa_g2g = 114.3 * wc.mm,  // note that grid plane must have
    // gap 4.7675mm for this number
    // to be consistent with above.
    // There's probably round-off
    // error in DocDB 203.

    // The "anode" cut off plane, here measured from APA
    // centerline, determines how close to the wires do we
    // consider any depo.  Anything closer will simply be
    // discarded, else it will either be drifted or "backed up" to
    // the response plane.  This is somewhat arbitrary choice.
    // Placing it w/in the response plane means any depos that are
    // "backed up" won't have proper field response.  But, the
    // tighter this is made, the less volume is simulated.
    local apa_plane = 0.5 * apa_g2g,  // pick it to be at the grid wires

    // The "response" plane is where the field response functions
    // start.  Garfield calcualtions start somewhere relative to
    // something, here's where that is made concrete.  This MUST
    // match what field response functions also used.
    response_plane: 10 * wc.cm,  // relative to collection wires
    local res_plane = 0.5 * apa_w2w + self.response_plane,

    // The cathode plane is like the anode cut off plane.  Any
    // depo not between the two is dropped prior to drifting.
    local cpa_plane = apa_cpa - 0.5 * cpa_thick,


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
        local sign = 2 * (n % 2) - 1,
        local centerline = sign * apa_cpa,
        wires: n,  // anode number
        name: 'apa%d' % n,
        faces:
          // top, front face is against cryo wall
          if sign > 0
          then [
            null,
            {
              anode: centerline - apa_plane,
              response: centerline - res_plane,
              cathode: centerline - cpa_plane,
            },
          ]
          // bottom, back face is against cryo wall
          else [
            {
              anode: centerline + apa_plane,
              response: centerline + res_plane,
              cathode: centerline + cpa_plane,
            },
            null,
          ],
      }
      for n in std.range(0, 3)
    ],

    // This describes some rough, overall bounding box.  It's not
    // directly needed but can be useful on the Jsonnet side, for
    // example when defining some simple kinematics.  It is
    // represented by a ray going from extreme corners of a
    // rectangular solid.  Again "wirecell-util wires-info" helps
    // to choose something.
    bounds: {
      tail: wc.point(-4.0, 0.0, 0.0, wc.m),
      head: wc.point(+4.0, 6.1, 7.0, wc.m),
    },
  },

  daq: super.daq {
    nticks: 6000,
  },

  adc: super.adc {
    resolution: 14,
    // reuse ProtoDUNE SP values
    baselines: [1003.4 * wc.millivolt, 1003.4 * wc.millivolt, 507.7 * wc.millivolt],
    fullscale: [0.2 * wc.volt, 1.6 * wc.volt],
  },

  elecs: [
    super.elec {
      // The FE amplifier gain in units of Voltage/Charge.
      gain: 7.8 * wc.mV / wc.fC,

      // The shaping (aka peaking) time of the amplifier shaper.
      shaping: 2.2 * wc.us,
    }
    for n in std.range(0, 3)
  ],

  elec: $.elecs[0],  // nominal

  sim: super.sim {

    // For running in LArSoft, the simulation must be in fixed time mode.
    fixed: true,

    // The "absolute" time (ie, in G4 time) that the lower edge of
    // of final readout tick #0 should correspond to.  This is a
    // "fixed" notion.
    local tick0_time = -250 * wc.us,

    // Open the ductor's gate a bit early.
    local response_time_offset = $.det.response_plane / $.lar.drift_speed,
    local response_nticks = wc.roundToInt(response_time_offset / $.daq.tick),

    ductor: {
      nticks: $.daq.nticks + response_nticks,
      readout_time: self.nticks * $.daq.tick,
      start_time: tick0_time - response_time_offset,
    },

    // To counter the enlarged duration of the ductor, a Reframer
    // chops off the little early, extra time.  Note, tags depend on how
    reframer: {
      tbin: response_nticks,
      nticks: $.daq.nticks,
    },

  },

  files: {
    wires: 'protodunehd-wires-larsoft-v1.json.bz2',

    fields: [
      'np04hd-garfield-6paths.json.bz2',
      'dune-garfield-1d565.json.bz2',
      'dune-garfield-1d565.json.bz2',
      'dune-garfield-1d565.json.bz2',
    ],

    // Noise models for different FE amplifier gains
    // Note: set gain value accordingly in the field of elecs
    // noise: 'protodunehd-noise-spectra-14mVfC-v1.json.bz2',
    noise: 'protodunehd-noise-spectra-7d8mVfC-v1.json.bz2',


    chresp: null,
  },
}
