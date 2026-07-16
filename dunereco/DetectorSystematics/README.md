# Wiremodifier module

This directory contains the WireModifier module for DUNE and a fhicl file to run it for HD-FD. 

This module allows to modify the waveforms of ROIs in simulated data based on their "true properties". The true properties are computed from the properties of the simulated energy deposits that can be matched to the ROI.

In the current implementation, the modification of the waveforms is meant to mimic a variation of a detector parameter from its nominal value (electron lifetime, diffusion coefficients, parameters of the modified box model for recombination etc.).
 
## Running Wiremod

The module can be run on simulated data containing recob:Wire, recob:Hits, original and shifted simulated energy deposits as will as sim channels. The label of these objects must be correctly given in the fhicl file. 
Wiremod can typically be run on the output of the detsim stage, by running the gaushit module beforehand.

By default, no modification of the waveform is applied. To apply a modification, the fhicl file must set one of the Apply...Var boolean to true. Examples:
- ApplyGainScale: for a uniform charge gain
- ApplyLifetimeVar: for an electron lifetime variation
- ApplyModBoxVar: for a variation of the modified box parameters values

The varied value of the corresponding parameters must be set as well with the correct units. Examples:
- ElectronicsGainScale: beware that this is not the value of the electronics gain scale, but the value of the uniform rescalign you want to apply
- LifetimeVar: value in microseconds of the new electron lifetime (nominal is 10400)
- ModBoxAlphaVar and ModBoxBetaVar: new values of the modified box model parameters. Alpha is unitless and nominal value is 0.93, Beta is in (g.kV)/(MeV.cm^2) and nominal value if 0.212.

To run wiremod, then just execute lar:
```lar -c wiremod.fcl -o output.root input.root```

### Interesting options

#### Edeps vs simChannels

The modifification of the waveforms requires access to true information. In the default configurayion, this is done by matching the simulated energy deposits to the subROIs and using the properties of the edeps to estimate the truth properties associated to the ROI. 
This matching requires a projection of the edeps positions in 3D space onto a TPC channel and can lead to ambiguities. 

Another option is to use simChannels and IDEs. The ROIs are matched to the simChannels without ambiguity dur to 3D geometry projection and the IDEs properties are used to estimate ROIs properties.

To use simchannels instead of edeps, set useSimChannels to true in the fcl file. 

 
#### Looser matching

The matching of edeps or IDEs to ROIs is quite strict. It requires strict projection of the simulated object inside the ROI (same channel and projected tick inside the ROI). 
You can loosen this criteria using the nTickTolerance and nNeighbourWires options.

nNeighbourWires: only useful when using edeps. It allows matching of edeps to wires neighbouring the one to chich the edep was projected. By setting it to > 0, you usually get 100% matching efficiency for ROIs leading an hit. Increasing it too much also increases the computing time, there is no need to set it larger than 2-3.

nTickTolerance: can be used for both edeps and simChannels. It sets the maximum distance allowed between the simulated deposit projected tick and the ROI border. Increasing this number typically increases the matching of low charge ROIs. 

NB: It is not clear yet what is causing the inefficiencies in matching ROIs to simulated object. Using these two options make things better but it is not sure that this is a correct solution yet. 

#### Multiplicative vs Additive modification of ROI

Now have two options to modify the ROI: multiply the waveform by the ratio of modified gaussian and old gaussian approximation or add the difference between the two. 

SBND has been using the multiplocative approach but ICARUS noticed that it could lead to explosion of the tails (due to low values of gaussian approximate at the tail, the ratio can explode). The additive approach prevents that. 

NB: The additive approach is being tested, not sure yet which approach better reproduces a detector variation. 

## Running full reco1 with Wiremod

To run the full reco1 stage with wiremod included, you need to make sure that the modules run after wiremod (typically gaushit, spsolve, hitfd) use as inputs the modified recob:Wire and recob:Hit.
The file dune_wiremod.fcl contains configuration to do so. First, a copy of gaushit called gauhitWireMod is created and its CalDataModuleLabel is modified. 
Then, the HitLabel of spsolve and hitfd are changed from "gaushit" to "gaushitWireMod".

## Making closure tests plots

### Simulated energy deposits matching plots

These plots allow to investigate the efficiency of matching a simulated energy deposits to  a signal ROI as a function of its energy or associated PDG code. 
The blue distributions represent energy deposits that were succesfully matched whereas the red distirbutions are for unmatched energy deposits. 

To produce this plot, set the boolean SaveEdepMatchingPlots to true in the fhicl file. 
The plots will be saved in a PDF file called Edep_matching_plots.pdf.

### Charge ratio plots

These plots serve as closure tests to make sure the charge modifcation of the ROI makes sense. 

If you want to produce this plot, set the boolean SaveChargeRatioPlot to true in the fhicl file. It can only be done for a Lifetime or Recombination Variation for now. 

The plot will contain ROI charge ratio (after/before recaling) vs dE/dx (for recombination) or drift distance (for attenuation) as well as the analytical expectation. The plot is save in a pdf file called ROI_charge_modification.pdf. 
