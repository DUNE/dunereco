////////////////////////////////////////////////////////////////////////
// Class:       TrackHitBacktracker
// Plugin Type: producer (art v2_03_00)
// File:        TrackHitBacktracker_module.cc
//
// Generated at Tue Sep 20 07:19:54 2016 by Matthew Thiesse using cetskelgen
// from cetlib version v1_19_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/RecoBaseArt/HitCreator.h"
#include "lardata/RecoBaseArt/TrackUtils.h"
#include "larcore/Geometry/Geometry.h"

#include <map>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <memory>

#include "TMath.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraph.h"

namespace dune {
  class TrackHitBacktracker;

  class HitStuff {
  public:
    HitStuff();
    raw::TDCtick_t startTick;
    raw::TDCtick_t endTick;
    float rms;
    float peak_time;
    float sigma_peak_time;
    float peak_amplitude;
    float sigma_peak_amplitude;
    float hit_integral;
    float hit_sigma_integral;
    float summedADC;
    float multiplicity;
    short int local_index;
    float goodness_of_fit;
    int dof;
    float preBaseline;
    float postBaseline;
    float preBaselineRMS;
    float postBaselineRMS;
  };

  class ChannelHits {
  public:
    ChannelHits();
    art::Ptr<recob::Wire> artWire;
    std::vector<raw::TDCtick_t> startTicks;
    std::vector<raw::TDCtick_t> endTicks;
    std::vector<float> peakTicks;
    std::vector<float> sigmaPeakTicks;
    std::vector<float> hitRMSs;
    std::vector<float> peakAmplitudes;
    std::vector<float> sigmaPeakAmplitudes;
    std::vector<float> sumADCs;
    std::vector<float> integrals;
    std::vector<float> sigmaIntegrals;
    std::vector<double> pitches;
    std::vector<double> dqdxs;
    unsigned int nhits;
    double chanz;
    raw::ChannelID_t chanid;
    unsigned int wireid;
    unsigned int tpcnum;
    int tpcrem2;
    bool assumedhits;
  };
}

dune::ChannelHits::ChannelHits()
  : nhits(0),
    chanz(-999),
    chanid(0),
    wireid(99999),
    tpcnum(99999),
    tpcrem2(0),
    assumedhits(true)
{}

dune::HitStuff::HitStuff()
  : startTick(0),
    endTick(0),
    rms(0),
    peak_time(0),
    sigma_peak_time(0),
    peak_amplitude(0),
    sigma_peak_amplitude(0),
    hit_integral(0),
    hit_sigma_integral(0),
    summedADC(0),
    multiplicity(0),
    local_index(0),
    goodness_of_fit(0),
    dof(0),
    preBaseline(0),
    postBaseline(0),
    preBaselineRMS(0),
    postBaselineRMS(0)
{}

class dune::TrackHitBacktracker : public art::EDProducer {
public:
  explicit TrackHitBacktracker(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrackHitBacktracker(TrackHitBacktracker const &) = delete;
  TrackHitBacktracker(TrackHitBacktracker &&) = delete;
  TrackHitBacktracker & operator = (TrackHitBacktracker const &) = delete;
  TrackHitBacktracker & operator = (TrackHitBacktracker &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;
  void beginJob() override;

private:
  void calculateHitStuff(dune::ChannelHits chits, std::vector<dune::HitStuff> & hitstuff);
  void getPulseStartEnd(float startTick, float endTick, float peakTick, int & start, int & end);

  std::string fTrackModuleLabel;
  std::string fWireModuleLabel;
  unsigned int fPreBaselineTicks;
  unsigned int fPostBaselineTicks;
  float fPulseStartWidthMult;
  float fPulseEndWidthMult;
  bool fMakeTree;

  TTree * fTree;
  int run;
  int event;
  int channel;
  int wire;
  int tpc;
  int signalsize;
  std::vector<float> signal;
  float prebaseline;
  float postbaseline;
  float prebaserms;
  float postbaserms;
  float integral;
  float sigmaintegral;
  float amplitude;
  float peaktick;
  int begintick;
  int endtick;
  float segmentlength;
  int trackid;
  int numtrajpts;
  double tracklength;
  bool isOnTrack;
  double dqdxatpt;
  int prevStartTick;
  int prevEndTick;
  float prevPeakTime;
  float prevSigmaPeakTime;
  float prevRMS;
  float prevPeakAmplitude;
  float prevSigmaPeakAmplitude;
  float prevSummedADC;
  float prevIntegral;
  float prevSigmaIntegral;
};


dune::TrackHitBacktracker::TrackHitBacktracker(fhicl::ParameterSet const & p)
{
  this->reconfigure(p);
  recob::HitCollectionCreator::declare_products(*this);
}

void dune::TrackHitBacktracker::produce(art::Event & e)
{
  run = e.run();
  event = e.event();
  
  art::ServiceHandle<geo::Geometry> fGeom;

  recob::HitCollectionCreator hcol(*this, e);

  art::Handle< std::vector< recob::Track> > trackh;
  e.getByLabel(fTrackModuleLabel,trackh);

  art::Handle< std::vector< recob::Wire> > wireh;
  e.getByLabel(fWireModuleLabel,wireh);

  art::FindManyP<recob::Hit> fh(trackh, e, fTrackModuleLabel);

  std::cout << trackh->size() << " tracks found" << std::endl;
  std::cout << wireh->size() << " wires found" << std::endl;

  if (trackh.isValid())
    {
      for (unsigned int i_trk = 0; i_trk < trackh->size(); ++i_trk)
	{
	  art::Ptr<recob::Track> ptrack(trackh, i_trk);
	  
	  trackid = ptrack->ID();
	  numtrajpts = ptrack->NumberTrajectoryPoints();
	  tracklength = ptrack->Length();

	  std::vector<art::Ptr<recob::Hit> > hitvec = fh.at(i_trk);

	  std::map<raw::ChannelID_t,dune::ChannelHits> chanhitmap;

	  for (unsigned int i_wire = 0; i_wire < wireh->size(); ++i_wire)
	    {
	      art::Ptr<recob::Wire> pwire(wireh,i_wire);
	      if (pwire->View() == geo::kZ)
		{
		  dune::ChannelHits chanhit;
		  chanhit.nhits = 0;
		  chanhit.artWire = pwire;
		  chanhit.chanid = pwire->Channel();
		  chanhit.wireid = fGeom->ChannelToWire(pwire->Channel())[0].Wire;
		  chanhit.tpcnum = fGeom->ChannelToWire(pwire->Channel())[0].TPC;
		  chanhit.tpcrem2 = chanhit.tpcnum % 2;
		  double wirexyz[3];
                  fGeom->Wire(*(fGeom->ChannelToWire(pwire->Channel()).begin())).GetCenter(wirexyz);
                  chanhit.chanz = static_cast<double>(wirexyz[2]);

		  for (size_t i_hit = 0; i_hit < hitvec.size(); ++i_hit)
		    {
		      art::Ptr<recob::Hit> hit = hitvec.at(i_hit);
		      if (hit->Channel() == chanhit.chanid)
			{
                          double wirePitch = fGeom->WirePitch(geo::kZ,chanhit.tpcnum,0);
                          double angleToVert = fGeom->WireAngleToVertical(geo::kZ,chanhit.tpcnum,0) - 0.5*::util::pi<>();
                          const TVector3& dir = ptrack->DirectionAtPoint(i_hit);
                          double cosgamma = std::abs(std::sin(angleToVert)*dir.Y() + std::cos(angleToVert)*dir.Z());
                          if (cosgamma < 1.e-5) continue;
			  
			  int start,end;
			  getPulseStartEnd(hit->StartTick(),hit->EndTick(),hit->PeakTime(),start,end);
			  chanhit.startTicks.push_back(start);
			  chanhit.endTicks.push_back(end);
			  chanhit.peakTicks.push_back(hit->PeakTime());
			  chanhit.hitRMSs.push_back(hit->RMS());
			  chanhit.sigmaPeakTicks.push_back(hit->SigmaPeakTime());
			  chanhit.peakAmplitudes.push_back(hit->PeakAmplitude());
			  chanhit.sigmaPeakAmplitudes.push_back(hit->SigmaPeakAmplitude());
			  chanhit.sumADCs.push_back(hit->SummedADC());
			  chanhit.integrals.push_back(hit->Integral());
			  chanhit.sigmaIntegrals.push_back(hit->SigmaIntegral());
			  chanhit.pitches.push_back(wirePitch/cosgamma);
			  chanhit.dqdxs.push_back(ptrack->DQdxAtPoint(i_hit,geo::kZ));
			  chanhit.nhits++;
			  chanhit.assumedhits = false;
			}
		    }
		  chanhitmap.insert(std::make_pair(pwire->Channel(),chanhit));
		}
	    }

	  float wirebeginz = 99999;
	  float wireendz = -99999;
	  raw::ChannelID_t wirebeginid = 9999;
	  raw::ChannelID_t wireendid = 0;

	  for (auto const & wireitr : chanhitmap)
	    {
	      if (wireitr.second.nhits > 0 && wireitr.second.chanz < wirebeginz)
		{
		  wirebeginz = wireitr.second.chanz;
		  wirebeginid = wireitr.first;
		}
	      if (wireitr.second.nhits > 0 && wireitr.second.chanz > wireendz)
		{
		  wireendz = wireitr.second.chanz;
		  wireendid = wireitr.first;
		}
	    }
	  std::cout << "Track " << i_trk << " starts on channel " << wirebeginid << " (wire=" << chanhitmap[wirebeginid].wireid << ",tpc=" << chanhitmap[wirebeginid].tpcnum << ") at z=" << wirebeginz << std::endl;
	  std::cout << "Track " << i_trk << " ends on channel " << wireendid << " (wire=" << chanhitmap[wireendid].wireid << ",tpc=" << chanhitmap[wireendid].tpcnum << ") at z=" << wireendz << std::endl;

	  if (wireendz < wirebeginz)
	    {
	      std::cout << "  No hits found!!! " << std::endl;
	      std::cout << "wirebegin Z=" << wirebeginz << "  wireend Z=" << wireendz << std::endl;
	      continue;
	    }

	  if (chanhitmap[wirebeginid].tpcrem2 != chanhitmap[wireendid].tpcrem2)
	    {
	      std::cout << "  Track crosses APA. Can't yet deal with this" << std::endl;
	      continue;
	    }

	  int driftdirection = chanhitmap[wirebeginid].tpcrem2;
	  for (auto it = chanhitmap.cbegin(); it != chanhitmap.cend() ; )
	    {
	      if ((*it).second.tpcrem2 != driftdirection) chanhitmap.erase(it++);
	      else ++it;
	    }

	  std::vector<std::pair<float,raw::ChannelID_t> > trackwires;
	  for (auto const & wireitr : chanhitmap)
	    {
	      if (wireitr.second.chanz >= wirebeginz && wireitr.second.chanz <= wireendz)
		{
		  trackwires.push_back(std::make_pair(wireitr.second.chanz,wireitr.second.chanid));
		  std::cout << wireitr.second.nhits << " reconstructed hit(s) found on channel " << wireitr.second.chanid << " at z=" << wireitr.second.chanz << std::endl;
		}
	    }
	  std::sort(trackwires.begin(),trackwires.end());

	  for (size_t i_tw = 0; i_tw < trackwires.size(); ++i_tw)
	    {	 
	      std::cout << "Getting chanhitmap[" << trackwires[i_tw].first << "[" << i_tw << "].second=" << trackwires[i_tw].second << "]" << std::endl;
	      dune::ChannelHits ch = chanhitmap[trackwires[i_tw].second];

	      std::vector<dune::HitStuff> hitstuffvec;
	      if (ch.nhits == 1)
		{
		  std::cout << "Hit found on channel " << ch.chanid << " (wire=" << ch.wireid << ",tpc=" << ch.tpcnum << ") at z=" << ch.chanz << " and time=(" << ch.startTicks[0] << "->" << ch.peakTicks[0] << "->" << ch.endTicks[0] << ") with rms=" << ch.hitRMSs[0] << std::endl;
		  calculateHitStuff(ch,hitstuffvec);
		}
	      else if (ch.nhits == 0 && ch.chanz >= wirebeginz && ch.chanz <= wireendz)
		{
		  int sub = 0;
		  dune::ChannelHits chnearlow = ch;
		  while (chnearlow.nhits != 1  && chnearlow.chanz >= wirebeginz)
		    {
		      ++sub;
		      chnearlow = chanhitmap[trackwires[i_tw-sub].second];
		    }

		  int add = 0;
		  dune::ChannelHits chnearhigh = ch;
		  while (chnearhigh.nhits != 1 && chnearhigh.chanz <= wireendz)
		    {
		      ++add;
		      chnearhigh = chanhitmap[trackwires[i_tw+add].second];
		    }

		  if (ch.tpcnum != chnearlow.tpcnum || ch.tpcnum != chnearhigh.tpcnum) continue;

		  std::cout << "Nearest channel to " << ch.chanid << " (wire=" << ch.wireid << ",tpc=" << ch.tpcnum << ") with a hit in -Z direction is " << ch.chanz - chnearlow.chanz << " cm away on wire " << chnearlow.chanid << " (wire=" << chnearlow.wireid << ",tpc=" << chnearlow.tpcnum << ")" << std::endl;
		  std::cout << "            time=(" << chnearlow.startTicks[0] << "->" << chnearlow.peakTicks[0] << "->" << chnearlow.endTicks[0] << ") with rms=" << chnearlow.hitRMSs[0] << std::endl;
		  std::cout << "Nearest channel to " << ch.chanid << " (wire=" << ch.wireid << ",tpc=" << ch.tpcnum << ") with a hit in +Z direction is " << chnearhigh.chanz - ch.chanz << " cm away on wire " << chnearhigh.chanid << " (wire=" << chnearhigh.wireid << ",tpc=" << chnearhigh.tpcnum << ")" << std::endl;
		  std::cout << "            time=(" << chnearhigh.startTicks[0] << "->" << chnearhigh.peakTicks[0] << "->" << chnearhigh.endTicks[0] << ") with rms=" << chnearhigh.hitRMSs[0] << std::endl;

		  int cls = chnearlow.startTicks[0];
		  int chs = chnearhigh.startTicks[0];
		  int cle = chnearlow.endTicks[0];
		  int che = chnearhigh.endTicks[0];
		  double clz = chnearlow.chanz;
		  double chz = chnearhigh.chanz;
		  std::cout << "            cls=" << cls << " chs=" << chs << " cle=" << cle << " che=" << che << " clz=" << clz << " chz=" << chz << std::endl;
		  std::cout << "            ch.chanz=" << ch.chanz << std::endl;

		  raw::TDCtick_t chstart = cls + (((chs - cls) / (chz - clz)) * (ch.chanz - clz));
		  raw::TDCtick_t chend   = cle + (((che - cle) / (chz - clz)) * (ch.chanz - clz));

		  double chpitch = (chnearlow.pitches[0] + chnearhigh.pitches[0]) / 2.0;

		  std::vector<float> waveform = ch.artWire->Signal();
		  raw::TDCtick_t peak = std::distance(waveform.begin(),std::max_element(waveform.begin()+chstart,waveform.begin()+chend));

		  ch.startTicks.push_back(chstart);
		  ch.endTicks.push_back(chend);
		  ch.peakTicks.push_back(peak);
		  ch.hitRMSs.push_back(sqrt(fabs(ch.startTicks[0]-ch.endTicks[0])));
		  ch.pitches.push_back(chpitch);
		  ch.nhits++;

		  std::cout << "ASSUMED Hit found on channel " << ch.chanid << " (wire=" << ch.wireid << ",tpc=" << ch.tpcnum << ") at z=" << ch.chanz << " and time=(" << ch.startTicks[0] << "->" << ch.peakTicks[0] << "->" << ch.endTicks[0] << ") with rms=" << ch.hitRMSs[0] << std::endl;


		  calculateHitStuff(ch,hitstuffvec);
		}
	      else
		{
		  continue;
		}

	      for (auto const & hitstuff : hitstuffvec)
		{
		  recob::HitCreator temphit(*(ch.artWire),
					    fGeom->ChannelToWire(ch.chanid)[0],
					    hitstuff.startTick,
					    hitstuff.endTick,
					    hitstuff.rms,
					    hitstuff.peak_time,
					    hitstuff.sigma_peak_time,
					    hitstuff.peak_amplitude,
					    hitstuff.sigma_peak_amplitude,
					    hitstuff.hit_integral,
					    hitstuff.hit_sigma_integral,
					    hitstuff.summedADC,
					    hitstuff.multiplicity,
					    hitstuff.local_index,
					    hitstuff.goodness_of_fit,
					    hitstuff.dof);
		  const recob::Hit hit(temphit.move());
		  hcol.emplace_back(std::move(hit),ch.artWire);

		  if (fMakeTree)
		    {
		      signal = ch.artWire->Signal();
		      signalsize = signal.size();
		      prebaseline = hitstuff.preBaseline;
		      postbaseline = hitstuff.postBaseline;
		      prebaserms = hitstuff.preBaselineRMS;
		      postbaserms = hitstuff.postBaselineRMS;
		      integral = hitstuff.hit_integral;
		      sigmaintegral = hitstuff.hit_sigma_integral;
		      amplitude = hitstuff.peak_amplitude;
		      peaktick = hitstuff.peak_time;
		      begintick = static_cast<int>(hitstuff.startTick);
		      endtick = static_cast<int>(hitstuff.endTick);
		      isOnTrack = !(ch.assumedhits);
		      channel = ch.chanid;
		      wire = ch.wireid;
		      tpc = ch.tpcnum;
		      segmentlength = ch.pitches[0];
		      if (isOnTrack)
			{
			  prevStartTick = ch.startTicks[0];
			  prevEndTick = ch.endTicks[0];
			  prevPeakTime = ch.peakTicks[0];
			  prevSigmaPeakTime = ch.sigmaPeakTicks[0];
			  prevRMS = ch.hitRMSs[0];
			  prevPeakAmplitude = ch.peakAmplitudes[0];
			  prevSigmaPeakAmplitude = ch.sigmaPeakAmplitudes[0];
			  prevSummedADC = ch.sumADCs[0];
			  prevIntegral = ch.integrals[0];
			  prevSigmaIntegral = ch.sigmaIntegrals[0];
			  dqdxatpt = ch.dqdxs[0];
			}
		      fTree->Fill();
		    }
		}
	    }	  
	}
    }
  std::cout << "Finished looking at tracks" << std::endl;
  hcol.put_into(e);
}

void dune::TrackHitBacktracker::reconfigure(fhicl::ParameterSet const & p)
{
  fWireModuleLabel = p.get<std::string>("WireModuleLabel");
  fTrackModuleLabel = p.get<std::string>("TrackModuleLabel");
  fPreBaselineTicks = p.get<unsigned int>("PreBaselineTicks",100);
  fPostBaselineTicks = p.get<unsigned int>("PostBaselineTicks",100);
  fPulseStartWidthMult = p.get<float>("PulseStartWidthMult",2.0);
  fPulseEndWidthMult = p.get<float>("PulseEndWidthMult",3.5);
  fMakeTree = p.get<bool>("MakeTree");
}

void dune::TrackHitBacktracker::beginJob()
{
  art::ServiceHandle<art::TFileService> fTfs;
  if (fMakeTree)
    {
      fTree = fTfs->make<TTree>("trackhit","trackhit");
      fTree->Branch("run",&run,"run/I");
      fTree->Branch("event",&event,"event/I");
      fTree->Branch("channel",&channel,"channel/I");
      fTree->Branch("wire",&wire,"wire/I");
      fTree->Branch("tpc",&tpc,"tpc/I");
      fTree->Branch("signalsize",&signalsize,"signalsize/I");
      fTree->Branch("signal",&signal);
      fTree->Branch("prebaseline",&prebaseline,"prebaseline/F");
      fTree->Branch("postbaseline",&postbaseline,"postbaseline/F");
      fTree->Branch("prebaserms",&prebaserms,"prebaserms/F");
      fTree->Branch("postbaserms",&postbaserms,"postbaserms/F");
      fTree->Branch("integral",&integral,"integral/F");
      fTree->Branch("sigmaintegral",&sigmaintegral,"sigmaintegral/F");
      fTree->Branch("amplitude",&amplitude,"amplitude/F");
      fTree->Branch("peaktick",&peaktick,"peaktick/F");
      fTree->Branch("begintick",&begintick,"begintick/I");
      fTree->Branch("endtick",&endtick,"endtick/I");
      fTree->Branch("segmentlength",&segmentlength,"segmentlength/F");
      fTree->Branch("trackid",&trackid,"trackid/I");
      fTree->Branch("numtrajpts",&numtrajpts,"numtrajpts/I");
      fTree->Branch("tracklength",&tracklength,"tracklength/D");
      fTree->Branch("isOnTrack",&isOnTrack,"isOnTrack/O");
      fTree->Branch("dqdxatpt",&dqdxatpt,"dqdxatpt/D");
      fTree->Branch("prevStartTick",&prevStartTick,"prevStartTick/I");
      fTree->Branch("prevEndTick",&prevEndTick,"prevEndTick/I");
      fTree->Branch("prevPeakTime",&prevPeakTime,"prevPeakTime/F");
      fTree->Branch("prevSigmaPeakTime",&prevSigmaPeakTime,"prevSigmaPeakTime/F");
      fTree->Branch("prevRMS",&prevRMS,"prevRMS/F");
      fTree->Branch("prevPeakAmplitude",&prevPeakAmplitude,"prevPeakAmplitude/F");
      fTree->Branch("prevSigmaPeakAmplitude",&prevSigmaPeakAmplitude,"prevSigmaPeakAmplitude/F");
      fTree->Branch("prevSummedADC",&prevSummedADC,"prevSummedADC/F");
      fTree->Branch("prevIntegral",&prevIntegral,"prevIntegral/F");
      fTree->Branch("prevSigmaIntegral",&prevSigmaIntegral,"prevSigmaIntegral/F");
    }
}

void dune::TrackHitBacktracker::getPulseStartEnd(float gausStartTick, float gausEndTick, float gausPeakTick, int & start, int & end)
{
  start = (int)(gausStartTick) - (int)(fPulseStartWidthMult*fabs(gausStartTick-gausPeakTick));
  end   = (int)(gausEndTick)   + (int)(fPulseEndWidthMult*fabs(gausEndTick-gausPeakTick));
}

void dune::TrackHitBacktracker::calculateHitStuff(dune::ChannelHits chits, std::vector<dune::HitStuff> & hitstuffvec)
{
  TF1 * gaus = new TF1("gaus","([0]/([2]*sqrt(2*3.1415926)))*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))+[3]+x*[4]",0,15000);
  gaus->SetNpx(30000);
  gaus->SetParLimits(2,1,14);

  for (unsigned int i = 0; i < chits.nhits; ++i)
    {
      std::vector<float> waveform = chits.artWire->Signal();
      int start = (int)(chits.startTicks[i]);
      int end   = (int)(chits.endTicks[i]);
      std::vector<float>::iterator pulseStart = (start < 0) ? waveform.begin() : waveform.begin()+start;
      std::vector<float>::iterator pulseEnd   = (end > (int)(waveform.size())) ? waveform.end() : waveform.begin()+end;
      std::vector<float>::iterator preBaselineStart = (start-(int)fPreBaselineTicks < 0) ? waveform.begin() : waveform.begin()+start-fPreBaselineTicks;
      std::vector<float>::iterator postBaselineEnd  = (end+(int)fPostBaselineTicks > (int)(waveform.size())) ? waveform.end() : waveform.begin()+end+fPostBaselineTicks;

      std::vector<float> pulseADCs(pulseStart,pulseEnd);
      std::vector<float> preBaseline(preBaselineStart,pulseStart);
      std::vector<float> postBaseline(pulseEnd,postBaselineEnd);
      
      if (preBaseline.size() < fPreBaselineTicks) continue;
      if (postBaseline.size() < fPostBaselineTicks) continue;

      float preBaselineMean = TMath::Mean(preBaseline.size(),preBaseline.data());
      float postBaselineMean = TMath::Mean(postBaseline.size(),postBaseline.data());

      float baselineSlope = (postBaselineMean-preBaselineMean)/(end-start);
      float baselineIntercept = preBaselineMean;

      for (size_t i = 0; i < preBaseline.size(); i++) preBaseline[i] -= preBaselineMean;
      for (size_t i = 0; i < postBaseline.size(); i++) postBaseline[i] -= postBaselineMean;
      for (size_t i = 0; i < pulseADCs.size(); i++) pulseADCs[i] -= (baselineIntercept + (i * baselineSlope));

      float preBaselineRMS  = TMath::RMS(preBaseline.size(),preBaseline.data());
      float postBaselineRMS  = TMath::RMS(postBaseline.size(),postBaseline.data());

      float pulseIntegral = std::accumulate(pulseADCs.begin(),pulseADCs.end(),0);

      //int peaktick = std::distance(waveform.begin(),std::max_element(pulseStart,pulseEnd));
      //gaus->SetParameter(1,peaktick);
      //gaus->SetParLimits(1,peaktick-5,peaktick+5);

      //TGraph * gr = new TGraph();
      //Int_t tick = 0;
      //for (auto adc : waveform)
      //	{
      //  gr->SetPoint(tick,(Double_t)tick,(Double_t)adc);
      //  ++tick;
      //}
      //gr->Fit(gaus,"BQ0");

      //std::cout << "Channel: " << chits.chanid << "  Hit: " << i << " PreBase: " << preBaselineMean << " (" << preBaselineRMS << ")  PostBase: " << postBaselineMean << " (" << postBaselineRMS << ")  PulseIntegral: " << pulseIntegral << std::endl;

      dune::HitStuff hs;
      hs.startTick = static_cast<raw::TDCtick_t>(start);
      hs.endTick = static_cast<raw::TDCtick_t>(end);
      hs.rms = sqrt((float)end - (float)start);
      //hs.rms = gaus->GetParameter(2);
      hs.peak_time = peaktick;
      hs.sigma_peak_time = sqrt(fabs(hs.endTick-hs.startTick));
      //hs.sigma_peak_time = gaus->GetParError(1);
      hs.peak_amplitude = *std::max_element(pulseADCs.begin(),pulseADCs.end());
      //hs.peak_amplitude = gaus->GetParameter(0)/(gaus->GetParameter(2)*sqrt(2*3.1415926));
      hs.sigma_peak_amplitude = sqrt(hs.peak_amplitude);
      hs.hit_integral = pulseIntegral;
      //hs.hit_integral = gaus->GetParameter(0);
      hs.hit_sigma_integral = sqrt(pulseADCs.size()) * TMath::RMS(pulseADCs.begin(),pulseADCs.end());
      //hs.hit_sigma_integral = gaus->GetParError(0);
      hs.summedADC = std::accumulate(waveform.begin()+start,waveform.end()+end,0);
      hs.multiplicity = 1;
      hs.local_index = -1;
      hs.goodness_of_fit = 1;
      hs.dof = (int)(end-start+1);
      hs.preBaseline = preBaselineMean;
      //hs.preBaseline = gaus->Eval(gaus->GetParameter(1)-5*gaus->GetParameter(2));
      hs.postBaseline = postBaselineMean;
      //hs.postBaseline = gaus->Eval(gaus->GetParameter(1)+5*gaus->GetParameter(2));
      hs.preBaselineRMS = preBaselineRMS;
      hs.postBaselineRMS = postBaselineRMS;
      hitstuffvec.push_back(hs);

      //delete gr;
    }
}

DEFINE_ART_MODULE(dune::TrackHitBacktracker)
