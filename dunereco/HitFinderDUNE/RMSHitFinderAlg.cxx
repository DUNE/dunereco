/****************************************

Searches waveforms for hits where the threshold
is determined on-the-fly from the actual channel noise,
rather than a set value for all channels.

October 2016
m.thiesse@sheffield.ac.uk

****************************************/

#include "RMSHitFinderAlg.h"

dune::RMSHitFinderAlg::RMSHitFinderAlg(fhicl::ParameterSet const & p)
{
  this->reconfigure(p);
  SetSearchTicks(-1,-1);
}

void dune::RMSHitFinderAlg::reconfigure(fhicl::ParameterSet const& p)
{
  fWindowWidth = p.get<int>("WindowWidth");
  fFilterWidth = p.get<float>("FilterWidth");
  fSigmaRiseThreshold = p.get<float>("SigmaRiseThreshold");
  fSigmaFallThreshold = p.get<float>("SigmaFallThreshold");
}

void dune::RMSHitFinderAlg::FindHits(dune::ChannelInformation & chan)
{
  if (fSearchTickStart < 0 || fSearchTickEnd < 0)
    {
      fSearchTickStart = 0; fSearchTickEnd = chan.signalSize;
    }
  std::vector<float> signalFilter(chan.signalFilterVec.begin()+fSearchTickStart,chan.signalFilterVec.begin()+fSearchTickEnd);
  //FilterWaveform(chan.signalVec,chan.signalFilterVec);
  //FilterWaveform(signal,signalFilter);
  //RobustRMSBase(chan.signalVec,chan.baseline,chan.rms);
  //RobustRMSBase(signalFilter,chan.baselineFilter,chan.rmsFilter);
  FindPulses(signalFilter,chan.baselineFilter,chan.rmsFilter,chan.pulse_ends);
  MergeHits(chan.pulse_ends);
}

void dune::RMSHitFinderAlg::FilterWaveform(std::vector<float> wf, std::vector<float> & fwf)
{
  fwf.clear();
  unsigned int wfs = wf.size();
  fwf.resize(wfs);
  std::vector<float> intermediate_wf(wfs);
  float filter_coef = static_cast<float>(TMath::Exp(static_cast<double>(-1.0/fFilterWidth)));
  float a0 = 1.0-filter_coef;
  float b1 = filter_coef;
  unsigned int order = 1;
  for (size_t i = 0; i < wfs; ++i)
    {
      if (i<order)
        {
          intermediate_wf[i] = wf[i];
        }
      else
        {
          intermediate_wf[i] = a0*wf[i]+b1*intermediate_wf[i-1];
        }
    }
  for (size_t i = 0; i < wfs; ++i)
    {
      if (i<order)
        {
          fwf[wfs-1-i] = intermediate_wf[wfs-1-i];
        }
      else
        {
          fwf[wfs-1-i] = a0*intermediate_wf[wfs-1-i]+b1*fwf[wfs-i];
        }
    }
}

void dune::RMSHitFinderAlg::RobustRMSBase(std::vector<float> wf, float & bl, float & r)
{
  unsigned int window_size = (unsigned int)(0.10*wf.size());
  std::vector<float> bl_collection;
  for (size_t i_wf = 0; i_wf < wf.size()-window_size; ++i_wf)
    {
      std::vector<float> partial_window(wf.begin()+i_wf,wf.begin()+i_wf+window_size);
      bl_collection.push_back(static_cast<float>(TMath::Mean(partial_window.size(),partial_window.data())));
    }
  bl = static_cast<float>(TMath::Median(bl_collection.size(),bl_collection.data()));

  std::vector<float> bl_sub_wf;
  for (size_t i_wf = 0; i_wf < wf.size(); ++i_wf) bl_sub_wf.push_back(wf[i_wf]-bl);

  std::vector<float> rms_collection;
  for (size_t i_wf = 0; i_wf < bl_sub_wf.size()-window_size; ++i_wf)
    {
      std::vector<float> partial_window(bl_sub_wf.begin()+i_wf,bl_sub_wf.begin()+i_wf+window_size);
      rms_collection.push_back(static_cast<float>(TMath::RMS(partial_window.size(),partial_window.data())));

    }
  r = static_cast<float>(TMath::Median(rms_collection.size(),rms_collection.data()));
}


void dune::RMSHitFinderAlg::FindPulses(std::vector<float> wf, float bl, float r, std::vector<std::pair<int,int> > & pulse_ends)
{
  pulse_ends.clear();
  int start = 0, end = 0;
  bool started = false;
  for (int i_wf = 0; i_wf < static_cast<int>(wf.size())-fWindowWidth; ++i_wf)
    {
      std::vector<float> window(wf.begin()+i_wf,wf.begin()+i_wf+fWindowWidth);
      float window_mean = static_cast<float>(TMath::Mean(window.size(),window.data()));
      if ((window_mean > bl+fSigmaRiseThreshold*r) && !started)
        {
          started = true;
          for (int i_wf_back = i_wf-1; i_wf_back >= 0; --i_wf_back)
            {
              std::vector<float> window_back(wf.begin()+i_wf_back,wf.begin()+i_wf_back+fWindowWidth);
              float window_mean_back = static_cast<float>(TMath::Mean(window_back.size(),window_back.data()));
              if (window_mean_back < bl+fSigmaFallThreshold*r)
                {
                  start = i_wf_back;
                  break;
                }
            }
          continue;
        }
      if ((window_mean < bl+fSigmaFallThreshold*r && window_mean > bl-fSigmaFallThreshold*r) && started)
        {
          started = false;
          end = i_wf+fWindowWidth;
          pulse_ends.push_back(std::make_pair(start+fSearchTickStart,end+fSearchTickStart));
          continue;
        }
    }
  if (started)
    {
      pulse_ends.push_back(std::make_pair(start+fSearchTickStart,wf.size()-1+fSearchTickStart));
    }
}

void dune::RMSHitFinderAlg::MergeHits(std::vector<std::pair<int,int> > & pulse_ends)
{
  std::vector<std::pair<int,int> > oldpulse_ends = std::move(pulse_ends);
  pulse_ends.clear();
  std::sort(oldpulse_ends.begin(),oldpulse_ends.end());
  int start = 0, end = 0;
  bool started = false;
  for (size_t i_p = 0; i_p < oldpulse_ends.size(); ++i_p)
    {
      if (!started) start = oldpulse_ends[i_p].first;
      end = oldpulse_ends[i_p].second;
      if (i_p < oldpulse_ends.size()-1 && oldpulse_ends[i_p+1].first < oldpulse_ends[i_p].second)
        {
          end = oldpulse_ends[i_p+1].second;
          started = true;
        }
      else
        {
          started = false;
          pulse_ends.push_back(std::make_pair(start,end));
        }
    }
}

