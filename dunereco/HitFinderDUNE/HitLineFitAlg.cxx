/********************************

Implementation of MLESAC algorithm for robust estimation
of model parameters in the presence of outliers.

October 2016
m.thiesse@sheffield.ac.uk

********************************/

#include "HitLineFitAlg.h"

dune::HitLineFitAlg::HitLineFitAlg(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
}

void dune::HitLineFitAlg::SetParameter(int i, double startValue, double minValue, double maxValue)
{
  fParIVal[i].start = startValue; fParIVal[i].min = minValue; fParIVal[i].max = maxValue;
}

bool dune::HitLineFitAlg::CheckModelParameters()
{
  fFitPolN = fParIVal.size()-1;
  int test = 0;
  for (auto & ipar : fParIVal)
    {
      if (ipar.first != test) return false;
      ++test;
    }
  if (fParIVal.size() < 2) return false;
  return true;
}

float dune::HitLineFitAlg::PointToLineDist(TVector3 ptloc, TVector3 linept1, TVector3 linept2)
{
  return (((ptloc-linept1).Cross(ptloc-linept2)).Mag()/(linept2-linept1).Mag());
}

void dune::HitLineFitAlg::DeterministicShuffle(std::vector<unsigned int> & vec)
{
  TRandom3 rand(fSeedValue);
  for (size_t i = vec.size()-1; i > 0; --i)
    {
      unsigned int randint =  (unsigned int)(rand.Uniform(i+1)); // gives random integer from [0,i) but we need [0,i], so do [0,i+1)...
      std::swap(vec[randint],vec[i]);
    }
}

void dune::HitLineFitAlg::SetHorizVertRanges(float hmin, float hmax, float vmin, float vmax)
{
  fHorizRangeMin = hmin; fHorizRangeMax = hmax;
  fVertRangeMin = vmin; fVertRangeMax = vmax;
}

int dune::HitLineFitAlg::FitLine(std::vector<HitLineFitData> & data, HitLineFitResults & bestfit)
{
  if (!CheckModelParameters()) 
    {
      throw cet::exception("HitLineFitAlg") << "Invalid fit parameters. Fix it!";
      return -9;
    }

  // define variables once
  size_t i;
  float horiz,vert,dist,ssr,diff,fiterr;
  TVector3 hitloc,lineloc1,lineloc2;
  int iterations;
  std::vector<float> distances;  

  // vectors of vector indices of points to request from the input data vector
  std::vector<unsigned int> points;
  std::vector<unsigned int> points_best;
  std::vector<unsigned int> bestkeys;
  std::vector<unsigned int> datakeys;

  // initialize values
  bestfit.fitsuccess = false;
  iterations = 0;
  fiterr = std::numeric_limits<float>::max();

  // construct TGraphAsymmErrors for doing the fits
  TGraphAsymmErrors * maybe = new TGraphAsymmErrors(); 
  TGraphAsymmErrors * maybebetter = new TGraphAsymmErrors();

  // define the fit model
  // currently uses a polynomial of N dimensions
  TF1 * model = new TF1("model",TString::Format("pol%u",fFitPolN),fHorizRangeMin,fHorizRangeMax);
  model->SetNpx(10000);
  model->SetParNames("constant","linear","quadratic","cubic","quartic","quintic");
  for (auto & ipar : fParIVal)
    {
      model->SetParLimits(ipar.first,ipar.second.min,ipar.second.max);
      model->SetParameter(ipar.first,ipar.second.start);
    }

  // steering parameters for the RANSAC algorithm
  // n (fMinStartPoints)       = minimum number of data points requred to fit model
  // k (fIterationsMultiplier) = maximum number of iterations allowed, defined as a multiple of the number of points in data set. 200 is usually more than enough
  // t (fInclusionThreshold)   = threshold value for model inclusion, in cm
  // d (fMinAlsoPoints)        = number of close data points requred to assert that a model fits well to data
  unsigned int n = std::max((unsigned int)fMinStartPoints,(unsigned int)(data.size()*0.01));
  if (n >= data.size()) return -1;
  int k = data.size()*fIterationsMultiplier;
  float t = fInclusionThreshold;
  unsigned int d = std::max(int(data.size()*0.05-n),int(fMinAlsoPoints));
  if (d < 2*n || n < 2) return -2;

  // make collection of indices to refer back to data vector
  for (i = 0; i < data.size(); ++i) datakeys.push_back(i);


  if (fLogLevel > 1) std::cout << "Minimum number of data points required to fit the model, n=" << n << "\n" 
			       << "Maximum number of iterations allowed, k=" << k << "\n"
			       << "Threshold value for model inclusion (cm), t=" << t << "\n"
			       << "Number of close data points required to assert a good fit, d=" << d << std::endl;


  // DO MAIN LOOP
  while (iterations < k)
    {
      ++iterations;

      // reset collections
      maybe->Clear();
      maybebetter->Clear();
      points.clear();
      points_best.clear();
      distances.clear();
      if (fLogLevel > 1) 
	{
	  if (iterations % 1000 == 0) std::cout << "Iteration # " << iterations << std::endl;
	}

      // Shuffle list of keys
      DeterministicShuffle(datakeys);

      // Randomly sample n points from the original data set
      for (i = 0; i < n; ++i)
        {
          points.push_back(datakeys[i]);
          HitLineFitData thisdata = data.at(datakeys[i]);
          maybe->SetPoint(i,thisdata.hitHoriz,thisdata.hitVert);
          maybe->SetPointError(i,thisdata.hitHorizErrLo,
			       thisdata.hitHorizErrHi,
			       thisdata.hitVertErrLo,
			       thisdata.hitVertErrHi);
        }

      // Do initial fit through these n points to the model
      TFitResultPtr r = maybe->Fit("model","SQCBR0");

      // Now, search through all data points and select those which are near the best fit model
      for (i = n; i < data.size(); ++i)
        {
          vert = data.at(datakeys[i]).hitVert;
          horiz = data.at(datakeys[i]).hitHoriz;
          hitloc.SetXYZ(horiz,vert,0.);
          double xloc1[1] = { horiz-1 };
          double xloc2[1] = { horiz+1 };
          lineloc1.SetXYZ(horiz-1,model->EvalPar(xloc1,r->GetParams()),0.);
          lineloc2.SetXYZ(horiz+1,model->EvalPar(xloc2,r->GetParams()),0.);
          dist = PointToLineDist(hitloc,lineloc1,lineloc2);
          if (dist <= t) points.push_back(datakeys[i]);
        }

      // If more inliers were found, then we've probably found a good track
      if (points.size() - n > d)
        {
	  // Make collection of point locations from improved set
          for (i = 0; i < points.size(); ++i)
            {
              HitLineFitData thisdata = data.at(points[i]);
              maybebetter->SetPoint(i,thisdata.hitHoriz,thisdata.hitVert);
              maybebetter->SetPointError(i,thisdata.hitHorizErrLo,thisdata.hitHorizErrHi,thisdata.hitVertErrLo,thisdata.hitVertErrHi);
            }

	  // Fit to improved data sample
          TFitResultPtr rb = maybebetter->Fit("model","SQCBR0");

	  // loop over one more time to find all nearby points, and calculate errors in the process
          ssr = 0;
          for (i = 0; i < points.size(); ++i)
            {
              vert = data.at(points[i]).hitVert;
              horiz = data.at(points[i]).hitHoriz;
              hitloc.SetXYZ(horiz,vert,0.);
              double xloc1[1] = { horiz-1 };
              double xloc2[1] = { horiz+1 };
              lineloc1.SetXYZ(horiz-1,model->EvalPar(xloc1,rb->GetParams()),0.);
              lineloc2.SetXYZ(horiz+1,model->EvalPar(xloc2,rb->GetParams()),0.);
              dist = PointToLineDist(hitloc,lineloc1,lineloc2);
              distances.push_back(fabs(dist));
              if (dist <= t)
                {
                  ssr += dist*dist;
                  points_best.push_back(points[i]);
                }
            }
	  if (points_best.size() < 2) continue;

	  // Computing the log likelihood for this model fit
          float sigma = TMath::Median(distances.size(),distances.data());
          float gamma = 0.5;
          float p_outlier_prob = 0;
          float v = 0.5;
          std::vector<float> p_inlier_prob(distances.size());
          for (int j = 0; j < 3; ++j)
            {
              for (i = 0; i < distances.size(); ++i)
                {
                  p_inlier_prob[i] = gamma*TMath::Exp(-(distances[i]*distances[i])/(2*sigma*sigma))/(TMath::Sqrt(2*TMath::Pi())*sigma);
                }
              p_outlier_prob = (1-gamma)/v;
              gamma = 0;
              for (i = 0; i < distances.size(); ++i)
                {
                  gamma += p_inlier_prob[i]/(p_inlier_prob[i]+p_outlier_prob);
                }
              if (distances.size() > 0) gamma /= distances.size();
            }
          float d_cur_penalty = 0;
          for (i = 0; i < distances.size(); ++i)
            {
              d_cur_penalty += TMath::Log(p_inlier_prob[i]+p_outlier_prob);
            }
          d_cur_penalty = (-d_cur_penalty);

	  // If a minimum -Log(L), then take this data set as "true"
	  // Also require that the slope is not zero
          if (d_cur_penalty < fiterr && fabs(model->GetParameter("linear")) > 0.0015)
            {
              diff = d_cur_penalty-fiterr;
              fiterr = d_cur_penalty;
	      for (auto & ipar : fParIVal)
		{
		  bestfit.bestVal[ipar.first] = model->GetParameter(ipar.first);
		  bestfit.bestValError[ipar.first] = model->GetParError(ipar.first);
		}
              bestfit.chi2 = model->GetChisquare();
              bestfit.ndf = model->GetNDF();
              bestfit.sum2resid = ssr;
              bestfit.mle = fiterr;
              bestfit.fitsuccess = true;

	      // Designate the "real" hits from the "fake" hits
              for (i = 0; i < data.size(); i++)
                {
                  if (std::find(points_best.begin(),points_best.end(),i) != points_best.end())
                    {
                      data.at(i).hitREAL = true;
                    }
                  else
                    {
                      data.at(i).hitREAL = false;
                    }
                }
	      if (fLogLevel > 1)
		{
		  std::cout << "-------------Found new minimum!-------------" << std::endl
			    << "FitError=" << fiterr << "  delta(fiterr)=" << diff << std::endl
			    << "Number of points included = " << points_best.size() << " out of " << data.size() << std::endl
			    << "--------------------------------------------" << std::endl;
		}
	    }
        }
    }
  if (bestfit.fitsuccess) return 1;
  return 0;
}

void dune::HitLineFitAlg::reconfigure(fhicl::ParameterSet const& p)
{
  fMinStartPoints = p.get<int>("MinStartPoints");
  fMinAlsoPoints = p.get<int>("MinAlsoPoints");
  fIterationsMultiplier = p.get<float>("IterationsMultiplier");
  fInclusionThreshold = p.get<float>("InclusionThreshold");
  fLogLevel = p.get<int>("LogLevel",1);
}
