
#include "HitLineFitAlg.h"

dune::HitLineFitAlg::HitLineFitAlg(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
}

void dune::HitLineFitAlg::SetCounterPositions(float c1vert, float c1horiz, float c2vert, float c2horiz)
{
  fC1Vert = c1vert; fC1Horiz = c1horiz; fC2Vert = c2vert; fC2Horiz = c2horiz;
}

float dune::HitLineFitAlg::hitGeomDist(TVector3 hitloc, TVector3 trigloc1, TVector3 trigloc2)
{
  return (((hitloc-trigloc1).Cross(hitloc-trigloc2)).Mag()/(trigloc2-trigloc1).Mag());
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

void dune::HitLineFitAlg::SetRanges(float hmin, float hmax, float vmin, float vmax)
{
  fHorizRangeMin = hmin; fHorizRangeMax = hmax;
  fVertRangeMin = vmin; fVertRangeMax = vmax;
}

int dune::HitLineFitAlg::FindTrack(std::vector<dune::HitInformation> & data, HitLineFitResults & bestfit)
{
  size_t i;
  float horiz,vert,dist,ssr,diff;
  std::vector<unsigned int> points;
  std::vector<unsigned int> points_best;
  TVector3 hitloc,lineloc1,lineloc2;
  bestfit.fitsuccess = false;
  TF1 * model = new TF1("model",TString::Format("pol%u",fFitPolN),fHorizRangeMin,fHorizRangeMax);
  model->SetNpx(10000);
  model->SetParNames("constant","linear","quadratic");
  model->SetParLimits(0,fVertRangeMin,fVertRangeMax);
  model->SetParLimits(1,((fC1Vert-fC2Vert)/(fC1Horiz-fC2Horiz))-0.15,((fC1Vert-fC2Vert)/(fC1Horiz-fC2Horiz))+0.15);
  model->SetParLimits(2,-0.0002,0.0002);
  unsigned int n = std::max((unsigned int)fMinStartPoints,(unsigned int)(data.size()*0.01));     //minimum number of data points required to fit the model
  if (n >= data.size()) return -1;
  int k = data.size()*fIterationsMultiplier;     //maximum number of iterations allowed
  float t = fInclusionThreshold;     //threshold value for model inclusion (cm)
  unsigned int d = std::max(int(data.size()*0.05-n),int(fMinAlsoPoints));     //number of close data values required to assert that a model fits well to data
  if (d >= data.size()) return -1;
  if (d < 2*n || n < 2) return -2;
  int iterations = 0;
  std::vector<unsigned int> bestkeys;
  std::vector<unsigned int> datakeys;
  for (i = 0; i < data.size(); ++i) datakeys.push_back(i);
  float fiterr = std::numeric_limits<float>::max();
  TGraphAsymmErrors * maybe = new TGraphAsymmErrors();
  TGraphAsymmErrors * maybebetter = new TGraphAsymmErrors();
  std::cout << "n=" << n << " k=" << k << " t=" << t << " d=" << d << std::endl;
  std::vector<float> distances;
  while (iterations < k)
    {
      maybe->Clear();
      maybebetter->Clear();
      points.clear();
      points_best.clear();
      distances.clear();
      if (iterations % 1000 == 0) std::cout << "iteration # " << iterations << std::endl;
      DeterministicShuffle(datakeys);
      for (i = 0; i < n; ++i)
        {
          points.push_back(datakeys[i]);
          dune::HitInformation thisdata = data.at(datakeys[i]);
          maybe->SetPoint(i,thisdata.hithoriz,thisdata.hitvert);
          maybe->SetPointError(i,thisdata.hithorizerrlo,thisdata.hithorizerrhi,thisdata.hitverterrlo,thisdata.hitverterrhi);
        }
      TFitResultPtr r = maybe->Fit("model","SQCBR0");
      for (i = n; i < data.size(); ++i)
        {
          vert = data.at(datakeys[i]).hitvert;
          horiz = data.at(datakeys[i]).hithoriz;
          hitloc.SetXYZ(horiz,vert,0.);
          double xloc1[1] = { horiz-1 };
          double xloc2[1] = { horiz+1 };
          lineloc1.SetXYZ(horiz-1,model->EvalPar(xloc1,r->GetParams()),0.);
          lineloc2.SetXYZ(horiz+1,model->EvalPar(xloc2,r->GetParams()),0.);
          dist = hitGeomDist(hitloc,lineloc1,lineloc2);
          if (dist <= t) points.push_back(datakeys[i]);
        }
      if (points.size() - n > d)
        {
          for (i = 0; i < points.size(); ++i)
            {
              dune::HitInformation thisdata = data.at(points[i]);
              maybebetter->SetPoint(i,thisdata.hithoriz,thisdata.hitvert);
              maybebetter->SetPointError(i,thisdata.hithorizerrlo,thisdata.hithorizerrhi,thisdata.hitverterrlo,thisdata.hitverterrhi);
            }
          TFitResultPtr rb = maybebetter->Fit("model","SQBR0");
          ssr = 0;
          for (i = 0; i < points.size(); ++i)
            {
              vert = data.at(points[i]).hitvert;
              horiz = data.at(points[i]).hithoriz;
              hitloc.SetXYZ(horiz,vert,0.);
              double xloc1[1] = { horiz-1 };
              double xloc2[1] = { horiz+1 };
              lineloc1.SetXYZ(horiz-1,model->EvalPar(xloc1,rb->GetParams()),0.);
              lineloc2.SetXYZ(horiz+1,model->EvalPar(xloc2,rb->GetParams()),0.);
              dist = hitGeomDist(hitloc,lineloc1,lineloc2);
              distances.push_back(fabs(dist));
              if (dist <= t)
                {
                  ssr += dist*dist;
                  points_best.push_back(points[i]);
                }
            }
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
          if (d_cur_penalty < fiterr && model->GetNpar() >= 2 && fabs(model->GetParameter("linear")) > 0.0015 && points_best.size() > 2)
            {
              diff = d_cur_penalty-fiterr;
              fiterr = d_cur_penalty;
              bestfit.fitconstant     = model->GetParameter("constant");
              bestfit.fitconstanterr  = model->GetParError(model->GetParNumber("constant"));
              bestfit.fitlinear       = model->GetParameter("linear");
              bestfit.fitlinearerr    = model->GetParError(model->GetParNumber("linear"));
              if (model->GetNpar() > 2)
                {
                  bestfit.fitquadratic    = model->GetParameter("quadratic");
                  bestfit.fitquadraticerr = model->GetParError(model->GetParNumber("quadratic"));
                }
              bestfit.fitchi2 = model->GetChisquare();
              bestfit.fitndf = model->GetNDF();
              bestfit.fitsumsqrresidual = ssr;
              bestfit.fitmle = fiterr;
              bestfit.fitsuccess = true;
              for (i = 0; i < data.size(); i++)
                {
                  if (std::find(points_best.begin(),points_best.end(),i) != points_best.end())
                    {
                      data.at(i).fitrealhit = true;
                    }
                  else
                    {
                      data.at(i).fitrealhit = false;
                    }
                }
              std::cout << "-------------Found new minimum!-------------" << std::endl;
              std::cout << "FitError=" << fiterr;
              if (fabs(diff) < 1e30) std::cout << "  delta(fiterr)=" << diff;
              std::cout << "\nNumber of points included = " << points_best.size() << " out of " << data.size() << std::endl;
            }
        }
      iterations++;
    }
  if (bestfit.fitsuccess) return 1;
  return 0;
}

void dune::HitLineFitAlg::reconfigure(fhicl::ParameterSet const& p)
{
  fFitPolN = p.get<int>("FitPolN");
  fMinStartPoints = p.get<int>("MinStartPoints");
  fMinAlsoPoints = p.get<int>("MinAlsoPoints");
  fIterationsMultiplier = p.get<float>("IterationsMultiplier");
  fInclusionThreshold = p.get<float>("InclusionThreshold");
}
