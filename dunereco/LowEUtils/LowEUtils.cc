#include "LowEUtils.h"

namespace solar
{
  LowEUtils::LowEUtils(fhicl::ParameterSet const &p)
      : fHitLabel(p.get<std::string>("HitLabel")),
        fGeometry(p.get<std::string>("Geometry")),
        fDetectorSizeX(p.get<double>("DetectorSizeX")),
        fClusterAlgoTime(p.get<double>("ClusterAlgoTime")),
        fClusterAlgoAdjChannel(p.get<int>("ClusterAlgoAdjChannel")),
        fClusterChargeVariable(p.get<std::string>("ClusterChargeVariable")),
        fClusterMatchNHit(p.get<int>("ClusterMatchNHit")),
        fClusterMatchCharge(p.get<double>("ClusterMatchCharge")),
        fClusterInd0MatchTime(p.get<double>("ClusterInd0MatchTime")),
        fClusterInd1MatchTime(p.get<double>("ClusterInd1MatchTime")),
        fClusterMatchTime(p.get<double>("ClusterMatchTime"))
  {
  }
  void LowEUtils::MakeClusterVector(std::vector<LowEClusterInfo> &ClusterVec, std::vector<std::vector<art::Ptr<recob::Hit>>> &Clusters, art::Event const &evt)
  {
    int ID = 0;
    for (std::vector<art::Ptr<recob::Hit>> Cluster : Clusters)
    {
      if (!Cluster.empty())
      {
        std::stable_sort(Cluster.begin(), Cluster.end(), [](art::Ptr<recob::Hit> a, art::Ptr<recob::Hit> b)
                         { return a->PeakTime() < b->PeakTime(); });
      }
      float StartWire = 0;
      float SigmaStartWire = 0;
      float StartTick = 1e6;
      float SigmaStartTick = 0;
      float StartCharge = 0;
      float StartAngle = 0;
      float StartOpeningAngle = 0;
      float EndWire = 0;
      float SigmaEndWire = 0;
      float EndTick = -1e6;
      float SigmaEndTick = 0;
      float EndCharge = 0;
      float EndAngle = 0;
      float EndOpeningAngle = 0;
      float Integral = 0;
      float IntegralStdDev = 0;
      float SummedADC = 0;
      float SummedADCstdDev = 0;
      int NHit = 0;
      float MultipleHitDensity = 0;
      float Width = 0;
      geo::View_t View;
      geo::PlaneID Plane;
      float WireCoord = 0;
      float SigmaWireCoord = 0;
      float TickCoord = 0;
      float SigmaTickCoord = 0;
      float IntegralAverage = 0;
      float SummedADCaverage = 0;
      float Charge = 0;
      float ChargeStdDev = 0;
      float ChargeAverage = 0;

      // For computing the average tick coordinate
      float StartIntegral = 0;
      float EndIntegral = 0;
      float TickCoordSum = 0;
      float PeakAmplitude = 0;
      float SummedAmplitude = 0;

      // Compute total number of PE and MaxPE.
      for (art::Ptr<recob::Hit> Hit : Cluster)
      {
        NHit++;
        SummedADC += Hit->ROISummedADC();
        SummedADCstdDev += Hit->ROISummedADC() * Hit->ROISummedADC();
        Integral += Hit->Integral();
        IntegralStdDev += Hit->Integral() * Hit->Integral();
        SummedAmplitude += Hit->PeakAmplitude();
        TickCoordSum += Hit->PeakTime() * Hit->PeakAmplitude();

        if (fClusterChargeVariable == "SummedADC")
        {
          Charge += SummedADC * Hit->RMS();
          ChargeStdDev += Charge * Charge;
        }

        if (Hit->PeakAmplitude() > PeakAmplitude)
        {
          PeakAmplitude = Hit->PeakAmplitude();
          WireCoord = Hit->WireID().Wire;
          View = Hit->View();
          Plane = Hit->WireID().planeID();
        }

        if (Hit->PeakTime() < StartTick)
        {
          StartTick = Hit->PeakTime();
          StartWire = Hit->WireID().Wire;
          SigmaStartTick = Hit->SigmaPeakTime();
          StartIntegral = Hit->Integral();
          if (fClusterChargeVariable == "SummedADC")
          {
            StartCharge = Hit->ROISummedADC() * Hit->RMS();
          }
        }

        if (Hit->PeakTime() > EndTick)
        {
          EndTick = Hit->PeakTime();
          EndWire = Hit->WireID().Wire;
          SigmaEndTick = Hit->SigmaPeakTime();
          EndIntegral = Hit->Integral();
          if (fClusterChargeVariable == "SummedADC")
          {
            EndCharge = Hit->ROISummedADC() * Hit->RMS();
          }
        }
      }

      TickCoord = TickCoordSum / SummedAmplitude;
      SummedADCaverage = SummedADC / NHit;
      SummedADCstdDev = sqrt(SummedADCstdDev / NHit - SummedADCaverage * SummedADCaverage);
      IntegralAverage = Integral / NHit;
      IntegralStdDev = sqrt(IntegralStdDev / NHit - IntegralAverage * IntegralAverage);
      Width = abs(EndTick - StartTick);

      if (fClusterChargeVariable == "Integral")
      {
        mf::LogInfo("LowEUtils") << "Charge variable set to 'Integral'";
        StartCharge = StartIntegral;
        EndCharge = EndIntegral;
        Charge = Integral;
        ChargeStdDev = IntegralStdDev;
        ChargeAverage = IntegralAverage;
      }
      if (fClusterChargeVariable == "SummedADC")
      {
        mf::LogInfo("LowEUtils") << "Charge variable set to 'SummedADC' (default is 'Integral')";
        ChargeStdDev = sqrt(ChargeStdDev / NHit - Charge * Charge);
        ChargeAverage = Charge / NHit;
      }
      else
      {
        mf::LogError("LowEUtils") << "Charge variable not recognized. Please choose between 'Integral' and 'SummedADC'";
        Charge = 0;
        ChargeStdDev = 0;
        ChargeAverage = 0;
      }

      ClusterVec.push_back(LowEClusterInfo{StartWire, SigmaStartWire, StartTick, SigmaStartTick, StartCharge, StartAngle, StartOpeningAngle, EndWire, SigmaEndWire, EndTick, SigmaEndTick, EndCharge, EndAngle, EndOpeningAngle, Integral, IntegralStdDev, SummedADC, SummedADCstdDev, NHit, MultipleHitDensity, Width, ID, View, Plane, WireCoord, SigmaWireCoord, TickCoord, SigmaTickCoord, IntegralAverage, SummedADCaverage, Charge, ChargeStdDev, ChargeAverage});
      ID++;
    }
    return;
  }
  void LowEUtils::CalcAdjHits(std::vector<recob::Hit> MyVec, std::vector<std::vector<recob::Hit>> &Clusters, TH1I *MyHist, TH1F *ADCIntHist, bool HeavDebug)
  /*
  Find adjacent hits in time and space:
  - MyVec is the vector of hits to be clustered
  - Clusters is the vector of clusters
  - MyHist is the histogram to be filled with the number of hits in each cluster
  - ADCIntHist is the histogram to be filled with the ADC integral of each cluster
  - HeavDebug is a boolean to turn on/off debugging statements
  */
  {
    const double TimeRange = fClusterAlgoTime;
    const int ChanRange = fClusterAlgoAdjChannel;
    unsigned int FilledHits = 0;
    unsigned int NumOriHits = MyVec.size();

    while (NumOriHits != FilledHits)
    {
      if (HeavDebug)
        std::cout << "\nStart of my while loop" << std::endl;
      std::vector<recob::Hit> AdjHitVec;
      AdjHitVec.push_back(MyVec[0]);
      MyVec.erase(MyVec.begin() + 0);
      int LastSize = 0;
      int NewSize = AdjHitVec.size();

      while (LastSize != NewSize)
      {
        std::vector<int> AddNow;
        for (size_t aL = 0; aL < AdjHitVec.size(); ++aL)
        {
          for (size_t nL = 0; nL < MyVec.size(); ++nL)
          {
            if (HeavDebug)
            {
              std::cout << "\t\tLooping though AdjVec " << aL << " and  MyVec " << nL
                        << " AdjHitVec - " << AdjHitVec[aL].Channel() << " & " << AdjHitVec[aL].PeakTime()
                        << " MVec - " << MyVec[nL].Channel() << " & " << MyVec[nL].PeakTime()
                        << " Channel " << abs((int)AdjHitVec[aL].Channel() - (int)MyVec[nL].Channel()) << " bool " << (bool)(abs((int)AdjHitVec[aL].Channel() - (int)MyVec[nL].Channel()) <= ChanRange)
                        << " Time " << abs(AdjHitVec[aL].PeakTime() - MyVec[nL].PeakTime()) << " bool " << (bool)(abs((double)AdjHitVec[aL].PeakTime() - (double)MyVec[nL].PeakTime()) <= TimeRange)
                        << std::endl;
            }

            if (abs((int)AdjHitVec[aL].Channel() - (int)MyVec[nL].Channel()) <= ChanRange &&
                abs((double)AdjHitVec[aL].PeakTime() - (double)MyVec[nL].PeakTime()) <= TimeRange)
            {

              if (HeavDebug)
                std::cout << "\t\t\tFound a new thing!!!" << std::endl;
              // --- Check that this element isn't already in AddNow.
              bool AlreadyPres = false;

              for (size_t zz = 0; zz < AddNow.size(); ++zz)
              {
                if (AddNow[zz] == (int)nL)
                  AlreadyPres = true;
              }

              if (!AlreadyPres)
                AddNow.push_back(nL);
            } // If this TPCHit is within the window around one of my other hits.
          } // Loop through my vector of colleciton plane hits.
        } // Loop through AdjHitVec

        // --- Now loop through AddNow and remove from Marley whilst adding to AdjHitVec
        std::sort(AddNow.begin(), AddNow.end());
        for (size_t aa = 0; aa < AddNow.size(); ++aa)
        {
          if (HeavDebug)
          {
            std::cout << "\tRemoving element " << AddNow.size() - 1 - aa << " from MyVec ===> "
                      << MyVec[AddNow[AddNow.size() - 1 - aa]].Channel() << " & " << MyVec[AddNow[AddNow.size() - 1 - aa]].PeakTime()
                      << std::endl;
          }

          AdjHitVec.push_back(MyVec[AddNow[AddNow.size() - 1 - aa]]);
          MyVec.erase(MyVec.begin() + AddNow[AddNow.size() - 1 - aa]); // This line creates segmentation fault
                                                                       // std::cout << "Erase works" << std::endl;
        }

        LastSize = NewSize;
        NewSize = AdjHitVec.size();
        if (HeavDebug)
        {
          std::cout << "\t---After that pass, AddNow was size " << AddNow.size() << " ==> LastSize is " << LastSize << ", and NewSize is " << NewSize
                    << "\nLets see what is in AdjHitVec...." << std::endl;
          for (size_t aL = 0; aL < AdjHitVec.size(); ++aL)
          {
            std::cout << "\tElement " << aL << " is ===> " << AdjHitVec[aL].Channel() << " & " << AdjHitVec[aL].PeakTime() << std::endl;
          }
        }
      } // while ( LastSize != NewSize )

      int NumAdjColHits = AdjHitVec.size();
      float SummedADCInt = 0;
      for (recob::Hit TPCHit : AdjHitVec)
        SummedADCInt += TPCHit.Integral();

      if (HeavDebug)
        std::cout << "After that loop, I had " << NumAdjColHits << " adjacent collection plane hits." << std::endl;

      MyHist->Fill(NumAdjColHits);
      ADCIntHist->Fill(SummedADCInt);
      FilledHits += NumAdjColHits;

      if (AdjHitVec.size() > 0)
        Clusters.push_back(AdjHitVec);
    }

    if (HeavDebug)
    {
      std::vector<double> avgChannel;
      std::vector<double> avgTick;
      std::vector<double> summedADCInt;

      for (std::vector<recob::Hit> hits : Clusters)
      {
        double adcInt = 0;
        double channel = 0;
        double tick = 0;

        for (recob::Hit TPCHit : hits)
        {
          tick += TPCHit.Integral() * TPCHit.PeakTime();
          channel += TPCHit.Integral() * TPCHit.Channel();
          adcInt += TPCHit.Integral();
        }
        if (adcInt != 0)
        {
          tick /= adcInt;
          channel /= adcInt;
        }
        summedADCInt.push_back(adcInt);
        avgTick.push_back(tick);
        avgChannel.push_back(channel);
      }

      for (int i = 0; i < int(avgTick.size() - 1); i++)
      {
        for (int j = i + 1; j < int(avgTick.size()); j++)
        {
          std::cout << avgChannel[i] << " " << avgChannel[j] << "  " << std::abs(avgChannel[i] - avgChannel[j]) << std::endl;
          std::cout << avgTick[i] << " " << avgTick[j] << "  " << std::abs(avgTick[i] - avgTick[j]) << std::endl;
          std::cout << summedADCInt[i] << " " << summedADCInt[j] << std::endl;
        }
      }
    }
    return;
  }

  void LowEUtils::CalcAdjHits(std::vector<art::Ptr<recob::Hit>> MyVec, std::vector<std::vector<art::Ptr<recob::Hit>>> &Clusters, std::vector<std::vector<int>> &ClusterIdx)
  /*
  Find adjacent hits in time and space:
  - MyVec is the vector of hits to be clustered
  - Clusters is the vector of clustered hits
  - ClusterIdx is the vector of indices of the hits in each cluster
  */
  {
    const double TimeRange = fClusterAlgoTime;
    const int ChanRange = fClusterAlgoAdjChannel;
    unsigned int FilledHits = 0;
    unsigned int NumOriHits = MyVec.size();
    std::vector<int> HitIdx;
    // Fill the hitidx vector with the index from 0 to MyVec.size()
    for (size_t i = 0; i < MyVec.size(); i++)
      HitIdx.push_back(i);

    while (NumOriHits != FilledHits)
    {
      std::vector<art::Ptr<recob::Hit>> AdjHitVec;
      std::vector<int> AdjHitIdx;

      AdjHitVec.push_back(MyVec[0]);
      AdjHitIdx.push_back(HitIdx[0]);

      MyVec.erase(MyVec.begin() + 0);
      HitIdx.erase(HitIdx.begin() + 0);

      int LastSize = 0;
      int NewSize = AdjHitVec.size();

      while (LastSize != NewSize)
      {
        std::vector<int> AddNow;
        for (size_t aL = 0; aL < AdjHitVec.size(); ++aL)
        {
          for (size_t nL = 0; nL < MyVec.size(); ++nL)
          {
            if (abs((int)AdjHitVec[aL]->Channel() - (int)MyVec[nL]->Channel()) <= ChanRange &&
                abs((double)AdjHitVec[aL]->PeakTime() - (double)MyVec[nL]->PeakTime()) <= TimeRange &&
                AdjHitVec[aL]->View() == MyVec[nL]->View() && AdjHitVec[aL]->SignalType() == MyVec[nL]->SignalType())
            {
              // --- Check that this element isn't already in AddNow.
              bool AlreadyPres = false;

              for (size_t zz = 0; zz < AddNow.size(); ++zz)
              {
                if (AddNow[zz] == (int)nL)
                  AlreadyPres = true;
              }

              if (!AlreadyPres)
                AddNow.push_back(nL);
            } // If this TPCHit is within the window around one of my other hits.
          } // Loop through my vector of colleciton plane hits.
        } // Loop through AdjHitVec

        // --- Now loop through AddNow and remove from Marley whilst adding to AdjHitVec
        std::sort(AddNow.begin(), AddNow.end());
        for (size_t aa = 0; aa < AddNow.size(); ++aa)
        {
          AdjHitVec.push_back(MyVec[AddNow[AddNow.size() - 1 - aa]]);
          MyVec.erase(MyVec.begin() + AddNow[AddNow.size() - 1 - aa]); // This line creates segmentation fault
          // Remove the corresponding index from the HitIdx vector
          AdjHitIdx.push_back(HitIdx[AddNow[AddNow.size() - 1 - aa]]);
          HitIdx.erase(HitIdx.begin() + AddNow[AddNow.size() - 1 - aa]);
        }

        LastSize = NewSize;
        NewSize = AdjHitVec.size();
      } // while ( LastSize != NewSize )

      int NumAdjColHits = AdjHitVec.size();
      FilledHits += NumAdjColHits;
      if (AdjHitVec.size() > 0)
      {
        Clusters.push_back(AdjHitVec);
        ClusterIdx.push_back(AdjHitIdx);
      }
    }
    return;
  }

  //......................................................
  std::vector<double> LowEUtils::ComputeRecoY(
      int Event,
      std::vector<int> &HIndTPC,
      std::vector<double> &Z,
      std::vector<double> &Time,
      std::vector<double> &IndZ,
      std::vector<double> &IndY,
      std::vector<double> &IndT,
      std::vector<double> &IndDir,
      bool HeavDebug)
  /*
   */
  {
    // Create the interpolator
    std::vector<double> RecoY = {};

    for (size_t i = 0; i < Z.size(); i++)
    {
      float ThisdT = -1e6;
      float ThisHZ = Z[i];
      float ThisHT = Time[i];
      float ThisHIndDir = IndDir[0];
      for (size_t j = 0; j < IndT.size(); j++)
      {
        float ThisRecoY = -1e6;
        if (abs(ThisHT - IndT[j]) < ThisdT)
        {
          ThisdT = abs(ThisHT - IndT[j]);
          float ThisHRefY = IndY[j];
          float ThisHRefZ = IndZ[j];
          ThisRecoY = ThisHRefY + (ThisHZ - ThisHRefZ) / (ThisHIndDir);
        }
        RecoY.push_back(ThisRecoY);
      }
    }
    return RecoY;
  }

  void LowEUtils::FillClusterVariables(std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
                                       std::vector<std::vector<int>> &ClNHits,
                                       std::vector<std::vector<float>> &ClT,
                                       std::vector<std::vector<float>> &ClCharge,
                                       bool HeavDebug)
  {
    for (size_t idx = 0; idx < Clusters.size(); idx++)
    {
      std::vector<std::vector<recob::Hit>> TheseClusters = Clusters[idx];
      for (size_t i = 0; i < TheseClusters.size(); i++)
      {
        float clustT = 0;
        float clustCharge = 0;
        std::vector<recob::Hit> ThisCluster = TheseClusters[i];
        for (recob::Hit hit : ThisCluster)
        {
          clustCharge += hit.Integral();
          clustT += hit.PeakTime();
        }
        ClT[idx].push_back(clustT / clustCharge);
        ClCharge[idx].push_back(clustCharge);
        ClNHits[idx].push_back(ThisCluster.size());
      } // End of loop over clusters
    } // End of loop over planes
  }

  void LowEUtils::MatchClusters(
      std::vector<std::vector<std::vector<recob::Hit>>> MatchedClusters,
      std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
      std::vector<std::vector<int>> &ClNHits,
      std::vector<std::vector<float>> &ClT,
      std::vector<std::vector<float>> &ClCharge,
      bool HeavDebug)
  {
    LowEUtils::FillClusterVariables(Clusters, ClNHits, ClT, ClCharge, HeavDebug);
    std::vector<std::vector<int>> MatchedClNHits = {{}, {}, {}};
    std::vector<std::vector<float>> MatchedClT = {{}, {}, {}};
    std::vector<std::vector<float>> MatchedClCharge = {{}, {}, {}};

    // --- Declare our variables to fill
    int MatchInd0Idx = -1, MatchInd1Idx = -1;
    double Ind0ClustdT = fClusterAlgoTime, Ind1ClustdT = fClusterAlgoTime;
    bool MatchInd0 = false, MatchInd1 = false;

    for (size_t ii = 0; ii < Clusters[2].size(); ii++)
    {
      if (Clusters[2][ii].empty())
      {
        continue;
      }
      // Reset variables for next match
      MatchInd0 = false;
      MatchInd1 = false;
      Ind0ClustdT = fClusterAlgoTime;
      Ind1ClustdT = fClusterAlgoTime;

      if (Clusters[0].empty())
      {
        continue;
      }
      for (int jj = 0; jj < int(Clusters[0].size()); jj++)
      {
        if (ClNHits[0][jj] < (1 - fClusterMatchNHit) * ClNHits[2][ii] || ClNHits[0][jj] > (1 + fClusterMatchNHit) * ClNHits[2][ii])
        {
          continue;
        } // Cut on number of hits of Ind0 cluster
        if (ClCharge[0][jj] < (1 - fClusterMatchCharge) * ClCharge[2][ii] || ClCharge[0][jj] > (1 + fClusterMatchCharge) * ClCharge[2][ii])
        {
          continue;
        } // Cut on charge of Ind0 cluster
        if (abs(ClT[2][ii] - ClT[0][jj]) < fClusterAlgoTime && abs(fClusterInd0MatchTime - abs(ClT[2][ii] - ClT[0][jj])) < abs(fClusterInd0MatchTime - Ind0ClustdT))
        {
          Ind0ClustdT = abs(ClT[2][ii] - ClT[0][jj]);
          MatchInd0 = true;
          MatchInd0Idx = jj;
        }
      }
      if (Clusters[1].empty())
      {
        continue;
      }
      for (int zz = 0; zz < int(Clusters[1].size()); zz++)
      {
        if (ClNHits[1][zz] < (1 - fClusterMatchNHit) * ClNHits[2][ii] || ClNHits[1][zz] > (1 + fClusterMatchNHit) * ClNHits[2][ii])
        {
          continue;
        } // Cut on number of hits of Ind1 cluster
        if (ClCharge[1][zz] < (1 - fClusterMatchCharge) * ClCharge[2][ii] || ClCharge[1][zz] > (1 + fClusterMatchCharge) * ClCharge[2][ii])
        {
          continue;
        } // Cut on charge of Ind1 cluster
        if (abs(ClT[2][ii] - ClT[1][zz]) < fClusterAlgoTime && abs(fClusterInd1MatchTime - abs(ClT[2][ii] - ClT[1][zz])) < abs(fClusterInd1MatchTime - Ind1ClustdT))
        {
          Ind1ClustdT = abs(ClT[2][ii] - ClT[1][zz]);
          MatchInd1 = true;
          MatchInd1Idx = zz;
        }
      } // Loop over ind1 clusters
      // Fill matched clusters according to the matching criteria
      if (MatchInd0 && MatchInd1)
      {
        MatchedClusters[0].push_back(Clusters[0][MatchInd0Idx]);
        MatchedClNHits[0].push_back(ClNHits[0][MatchInd0Idx]);
        MatchedClT[0].push_back(ClT[0][MatchInd0Idx]);
        MatchedClCharge[0].push_back(ClCharge[0][MatchInd0Idx]);

        MatchedClusters[1].push_back(Clusters[1][MatchInd1Idx]);
        MatchedClNHits[1].push_back(ClNHits[1][MatchInd1Idx]);
        MatchedClT[1].push_back(ClT[1][MatchInd1Idx]);
        MatchedClCharge[1].push_back(ClCharge[1][MatchInd1Idx]);

        MatchedClusters[2].push_back(Clusters[2][ii]);
        MatchedClNHits[2].push_back(ClNHits[2][ii]);
        MatchedClT[2].push_back(ClT[2][ii]);
        MatchedClCharge[2].push_back(ClCharge[2][ii]);
      }
      else if (MatchInd0 && !MatchInd1)
      {
        MatchedClusters[0].push_back(Clusters[0][MatchInd0Idx]);
        MatchedClNHits[0].push_back(ClNHits[0][MatchInd0Idx]);
        MatchedClT[0].push_back(ClT[0][MatchInd0Idx]);
        MatchedClCharge[0].push_back(ClCharge[0][MatchInd0Idx]);

        MatchedClusters[2].push_back(Clusters[2][ii]);
        MatchedClNHits[2].push_back(ClNHits[2][ii]);
        MatchedClT[2].push_back(ClT[2][ii]);
        MatchedClCharge[2].push_back(ClCharge[2][ii]);
        // Fill missing cluster with empty vector
        MatchedClusters[1].push_back({});
        MatchedClNHits[1].push_back(0);
        MatchedClT[1].push_back(0);
        MatchedClCharge[1].push_back(0);
      }
      else if (!MatchInd0 && MatchInd1)
      {
        MatchedClusters[1].push_back(Clusters[1][MatchInd1Idx]);
        MatchedClNHits[1].push_back(ClNHits[1][MatchInd1Idx]);
        MatchedClT[1].push_back(ClT[1][MatchInd1Idx]);
        MatchedClCharge[1].push_back(ClCharge[1][MatchInd1Idx]);

        MatchedClusters[2].push_back(Clusters[2][ii]);
        MatchedClNHits[2].push_back(ClNHits[2][ii]);
        MatchedClT[2].push_back(ClT[2][ii]);
        MatchedClCharge[2].push_back(ClCharge[2][ii]);
        // Fill missing cluster with empty vector
        MatchedClusters[0].push_back({});
        MatchedClNHits[0].push_back(0);
        MatchedClT[0].push_back(0);
        MatchedClCharge[0].push_back(0);
      }
    }
    ClNHits = MatchedClNHits;
    ClT = MatchedClT;
    ClCharge = MatchedClCharge;
    return;
  }
} // namespace solar