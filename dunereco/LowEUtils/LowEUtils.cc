#include "LowEUtils.h"

using namespace producer;

namespace lowe
{
  LowEUtils::LowEUtils(fhicl::ParameterSet const &p)
  : fHitLabel(p.get<std::string>("HitLabel", "hitfd")),
    fClusterAlgoTime(p.get<double>("ClusterAlgoTime", 25.0)), // Time threshold for clustering [ticks]
    fClusterAlgoAdjChannel(p.get<int>("ClusterAlgoAdjChannel", 3)), // Channel threshold for clustering
    fClusterChargeVariable(p.get<std::string>("ClusterChargeVariable", "Integral")), // Variable to use for charge calculation
    fClusterMatchNHit(p.get<int>("ClusterMatchNHit", 2)), // NHit fraction to match clusters. abs(NHitsCol - NHitsInd) / NHitsCol < ClusterMatchNHit.
    fClusterMatchCharge(p.get<double>("ClusterMatchCharge", 0.6)), // Charge fraction to match clusters. abs(ChargeCol - ChargeInd) / ChargeCol < ClusterMatchCharge.
    fClusterInd0MatchTime(p.get<double>("ClusterInd0MatchTime", 0.0)), // Goal time difference to match clusters. abs(TimeCol - TimeInd) < ClusterInd0MatchTime. [ticks]
    fClusterInd1MatchTime(p.get<double>("ClusterInd1MatchTime", 0.0)), // Goal time difference to match clusters. abs(TimeCol - TimeInd) < ClusterInd1MatchTime. [ticks]
    fClusterMatchTime(p.get<double>("ClusterMatchTime", 20.0)), // Max time difference to match clusters. abs(TimeCol - TimeInd) < ClusterMatchTime. [ticks]
    fClusterPreselectionSignal(p.get<bool>("ClusterPreselectionSignal", false)), // Whether to apply signal (purity > 0) preselection
    fClusterPreselectionNHits(p.get<int>("ClusterPreselectionNHits", 3)), // Minimum number of hits in a cluster to consider it for further processing
    fAdjClusterRad(p.get<double>("AdjClusterRad", 100)), // Radius for adjacent cluster search [cm]
    fAdjClusterSingleMatch(p.get<bool>("AdjClusterSingleMatch", false)), // Whether to match adjacent clusters to the first primary cluster only
    fAdjOpFlashMinNHitCut(p.get<int>("AdjOpFlashMinNHitCut", 3)),
    fAdjOpFlashX(p.get<double>("AdjOpFlashX", 100.0)), // X coordinate for flash projection [cm]
    fAdjOpFlashY(p.get<double>("AdjOpFlashY", 100.0)), // Y coordinate for flash projection [cm]
    fAdjOpFlashZ(p.get<double>("AdjOpFlashZ", 100.0)), // Z coordinate for flash projection [cm]
    fAdjOpFlashMinPECut(p.get<double>("AdjOpFlashMinPECut", 20.0)), // Minimum PE for flash selection
    fAdjOpFlashMaxPERatioCut(p.get<double>("AdjOpFlashMaxPERatioCut", 1)), // Maximum PE ratio for flash selection
    fAdjOpFlashMembraneProjection(p.get<bool>("AdjOpFlashMembraneProjection", false)), // Whether to project flashes onto the membrane
    fAdjOpFlashEndCapProjection(p.get<bool>("AdjOpFlashEndCapProjection", false)), // Whether to project flashes onto the end cap
    producer(new ProducerUtils(p))
  {
    // Initialize the LowEUtils instance
    producer->PrintInColor("LowEUtils initialized with parameters from FHiCL configuration.", ProducerUtils::GetColor("green"), "Debug");
  }

  void LowEUtils::MakeClusterVector(std::vector<RawPerPlaneCluster> &ClusterVec, std::vector<std::vector<art::Ptr<recob::Hit>>> &Clusters, art::Event const &evt)
  {
    mf::LogDebug("LowEUtils") << "Charge variable set to " << fClusterChargeVariable;
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
        StartCharge = StartIntegral;
        EndCharge = EndIntegral;
        Charge = Integral;
        ChargeStdDev = IntegralStdDev;
        ChargeAverage = IntegralAverage;
      }
      else if (fClusterChargeVariable == "SummedADC")
      {
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

      ClusterVec.push_back(RawPerPlaneCluster{StartWire, SigmaStartWire, StartTick, SigmaStartTick, StartCharge, StartAngle, StartOpeningAngle, EndWire, SigmaEndWire, EndTick, SigmaEndTick, EndCharge, EndAngle, EndOpeningAngle, Integral, IntegralStdDev, SummedADC, SummedADCstdDev, NHit, MultipleHitDensity, Width, ID, View, Plane, WireCoord, SigmaWireCoord, TickCoord, SigmaTickCoord, IntegralAverage, SummedADCaverage, Charge, ChargeStdDev, ChargeAverage});
      ID++;
    }
    return;
  }

  void LowEUtils::CalcAdjHits(std::vector<recob::Hit> MyVec, std::vector<std::vector<recob::Hit>> &Clusters, TH1I *MyHist, TH1F *ADCIntHist, bool debug)
  /*
  Find adjacent hits in time and space:
  - MyVec is the vector of hits to be clustered
  - Clusters is the vector of clusters
  - MyHist is the histogram to be filled with the number of hits in each cluster
  - ADCIntHist is the histogram to be filled with the ADC integral of each cluster
  - debug is a boolean to turn on/off debugging statements
  */
  {
    const double TimeRange = fClusterAlgoTime;
    const int ChanRange = fClusterAlgoAdjChannel;
    unsigned int FilledHits = 0;
    unsigned int NumOriHits = MyVec.size();

    while (NumOriHits != FilledHits)
    {
      if (debug)
        std::cout << "\nStart of my while loop" << std::endl;
      std::vector<recob::Hit> AdjClusterVec;
      AdjClusterVec.push_back(MyVec[0]);
      MyVec.erase(MyVec.begin() + 0);
      int LastSize = 0;
      int NewSize = AdjClusterVec.size();

      while (LastSize != NewSize)
      {
        std::vector<int> AddNow;
        for (size_t aL = 0; aL < AdjClusterVec.size(); ++aL)
        {
          for (size_t nL = 0; nL < MyVec.size(); ++nL)
          {
            if (debug)
            {
              std::cout << "\t\tLooping though AdjVec " << aL << " and  MyVec " << nL
                        << " AdjClusterVec - " << AdjClusterVec[aL].Channel() << " & " << AdjClusterVec[aL].PeakTime()
                        << " MVec - " << MyVec[nL].Channel() << " & " << MyVec[nL].PeakTime()
                        << " Channel " << abs((int)AdjClusterVec[aL].Channel() - (int)MyVec[nL].Channel()) << " bool " << (bool)(abs((int)AdjClusterVec[aL].Channel() - (int)MyVec[nL].Channel()) <= ChanRange)
                        << " Time " << abs(AdjClusterVec[aL].PeakTime() - MyVec[nL].PeakTime()) << " bool " << (bool)(abs((double)AdjClusterVec[aL].PeakTime() - (double)MyVec[nL].PeakTime()) <= TimeRange)
                        << std::endl;
            }

            if (abs((int)AdjClusterVec[aL].Channel() - (int)MyVec[nL].Channel()) <= ChanRange &&
                abs((double)AdjClusterVec[aL].PeakTime() - (double)MyVec[nL].PeakTime()) <= TimeRange)
            {

              if (debug)
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
        } // Loop through AdjClusterVec

        // --- Now loop through AddNow and remove from Marley whilst adding to AdjClusterVec
        std::sort(AddNow.begin(), AddNow.end());
        for (size_t aa = 0; aa < AddNow.size(); ++aa)
        {
          if (debug)
          {
            std::cout << "\tRemoving element " << AddNow.size() - 1 - aa << " from MyVec ===> "
                      << MyVec[AddNow[AddNow.size() - 1 - aa]].Channel() << " & " << MyVec[AddNow[AddNow.size() - 1 - aa]].PeakTime()
                      << std::endl;
          }

          AdjClusterVec.push_back(MyVec[AddNow[AddNow.size() - 1 - aa]]);
          MyVec.erase(MyVec.begin() + AddNow[AddNow.size() - 1 - aa]); // This line creates segmentation fault
                                                                       // std::cout << "Erase works" << std::endl;
        }

        LastSize = NewSize;
        NewSize = AdjClusterVec.size();
        if (debug)
        {
          std::cout << "\t---After that pass, AddNow was size " << AddNow.size() << " ==> LastSize is " << LastSize << ", and NewSize is " << NewSize
                    << "\nLets see what is in AdjClusterVec...." << std::endl;
          for (size_t aL = 0; aL < AdjClusterVec.size(); ++aL)
          {
            std::cout << "\tElement " << aL << " is ===> " << AdjClusterVec[aL].Channel() << " & " << AdjClusterVec[aL].PeakTime() << std::endl;
          }
        }
      } // while ( LastSize != NewSize )

      int NumAdjColHits = AdjClusterVec.size();
      float SummedADCInt = 0;
      for (recob::Hit TPCHit : AdjClusterVec)
        SummedADCInt += TPCHit.Integral();

      if (debug)
        std::cout << "After that loop, I had " << NumAdjColHits << " adjacent collection plane hits." << std::endl;

      MyHist->Fill(NumAdjColHits);
      ADCIntHist->Fill(SummedADCInt);
      FilledHits += NumAdjColHits;

      if (AdjClusterVec.size() > 0)
        Clusters.push_back(AdjClusterVec);
    }

    if (debug)
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
      std::vector<art::Ptr<recob::Hit>> AdjClusterVec;
      std::vector<int> AdjClusterIdx;

      AdjClusterVec.push_back(MyVec[0]);
      AdjClusterIdx.push_back(HitIdx[0]);

      MyVec.erase(MyVec.begin() + 0);
      HitIdx.erase(HitIdx.begin() + 0);

      int LastSize = 0;
      int NewSize = AdjClusterVec.size();

      while (LastSize != NewSize)
      {
        std::vector<int> AddNow;
        for (size_t aL = 0; aL < AdjClusterVec.size(); ++aL)
        {
          for (size_t nL = 0; nL < MyVec.size(); ++nL)
          {
            if (abs((int)AdjClusterVec[aL]->Channel() - (int)MyVec[nL]->Channel()) <= ChanRange &&
                abs((double)AdjClusterVec[aL]->PeakTime() - (double)MyVec[nL]->PeakTime()) <= TimeRange &&
                AdjClusterVec[aL]->View() == MyVec[nL]->View() && 
                AdjClusterVec[aL]->SignalType() == MyVec[nL]->SignalType())
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
        } // Loop through AdjClusterVec

        // --- Now loop through AddNow and remove from Marley whilst adding to AdjClusterVec
        std::sort(AddNow.begin(), AddNow.end());
        for (size_t aa = 0; aa < AddNow.size(); ++aa)
        {
          AdjClusterVec.push_back(MyVec[AddNow[AddNow.size() - 1 - aa]]);
          MyVec.erase(MyVec.begin() + AddNow[AddNow.size() - 1 - aa]); // This line creates segmentation fault
          // Remove the corresponding index from the HitIdx vector
          AdjClusterIdx.push_back(HitIdx[AddNow[AddNow.size() - 1 - aa]]);
          HitIdx.erase(HitIdx.begin() + AddNow[AddNow.size() - 1 - aa]);
        }

        LastSize = NewSize;
        NewSize = AdjClusterVec.size();
      } // while ( LastSize != NewSize )

      int NumAdjColHits = AdjClusterVec.size();
      FilledHits += NumAdjColHits;
      if (AdjClusterVec.size() > 0)
      {
        Clusters.push_back(AdjClusterVec);
        ClusterIdx.push_back(AdjClusterIdx);
      }
    }
    return;
  }

  //......................................................
  void LowEUtils::FillClusterHitVectors(
    std::vector<recob::Hit> Cluster,
    std::vector<int> &TPC,
    std::vector<int> &Channel,
    std::vector<double> &Charge,
    std::vector<double> &Time,
    double &SigmaTime,
    double &DeltaTime,
    std::vector<double> &Y,
    std::vector<double> &Z,
    std::vector<double> &Dir,
    detinfo::DetectorClocksData ClockData)
  /*
   */
  {
    // --- Clear the vectors
    TPC = {};
    Channel = {};
    Charge = {};
    Time = {};
    Y = {};
    Z = {};
    Dir = {};
    // --- Get the wire readout object
    geo::WireReadoutGeom const &wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    double MaxTime = -1e6;
    double MinTime = 1e6;
    // --- Declare our vectors to fill
    std::vector<int> TrackID;
    TrackID = {};

    for (recob::Hit ThisHit : Cluster)
    {
      TPC.push_back(ThisHit.WireID().TPC);
      Channel.push_back(ThisHit.Channel());
      Charge.push_back(ThisHit.Integral());
      Time.push_back(ThisHit.PeakTime());
      if (ThisHit.PeakTime() > MaxTime)
      {
        MaxTime = ThisHit.PeakTime();
      }
      if (ThisHit.PeakTime() < MinTime)
      {
        MinTime = ThisHit.PeakTime();
      }
      const geo::WireGeo *ThisWire = wireReadout.WirePtr(ThisHit.WireID());
      geo::Point_t hXYZ = ThisWire->GetCenter();
      geo::Point_t sXYZ = ThisWire->GetStart();
      geo::Point_t eXYZ = ThisWire->GetEnd();
      Y.push_back(hXYZ.Y());
      Z.push_back(hXYZ.Z());
      geo::Vector_t Direction = eXYZ - sXYZ;
      float dyds = Direction.Y();
      float dzds = Direction.Z();
      Dir.push_back(dzds / dyds);
    }
    DeltaTime = MaxTime - MinTime;
    SigmaTime = LowEUtils::STD(Time);
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
      bool debug)
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
                                       bool debug)
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
  void LowEUtils::FillClusterVariables(
    std::set<int> SignalTrackIDs,
    std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
    std::vector<std::vector<int>> &ClMainID,
    std::vector<std::vector<int>> &ClNHits,
    std::vector<std::vector<int>> &ClChannel,
    std::vector<std::vector<float>> &ClT,
    std::vector<std::vector<float>> &ClY,
    std::vector<std::vector<float>> &ClZ,
    std::vector<std::vector<float>> &ClDir,
    std::vector<std::vector<float>> &ClCharge,
    std::vector<std::vector<float>> &ClPurity,
    std::vector<std::vector<float>> &ClCompleteness,
    detinfo::DetectorClocksData const &clockData,
    bool debug)
  {
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    geo::WireReadoutGeom const &wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    std::vector<float> globalSignalCharge = {0, 0, 0};
    for (int idx = 0; idx < 3; idx++)
    {
      std::vector<std::vector<recob::Hit>> TheseClusters = Clusters[idx];
      for (size_t i = 0; i < TheseClusters.size(); i++)
      {
        int mainTrID = -1;
        int mainChannel = -1;
        double mainCharge = 0;
        float clustT = 0, clustY = 0, clustZ = 0, clustDir = 0;
        float clustCharge = 0, clustPurity = 0;
        std::vector<recob::Hit> ThisCluster = TheseClusters[i];
        for (recob::Hit hit : ThisCluster)
        {
          int mainHitTrID = -1;
          double mainEFrac = 0;
          double hitCharge = hit.Integral();
          std::vector<sim::TrackIDE> ThisHitIDE = bt_serv->HitToTrackIDEs(clockData, hit);

          for (size_t ideL = 0; ideL < ThisHitIDE.size(); ++ideL)
          {
            if (ThisHitIDE[ideL].energyFrac > mainEFrac)
            {
              mainEFrac = ThisHitIDE[ideL].energyFrac;
              mainHitTrID = abs(ThisHitIDE[ideL].trackID);
            }
          }
          const geo::WireGeo *ThisWire = wireReadout.WirePtr(hit.WireID());
          geo::Point_t hXYZ = ThisWire->GetCenter();
          geo::Point_t sXYZ = ThisWire->GetStart();
          geo::Point_t eXYZ = ThisWire->GetEnd();
          clustCharge += hitCharge;
          if (clustCharge > mainCharge)
          {
            mainTrID = mainHitTrID;
            mainCharge = clustCharge;
            mainChannel = hit.Channel();
          };
          clustT += hit.PeakTime() * hitCharge;
          clustY += hXYZ.Y() * hitCharge;
          clustZ += hXYZ.Z() * hitCharge;
          geo::Vector_t Direction = eXYZ - sXYZ;
          clustDir += hitCharge * Direction.Z() / Direction.Y();
          if (SignalTrackIDs.find(mainHitTrID) != SignalTrackIDs.end())
          {
            globalSignalCharge[idx] += hitCharge;
            clustPurity += hitCharge;
          }
        }
        ClMainID[idx].push_back(mainTrID);
        ClNHits[idx].push_back(ThisCluster.size());
        ClChannel[idx].push_back(mainChannel);
        ClCharge[idx].push_back(clustCharge);
        ClT[idx].push_back(clustT / clustCharge);
        ClY[idx].push_back(clustY / clustCharge);
        ClZ[idx].push_back(clustZ / clustCharge);
        ClDir[idx].push_back(clustDir / clustCharge);
        ClPurity[idx].push_back(clustPurity / clustCharge);
      } // End of loop over clusters
      for (size_t i = 0; i < TheseClusters.size(); i++)
      {
        ClCompleteness[idx].push_back(ClPurity[idx][i] * ClCharge[idx][i] / globalSignalCharge[idx]);
      }
    } // End of loop over planes
  }

  void LowEUtils::MatchClusters(
      std::set<int> SignalTrackIDs,
      std::vector<std::vector<int>> &MatchedClustersIdx,
      std::vector<std::vector<std::vector<recob::Hit>>> &MatchedClusters,
      std::vector<std::vector<int>> ClustersIdx,
      std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
      std::vector<std::vector<int>> &ClMainID,
      std::vector<std::vector<int>> &ClNHits,
      std::vector<std::vector<int>> &ClChannel,
      std::vector<std::vector<float>> &ClT,
      std::vector<std::vector<float>> &ClY,
      std::vector<std::vector<float>> &ClZ,
      std::vector<std::vector<float>> &ClDir,
      std::vector<std::vector<float>> &ClCharge,
      std::vector<std::vector<float>> &ClPurity,
      std::vector<std::vector<float>> &ClCompleteness,
      detinfo::DetectorClocksData const &clockData,
      bool debug)
  {
    LowEUtils::FillClusterVariables(SignalTrackIDs, Clusters, ClMainID, ClNHits, ClChannel, ClT, ClY, ClZ, ClDir, ClCharge, ClPurity, ClCompleteness, clockData, debug);
    std::vector<std::vector<int>> MatchedClMainID = {{}, {}, {}}, MatchedClNHits = {{}, {}, {}}, MatchedClChannel = {{}, {}, {}};
    std::vector<std::vector<float>> MatchedClT = {{}, {}, {}}, MatchedClY = {{}, {}, {}}, MatchedClZ = {{}, {}, {}}, MatchedClDir = {{}, {}, {}};
    std::vector<std::vector<float>> MatchedClCompleteness = {{}, {}, {}}, MatchedClCharge = {{}, {}, {}}, MatchedClPurity = {{}, {}, {}};

    // --- Declare our variables to fill
    int MatchInd0Idx = -1, MatchInd1Idx = -1;
    double Ind0ClustdT = fClusterAlgoTime, Ind1ClustdT = fClusterAlgoTime;
    bool MatchInd0 = false, MatchInd1 = false;

    // Create an index vector with the time sorted indices for each plane
    std::vector<std::vector<int>> ClustersIdxSorted = {{}, {}, {}};
    for (int i = 0; i < 3; i++)
    {
      std::vector<std::pair<float, int>> ClTIdx;
      for (size_t j = 0; j < ClT[i].size(); j++)
      {
        ClTIdx.push_back(std::make_pair(ClT[i][j], j));
      }
      std::sort(ClTIdx.begin(), ClTIdx.end());
      for (auto &pair : ClTIdx)
      {
        ClustersIdxSorted[i].push_back(pair.second);
      }
    }

    for (size_t ii = 0; ii < ClustersIdxSorted[2].size(); ii++)
    {
      int index = ClustersIdxSorted[2][ii];
      if (Clusters[2][index].empty())
      {
        continue;
      }
      // Reset variables for next match
      MatchInd0 = false;
      MatchInd1 = false;
      Ind0ClustdT = fClusterMatchTime;
      Ind1ClustdT = fClusterMatchTime;
      MatchInd0Idx = -1;
      MatchInd1Idx = -1;
      // Initialize the last sorted indices for Ind0 and Ind1
      int LastSortedInd0 = 0, LastSortedInd1 = 0;
      // std::cout << " - Matching cluster " << ii << " with index " << index << " and time " << ClT[2][index] << std::endl;
      for (int jj = LastSortedInd0; jj < int(ClustersIdxSorted[0].size()); jj++)
      {
        int index0 = ClustersIdxSorted[0][jj];
        if (ClT[2][index] - ClT[0][index0] > fClusterMatchTime ||
            abs(ClNHits[0][index0] - ClNHits[2][index]) / ClNHits[2][index] > fClusterMatchNHit || 
            abs(ClCharge[0][index0] - ClCharge[2][index]) / ClCharge[2][index] > fClusterMatchCharge)
        {
          continue;
        } // Cut on number of hits of Ind0 cluster
        if (abs(ClT[2][index] - ClT[0][index0]) <= fClusterMatchTime)
        {
          // std::cout << "    Checking Ind0 cluster " << index0 << " with index " << jj << " and time " << ClT[0][index0] << std::endl;
          if (abs(fClusterInd0MatchTime - abs(ClT[2][index] - ClT[0][index0])) < abs(fClusterInd0MatchTime - Ind0ClustdT))
          {
            MatchInd0 = true;
            Ind0ClustdT = abs(ClT[2][index] - ClT[0][index0]);
            MatchInd0Idx = index0;
            LastSortedInd0 = jj; // Store the last sorted index for Ind0
          }
        }
        else if ((ClT[0][index0] - ClT[2][index]) > fClusterMatchTime)
        {
          break; // Exit loop if time difference exceeds threshold
        }
      }

      for (int zz = LastSortedInd1; zz < int(ClustersIdxSorted[1].size()); zz++)
      {
        int index1 = ClustersIdxSorted[1][zz];
        if (ClT[1][index1] - ClT[2][index] > fClusterMatchTime ||
            abs(ClNHits[1][index1] - ClNHits[2][index]) / ClNHits[2][index] > fClusterMatchNHit || 
            abs(ClCharge[1][index1] - ClCharge[2][index]) / ClCharge[2][index] > fClusterMatchCharge)
        {
          continue;
        } // Cut on number of hits of Ind1 cluster
        if (abs(ClT[2][index] - ClT[1][index1]) <= fClusterMatchTime)
        {
          // std::cout << "    Checking Ind1 cluster " << index1 << " with index " << zz << " and time " << ClT[1][index1] << std::endl;
          if (abs(fClusterInd1MatchTime - abs(ClT[2][index] - ClT[1][index1])) < abs(fClusterInd1MatchTime - Ind1ClustdT))
          {
            MatchInd1 = true;
            Ind1ClustdT = abs(ClT[2][index] - ClT[1][index1]);
            MatchInd1Idx = index1;
            LastSortedInd1 = zz; // Store the last sorted index for Ind1
          }
        }
        else if ((ClT[1][index1] - ClT[2][index]) > fClusterMatchTime)
        {
          break; // Exit loop if time difference exceeds threshold
        }
      } // Loop over ind1 clusters
      // Fill matched clusters according to the matching criteria
      float ClY1, ClY0, ClY2;
      
      if (!MatchInd0 && !MatchInd1)
      {
        // std::cout << "    No match found for index " << index << ". Skipping..." << std::endl;
        continue;
      } // No match found, skip to next index
      
      else if (MatchInd0 && MatchInd1)
      {
        ClY0 = ClY[0][MatchInd0Idx] + (ClZ[2][index] - ClZ[0][MatchInd0Idx]) / (ClDir[0][MatchInd0Idx]);
        ClY1 = ClY[1][MatchInd1Idx] + (ClZ[2][index] - ClZ[1][MatchInd1Idx]) / (ClDir[1][MatchInd1Idx]);
        ClY2 = (ClY0 + ClY1) / 2;
        // std::cout << "\tMatched all three planes" << std::endl;
        MatchedClustersIdx[0].push_back(ClustersIdx[0][MatchInd0Idx]);
        MatchedClusters[0].push_back(Clusters[0][MatchInd0Idx]);
        MatchedClMainID[0].push_back(ClMainID[0][MatchInd0Idx]);
        MatchedClNHits[0].push_back(ClNHits[0][MatchInd0Idx]);
        MatchedClChannel[0].push_back(ClChannel[0][MatchInd0Idx]);
        MatchedClT[0].push_back(ClT[0][MatchInd0Idx]);
        MatchedClY[0].push_back(ClY0);
        MatchedClZ[0].push_back(ClZ[0][MatchInd0Idx]);
        MatchedClCharge[0].push_back(ClCharge[0][MatchInd0Idx]);
        MatchedClPurity[0].push_back(ClPurity[0][MatchInd0Idx]);
        MatchedClCompleteness[0].push_back(ClCompleteness[0][MatchInd0Idx]);

        MatchedClustersIdx[1].push_back(ClustersIdx[1][MatchInd1Idx]);
        MatchedClusters[1].push_back(Clusters[1][MatchInd1Idx]);
        MatchedClMainID[1].push_back(ClMainID[1][MatchInd1Idx]);
        MatchedClNHits[1].push_back(ClNHits[1][MatchInd1Idx]);
        MatchedClChannel[1].push_back(ClChannel[1][MatchInd1Idx]);
        MatchedClT[1].push_back(ClT[1][MatchInd1Idx]);
        MatchedClY[1].push_back(ClY1);
        MatchedClZ[1].push_back(ClZ[1][MatchInd1Idx]);
        MatchedClCharge[1].push_back(ClCharge[1][MatchInd1Idx]);
        MatchedClPurity[1].push_back(ClPurity[1][MatchInd1Idx]);
        MatchedClCompleteness[1].push_back(ClCompleteness[1][MatchInd1Idx]);
        // ThisRecoY = ThisHRefY + (ThisHZ - ThisHRefZ) / (ThisHIndDir);
      }
      
      else if (MatchInd0 && !MatchInd1)
      {
        ClY0 = ClY[0][MatchInd0Idx] + (ClZ[2][index] - ClZ[0][MatchInd0Idx]) / (ClDir[0][MatchInd0Idx]);
        ClY1 = -1e6;
        ClY2 = ClY0;
        // std::cout << "\tMatched only Ind0 and Col" << std::endl;
        MatchedClustersIdx[0].push_back(ClustersIdx[0][MatchInd0Idx]);
        MatchedClusters[0].push_back(Clusters[0][MatchInd0Idx]);
        MatchedClMainID[0].push_back(ClMainID[0][MatchInd0Idx]);
        MatchedClNHits[0].push_back(ClNHits[0][MatchInd0Idx]);
        MatchedClChannel[0].push_back(ClChannel[0][MatchInd0Idx]);
        MatchedClT[0].push_back(ClT[0][MatchInd0Idx]);
        MatchedClY[0].push_back(ClY0);
        MatchedClZ[0].push_back(ClZ[0][MatchInd0Idx]);
        MatchedClCharge[0].push_back(ClCharge[0][MatchInd0Idx]);
        MatchedClPurity[0].push_back(ClPurity[0][MatchInd0Idx]);
        MatchedClCompleteness[0].push_back(ClCompleteness[0][MatchInd0Idx]);
        // Fill missing cluster with empty vector
        MatchedClustersIdx[1].push_back({});
        MatchedClusters[1].push_back({});
        MatchedClMainID[1].push_back(-1);
        MatchedClNHits[1].push_back(0);
        MatchedClChannel[1].push_back(0);
        MatchedClT[1].push_back(-1e6);
        MatchedClY[1].push_back(-1e6);
        MatchedClZ[1].push_back(-1e6);
        MatchedClCharge[1].push_back(-1e6);
        MatchedClPurity[1].push_back(-1);
        MatchedClCompleteness[1].push_back(-1);

      }
      
      else if (!MatchInd0 && MatchInd1)
      {
        ClY0 = -1e6;
        ClY1 = ClY[1][MatchInd1Idx] + (ClZ[2][index] - ClZ[1][MatchInd1Idx]) / (ClDir[1][MatchInd1Idx]);
        ClY2 = ClY1;
        // std::cout << "\tMatched only Ind1 and Col" << std::endl;
        MatchedClustersIdx[1].push_back(ClustersIdx[1][MatchInd1Idx]);
        MatchedClusters[1].push_back(Clusters[1][MatchInd1Idx]);
        MatchedClMainID[1].push_back(ClMainID[1][MatchInd1Idx]);
        MatchedClNHits[1].push_back(ClNHits[1][MatchInd1Idx]);
        MatchedClChannel[1].push_back(ClChannel[1][MatchInd1Idx]);
        MatchedClT[1].push_back(ClT[1][MatchInd1Idx]);
        MatchedClY[1].push_back(ClY1);
        MatchedClZ[1].push_back(ClZ[1][MatchInd1Idx]);
        MatchedClCharge[1].push_back(ClCharge[1][MatchInd1Idx]);
        MatchedClPurity[1].push_back(ClPurity[1][MatchInd1Idx]);
        MatchedClCompleteness[1].push_back(ClCompleteness[1][MatchInd1Idx]);
        // Fill missing cluster with empty vector
        MatchedClustersIdx[0].push_back({});
        MatchedClusters[0].push_back({});
        MatchedClMainID[0].push_back(-1);
        MatchedClNHits[0].push_back(0);
        MatchedClChannel[0].push_back(0);
        MatchedClT[0].push_back(-1e6);
        MatchedClY[0].push_back(-1e6);
        MatchedClZ[0].push_back(-1e6);
        MatchedClCharge[0].push_back(-1e6);
        MatchedClPurity[0].push_back(-1);
        MatchedClCompleteness[0].push_back(-1);
      }

      // std::cout << "    ***Matched cluster " << ii << " with index " << index << " to Ind0: " << MatchInd0Idx << ", Ind1: " << MatchInd1Idx
                // << ", ClY2: " << ClY2 << std::endl;
      // Fill the matched collection plane cluster
      MatchedClustersIdx[2].push_back(ClustersIdx[2][index]);
      MatchedClusters[2].push_back(Clusters[2][index]);
      MatchedClMainID[2].push_back(ClMainID[2][index]);
      MatchedClNHits[2].push_back(ClNHits[2][index]);
      MatchedClChannel[2].push_back(ClChannel[2][index]);
      MatchedClT[2].push_back(ClT[2][index]);
      MatchedClY[2].push_back(ClY2);
      MatchedClZ[2].push_back(ClZ[2][index]);
      MatchedClCharge[2].push_back(ClCharge[2][index]);
      MatchedClPurity[2].push_back(ClPurity[2][index]);
      MatchedClCompleteness[2].push_back(ClCompleteness[2][index]);
    }
    ClMainID = MatchedClMainID;
    ClNHits = MatchedClNHits;
    ClChannel = MatchedClChannel;
    ClT = MatchedClT;
    ClY = MatchedClY;
    ClZ = MatchedClZ;
    ClDir = MatchedClDir;
    ClCharge = MatchedClCharge;
    ClPurity = MatchedClPurity;
    ClCompleteness = MatchedClCompleteness;
    return;
  }

  void LowEUtils::MatchClusters(
      std::vector<std::vector<std::vector<recob::Hit>>> &MatchedClusters,
      std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
      std::vector<std::vector<int>> &ClNHits,
      std::vector<std::vector<float>> &ClT,
      std::vector<std::vector<float>> &ClCharge,
      bool debug)
  {
    LowEUtils::FillClusterVariables(Clusters, ClNHits, ClT, ClCharge, debug);
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
        MatchedClT[1].push_back(-1e6);
        MatchedClCharge[1].push_back(-1e6);
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
        MatchedClT[0].push_back(-1e6);
        MatchedClCharge[0].push_back(-1e6);
      }
    }
    ClNHits = MatchedClNHits;
    ClT = MatchedClT;
    ClCharge = MatchedClCharge;
    return;
  }

  //......................................................
  double LowEUtils::STD(const std::vector<double>& Vec)
  {
    if (Vec.size() == 0)
    {
      return 0.0;
    }
    double mean = std::accumulate(Vec.begin(), Vec.end(), 0.0) / Vec.size();
    double variance = 0.0;

    for (double value : Vec) {
        variance += (value - mean) * (value - mean);
    }

    variance /= Vec.size();
    return std::sqrt(variance);
  }

  //......................................................
  void LowEUtils::FindPrimaryClusters(const std::vector<art::Ptr<solar::LowECluster>> &SolarClusterVector, std::vector<bool> &EventCandidateFound, std::vector<std::vector<art::Ptr<solar::LowECluster>>> &EventCandidateVector, std::vector<std::vector<int>> &EventCandidateIdx, const detinfo::DetectorClocksData &clockData, const art::Event &evt)
  {
    // This is the low energy primary cluster algorithm. It groups all input clusters into event candidates by finding the primary clusters (charge > adjacent clusters up to distance fAdjClusterRad).
    // The algorithm outputs vectors of clusters where the first entry is the primary cluster and the rest are the corresponding adjacent clusters.

    // Initialize the vector of EventCandidateVector and the vector of indices.
    EventCandidateFound.clear();
    EventCandidateVector.clear();
    EventCandidateIdx.clear();

    // Find correct drift time and distance relation from geometry service
    art::ServiceHandle<geo::Geometry> geom;
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    geo::CryostatID c(0);
    const geo::CryostatGeo& cryostat = geom->Cryostat(c);
    const geo::TPCGeo& tpcg = cryostat.TPC(0);
    const double driftLength = tpcg.DriftDistance(); // in cm
    const double driftTime = driftLength / detProp.DriftVelocity(); // in microseconds

    // Create a sorting index based on the input vector.
    std::vector<int> sortingIndex(SolarClusterVector.size());
    std::iota(sortingIndex.begin(), sortingIndex.end(), 0);
    std::stable_sort(sortingIndex.begin(), sortingIndex.end(), [&](int a, int b)
                     { return SolarClusterVector[a]->getAverageTime() < SolarClusterVector[b]->getAverageTime(); });

    // Create a vector to track if a cluster has been evaluated (found primary or atributted to one).
    std::vector<bool> EvaluatedCluster(SolarClusterVector.size(), false);

    std::string sEventCandidateFinding = "LowEUtils::FindPrimaryClusters " + ProducerUtils::str(SolarClusterVector.size()) + " clusters found in the event\n";
    for (auto it = sortingIndex.begin(); it != sortingIndex.end(); ++it)
    {
      if (EvaluatedCluster[*it]) {
        EventCandidateFound.push_back(false);
        continue; // Skip if the cluster has already been evaluated
      }
      
      std::vector<int> AdjClusterIdx = {};
      std::vector<art::Ptr<solar::LowECluster>> AdjClusterVec = {};
      std::vector<bool> EvaluatedAdjCluster(SolarClusterVector.size(), false);

      const auto &cluster = SolarClusterVector[*it];
      if (fClusterPreselectionSignal && cluster->getPurity() == 0) {
        EventCandidateFound.push_back(false);
        continue;
      }

      if (cluster->getNHits() <= fClusterPreselectionNHits) {
        EventCandidateFound.push_back(false);
        continue;
      }

      bool PreselectionCluster = true;

      // If a trigger hit is found, start a new cluster with the hits around it that are within the time and radius range
      EvaluatedCluster[*it] = true;
      AdjClusterIdx.push_back(*it);
      AdjClusterVec.push_back(cluster);
      sEventCandidateFinding += "Trigger cluster found: NHits " + ProducerUtils::str(cluster->getNHits()) + " Channel " + ProducerUtils::str(cluster->getMainChannel()) + " Time " + ProducerUtils::str(cluster->getAverageTime()) + " Charge " + ProducerUtils::str(cluster->getTotalCharge()) + " Purity " + ProducerUtils::str(cluster->getPurity()) + "\n";
      
      // Make use of the fact that the clusters are sorted in time to only consider the clusters that are adjacent in the vector up to a certain time range
      for (auto it2 = it; it2 != sortingIndex.end(); ++it2) // Start from the next element
      {
        // make sure we don't go out of bounds and the pointer is valid
        if (it == sortingIndex.end()) {
          sEventCandidateFinding += "\tBreaking time loop at end of vector\n";
          break;
        }
        if (it2 == it) {
          continue; // Skip the current cluster itself
        }

        auto &adjcluster = SolarClusterVector[*it2]; // Update adjcluster here

        float dTcluster1 = adjcluster->getAverageTime() - cluster->getAverageTime();
        float dXcluster1 = dTcluster1 * driftLength / driftTime;
        if (std::abs(dTcluster1) > fAdjClusterRad * driftTime / driftLength) {
          sEventCandidateFinding += "\tBreaking time loop at dT " + ProducerUtils::str(dTcluster1) + " us\n";
          break;
        }

        // If AdjClusterSingleMatch is true and adjcluster has already been evaluated, skip
        if (fAdjClusterSingleMatch && EvaluatedCluster[*it2] == true) {
          sEventCandidateFinding += "\tSkipping already evaluated adjacent cluster: NHits " + ProducerUtils::str(adjcluster->getNHits()) + " Channel " + ProducerUtils::str(adjcluster->getMainChannel()) + " Time " + ProducerUtils::str(adjcluster->getAverageTime()) + " Charge " + ProducerUtils::str(adjcluster->getTotalCharge()) + " Purity " + ProducerUtils::str(adjcluster->getPurity()) + "\n";
          continue;
        }
        
        // If cluster has already been clustered, skip
        if (EvaluatedAdjCluster[*it2] == true) {
          sEventCandidateFinding += "\tSkipping already evaluated cluster: NHits " + ProducerUtils::str(adjcluster->getNHits()) + " Channel " + ProducerUtils::str(adjcluster->getMainChannel()) + " Time " + ProducerUtils::str(adjcluster->getAverageTime()) + " Charge " + ProducerUtils::str(adjcluster->getTotalCharge()) + " Purity " + ProducerUtils::str(adjcluster->getPurity()) + "\n";
          continue;
        }

        auto ref4 = TVector3(0, cluster->getY(), cluster->getZ()) - TVector3(dXcluster1, adjcluster->getY(), adjcluster->getZ());
        if (ref4.Mag() < fAdjClusterRad)
        {
          // sEventCandidateFinding += "\tFound adjacent cluster: NHits " + ProducerUtils::str(adjcluster->getNHits()) + " Channel " + ProducerUtils::str(adjcluster->getMainChannel()) + " Time " + ProducerUtils::str(adjcluster->getAverageTime()) + " Charge " + ProducerUtils::str(adjcluster->getTotalCharge()) + "\n";
          if (adjcluster->getTotalCharge() > cluster->getTotalCharge())
          {
            sEventCandidateFinding += "¡¡¡Found bigger cluster: NHits " + ProducerUtils::str(adjcluster->getNHits()) + " Channel " + ProducerUtils::str(adjcluster->getMainChannel()) + " Time " + ProducerUtils::str(adjcluster->getAverageTime()) + " Charge " + ProducerUtils::str(adjcluster->getTotalCharge()) + " Purity " + ProducerUtils::str(adjcluster->getPurity()) + "\n";
            EvaluatedCluster[*it] = false;
            PreselectionCluster = false;

            // Reset the EvaluatedCluster values for the clusters that have been added to the cluster
            for (auto it3 = AdjClusterIdx.begin(); it3 != AdjClusterIdx.end(); ++it3)
            {
              sEventCandidateFinding += "---Removing cluster: NHits " + ProducerUtils::str(SolarClusterVector[*it3]->getNHits()) + " Channel " + ProducerUtils::str(SolarClusterVector[*it3]->getMainChannel()) + " Time " + ProducerUtils::str(SolarClusterVector[*it3]->getAverageTime()) + " Charge " + ProducerUtils::str(SolarClusterVector[*it3]->getTotalCharge()) + " Purity " + ProducerUtils::str(SolarClusterVector[*it3]->getPurity()) + "\n";
              EvaluatedCluster[*it3] = false;
              EvaluatedAdjCluster[*it3] = false; // Reset the evaluated status for the adjacent clusters
            }
            break;
          }
          else{
            sEventCandidateFinding += "+++Adding cluster: NHits " + ProducerUtils::str(adjcluster->getNHits()) + " Channel " + ProducerUtils::str(adjcluster->getMainChannel()) + " Time " + ProducerUtils::str(adjcluster->getAverageTime()) + " Charge " + ProducerUtils::str(adjcluster->getTotalCharge()) + " Purity " + ProducerUtils::str(adjcluster->getPurity()) + "\n";
            AdjClusterVec.push_back(adjcluster);
            AdjClusterIdx.push_back(*it2);
            EvaluatedCluster[*it2] = true;
            EvaluatedAdjCluster[*it2] = true; // Mark this adjacent cluster as evaluated
          }
        }
      }

      if (PreselectionCluster == false) {
        sEventCandidateFinding += "\tSkipping non-preselected cluster\n";
        EvaluatedCluster[*it] = false;
        EventCandidateFound.push_back(false);
        continue; // Skip if the cluster has been reset
      }

      for (auto it4 = it; it4 != sortingIndex.begin() - 1 ; --it4) // Start from the previous element
      {
        // make sure we don't go out of bounds and the pointer is valid
        if (it == sortingIndex.begin()) {
          sEventCandidateFinding += "\tBreaking time loop at beginning of vector\n";
          break;
        }
        if (it4 == it) {
          continue; // Skip the current cluster itself
        }
        
        auto &adjcluster = SolarClusterVector[*it4];

        float dTcluster2 = cluster->getAverageTime() - adjcluster->getAverageTime();
        float dXcluster2 = dTcluster2 * driftLength / driftTime;
        if (std::abs(dTcluster2) > fAdjClusterRad * driftTime / driftLength) {
          sEventCandidateFinding += "\tBreaking time loop at dT " + ProducerUtils::str(dTcluster2) + " us\n";
          break;
        }

        // If AdjClusterSingleMatch is true and adjcluster has already been evaluated, skip
        if (fAdjClusterSingleMatch && EvaluatedCluster[*it4] == true) {
          sEventCandidateFinding += "\tSkipping already evaluated adjacent cluster: NHits " + ProducerUtils::str(adjcluster->getNHits()) + " Channel " + ProducerUtils::str(adjcluster->getMainChannel()) + " Time " + ProducerUtils::str(adjcluster->getAverageTime()) + " Charge " + ProducerUtils::str(adjcluster->getTotalCharge()) + " Purity " + ProducerUtils::str(adjcluster->getPurity()) + "\n";
          continue;
        }

        // if cluster has already been clustered, skip
        if (EvaluatedAdjCluster[*it4] == true){
          sEventCandidateFinding += "\tSkipping already evaluated cluster: NHits " + ProducerUtils::str(adjcluster->getNHits()) + " Channel " + ProducerUtils::str(adjcluster->getMainChannel()) + " Time " + ProducerUtils::str(adjcluster->getAverageTime()) + " Charge " + ProducerUtils::str(adjcluster->getTotalCharge()) + " Purity " + ProducerUtils::str(adjcluster->getPurity()) + "\n";
          continue;
        }

        auto ref4 = TVector3(0, cluster->getY(), cluster->getZ()) - TVector3(dXcluster2, adjcluster->getY(), adjcluster->getZ());
        if (ref4.Mag() < fAdjClusterRad)
        {
          // sEventCandidateFinding += "\tFound adjacent cluster: NHits " + ProducerUtils::str(adjcluster->getNHits()) + " Channel " + ProducerUtils::str(adjcluster->getMainChannel()) + " Time " + ProducerUtils::str(adjcluster->getAverageTime()) + " Charge " + ProducerUtils::str(adjcluster->getTotalCharge()) + "\n";
          if (adjcluster->getTotalCharge() > cluster->getTotalCharge())
          {
            sEventCandidateFinding += "¡¡¡Found bigger cluster: NHits " + ProducerUtils::str(adjcluster->getNHits()) + " Channel " + ProducerUtils::str(adjcluster->getMainChannel()) + " Time " + ProducerUtils::str(adjcluster->getAverageTime()) + " Charge " + ProducerUtils::str(adjcluster->getTotalCharge()) + " Purity " + ProducerUtils::str(adjcluster->getPurity()) + "\n";
            EvaluatedCluster[*it] = false;
            PreselectionCluster = false;

            // Reset the EvaluatedCluster values for the clusters that have been added to the cluster
            for (auto it5 = AdjClusterIdx.begin(); it5 != AdjClusterIdx.end(); ++it5)
            {
              sEventCandidateFinding += "---Removing cluster: NHits " + ProducerUtils::str(SolarClusterVector[*it5]->getNHits()) + " Channel " + ProducerUtils::str(SolarClusterVector[*it5]->getMainChannel()) + " Time " + ProducerUtils::str(SolarClusterVector[*it5]->getAverageTime()) + " Charge " + ProducerUtils::str(SolarClusterVector[*it5]->getTotalCharge()) + " Purity " + ProducerUtils::str(SolarClusterVector[*it5]->getPurity()) + "\n";
              EvaluatedCluster[*it5] = false;
              EvaluatedAdjCluster[*it5] = false; // Reset the evaluated status for the adjacent clusters
            }
            break;
          }
          else {
            sEventCandidateFinding += "+++Adding cluster: NHits " + ProducerUtils::str(adjcluster->getNHits()) + " Channel " + ProducerUtils::str(adjcluster->getMainChannel()) + " Time " + ProducerUtils::str(adjcluster->getAverageTime()) + " Charge " + ProducerUtils::str(adjcluster->getTotalCharge()) + " Purity " + ProducerUtils::str(adjcluster->getPurity()) + "\n";
            AdjClusterVec.push_back(adjcluster);
            AdjClusterIdx.push_back(*it4);
            EvaluatedCluster[*it4] = true;
            EvaluatedAdjCluster[*it4] = true; // Mark this adjacent cluster as evaluated
          }
        }
      }

      if (PreselectionCluster == false)
      {
        sEventCandidateFinding += "\tSkipping non-preselected cluster\n";
        EvaluatedCluster[*it] = false;
        EventCandidateFound.push_back(false);
        continue; // Skip if the cluster has been reset
      }

      // Store the original indices of the clustered clusters
      sEventCandidateFinding += "***Preselection cluster: NHits " + ProducerUtils::str(cluster->getNHits()) + " Channel " + ProducerUtils::str(cluster->getMainChannel()) + " Time " + ProducerUtils::str(cluster->getAverageTime()) + " Charge " + ProducerUtils::str(cluster->getTotalCharge()) + " AdjClNum " + ProducerUtils::str(AdjClusterVec.size()-1) + "\n";
      EventCandidateVector.push_back(std::move(AdjClusterVec));
      EventCandidateIdx.push_back(std::move(AdjClusterIdx));
      EventCandidateFound.push_back(true);
    }

    // Check that the EventCandidateFound vector is the same size as the SolarClusterVector
    if (EventCandidateFound.size() != SolarClusterVector.size()) {
      producer->PrintInColor("Error: EventCandidateFound vector size does not match SolarClusterVector size", ProducerUtils::GetColor("red"), "Error");
    }

    sEventCandidateFinding += "Total number of event candidates found: " + ProducerUtils::str(EventCandidateVector.size()) + "\n";
    ProducerUtils::PrintInColor(sEventCandidateFinding, ProducerUtils::GetColor("blue"), "Debug");
    return;
  }

  int LowEUtils::MatchPDSFlash(
    const std::vector<art::Ptr<solar::LowECluster>> &SolarClusterVector,
    const std::vector<art::Ptr<recob::OpFlash>> &PDSFlashes,
    const detinfo::DetectorClocksData &clockData,
    const art::Event &evt,
    bool debug)
  {
    std::string sFlashMatching = "LowEUtils::MatchPDSFlash " + ProducerUtils::str(SolarClusterVector.size()) + " clusters and " + ProducerUtils::str(PDSFlashes.size()) + " flashes found in the event\n";
    // get drift properties
    art::ServiceHandle<geo::Geometry> geom;
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    double driftLength = 0;
    double driftTime = 0; // in cm/us, but we will use it as a time variable
    
    // Get geometry from the geometry service
    std::string fGeometry = "";
    std::string geoName = geom->DetectorName();
    // Get the drift length from the geometry service
    geo::CryostatID c(0);
    const geo::CryostatGeo &cryostat = geom->Cryostat(c);
    const geo::TPCGeo &tpcg = cryostat.TPC(0);
    
    if (geoName.find("dune10kt") != std::string::npos)
    {
        driftLength = tpcg.DriftDistance();
        driftTime = driftLength / detProp.DriftVelocity(); // in us
        fGeometry = "HD";
    }
    else if(geoName.find("dunevd10kt") != std::string::npos)
    {
        driftLength = tpcg.DriftDistance();
        driftTime = driftLength / detProp.DriftVelocity();
        fGeometry = "VD";
    }
    else
    {
        producer->PrintInColor("Unknown geometry: " + geoName, ProducerUtils::GetColor("red"));
        return -1; // Unknown geometry
    }

    sFlashMatching += "Drift length: " + ProducerUtils::str(driftLength) + " cm, Drift time: " + ProducerUtils::str(driftTime) + " us, Geometry: " + fGeometry + "\n";
    
    // Match according to the main cluster time and position
    if (SolarClusterVector.empty() || PDSFlashes.empty()) {
        return -1; // No clusters or flashes to match
    }
    // Assuming the first cluster is the main one
    const auto &mainCluster = SolarClusterVector.front();
    // Get the time and position of the main cluster
    double clusterTime = mainCluster->getAverageTime();
    double clusterX = -1e6;    
    double clusterY = mainCluster->getY();
    double clusterZ = mainCluster->getZ();
    
    int matchedFlashIndex = -1;
    double matchedFlashPE = 0.0;
    double matchedRecoX = -1e6;
    // Loop through the flashes to find the best match
    for (int i = 0; i < int(PDSFlashes.size()); i++)
    {
        // Reset cluster X for each flash
        const auto &flash = PDSFlashes[i];
        double flashPE = flash->TotalPE();
        // double flashR = -1e6;
        // Check if the flash time is within the acceptable range
        double flashTime = flash->Time();
        double dT = clusterTime - flashTime; // Time difference between the cluster and the flash
        // If dT is bigger that drift time, skip this flash
        if (std::abs(dT) > driftTime){
            continue; // Skip this flash if it's too far in time
        }
        
        if (int(flash->PEs().size()) < fAdjOpFlashMinNHitCut || flash->TotalPE() < fAdjOpFlashMinPECut || *std::max_element(flash->PEs().begin(), flash->PEs().end()) / flash->TotalPE() > fAdjOpFlashMaxPERatioCut)
        {
            continue; // Skip flashes with insufficient charge or hits
        }

        // Calculate the distance in space
        double flashY = flash->YCenter();
        double flashZ = flash->ZCenter();
        // if (!producer) {
        //     fhicl::ParameterSet pset; // Create a parameter set (you may need to initialize it properly)
        //     producer = std::make_unique<producer::ProducerUtils>(pset);
        // }
        producer->ComputeDistanceX(clusterX, clusterTime, flashTime, driftLength, driftTime);
        if (fGeometry == "HD")
        {
            // For DUNE 10kt geometry, we have different projections based on the plane
            // Change the sign of clusterX for the collection plane
            if (flash->Frame() == 0) // Collection plane
            {
                clusterX = -clusterX; // Convert to the collection plane coordinate system
            }
            if (pow(clusterY - flashY, 2) / pow(fAdjOpFlashY, 2) + pow(clusterZ - flashZ, 2) / pow(fAdjOpFlashZ, 2) > 1)
            {
                continue;
            }
            // flashR = sqrt(pow(clusterY - flashY, 2) + pow(clusterZ - flashZ, 2));
        }
        else if (fGeometry == "VD")
        {
          // Convert clusterX to the VD geometry [-driftLength/2, driftLength/2]
          if (flash->Frame() == 0) // Cathode flashes
          {
              if (pow(clusterY - flashY, 2) / pow(fAdjOpFlashY, 2) + pow(clusterZ - flashZ, 2) / pow(fAdjOpFlashZ, 2) > 1)
              {
                  continue;
              }
              // flashR = sqrt(pow(clusterY - flashY, 2) + pow(clusterZ - flashZ, 2));
          }
          else if (flash->Frame() == 1 || flash->Frame() == 2) // Membrane flashes
          {
              if (fAdjOpFlashMembraneProjection){
                  if (pow(clusterX, 2) / pow(fAdjOpFlashX, 2) + pow(clusterZ - flashZ, 2) / pow(fAdjOpFlashZ, 2) > 1){
                      continue;
                  }
              }
              else{
                  if (pow(clusterX, 2) / pow(fAdjOpFlashX, 2) + pow(clusterY - flashY, 2) / pow(fAdjOpFlashY, 2) + pow(clusterZ - flashZ, 2) / pow(fAdjOpFlashZ, 2) > 1){
                      continue;
                  }
              }
              // flashR = sqrt(pow(clusterX, 2) + pow(clusterZ - flashZ, 2));
          } 
          else if (flash->Frame() == 3 || flash->Frame() == 4) // End-Cap flashes
          {
              if (fAdjOpFlashEndCapProjection){
                  if (pow(clusterX, 2) / pow(fAdjOpFlashX, 2) + pow(clusterY - flashY, 2) / pow(fAdjOpFlashY, 2) > 1){
                      continue;
                  }
              }
              else{
                  if (pow(clusterX, 2) / pow(fAdjOpFlashX, 2) + pow(clusterY - flashY, 2) / pow(fAdjOpFlashY, 2) + pow(clusterZ - flashZ, 2) / pow(fAdjOpFlashZ, 2) > 1){
                      continue;
                  }
              }
              // flashR = sqrt(pow(clusterX, 2) + pow(clusterY - flashY, 2) + pow(clusterZ - flashZ, 2));
          } 
          clusterX = driftLength / 2 - clusterX; // Convert to the VD geometry coordinate system
        }
        else
        {
            continue; // Unknown geometry, skip this matching
        }

        if (flashPE > matchedFlashPE || matchedFlashIndex == -1) {
            // If this flash has more PE than the previous best match, update the match
            matchedFlashPE = flashPE;
            matchedFlashIndex = i;
            matchedRecoX = clusterX;
        }
    }
    if (debug)
    {
        if (matchedFlashIndex != -1)
        {
            sFlashMatching += "Matched flash: Index " + ProducerUtils::str(matchedFlashIndex) + " PE " + ProducerUtils::str(matchedFlashPE) + " Time " + ProducerUtils::str(PDSFlashes[matchedFlashIndex]->Time()) + " Y " + ProducerUtils::str(PDSFlashes[matchedFlashIndex]->YCenter()) + " Z " + ProducerUtils::str(PDSFlashes[matchedFlashIndex]->ZCenter()) + "\n";
            // Print recon info of the cluster
            sFlashMatching += "Cluster: Time " + ProducerUtils::str(mainCluster->getAverageTime()) + " Vertex " + ProducerUtils::str(matchedRecoX) + ", " + ProducerUtils::str(mainCluster->getY()) + ", " + ProducerUtils::str(mainCluster->getZ()) + "\n";
        }
        else
        {
            sFlashMatching += "No flash matched to cluster: NHits " + ProducerUtils::str(mainCluster->getNHits()) + " Channel " + ProducerUtils::str(mainCluster->getMainChannel()) + " Time " + ProducerUtils::str(mainCluster->getAverageTime()) + " Charge " + ProducerUtils::str(mainCluster->getTotalCharge()) + "\n";
        }
        ProducerUtils::PrintInColor(sFlashMatching, ProducerUtils::GetColor("blue"), "Debug");
    }
    return matchedFlashIndex; // Return the index of the matched flash or -1 if no match found
  } // MatchPDSFlash
} // namespace lowe