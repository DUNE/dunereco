/**
 *  @file   SpacePointBuilderDUNE_module.cc
 *
 *          Class:       SpacePointBuilderDUNE
 *          Module Type: Producer
 *
 *  @brief  Producer module that builds 3D space points from recob::Hit objects via an
 *          IHit3DBuilder tool, with no clustering/pattern-recognition step.
 *
 *          Cluster3DDUNE_module.cc was originally a full 3D pattern-recognition package
 *          (clustering, merging, path finding, PCA, seeds, PFParticles, ...), and is now
 *          mostly run in its "MakeSpacePointsOnly" mode, where all of that is skipped.
 *          This module is that mode, pulled out on its own so the dead clustering code
 *          (and the output products that go with it) don't have to be carried along.
 *
 *          Every reco::ClusterHit3D the Hit3DBuilder tool makes is saved unconditionally as
 *          a recob::SpacePoint -- there is no clustering step whose completion needs to be
 *          tracked with a status bit here, unlike Cluster3DDUNE_module.cc's free-points loop
 *          (which decides what to save by checking whether MADESPACEPOINT is set, a bit that
 *          is *supposed* to mean "already saved via a formed cluster" but is also set by some
 *          Hit3DBuilder tools at triplet-creation time for an unrelated bookkeeping purpose --
 *          in MakeSpacePointsOnly mode, where clustering never runs to give that bit its
 *          intended meaning, that silently drops every genuine triplet those tools make).
 *
 *          Each saved point is separately classified as a 3-hit (triplet) or 2-hit
 *          (orphan/mythical/bad-channel) point by counting non-null entries in
 *          ClusterHit3D::getHits() -- this works for any IHit3DBuilder implementation,
 *          without relying on how it manages the MADESPACEPOINT bit internally.
 *
 *  @author usher@slac.stanford.edu
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib/cpu_timer.h"

// LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larreco/RecoAlg/Cluster3DAlgs/Cluster3D.h"
#include "larreco/RecoAlg/Cluster3DAlgs/IHit3DBuilder.h"

// ROOT includes
#include "TTree.h"

// std includes
#include <memory>
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_cluster3d {

  using RecobHitVector = art::PtrVector<recob::Hit>;

  /**
   *  @brief  Definition of the SpacePointBuilderDUNE class
   */
  class SpacePointBuilderDUNE : public art::EDProducer {
  public:
    explicit SpacePointBuilderDUNE(fhicl::ParameterSet const& pset);

  private:
    void beginJob() override;
    void produce(art::Event& evt) override;

    /**
     *  @brief Initialize the internal monitoring
     */
    void InitializeMonitoring();

    /**
     *   Algorithm parameters
     */
    bool m_enableMonitoring; ///< Turn on monitoring of this algorithm

    /**
     *   Tree variables for output
     */
    TTree* m_pRecoTree;      ///<
    int m_run;               ///<
    int m_event;             ///<
    int m_hits;              ///< Number of 2D hits used to build space points
    int m_hits3D;            ///< Number of 3D space points made
    int m_hits3DTriplet;     ///< Number of those which are genuine 3-plane triplets
    int m_hits3DPair;        ///< Number of those which are 2-hit orphan/mythical/bad-channel points
    float m_artHitsTime;     ///< Time to recover/organize the input 2D hits
    float m_makeHitsTime;    ///< Time to build the 3D space points
    float m_saveTime;        ///< Time to convert 3D hits into recob::SpacePoint output

    // Algorithm
    std::unique_ptr<lar_cluster3d::IHit3DBuilder> m_hit3DBuilderAlg; ///<  Builds the 3D hits to operate on
  };

  DEFINE_ART_MODULE(SpacePointBuilderDUNE)

} // namespace lar_cluster3d

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_cluster3d {

  SpacePointBuilderDUNE::SpacePointBuilderDUNE(fhicl::ParameterSet const& pset) : EDProducer{pset}
  {
    m_enableMonitoring = pset.get<bool>("EnableMonitoring", false);

    m_hit3DBuilderAlg = art::make_tool<lar_cluster3d::IHit3DBuilder>(
      pset.get<fhicl::ParameterSet>("Hit3DBuilderAlg"));

    // Handle special case of Space Point building outputting a new hit collection
    m_hit3DBuilderAlg->produces(producesCollector());

    produces<std::vector<recob::SpacePoint>>();
    produces<art::Assns<recob::Hit, recob::SpacePoint>>();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void SpacePointBuilderDUNE::beginJob()
  {
    if (m_enableMonitoring) this->InitializeMonitoring();
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void SpacePointBuilderDUNE::produce(art::Event& evt)
  {
    mf::LogInfo("SpacePointBuilderDUNE") << " *** SpacePointBuilderDUNE::produce(...)  [Run="
                                         << evt.run() << ", Event=" << evt.id().event()
                                         << "] Starting Now! *** ";

    cet::cpu_timer theClockTotal;
    cet::cpu_timer theClockSave;

    if (m_enableMonitoring) theClockTotal.start();

    IHit3DBuilder::RecobHitToPtrMap hitToPtrMap;
    std::unique_ptr<reco::HitPairList> hitPairList(
      new reco::HitPairList); // Potentially lots of hits, use heap instead of stack

    // Call the algorithm that builds 3D hits and stores the (filtered) 2D hit collection
    m_hit3DBuilderAlg->Hit3DBuilder(evt, *hitPairList, hitToPtrMap);

    if (m_enableMonitoring) theClockSave.start();

    auto spacePointVec = std::make_unique<std::vector<recob::SpacePoint>>();
    auto spHitAssns = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>();

    spacePointVec->reserve(hitPairList->size());

    int nTriplet(0);
    int nPair(0);

    // Every 3D hit the builder made gets saved here -- see the file-level comment for why this
    // doesn't (and shouldn't) check MADESPACEPOINT to decide what to save.
    for (auto& hitPair : *hitPairList) {
      // Kept for parity with Cluster3DDUNE_module.cc's MakeAndSaveSpacePoints -- currently a
      // no-op since no Hit3DBuilder sets this, but cheap to honor if one ever does.
      if (hitPair.bitsAreSet(reco::ClusterHit3D::REJECTEDHIT)) continue;

      double spacePointPos[] = {
        hitPair.getPosition()[0], hitPair.getPosition()[1], hitPair.getPosition()[2]};
      double spacePointErr[] = {1., 0., 0., 1., 0., 1.};
      double chisq = hitPair.getHitChiSquare();

      RecobHitVector recobHits;
      int nHits(0);

      for (const auto hit : hitPair.getHits()) {
        if (!hit) continue;

        nHits++;
        recobHits.push_back(hitToPtrMap[hit->getHit()]);
      }

      if (nHits > 2)
        nTriplet++;
      else
        nPair++;

      spacePointErr[1] = hitPair.getTotalCharge();
      spacePointErr[3] = hitPair.getChargeAsymmetry();

      spacePointVec->push_back(
        recob::SpacePoint(spacePointPos, spacePointErr, chisq, spacePointVec->size()));

      if (!recobHits.empty()) {
        for (auto& hit : recobHits)
          util::CreateAssn(evt, *spacePointVec, hit, *spHitAssns);
      }
    }

    mf::LogInfo("SpacePointBuilderDUNE") << "total num space points: " << spacePointVec->size()
                                         << ", 3-hit: " << nTriplet << ", 2-hit: " << nPair;

    if (m_enableMonitoring) {
      theClockSave.stop();
      theClockTotal.stop();

      m_run = evt.run();
      m_event = evt.id().event();
      m_hits = static_cast<int>(hitToPtrMap.size());
      m_hits3D = static_cast<int>(spacePointVec->size());
      m_hits3DTriplet = nTriplet;
      m_hits3DPair = nPair;
      m_artHitsTime = m_hit3DBuilderAlg->getTimeToExecute(IHit3DBuilder::COLLECTARTHITS);
      m_makeHitsTime = m_hit3DBuilderAlg->getTimeToExecute(IHit3DBuilder::BUILDTHREEDHITS);
      m_saveTime = theClockSave.accumulated_real_time();
      m_pRecoTree->Fill();
    }

    evt.put(std::move(spacePointVec));
    evt.put(std::move(spHitAssns));

    return;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  void SpacePointBuilderDUNE::InitializeMonitoring()
  {
    art::ServiceHandle<art::TFileService> tfs;
    m_pRecoTree = tfs->make<TTree>("monitoring", "SpacePointBuilderDUNE");
    m_pRecoTree->Branch("run", &m_run, "run/I");
    m_pRecoTree->Branch("event", &m_event, "event/I");
    m_pRecoTree->Branch("hits", &m_hits, "hits/I");
    m_pRecoTree->Branch("hits3D", &m_hits3D, "hits3D/I");
    m_pRecoTree->Branch("hits3DTriplet", &m_hits3DTriplet, "hits3DTriplet/I");
    m_pRecoTree->Branch("hits3DPair", &m_hits3DPair, "hits3DPair/I");
    m_pRecoTree->Branch("artHitsTime", &m_artHitsTime, "artHitsTime/F");
    m_pRecoTree->Branch("makeHitsTime", &m_makeHitsTime, "makeHitsTime/F");
    m_pRecoTree->Branch("saveTime", &m_saveTime, "saveTime/F");
  }

} // namespace lar_cluster3d
