/**
 *  @file   larpandora/LArPandoraInterface/PFParticleStitcher.h
 *
 *  @brief  Producer module for stitching together recob::PFParticles
 *
 */
#ifndef PFPARTICLE_STITCHER_H
#define PFPARTICLE_STITCHER_H 1

// Framework Includes
#include "art/Framework/Core/EDProducer.h"

// Pandora includes
#include "Objects/CartesianVector.h"
#include "larpandoracontent/LArObjects/LArTrackPfo.h"

// Local LArPandora includes
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// ROOT includes
#include "TTree.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  PFParticleStitcher class
 */
class PFParticleStitcher : public art::EDProducer
{
public:    

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    PFParticleStitcher(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~PFParticleStitcher();

    virtual void beginJob();
    virtual void endJob();
    virtual void produce(art::Event &evt);

protected:

    /**
     *  @brief  PFParticleTrack class
     */
    class PFParticleTrack
    {
    public:

        /**
         *  @brief  Constructor
         *
         *  @param  SpacePointVector
         */
        PFParticleTrack(const art::Ptr<recob::Track> track);

        /**
         *  @brief  Destructor
         */
        ~PFParticleTrack();

        /**
         *  @brief  Return vertex position 
         */
        pandora::CartesianVector GetInnerPosition() const;

        /**
         *  @brief  Return vertex direction 
         */
        pandora::CartesianVector GetInnerDirection() const;

        /**
         *  @brief  Return end position
         */
        pandora::CartesianVector GetOuterPosition()  const;

        /**
         *  @brief  Return end direction
         */
        pandora::CartesianVector GetOuterDirection() const;

    private:

        pandora::CartesianVector       m_innerPosition;    ///<
        pandora::CartesianVector       m_innerDirection;   ///<
        pandora::CartesianVector       m_outerPosition;    ///<
        pandora::CartesianVector       m_outerDirection;   ///<
    };

    /**
     *  @brief  ParticleAssociation class
     */
    class ParticleAssociation
    {
    public:
        /**
         *  @brief  Vertex enumeration
         */
        enum VertexType
        {
            UNDEFINED = 0,
            INNER     = 1,
            OUTER     = 2
        };

        /**
         *  @brief  Constructor
         *
         *  @param  parent
         *  @param  daughter
         *  @param  fom
         */
        ParticleAssociation(const VertexType parent, const VertexType daughter, const float fom);

        /**
         *  @brief  Get parent
         *
         *  @return the parent
         */
        VertexType GetParent() const;

        /**
         *  @brief  Get daughter
         *
         *  @return the daughter
         */
        VertexType GetDaughter() const;

        /**
         *  @brief  Get figure of merit
         *
         *  @return the figure of merit
         */
        float GetFigureOfMerit() const;

    private:
        VertexType      m_parent;           ///<
        VertexType      m_daughter;         ///<
        float           m_fom;              ///<
    };

    typedef std::map<const art::Ptr<recob::PFParticle>, ParticleAssociation> ParticleAssociationMap;
    typedef std::map<const art::Ptr<recob::PFParticle>, ParticleAssociationMap> ParticleAssociationMatrix;

    typedef std::map< art::Ptr<recob::PFParticle>, const PFParticleTrack > PFParticleTrackMap;
    typedef std::map< art::Ptr<recob::PFParticle>, const unsigned int >   PFParticleVolumeMap;

    typedef std::set< art::Ptr<recob::PFParticle> >                       PFParticleList;
    typedef std::map< art::Ptr<recob::PFParticle>, PFParticleList >       PFParticleMergeMap;

    typedef std::map< art::Ptr<recob::PFParticle>, lar_content::LArTrackStateVector > PFParticleTrajectoryMap;

    /**
     *  @brief Configure this module
     *  
     *  @param pset  the set of input parameters
     */
    void reconfigure(fhicl::ParameterSet const &pset);

    /**
     *  @brief Initialize the internal monitoring
     */
    void InitializeMonitoring();

    /**
     *  @brief Produce the ART output
     *  
     *  @param evt  the ART event
     *  @param particleMap  mapping between the old PFParticles and their particle ID codes
     *  @param particleMerges  map of PFParticles to be merged together in new PFParticle output
     *  @param particlesTrajectories  map of Trajectory Points to be merged together in new PFParticle output 
     *  @param particlesToVertices  mapping between old PFParticles to their vertices
     *  @param particlesToClusters  mapping between old PFParticles and their clusters
     *  @param particlesToSpacePoints  mapping between old PFParticles and their space points
     *  @param particlesToHits  mapping between old PFParticles and their hits
     */
    void ProduceArtOutput(art::Event &evt, const PFParticleMap &particleMap, const PFParticleMergeMap &particleMerges,
	const PFParticleTrajectoryMap &particleTrajectories, const PFParticlesToVertices &particlesToVertices,
        const PFParticlesToClusters &particlesToClusters, const PFParticlesToSpacePoints &particlesToSpacePoints, 
        const PFParticlesToHits &particlesToHits);
    
    /**
     *  @brief Get drift volume ID code from cryostat and TPC number
     *  
     *  @param cstat  the cryostat number
     *  @param tpc  the TPC number
     */
    virtual unsigned int GetVolumeID(const unsigned int cstat, const unsigned int tpc) const = 0;

    /**
     *  @brief  Get drift volume ID code for a set of hts
     *  
     *  @param hitVector  the input vector of hits
     */
    unsigned int GetVolumeID(const HitVector &hitVector) const;

    /**
     *  @brief  Get drift volume center for a set of hits
     *  
     *  @param hitVector  the input vector of hits
     */
    TVector3 GetVolumeCenter(const HitVector &hitVector) const;

    /**
     *  @brief Match PFParticles and populate the matrix of particle associations
     *  
     *  @param particleVector  the input vector of PFParticles 
     *  @param particleVolumeMap  mapping between PFParticles and their drift volume ID codes
     *  @param particleTrackMap  mapping betwen PFParticles and their particle track objects
     *  @param associationMatrix  matrix of associations between pairs of PFParticles
     */
    void CreateParticleMatches(const PFParticleVector &particleVector, const PFParticleVolumeMap &particleVolumeMap,
        const PFParticleTrackMap &particleTrackMap, ParticleAssociationMatrix &associationMatrix) const;

    /**
     *  @brief Match PFParticles and populate the matrix of particle associations
     *  
     *  @param particle1  the first PFParticle
     *  @param particle2  the second PFParticle
     *  @param particleTrackMap  the mapping between PFParticles and their particle track objects
     *  @param particleAssociationMatrix  matrix of associations between pairs of PFParticles
     */
    void CreateParticleMatches(const art::Ptr<recob::PFParticle> particle1, const art::Ptr<recob::PFParticle> particle2,
        const PFParticleTrackMap &particleTrackMap, ParticleAssociationMatrix &particleAssociationMatrix) const;

    /**
     *  @brief Select PFParticle associations that will become merges
     *  
     *  @param particleMap  mapping between PFParticles and their particle ID codes
     *  @param particleAssociationMatrix  matrix of associations between pairs of PFParticles
     *  @param particleMergeMap  the output map of PFParticle merges
     */
    void SelectParticleMatches(const PFParticleMap &particleMap, const ParticleAssociationMatrix &particleAssociationMatrix,
        PFParticleMergeMap &particleMergeMap) const;

    /**
     *  @brief Write the selected particles matches to a ROOT file
     *
     *  @param particleMergeMap  mapping between matched particles
     *  @param particleTrackMap  mapping between particles and particle tracks
     */
    void WriteParticleMatches(const PFParticleMergeMap &particleMergeMap, const PFParticleTrackMap &particleTrackMap);

    /**
     *  @brief Collate the PFParticle merges
     *
     *  @param particleVector  the input vector of PFParticle objects
     *  @param inputMergeMap  the input matrix of PFParticle merges
     *  @param outputMergeMap  the collated map of PFParticle merges
     */
    void SelectParticleMerges(const PFParticleVector &particleVector, const PFParticleMergeMap &inputMergeMap,
        PFParticleMergeMap &outputMergeMap) const;

    /**
     *  @brief Order the PFParticle merges
     *
     *  @param particleMerges  the collated map of PFParticle merges
     *  @param particlesToTracks  the mapping between PFParticles and their tracks
     *  @param particlesToHits  the mapping between PFParticles and their hits
     */
    void OrderParticleMerges(const PFParticleMergeMap &particleMerges, const PFParticlesToTracks &particlesToTracks, 
        const PFParticlesToHits &particlesToHits, PFParticleTrajectoryMap &particleTrajectories) const;

    /**
     *  @brief Collect connected sets of PFParticle merges from a map of associations
     *
     *  @param trackParticle  the first PFParticle in the set
     *  @param currentParticle  the current PFParticle in the set
     *  @param particleMergeMap  the input map of PFParticle associations
     *  @param vetoList  the list of PFParticles already collected
     *  @param associatedParticleList  the running list of PFParticles associated to the first PFParticle via the association map
     */
    void CollectAssociatedParticles(art::Ptr<recob::PFParticle> trackParticle, art::Ptr<recob::PFParticle> currentParticle, 
        const PFParticleMergeMap &particleMergeMap, const PFParticleList &vetoList, PFParticleList &associatedParticleList) const;

    /**
     *  @brief Get closest pair of vertices for two track particles
     *
     *  @param track1  the first track particle
     *  @param track2  the second track particle
     *  @param vertexType1  the closest vertex from the first track particle
     *  @param vertexType2  the closest vertex from the second track particle
     */
    void GetClosestVertices(const PFParticleTrack& track1, const PFParticleTrack &track2, 
        ParticleAssociation::VertexType &vertexType1, ParticleAssociation::VertexType &vertexType2) const; 

    /**
     *  @brief Calculate impact parameter between two trajectories (3D method, uses x corodinate)
     *
     *  @param initialPosition  the vertex position of the first trajectory
     *  @param initialDirection  the vertex direction of the first trajectory
     *  @param targetPosition  the vertex postition of the second trajectory
     *  @param longitudinal  the longitudinal projection from the first to the second trajectory
     *  @param transverse  the transverse projection from the first to the second trajectory
     */
    void GetImpactParameters3D(const pandora::CartesianVector &initialPosition, const pandora::CartesianVector &initialDirection, 
        const pandora::CartesianVector &targetPosition, float &longitudinal, float &transverse) const;

    /**
     *  @brief Calculate impact parameter between two trajectories (2D method, doesn't use x coordinate)
     *
     *  @param initialPosition  the vertex position of the first trajectory
     *  @param initialDirection  the vertex direction of the first trajectory
     *  @param targetPosition  the vertex position of the second trajectory
     *  @param longitudinal  the longitudinal projection from the first to the second trajectory
     *  @param transverse  the transverse projection from the first to the second trajectory
     */
    void GetImpactParameters2D(const pandora::CartesianVector &initialPosition, const pandora::CartesianVector &initialDirection, 
        const pandora::CartesianVector &targetPosition, float &longitudinal, float &transverse) const;

    /**
     *  @brief Calculate longitudinal displacement cut from maximum allowed displacement in X
     *
     *  @param direction  the input trajectory
     */
    float GetLongitudinalDisplacementCut3D(const pandora::CartesianVector &direction) const; 

    /**
     *  @brief Calculate longitudinal displacement cut from maximum allowed displacement in X
     *
     *  @param direction  the input trajectory
     */
    float GetLongitudinalDisplacementCut2D(const pandora::CartesianVector &direction) const;

    /**
     *  @brief Calculate displacement in x coordinate between two trajectories
     *
     *  @param initialPosition  the vertex position of the first trajectory
     *  @param initialDirection  the vertex direction of the first trajectory
     *  @param targetPosition  the vertex position of the second trajectory
     */
    float GetDeltaX(const pandora::CartesianVector &initialPosition, const pandora::CartesianVector &initialDirection, 
        const pandora::CartesianVector &targetPosition) const;

    bool          m_enableMonitoring;                 ///<
    std::string   m_particleLabel;                    ///<
    std::string   m_trackLabel;                       ///<

    bool          m_enableStitching;                  ///<
    bool          m_useXcoordinate;                   ///<
    float         m_minCosRelativeAngle;              ///<
    float         m_maxLongitudinalDisplacementX;     ///<
    float         m_maxTransverseDisplacement;        ///<
    float         m_relaxCosRelativeAngle;            ///<
    float         m_relaxTransverseDisplacement;      ///<

    TTree        *m_pRecoTree;                        ///<
    int           m_run;                              ///<
    int           m_event;                            ///<
    int           m_particle1;                        ///<
    int           m_particle2;                        ///<
    float         m_cosRelativeAngle;                 ///<
    float         m_longitudinalDisplacementCut;      ///<
    float         m_longitudinalDisplacement;         ///<
    float         m_transverseDisplacement;           ///<
    float         m_deltaX;                           ///<
};

} // namespace lar_pandora

#endif // #ifndef PFPARTICLE_STITCHER_H
