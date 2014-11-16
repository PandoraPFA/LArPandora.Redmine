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

// Local includes
#include "LArPandoraCollector.h"
#include "PFParticleSeed.h"

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

    typedef std::map< art::Ptr<recob::PFParticle>, const PFParticleSeed > PFParticleSeedMap;
    typedef std::map< art::Ptr<recob::PFParticle>, const unsigned int >   PFParticleVolumeMap;

    typedef std::set< art::Ptr<recob::PFParticle> >                       PFParticleList;
    typedef std::map< art::Ptr<recob::PFParticle>, PFParticleList >       PFParticleMergeMap;

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
     *  @param particleMerges  map of PFParticles to be merged together in the new PFParticle output
     *  @param particlesToClusters  mapping between old PFParticles and their clusters
     *  @param particlesToSpacePoints  mapping between old PFParticles and their space points
     */
    void ProduceArtOutput(art::Event &evt, const PFParticleMap &particleMap, const PFParticleMergeMap &particleMerges, 
        const PFParticlesToClusters &particlesToClusters, const PFParticlesToSpacePoints &particlesToSpacePoints);
    
    /**
     *  @brief Get drift volume ID code from cryostat and TPC number
     *  
     *  @param cstat  the cryostat number
     *  @param tpc  the TPC number
     */
    virtual unsigned int GetVolumeID(const unsigned int cstat, const unsigned int tpc) const = 0;

    /**
     *  @brief  Get drift volume ID code for a set of space points
     *  
     *  @param spacePoints  the input vector of space points
     *  @param spacePointsToHits  mapping between space points and 2D hits (which hold the wire ID)
     */
    unsigned int GetVolumeID(const SpacePointVector &spacePoints, const SpacePointsToHits &spacePointsToHits) const;

    /**
     *  @brief Match PFParticles and populate the matrix of particle associations
     *  
     *  @param particleVector  the input vector of PFParticles 
     *  @param particleVolumeMap  mapping between PFParticles and their drift volume ID codes
     *  @param particleSeedMap  mapping betwen PFParticles and their particle seed objects
     *  @param associationMatrix  matrix of associations between pairs of PFParticles
     */
    void CreateParticleMatches(const PFParticleVector &particleVector, const PFParticleVolumeMap &particleVolumeMap,
        const PFParticleSeedMap &particleSeedMap, ParticleAssociationMatrix &associationMatrix) const;

    /**
     *  @brief Match PFParticles and populate the matrix of particle associations
     *  
     *  @param particle1  the first PFParticle
     *  @param particle2  the second PFParticle
     *  @param particleSeedMap  the mapping between PFParticles and their particle seed objects
     *  @param particleAssociationMatrix  matrix of associations between pairs of PFParticles
     */
    void CreateParticleMatches(const art::Ptr<recob::PFParticle> particle1, const art::Ptr<recob::PFParticle> particle2,
        const PFParticleSeedMap &particleSeedMap, ParticleAssociationMatrix &particleAssociationMatrix) const;

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
     *  @param particleSeedMap  mapping between particles and particle seeds
     */
    void WriteParticleMatches(const PFParticleMergeMap &particleMergeMap, const PFParticleSeedMap &particleSeedMap);

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
     *  @brief Collect connected sets of PFParticle merges from a map of associations
     *
     *  @param seedParticle  the first PFParticle in the set
     *  @param currentParticle  the current PFParticle in the set
     *  @param particleMergeMap  the input map of PFParticle associations
     *  @param vetoList  the list of PFParticles already collected
     *  @param associatedParticleList  the running list of PFParticles associated to the first PFParticle via the association map
     */
    void CollectAssociatedParticles(art::Ptr<recob::PFParticle> seedParticle, art::Ptr<recob::PFParticle> currentParticle, 
        const PFParticleMergeMap &particleMergeMap, const PFParticleList &vetoList, PFParticleList &associatedParticleList) const;

    /**
     *  @brief Get closest pair of vertices for two seed particles
     *
     *  @param seed1  the first seed particle
     *  @param seed2  the second seed particle
     *  @param vertexType1  the closest vertex from the first seed particle
     *  @param vertexType2  the closest vertex from the second seed particle
     */
    void GetClosestVertices(const PFParticleSeed& seed1, const PFParticleSeed &seed2, 
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

    bool          m_useXcoordinate;                   ///<
    float         m_minCosRelativeAngle;              ///<
    float         m_maxLongitudinalDisplacementX;     ///<
    float         m_maxTransverseDisplacement;        ///<

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
