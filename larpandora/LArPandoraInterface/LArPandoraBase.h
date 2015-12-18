/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraBase.h
 *
 *  @brief  Base producer module for reconstructing recob::PFParticles from recob::Hits
 *
 */

#ifndef LAR_PANDORA_BASE_H
#define LAR_PANDORA_BASE_H 1

// Framework Includes
#include "art/Framework/Core/EDProducer.h"

// Pandora includes
#include "Api/PandoraApi.h"

// Local LArPandora includes
#include "LArPandoraInterface/LArPandoraCollector.h"

// std includes
#include <string>

namespace pandora {class Pandora;}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  LArPandoraBase class
 */
class LArPandoraBase : public art::EDProducer
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    LArPandoraBase(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~LArPandoraBase();

    virtual void beginJob();

protected:

    typedef std::map< int, const pandora::Pandora* > PandoraInstanceMap;
    typedef std::map< const pandora::Pandora*, std::vector<int> > PandoraAddressList;
  
    /**
     *  @brief Create Pandora instances and register Pandora geometry
     */
    virtual void ConfigurePandoraGeometry() const = 0;

    /**
     *  @brief Return an identification number for the drift volume
     *
     *  @param cstat the ID number of the cryostat
     *  @param tpc the ID number of the TPC
     */
    virtual unsigned int GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const = 0;

    /**
     *  @brief Register the Pandora algorithms, helper functions and geometry
     */
    void InitializePandora();

    /**
     *  @brief Create the Pandora 2D hits from the ART hits
     *
     *  @param hits  the input list of ART hits for this event
     *  @param hitMap  mapping from Pandora hit addresses to ART hits
     */
    void CreatePandoraHits2D(const HitVector &hitVector, HitMap &hitMap) const;

    /**
     *  @brief Create the Pandora 3D hits from the ART space points
     *
     *  @param spacepoints  the input list of ART spacepoints for this event
     *  @param spacePointsToHits  the mapping between ART space points and ART hits
     *  @param spacePointMap  mapping from Pandora hit addresses to ART space points
     */
    void CreatePandoraHits3D(const SpacePointVector &spacePointVector, const SpacePointsToHits &spacePointsToHits,
        SpacePointMap &spacePointMap) const;

    /**
     *  @brief Create the Pandora MC particles from the MC particles
     *
     *  @param truthToParticles  mapping from MC truth to MC particles
     *  @param particlesToTruth  mapping from MC particles to MC truth
     */
    void CreatePandoraParticles(const MCTruthToMCParticles &truthToParticles, 
        const MCParticlesToMCTruth &particlesToTruth) const;

    /**
     *  @brief Create 2D projections of the Pandora MC particles
     *
     *  @param particleVector  the input vector of MC particles
     */
    void CreatePandoraParticles2D(const MCParticleVector &particleVector) const;

    /**
     *  @brief Create links between the 2D hits and Pandora MC particles
     *
     *  @param hitMap  mapping from Pandora hit addresses to ART hits
     *  @param hitToParticleMap  mapping from each ART hit to its underlying G4 track ID
     */
    void CreatePandoraLinks2D(const HitMap &hitMap, const HitsToTrackIDEs &hitToParticleMap) const;

    /**
     *  @brief Create the Pandora instances
     */
    void CreatePandoraInstances();

    /**
     *  @brief Process the event using each Pandora instance
     */
    void RunPandoraInstances() const;

    /**
     *  @brief Reset each Pandora instance
     */
    void ResetPandoraInstances() const;

    /**
     *  @brief Delete each Pandora instance
     */
    void DeletePandoraInstances();

    /**
     *  @brief Loop over MC trajectory points and identify start and end points within the detector
     *
     *  @param volumeID  the drift volume
     *  @param particle  the true particle
     *  @param startT  the first trajectory point in the detector
     *  @param endT  the last trajectory point in the detector
     */
    void GetTrueStartAndEndPoints(const unsigned int volumeID, const art::Ptr<simb::MCParticle> &particle,
        int &startT, int& endT) const;

    /**
     *  @brief Loop over MC trajectory points and identify start and end points within a given cryostat and TPC
     *
     *  @param cstat  the cryostat
     *  @param tpc  the TPC
     *  @param particle  the true particle
     *  @param startT  the first trajectory point in the detector
     *  @param endT  the last trajectory point in the detector
     */
    void GetTrueStartAndEndPoints(const unsigned int cstat, const unsigned int tpc,
        const art::Ptr<simb::MCParticle> &particle, int &startT, int& endT) const;

    /**
     *  @brief Use detector and time services to get a true X offset for a given trajectory point
     *
     *  @param particle  the true particle
     *  @param nT  the trajectory point
     */
    float GetTrueX0(const art::Ptr<simb::MCParticle> &particle, const int nT) const;

    /**
     *  @brief Convert charge in ADCs to approximate MIPs
     *
     *  @param hit_Charge  the input charge
     *  @param hit_View  the input view number
     */
    double GetMips(const double hit_Charge, const geo::View_t hit_View) const;


    bool                 m_enableProduction;      ///<
    bool                 m_enableMCParticles;     ///<
    bool                 m_enableMonitoring;      ///<

    std::string          m_configFile;            ///<
    std::string          m_geantModuleLabel;      ///<
    std::string          m_hitfinderModuleLabel;  ///<
    std::string          m_spacepointModuleLabel; ///<
    std::string          m_pandoraModuleLabel;    ///<

    bool                 m_useHitWidths;          ///<
    int                  m_uidOffset;             ///<

    double               m_dx_cm;                 ///<
    double               m_int_cm;                ///<
    double               m_rad_cm;                ///< 
    double               m_dEdX_max;              ///< 
    double               m_dEdX_mip;              ///< 
    double               m_mips_to_gev;           ///< 
    double               m_recombination_factor;  ///<

    PandoraInstanceMap   m_pandoraInstanceMap;    ///<
};

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_BASE_H
