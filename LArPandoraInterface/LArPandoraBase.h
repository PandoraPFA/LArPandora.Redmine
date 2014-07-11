/**
 *  @file   LArPandora/LArPandoraBase.cc
 * 
 *  @brief  lar pandora base class for producer module
 * 
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"

// LArSoft includes
#include "MCCheater/BackTracker.h"
#include "RecoBase/Hit.h"

// Pandora includes
#include "Api/PandoraApi.h"

// ROOT includes
#include "TTree.h"

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

    void beginJob();
    void endJob();
    void produce(art::Event &evt);

    void reconfigure(fhicl::ParameterSet const &pset);

protected:
    typedef std::vector< art::Ptr<simb::MCParticle> > ParticleVector;
    typedef std::vector< art::Ptr<recob::Hit> >  HitVector;
    typedef std::map< art::Ptr<recob::Hit>, std::vector<cheat::TrackIDE> > HitToParticleMap;
    typedef std::map< art::Ptr<simb::MCTruth>, std::vector<int> > TruthToParticleMap;
    typedef std::map< int, art::Ptr<simb::MCParticle> > ParticleMap;
    typedef std::map< int, art::Ptr<recob::Hit> > HitMap;

    /**
     *  @brief  Event Preparation 
     * 
     *  @param  evt  the ART event 
     */
    void PrepareEvent(const art::Event &evt);  

    /**
     *  @brief Register the Pandora algorithms, helper functions and geometry
     */
    void InitializePandora();

    /**
     *  @brief  Extract the ART hits and the ART hit-particle relationships
     * 
     *  @param  evt  the ART event 
     *  @param  hits  the ART hits for this event
     *  @param  hitToParticleMap  mapping from each ART hit to its underlying G4 track ID
     */
    void CollectArtHits(const art::Event &evt, HitVector &hits, HitToParticleMap &hitToParticleMap) const;

    /**
     *  @brief Extract the ART MC particles
     * 
     *  @param evt  the ART event 
     *  @param particleMap  mapping from each G4 track ID to its corresponding ART MC particle
     *  @param truthMap  mapping from each G4 track ID to its corresponding ART MC primary particle
     */
    void CollectArtParticles(const art::Event &evt, ParticleMap &particleMap, TruthToParticleMap &truthToParticleMap) const;

    /**
     *  @brief Register the Pandora geometry
     */
    virtual void CreatePandoraGeometry() = 0;

    /**
     *  @brief Create the Pandora hits from the ART hits
     * 
     *  @param hits  the ART hits for this event
     *  @param hitMap  mapping from Pandora hit ID to ART hit 
     */
    virtual void CreatePandoraHits(const HitVector &hits, HitMap &hitMap) const = 0;

    /**
     *  @brief Create the Pandora MC particles from the MC particles
     * 
     *  @param particleMap  mapping from each G4 track ID to its corresponding ART MC particle 
     *  @param truthMap  mapping from each G4 track ID to its corresponding ART MC primary particle
     */
    virtual void CreatePandoraParticles(const ParticleMap &particleMap, const TruthToParticleMap &truthToParticleMap) const = 0;

    /**
     *  @brief Create the Pandora hit-particle links
     * 
     *  @param hitMap  mapping from Pandora hit ID to ART hit  
     *  @param hitToParticleMap  mapping from each ART hit to its underlying G4 track ID
     */
    virtual void CreatePandoraLinks(const HitMap &hitMap, const HitToParticleMap &hitToParticleMap) const = 0;

    /**
     *  @brief Convert the Pandora PFOs into ART clusters and write into ART event
     * 
     *  @param evt  the ART event  
     *  @param hitMap  mapping from Pandora hit ID to ART hit   
     */
    virtual void ProduceArtOutput(art::Event &evt, const HitMap &hitMap) const = 0;

    /**
     *  @brief Process the event using Pandora
     */
    void RunPandora() const;

    /**
     *  @brief Reset Pandora
     */
    void ResetPandora() const; 

    /**
     *  @brief Initialize the internal monitoring
     */
    void InitializeMonitoring();

    /**
     *  @brief Loop over MC trajectory points and identify start and end points within detector
     *
     *  @param particle  the true particle
     *  @param startT  the first trajectory point in the detector
     *  @param endT  the last trajectory point in the detector
     */
    void GetStartAndEndPoints(const art::Ptr<simb::MCParticle> &particle, int &startT, int& endT) const;



    bool               m_enableProduction;      ///< 
    bool               m_enableMCParticles;     ///< 
    bool               m_enableMonitoring;      ///<
    std::string        m_configFile;            ///< 
    std::string        m_geantModuleLabel;      ///< 
    std::string        m_hitfinderModuleLabel;  ///< 
    pandora::Pandora  *m_pPandora;              ///< 

    TTree             *m_pRecoTree;             ///< 
    int                m_run;                   ///< 
    int                m_event;                 ///< 
    int                m_hits;                  ///<
    float              m_time;                  ///< 

};

} // namespace lar_pandora
