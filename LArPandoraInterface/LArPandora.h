/**
 *  @file   LArPandora/LArPandora.h
 * 
 *  @brief  Header file for the lar pandora producer.
 * 
 *  $Log: $
 */

#ifndef LAR_PANDORA_H
#define LAR_PANDORA_H

// Framework Includes
#include "art/Framework/Core/EDProducer.h"

// LArSoft includes
#include "MCCheater/BackTracker.h"

// std includes
#include <string>

namespace pandora {class Pandora;}
namespace recob {class Hit;}

class TTree;

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  LArPandora class
 */
class LArPandora : public art::EDProducer
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    LArPandora(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~LArPandora();

    void beginJob();
    void endJob();
    void produce(art::Event &evt);
    void reconfigure(fhicl::ParameterSet const &pset);

private:
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
    void InitializePandora() const;

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
     *  @brief Create the Pandora hits from the ART hits
     * 
     *  @param hits  the ART hits for this event
     *  @param hitMap  mapping from Pandora hit ID to ART hit 
     */
    void CreatePandoraHits(const HitVector &hits, HitMap &hitMap) const;

    /**
     *  @brief Create the Pandora MC particles from the MC particles
     * 
     *  @param particleMap  mapping from each G4 track ID to its corresponding ART MC particle 
     *  @param truthMap  mapping from each G4 track ID to its corresponding ART MC primary particle
     */
    void CreatePandoraParticles(const ParticleMap &particleMap, const TruthToParticleMap &truthToParticleMap) const;

    /**
     *  @brief Create the Pandora hit-particle links
     * 
     *  @param hitMap  mapping from Pandora hit ID to ART hit  
     *  @param hitToParticleMap  mapping from each ART hit to its underlying G4 track ID
     */
    void CreatePandoraLinks(const HitMap &hitMap, const HitToParticleMap &hitToParticleMap) const;

    /**
     *  @brief Process the event using Pandora
     */
    void RunPandora() const;

    /**
     *  @brief Convert the Pandora PFOs into ART clusters and write into ART event
     * 
     *  @param evt  the ART event  
     *  @param hitMap  mapping from Pandora hit ID to ART hit   
     */
    void ProduceArtClusters(art::Event &evt, const HitMap &hitMap) const;

    /**
     *  @brief Reset Pandora
     */
    void ResetPandora() const; 

    /**
     *  @brief Initialize the internal monitoring
     */
    void InitializeMonitoring();


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

#endif  // LAR_PANDORA_H

