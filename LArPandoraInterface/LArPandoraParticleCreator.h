/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraParticleCreator.h
 *
 *  @brief  Base producer module for reconstructing recob::PFParticles from recob::Hits
 *
 */

#ifndef LAR_PANDORA_PARTICLE_CREATOR_H
#define LAR_PANDORA_PARTICLE_CREATOR_H 1

// Local includes
#include "LArPandoraBase.h"

// ROOT includes
#include "TTree.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  LArPandoraParticleCreator class
 */
class LArPandoraParticleCreator : public LArPandoraBase
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    LArPandoraParticleCreator(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~LArPandoraParticleCreator();

    void beginJob();
    void produce(art::Event &evt);

protected:

    /**
     *  @brief Convert the Pandora PFOs into ART clusters and write into ART event
     *
     *  @param evt  the ART event
     *  @param hitMap  mapping from Pandora hit ID to ART hit
     */
    void ProduceArtOutput(art::Event &evt, const HitMap &hitMap) const;

    /**
     *  @brief Initialize the internal monitoring
     */
    void InitializeMonitoring();

    bool      m_buildTracks;           ///<
    bool      m_buildShowers;          ///<

    TTree    *m_pRecoTree;             ///<
    int       m_run;                   ///<
    int       m_event;                 ///<
    int       m_hits;                  ///<
    float     m_collectionTime;        ///<
    float     m_inputTime;             ///<
    float     m_processTime;           ///<
    float     m_outputTime;            ///<
};

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_PARTICLE_CREATOR_H
