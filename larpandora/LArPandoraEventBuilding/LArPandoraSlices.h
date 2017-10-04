/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraSlices.h
 *
 *  @brief  A mapping between two different reconstruction hypotheses on the same input hit collection
 */

#ifndef LAR_PANDORA_SLICES_H
#define LAR_PANDORA_SLICES_H 1

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"

#include <memory>
#include <algorithm>
#include <map>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @breif LArPandoraSlices class
 */
class LArPandoraSlices
{
public:

    /**
     *  @breif  Constructor from an art::Event
     */
    LArPandoraSlices( art::EDProducer *   pProducer,
                      art::Event *        pEvent, 
                      const std::string & crRecoProducerLabel,
                      const std::string & nuRecoProducerLabel,
                      const std::string & hitProducerLabel );

    typedef unsigned int SliceId;

    /**
     *  @brief  Get the vector of slice Id
     */
    std::vector< SliceId > GetSlices();

    /**
     *  @brief  Get the vector of PFParticles in a given slice reconstructed as a cosmic ray
     */
    std::vector< art::Ptr< recob::PFParticle > > GetSliceAsCR( const SliceId & id );

    /**
     *  @brief  Get the vector of PFParticles in a given slice reconstructed as a neutrino
     */
    std::vector< art::Ptr< recob::PFParticle > > GetSliceAsNu( const SliceId & id );

    /**
     *  @brief  Identify a given slice as being the neutrino event
     */
    void IdSliceAsNu( const SliceId & id );

    /**
     *  @brief  Produce the CR tags for the slices, and write to the event
     */
    void WriteTags();

private:
    
    // Meta data
    art::EDProducer *  m_pProducer;            ///<
    art::Event *       m_pEvent;               ///<

    // Labels
    const std::string  m_crRecoProducerLabel;  ///<
    const std::string  m_nuRecoProducerLabel;  ///<
    const std::string  m_hitProducerLabel;     ///<

    // Top-level PFParticles
    std::map< SliceId, std::vector< art::Ptr< recob::PFParticle > > > m_crSlicePFParticles;  ///<
    std::map< SliceId, std::vector< art::Ptr< recob::PFParticle > > > m_nuSlicePFParticles;  ///<

    // Neutrino Id decisions
    bool    m_doesEventContainNeutrino;  ///<
    SliceId m_nuSliceId;                 ///<

    /**
     *  @brief  Use the associations between PFParticles and hits to identify the slices
     */
    void IdentifySlices();

    /**
     *  @brief  Tag the PFParticles in the supplied vector as specified 
     */
    void WriteTag( const bool &                                          shouldTagAsNeutrino, 
                   const std::vector< art::Ptr< recob::PFParticle > > &  pfParticleVector );

    /**
     *  @brief  Get a vector of top-level PFParticles from an art handle 
     */
    void GetPrimaryPFParticles( const art::Handle< std::vector< recob::PFParticle > > &  pfParticleHandle, 
                                std::vector< art::Ptr< recob::PFParticle > > &           primaryPFParticles );

    
    /**
     *  @brief  Get mapping between PFParticles and their id
     */
    void GetIdToPFParticleMap( const art::Handle< std::vector< recob::PFParticle > > &  pfParticleHandle,
                               std::map< size_t, art::Ptr< recob::PFParticle > > &      idToPFParticleMap );

    
    void GetDownstreamPFParticles( const std::vector< art::Ptr< recob::PFParticle > > &                                      primaryPFParticles, 
                                   const std::map< size_t, art::Ptr< recob::PFParticle > > &                                 idToPFParticleMap, 
                                   std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PFParticle > > > & primaryToDownstreamPFParticles );

    void GetDownstreamPFParticles( const art::Ptr< recob::PFParticle > &                      part, 
                                   const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                   std::vector< art::Ptr< recob::PFParticle > > &             downstreamPFParticles );

    void GetPFParticleToClusterMap( const art::Handle< std::vector< recob::PFParticle > > &                                pfParticleHandle,
                                    const std::string &                                                                    producerLabel, 
                                    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Cluster > > > & pfParticleToClusterMap );
};

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_SLICES_H
