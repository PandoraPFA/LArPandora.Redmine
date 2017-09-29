/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraEvent.h
 *
 *  @brief  A description of all outputs from an instance of pandora with functionality to filter and merge multiple output
 */

#ifndef LAR_PANDORA_EVENT_H
#define LAR_PANDORA_EVENT_H 1

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
 *  @breif LArPandoraEvent class
 */
class LArPandoraEvent
{
public:

    /**
     *  @breif  Constructor from an art::Event
     */
    LArPandoraEvent( art::EDProducer * pProducer,
                     art::Event *      pEvent, 
                     std::string       inputProducerLabel,
                     std::string       hitProducerLabel );

    /**
     *  @breif  Produce a copy of the event keeping only the collections that are associated with a top-level particle whose Pdg code
     *          is a neutrino (non-neutrino) if shouldProduceNeutrinos is set to true (false)
     *
     *  @param  shouldProduceNeutrinos  if the returned event should contain neutrinos (or non-neutrinos)
     */
    LArPandoraEvent FilterByPdgCode( const bool shouldProduceNeutrinos );

    /**
     *  @breif  Produce a copy of the event keeping only the collections that are associated with a top-level particle that is not 
     *          tagged as a neutrino (non-neutrino) if shouldProduceNeutrinos is set to true (false)
     *
     *  @param  shouldProduceNeutrinos  if the returned event should contain neutrinos (or non-neutrinos)
     *  @param  tagProducerLabel        label for the producer of the CRTags
     */
    LArPandoraEvent FilterByCRTag( const bool          shouldProduceNeutrinos,
                                   const std::string & tagProducerLabel );

    /**
     *  @breif  Write (put) the collections in this LArPandoraEvent to the art::Event
     */
    void WriteToEvent();

    /**
     *  @breif  Merge collections from two events into one
     */
    LArPandoraEvent Merge( LArPandoraEvent & other );


    /**
     *  @breif  Class to handle the required producer labels 
     */
    class Labels
    {
    public:

        /**
         *  @breif  Minimal parametrised constructor.
         *          Sets all collection labels to be the same as the PFParticle producer label 
         */
        Labels( std::string pfParticleProducerLabel,
                std::string hitProducerLabel );

        // Setter functions

        void SetSpacePointProducerLabel( const std::string & label );
        void SetClusterProducerLabel( const std::string & label );
        void SetVertexProducerLabel( const std::string & label );
        void SetTrackProducerLabel( const std::string & label );
        void SetShowerProducerLabel( const std::string & label );
        void SetSeedProducerLabel( const std::string & label );
        void SetPCAxisProducerLabel( const std::string & label );

    private:

        std::string  m_pfParticleLabel;
        std::string  m_spacePointLabel;
        std::string  m_clusterLabel;
        std::string  m_vertexLabel;
        std::string  m_trackLabel;
        std::string  m_showerLabel;
        std::string  m_seedLabel;
        std::string  m_pcAxisLabel;
        std::string  m_hitLabel;
    };

private:
    
    // Meta data
    art::EDProducer * m_pProducer;            ///<
    art::Event *      m_pEvent;               ///<
    std::string       m_inputProducerLabel;   ///<
    std::string       m_hitProducerLabel;     ///<

    // Collections
    std::vector< art::Ptr< recob::PFParticle > > m_pfParticles;    ///<
    std::vector< art::Ptr< recob::SpacePoint > > m_spacePoints;    ///<
    std::vector< art::Ptr< recob::Cluster > >    m_clusters;       ///<
    std::vector< art::Ptr< recob::Vertex > >     m_vertices;       ///<
    std::vector< art::Ptr< recob::Track > >      m_tracks;         ///<
    std::vector< art::Ptr< recob::Shower > >     m_showers;        ///<
    std::vector< art::Ptr< recob::Seed > >       m_seeds;          ///<
    std::vector< art::Ptr< recob::PCAxis > >     m_pcAxes;         ///<
    std::vector< art::Ptr< recob::Hit> >         m_hits;           ///<

    // Association maps
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::SpacePoint > > >    m_pfParticleSpacePointMap;    ///<
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Cluster > > >       m_pfParticleClusterMap;       ///<
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Vertex > > >        m_pfParticleVertexMap;        ///<
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Track > > >         m_pfParticleTrackMap;         ///<
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Shower > > >        m_pfParticleShowerMap;        ///<
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Seed > > >          m_pfParticleSeedMap;          ///<
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PCAxis > > >        m_pfParticlePCAxisMap;        ///<

    std::map< art::Ptr< recob::SpacePoint >, std::vector< art::Ptr< recob::Hit > > >           m_spacePointHitMap;           ///<
    std::map< art::Ptr< recob::Cluster >   , std::vector< art::Ptr< recob::Hit > > >           m_clusterHitMap;              ///<
    std::map< art::Ptr< recob::Track >     , std::vector< art::Ptr< recob::Hit > > >           m_trackHitMap;                ///<
    std::map< art::Ptr< recob::Shower >    , std::vector< art::Ptr< recob::Hit > > >           m_showerHitMap;               ///<
    std::map< art::Ptr< recob::Seed >      , std::vector< art::Ptr< recob::Hit > > >           m_seedHitMap;                 ///<

    std::map< art::Ptr< recob::Shower >    , std::vector< art::Ptr< recob::PCAxis > > >        m_showerPCAxisMap;            ///<
    
    // -------------------------------------------------------------------------------------------------------------------------------------
    
    /**
     *  @breif  Get the collections and associations from m_pEvent with the required labels
     */
    void GetCollections();

    /**
     *  @breif  Gets a given collection from m_pEvent with the label supplied
     *
     *  @param  inputLabel        a label for the producer of the collection required
     *  @param  outputHandle      the output art Handle to required collection
     *  @param  outputCollection  the required collection
     */
    template < class T >
    void GetCollection( const std::string &                inputLabel, 
                        art::Handle< std::vector< T > > &  outputHandle, 
                        std::vector< art::Ptr< T > > &     outputCollection );

    /**
     *  @breif  Get the mapping between two collections using the specified label
     *
     *  @param  inputLabel            a label for the producer of the association required
     *  @param  inputHandleT          the input art Handle to the first collection
     *  @param  inputHandleU          the input art Handle to the second collection
     *  @param  outputAssociationMap  output mapping between the two data types supplied (T -> U)
     */
    template < class T, class U >
    void GetAssociationMap( const std::string &                                        inputLabel, 
                            art::Handle< std::vector< T > > &                          inputHandleT, 
                            std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > &  outputAssociationMap );

    /**
     *  @brief  Filters primary PFParticles from the m_pfParticles
     *
     *  @param  primaryPFParticles  output vector of all primary PFParticles in the input vector
     */
    void GetPrimaryPFParticles( std::vector< art::Ptr< recob::PFParticle > > &  primaryPFParticles );

    /**
     *  @brief  Filters PFParticles based on their Pdg from the inputPFParticles
     *
     *  @param  shouldProduceNeutrinos  if the filtered particle vector should contain neutrinos (or non-neutrinos)
     *  @param  inputPFParticles        input vector of PFParticles
     *  @param  filteredPFParticles     output vector of filtered PFParticles
     */
    void GetFilteredParticlesByPdgCode( const bool                                            shouldProduceNeutrinos, 
                                        const std::vector< art::Ptr< recob::PFParticle > > &  inputPFParticles, 
                                        std::vector< art::Ptr< recob::PFParticle > > &        outputPFParticles );

    /**
     *  @brief  Filters PFParticles based on their Pdg from the inputPFParticles
     *
     *  @param  shouldProduceNeutrinos  if the filtered particle vector should contain neutrinos (or non-neutrinos)
     *  @param  tagProducerLabel        the label for the producer of the CR tags
     *  @param  inputPFParticles        input vector of PFParticles
     *  @param  filteredPFParticles     output vector of filtered PFParticles
     */
    void GetFilteredParticlesByCRTag(   const bool                                            shouldProduceNeutrinos, 
                                        const std::string &                                   tagProducerLabel,
                                        const std::vector< art::Ptr< recob::PFParticle > > &  inputPFParticles, 
                                        std::vector< art::Ptr< recob::PFParticle > > &        outputPFParticles );

    /**
     *  @brief  Produce a mapping between PFParticles and their ID
     *
     *  @param  idToPFParticleMap   output mapping between PFParticles and their IDs
     */
    void GetIdToPFParticleMap( std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap );

    /**
     *  \brief  Collects all PFParticles downstream (children, grandchildren, ...) of a given PFParticle
     *
     *  \param  part                         input PFParticle
     *  \param  idToPFParticleMap            input mapping between PFParticles and their IDs
     *  \param  idToDownstreamPFParticleMap  output mapping between of all downstream PFParticles and their IDs
     */
    void GetDownstreamPFParticles( art::Ptr< recob::PFParticle >                              part, 
                                   const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                   std::map< size_t, art::Ptr< recob::PFParticle > > &        idToDownstreamPFParticleMap );

    /**
     *  @brief  Collects all PFParticles downstream (children, grandchildren, ...) of a given vector of PFParticle
     *
     *  @param  inputPFParticles             input vector of PFParticles
     *  @param  idToPFParticleMap            input mapping between PFParticles and their IDs
     *  @param  idToDownstreamPFParticleMap  output mapping between of all downstream PFParticles and their IDs
     */
    void GetDownstreamPFParticles( const std::vector< art::Ptr< recob::PFParticle > > &       inputPFParticles, 
                                   const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                   std::map< size_t, art::Ptr< recob::PFParticle > > &        idToDownstreamPFParticleMap );
    
    /**
     *  @brief  Converts a map from id to particles to a vector of particles
     * 
     *  @param  idToPFParticleMap  input map of PFParticles
     *  @param  pfParticleVector   output vector of PFParticles
     */
    void GetParticleVector( const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                            std::vector< art::Ptr< recob::PFParticle > > &             pfParticleVector );
    
    /**
     *  @brief  Replaces the collections and associations in the supplied event by any objects in this's collections that are associated 
     *          with a PFParticle in the selection supplied.
     * 
     *  @param  filteredEvent      input event with which to give the filtered collections
     *  @param  pfParticleVector   input vector of selected particles 
     */
    void FilterByParticleSelection( LArPandoraEvent &                               filteredEvent, 
                                    std::vector< art::Ptr< recob::PFParticle > > &  selectedPFParticles );


    /**
     *  @brief  Collects all objects of type U associated to a given object of type T
     *
     *  @param  anObject         an input object of type T with which we want to collect associated objects of type U
     *  @param  associationTtoU  the general input association between objects of type U and T
     *  @param  associatedU      output vector of objects of type U associated with anObject
     */
    template < class T, class U >
    void CollectAssociated( const art::Ptr< T > &                                            anObject, 
                            std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > &        associationTtoU, 
                            std::vector< art::Ptr< U > > &                                   associatedU );

    /**
     *   @brief  Gets the mapping between two filtered collections 
     *
     *   @param  collectionT            a first filtered collection
     *   @param  collectionU            a second filtered collection
     *   @param  inputAssociationTtoU   mapping between the two unfiltered collections
     *   @param  outputAssociationTtoU  mapping between the two filtered collections
     *
     *   @return mapping between the filtered collections
     */
    template < class T, class U >
    void GetFilteredAssociationMap( const std::vector< art::Ptr< T > > &                      collectionT, 
                                    const std::vector< art::Ptr< U > > &                      collectionU, 
                                    std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > & inputAssociationTtoU,
                                    std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > & outputAssociationTtoU );

    /**
     *  @brief  Write a given collection to the event
     */
    template < class T >
    void WriteCollection( const std::vector< art::Ptr< T > > & collection );

    /**
     *  @brief  Write a given association to the event
     */
    template < class T, class U >
    void WriteAssociation( const std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > & associationMap );

    /**
     *  @brief  Append a collection onto an other collection
     *
     *  @param  collectionToMerge  the collection to accept
     *  @param  collection         the collection to append
     */
    template < class T >
    void MergeCollection( std::vector< art::Ptr< T > > &  collectionToMerge, 
                          std::vector< art::Ptr< T > > &  collection );

    /**
     *  @brief  Append an association to another association
     *
     *  @param  associationToMerge  the association to accept
     *  @param  association         the association to append
     */
    template < class T, class U >
    void MergeAssociation( std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > &  associationToMerge, 
                           std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > &  association );

    // Useful PDG codes for readability
    enum Pdg {
        nue   = 12,
        numu  = 14,
        nutau = 16
    };
};

//------------------------------------------------------------------------------------------------------------------------------------------

template < class T >
inline void LArPandoraEvent::GetCollection( const std::string &                inputLabel, 
                                            art::Handle< std::vector< T > > &  outputHandle, 
                                            std::vector< art::Ptr< T > > &     outputCollection )
{
    m_pEvent->getByLabel( inputLabel, outputHandle);   

    for( unsigned int i = 0; i != outputHandle->size(); i++ ) {
        art::Ptr< T > object( outputHandle, i );
        outputCollection.push_back( object );
    } 
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
template < class T, class U >
inline void LArPandoraEvent::GetAssociationMap( const std::string &                                        inputLabel, 
                                                art::Handle< std::vector< T > > &                          inputHandleT, 
                                                std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > &  outputAssociationMap )
{
    art::FindManyP< U > assoc( inputHandleT, (*m_pEvent), inputLabel );

    for ( unsigned int iT = 0; iT < inputHandleT->size(); iT++ ) {

        art::Ptr< T > objectT( inputHandleT, iT );

        if ( outputAssociationMap.find( objectT ) == outputAssociationMap.end() ) {
            std::vector< art::Ptr< U > > emptyVect;
            outputAssociationMap.insert( typename std::map< art::Ptr< T >, std::vector< art::Ptr< U > > >::value_type( objectT, emptyVect ) );
        }

        for ( art::Ptr< U > objectU : assoc.at( objectT.key() ) ) { 
            outputAssociationMap[ objectT ].push_back( objectU );        
        }
    } 
}

//------------------------------------------------------------------------------------------------------------------------------------------

template < class T, class U >
inline void LArPandoraEvent::CollectAssociated( const art::Ptr< T > &                                            anObject, 
                                                std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > &        associationTtoU, 
                                                std::vector< art::Ptr< U > > &                                   associatedU )
{
    std::vector< art::Ptr< U > > associatedObjects = associationTtoU[ anObject ];
    associatedU.insert( associatedU.end(), associatedObjects.begin(), associatedObjects.end() );
}

//------------------------------------------------------------------------------------------------------------------------------------------

template < class T, class U >
inline void LArPandoraEvent::GetFilteredAssociationMap( const std::vector< art::Ptr< T > > &                      collectionT, 
                                       const std::vector< art::Ptr< U > > &                      collectionU, 
                                       std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > & inputAssociationTtoU,
                                       std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > & outputAssociationTtoU )
{

    for ( art::Ptr< T > objectT : collectionT ) {
        
        std::vector< art::Ptr< U > > emptyVector;
        outputAssociationTtoU.insert( typename std::map< art::Ptr< T >, std::vector< art::Ptr< U > > >::value_type( objectT, emptyVector ) );

        for ( art::Ptr< U > objectU : inputAssociationTtoU[ objectT ] ) {

            // Check that the objectU is in collectionU
            typename std::vector< art::Ptr< U > >::const_iterator associatedObjectIter = std::find( collectionU.begin(), collectionU.end(), objectU );
            if ( associatedObjectIter == collectionU.end() ) continue;
       
            outputAssociationTtoU[ objectT ].push_back( objectU );
        }
    }
}
        
//------------------------------------------------------------------------------------------------------------------------------------------

template < class T >
inline void LArPandoraEvent::WriteCollection( const std::vector< art::Ptr< T > > & collection )
{
  std::unique_ptr< std::vector< T > > output( new std::vector< T > );

  for ( art::Ptr< T > object : collection )
    output->push_back( *object );
  
  m_pEvent->put( std::move( output ) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

template < class T, class U >
inline void LArPandoraEvent::WriteAssociation( const std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > & associationMap )
{
  
  std::unique_ptr< art::Assns< T, U > > outputAssn( new art::Assns< T, U > );

  for ( typename std::map< art::Ptr< T >, std::vector< art::Ptr< U > > >::const_iterator it=associationMap.begin(); it != associationMap.end(); ++it ) {
    art::Ptr< T > objectT = it->first;
    for ( art::Ptr< U > objectU : it->second ) {
        util::CreateAssn( *m_pProducer, *m_pEvent, objectU , objectT, *outputAssn );
    }
  }

  m_pEvent->put( std::move( outputAssn ) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

template < class T >
inline void LArPandoraEvent::MergeCollection( std::vector< art::Ptr< T > > &  collectionToMerge, 
                                              std::vector< art::Ptr< T > > &  collection )
{
    collectionToMerge.insert( collectionToMerge.end(), collection.begin(), collection.end() );
}

//------------------------------------------------------------------------------------------------------------------------------------------

template < class T, class U >
inline void LArPandoraEvent::MergeAssociation( std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > &  associationToMerge, 
                              std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > &  association )
{
    associationToMerge.insert( association.begin(), association.end() );
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_EVENT_H
