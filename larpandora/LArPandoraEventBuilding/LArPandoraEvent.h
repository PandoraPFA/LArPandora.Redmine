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
    void GetCollection( std::string                        inputLabel, 
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
    void GetAssociationMap( std::string                                                inputLabel, 
                            art::Handle< std::vector< T > > &                          inputHandleT, 
                            art::Handle< std::vector< U > > &                          inputHandleU, 
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
    void GetFilteredParticlesByPdgCode( bool                                                  shouldProduceNeutrinos, 
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
    void GetFilteredParticlesByCRTag(   bool                                                  shouldProduceNeutrinos, 
                                        std::string                                           tagProducerLabel,
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

public:

    /**
     *  @breif  Constructor from an art::Event
     */
    LArPandoraEvent( art::EDProducer * pProducer,
                     art::Event *      pEvent, 
                     std::string       inputProducerLabel,
                     std::string       hitProducerLabel );

    /**
     *  @breif  Copy constructor
     */
    LArPandoraEvent( const LArPandoraEvent & other );


    /**
     *  @breif  Produce a copy of the event keeping only the collections that are associated with a top-level particle whose Pdg code
     *          is a neutrino (non-neutrino) if shouldProduceNeutrinos is set to true (false)
     *
     *  @param  shouldProduceNeutrinos  if the returned event should contain neutrinos (or non-neutrinos)
     */
    LArPandoraEvent FilterByPdgCode( bool shouldProduceNeutrinos );

    /**
     *  @breif  Produce a copy of the event keeping only the collections that are associated with a top-level particle that is not 
     *          tagged as a neutrino (non-neutrino) if shouldProduceNeutrinos is set to true (false)
     *
     *  @param  shouldProduceNeutrinos  if the returned event should contain neutrinos (or non-neutrinos)
     *  @param  tagProducerLabel        label for the producer of the CRTags
     */
    LArPandoraEvent FilterByCRTag( bool         shouldProduceNeutrinos,
                                   std::string  tagProducerLabel );

    /**
     *  @breif  Write (put) the collections in this LArPandoraEvent to the art::Event
     */
    void WriteToEvent();

    /**
     *  @breif  Merge collections from two events into one
     */
    LArPandoraEvent Merge( LArPandoraEvent & other );
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArPandoraEvent::LArPandoraEvent( art::EDProducer * pProducer,
                                         art::Event *      pEvent,
                                         std::string       inputProducerLabel,
                                         std::string       hitProducerLabel ) :
    m_pProducer( pProducer ),
    m_pEvent( pEvent ), 
    m_inputProducerLabel( inputProducerLabel ),
    m_hitProducerLabel( hitProducerLabel )
{

    this->GetCollections();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArPandoraEvent::LArPandoraEvent( const LArPandoraEvent & other ) :
    m_pProducer( other.m_pProducer ),
    m_pEvent( other.m_pEvent),
    m_inputProducerLabel( other.m_inputProducerLabel ),
    m_hitProducerLabel( other.m_hitProducerLabel),
    m_pfParticles( other.m_pfParticles ),
    m_spacePoints( other.m_spacePoints ),
    m_clusters( other.m_clusters ),   
    m_vertices( other.m_vertices ),   
    m_tracks( other.m_tracks ),     
    m_showers( other.m_showers ),    
    m_seeds( other.m_seeds ),      
    m_pcAxes( other.m_pcAxes ),     
    m_hits( other.m_hits ),       
    m_pfParticleSpacePointMap( other.m_pfParticleSpacePointMap ),  
    m_pfParticleClusterMap( other.m_pfParticleClusterMap ),     
    m_pfParticleVertexMap( other.m_pfParticleVertexMap ),      
    m_pfParticleTrackMap( other.m_pfParticleTrackMap ),       
    m_pfParticleShowerMap( other.m_pfParticleShowerMap ),      
    m_pfParticleSeedMap( other.m_pfParticleSeedMap ),        
    m_pfParticlePCAxisMap( other.m_pfParticlePCAxisMap ),      
    m_spacePointHitMap( other.m_spacePointHitMap ),         
    m_clusterHitMap( other.m_clusterHitMap ),            
    m_trackHitMap( other.m_trackHitMap ),              
    m_showerHitMap( other.m_showerHitMap ),             
    m_seedHitMap( other.m_seedHitMap ),               
    m_showerPCAxisMap( other.m_showerPCAxisMap ) 
{}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArPandoraEvent::GetCollections()
{
    // Get each of the collections in turn
    art::Handle< std::vector< recob::PFParticle > > pfParticleHandle;
    this->GetCollection( m_inputProducerLabel, pfParticleHandle, m_pfParticles );

    art::Handle< std::vector< recob::SpacePoint > > spacePointHandle;
    this->GetCollection( m_inputProducerLabel, spacePointHandle, m_spacePoints ); 

    art::Handle< std::vector< recob::Cluster > > clusterHandle;
    this->GetCollection( m_inputProducerLabel, clusterHandle, m_clusters ); 

    art::Handle< std::vector< recob::Vertex > > vertexHandle;
    this->GetCollection( m_inputProducerLabel, vertexHandle, m_vertices ); 

    art::Handle< std::vector< recob::Track > > trackHandle;
    this->GetCollection( m_inputProducerLabel, trackHandle, m_tracks ); 

    art::Handle< std::vector< recob::Shower > > showerHandle;
    this->GetCollection( m_inputProducerLabel, showerHandle, m_showers ); 

    art::Handle< std::vector< recob::Seed > > seedHandle;
    this->GetCollection( m_inputProducerLabel, seedHandle, m_seeds ); 

    art::Handle< std::vector< recob::PCAxis > > pcAxisHandle;
    this->GetCollection( m_inputProducerLabel, pcAxisHandle, m_pcAxes ); 

    art::Handle< std::vector< recob::Hit > > hitHandle;
    this->GetCollection( m_hitProducerLabel, hitHandle, m_hits ); 

    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, spacePointHandle, m_pfParticleSpacePointMap );
    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, clusterHandle   , m_pfParticleClusterMap    );
    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, vertexHandle    , m_pfParticleVertexMap     );
    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, trackHandle     , m_pfParticleTrackMap      );
    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, showerHandle    , m_pfParticleShowerMap     );
    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, seedHandle      , m_pfParticleSeedMap       );
    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, pcAxisHandle    , m_pfParticlePCAxisMap     );

    this->GetAssociationMap( m_inputProducerLabel, spacePointHandle, hitHandle, m_spacePointHitMap );
    this->GetAssociationMap( m_inputProducerLabel, clusterHandle   , hitHandle, m_clusterHitMap    );
    this->GetAssociationMap( m_inputProducerLabel, trackHandle     , hitHandle, m_trackHitMap      );
    this->GetAssociationMap( m_inputProducerLabel, showerHandle    , hitHandle, m_showerHitMap     );
    this->GetAssociationMap( m_inputProducerLabel, seedHandle      , hitHandle, m_seedHitMap       );

    this->GetAssociationMap( m_inputProducerLabel, showerHandle, pcAxisHandle, m_showerPCAxisMap );
}

//------------------------------------------------------------------------------------------------------------------------------------------

template < class T >
inline void LArPandoraEvent::GetCollection( std::string                        inputLabel, 
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
inline void LArPandoraEvent::GetAssociationMap( std::string                                                inputLabel, 
                                                art::Handle< std::vector< T > > &                          inputHandleT, 
                                                art::Handle< std::vector< U > > &                          inputHandleU, 
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

inline LArPandoraEvent LArPandoraEvent::FilterByPdgCode( bool shouldProduceNeutrinos )
{
    // Get the primary top-level PFParticles that we want to keep
    std::vector< art::Ptr< recob::PFParticle > > primaryPFParticles;
    this->GetPrimaryPFParticles( primaryPFParticles );
  
    std::vector< art::Ptr< recob::PFParticle > > filteredPFParticles;
    this->GetFilteredParticlesByPdgCode( shouldProduceNeutrinos, primaryPFParticles, filteredPFParticles );

    // Collect all daughter PFParticles to produce the final list of selected PFParticles to persist
    std::map< size_t, art::Ptr< recob::PFParticle > > idToPFParticleMap;
    this->GetIdToPFParticleMap( idToPFParticleMap );

    std::map< size_t, art::Ptr< recob::PFParticle > > idToSelectedPFParticleMap;
    this->GetDownstreamPFParticles( filteredPFParticles, idToPFParticleMap, idToSelectedPFParticleMap);

    std::vector< art::Ptr< recob::PFParticle > > selectedPFParticles;
    this->GetParticleVector( idToSelectedPFParticleMap, selectedPFParticles );

    // Collect all related collections
    LArPandoraEvent filteredEvent( *this );
    this->FilterByParticleSelection( filteredEvent, selectedPFParticles );

    return filteredEvent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArPandoraEvent LArPandoraEvent::FilterByCRTag( bool         shouldProduceNeutrinos, 
                                                       std::string  tagProducerLabel )
{
    // Get the primary top-level PFParticles that we want to keep
    std::vector< art::Ptr< recob::PFParticle > > primaryPFParticles;
    this->GetPrimaryPFParticles( primaryPFParticles );

    std::vector< art::Ptr< recob::PFParticle > > filteredPFParticles;
    this->GetFilteredParticlesByCRTag( shouldProduceNeutrinos, tagProducerLabel, primaryPFParticles, filteredPFParticles );

    // Collect all daughter PFParticles to produce the final list of selected PFParticles to persist
    std::map< size_t, art::Ptr< recob::PFParticle > > idToPFParticleMap;
    this->GetIdToPFParticleMap( idToPFParticleMap );

    std::map< size_t, art::Ptr< recob::PFParticle > > idToSelectedPFParticleMap;
    this->GetDownstreamPFParticles( filteredPFParticles, idToPFParticleMap, idToSelectedPFParticleMap);

    std::vector< art::Ptr< recob::PFParticle > > selectedPFParticles;
    this->GetParticleVector( idToSelectedPFParticleMap, selectedPFParticles );

    // Collect all related collections
    LArPandoraEvent filteredEvent( *this );
    this->FilterByParticleSelection( filteredEvent, selectedPFParticles );

    return filteredEvent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArPandoraEvent::GetPrimaryPFParticles( std::vector< art::Ptr< recob::PFParticle > > & primaryPFParticles )
{
    for ( art::Ptr< recob::PFParticle > part : m_pfParticles ) 
        if ( part->IsPrimary() )
            primaryPFParticles.push_back( part );
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArPandoraEvent::GetFilteredParticlesByPdgCode( bool                                                  shouldProduceNeutrinos, 
                                                            const std::vector< art::Ptr< recob::PFParticle > > &  inputPFParticles, 
                                                            std::vector< art::Ptr< recob::PFParticle > > &        outputPFParticles )
{
    for ( art::Ptr< recob::PFParticle > part : inputPFParticles ) {
        unsigned int pdg = std::abs( part->PdgCode() );
        bool isNeutrino = ( pdg == nue || pdg == numu || pdg == nutau );

        if ( ( shouldProduceNeutrinos && isNeutrino ) || ( !shouldProduceNeutrinos && !isNeutrino ) ) {
            outputPFParticles.push_back( part );
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArPandoraEvent::GetFilteredParticlesByCRTag(   bool                                                  shouldProduceNeutrinos, 
                                                            std::string                                           tagProducerLabel,
                                                            const std::vector< art::Ptr< recob::PFParticle > > &  inputPFParticles, 
                                                            std::vector< art::Ptr< recob::PFParticle > > &        outputPFParticles )
{

    art::Handle< std::vector< recob::PFParticle > > pfParticleHandle;
    m_pEvent->getByLabel(m_inputProducerLabel, pfParticleHandle);

    art::FindManyP< anab::CosmicTag > pfParticleTagAssoc( pfParticleHandle, *m_pEvent, tagProducerLabel );
    
    for ( art::Ptr< recob::PFParticle > part : inputPFParticles ) {
        const std::vector< art::Ptr< anab::CosmicTag > > cosmicTags = pfParticleTagAssoc.at( part.key() );

        /// @todo error gracefully
        if (cosmicTags.size() != 1) exit(1);

        art::Ptr< anab::CosmicTag > cosmicTag = cosmicTags.front();
        bool isNeutrino = ( cosmicTag->CosmicType() == anab::kNotTagged );

        if ( ( shouldProduceNeutrinos && isNeutrino ) || ( !shouldProduceNeutrinos && !isNeutrino ) ) {
            outputPFParticles.push_back( part );
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArPandoraEvent::GetIdToPFParticleMap( std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap ) 
{
    for ( art::Ptr< recob::PFParticle > part : m_pfParticles ) 
      idToPFParticleMap.insert( std::map< size_t, art::Ptr< recob::PFParticle > >::value_type( part->Self(), part ) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArPandoraEvent::GetDownstreamPFParticles( const std::vector< art::Ptr< recob::PFParticle > > &       inputPFParticles, 
                                                       const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                                       std::map< size_t, art::Ptr< recob::PFParticle > > &        idToDownstreamPFParticleMap )
{
  for ( art::Ptr< recob::PFParticle > part : inputPFParticles )
    this->GetDownstreamPFParticles( part, idToPFParticleMap, idToDownstreamPFParticleMap );
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArPandoraEvent::GetDownstreamPFParticles( art::Ptr< recob::PFParticle >             part, 
                                      const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                      std::map< size_t, art::Ptr< recob::PFParticle > > &        idToDownstreamPFParticleMap )
{
    if ( idToDownstreamPFParticleMap.find( part->Self() ) == idToDownstreamPFParticleMap.end() )
        idToDownstreamPFParticleMap.insert( std::map< size_t, art::Ptr< recob::PFParticle > >::value_type( part->Self(), part ) );

    for ( size_t daughterId : part->Daughters() ) {
        std::map< size_t, art::Ptr< recob::PFParticle > >::const_iterator daughterIt = idToPFParticleMap.find( daughterId );

        /// @todo handle error gracefully!
        if ( daughterIt == idToPFParticleMap.end() ) std::exit(1);

        this->GetDownstreamPFParticles( daughterIt->second, idToPFParticleMap, idToDownstreamPFParticleMap );
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArPandoraEvent::GetParticleVector( const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                                std::vector< art::Ptr< recob::PFParticle > > &             pfParticleVector )
{
    for ( std::map< size_t, art::Ptr< recob::PFParticle > >::const_iterator it=idToPFParticleMap.begin(); it != idToPFParticleMap.end(); ++it ) {
        art::Ptr< recob::PFParticle > part = it->second;
        pfParticleVector.push_back( part );
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArPandoraEvent::FilterByParticleSelection( LArPandoraEvent &                               filteredEvent, 
                                                        std::vector< art::Ptr< recob::PFParticle > > &  selectedPFParticles )
{
    filteredEvent.m_pfParticles = selectedPFParticles;

    std::vector< art::Ptr< recob::PFParticle> > selectedParticles;
    std::vector< art::Ptr< recob::SpacePoint> > selectedSpacePoints;
    std::vector< art::Ptr< recob::Cluster>    > selectedClusters;
    std::vector< art::Ptr< recob::Seed>       > selectedSeeds;
    std::vector< art::Ptr< recob::Vertex>     > selectedVertices;
    std::vector< art::Ptr< recob::Track>      > selectedTracks;
    std::vector< art::Ptr< recob::Shower>     > selectedShowers; 
    std::vector< art::Ptr< recob::PCAxis>     > selectedPCAxes;

    for ( art::Ptr< recob::PFParticle > part : selectedPFParticles ) {
        this->CollectAssociated( part, m_pfParticleSpacePointMap, selectedSpacePoints );
        this->CollectAssociated( part, m_pfParticleClusterMap   , selectedClusters    );
        this->CollectAssociated( part, m_pfParticleSeedMap      , selectedSeeds       );
        this->CollectAssociated( part, m_pfParticleVertexMap    , selectedVertices    );
        this->CollectAssociated( part, m_pfParticleTrackMap     , selectedTracks      );
        this->CollectAssociated( part, m_pfParticleShowerMap    , selectedShowers     );
        this->CollectAssociated( part, m_pfParticlePCAxisMap    , selectedPCAxes      );
    }

    filteredEvent.m_spacePoints = selectedSpacePoints;
    filteredEvent.m_clusters    = selectedClusters;
    filteredEvent.m_seeds       = selectedSeeds;
    filteredEvent.m_vertices    = selectedVertices;
    filteredEvent.m_tracks      = selectedTracks;
    filteredEvent.m_showers     = selectedShowers;
    filteredEvent.m_pcAxes      = selectedPCAxes;

    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::SpacePoint > > >    selectedPFParticleSpacePointMap;
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Cluster > > >       selectedPFParticleClusterMap;  
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Vertex > > >        selectedPFParticleVertexMap;  
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Track > > >         selectedPFParticleTrackMap;  
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Shower > > >        selectedPFParticleShowerMap;
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Seed > > >          selectedPFParticleSeedMap; 
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PCAxis > > >        selectedPFParticlePCAxisMap;
    std::map< art::Ptr< recob::SpacePoint >, std::vector< art::Ptr< recob::Hit > > >           selectedSpacePointHitMap;  
    std::map< art::Ptr< recob::Cluster >   , std::vector< art::Ptr< recob::Hit > > >           selectedClusterHitMap;    
    std::map< art::Ptr< recob::Track >     , std::vector< art::Ptr< recob::Hit > > >           selectedTrackHitMap;     
    std::map< art::Ptr< recob::Shower >    , std::vector< art::Ptr< recob::Hit > > >           selectedShowerHitMap;
    std::map< art::Ptr< recob::Seed >      , std::vector< art::Ptr< recob::Hit > > >           selectedSeedHitMap; 
    std::map< art::Ptr< recob::Shower >    , std::vector< art::Ptr< recob::PCAxis > > >        selectedShowerPCAxisMap;

    this->GetFilteredAssociationMap( selectedPFParticles, selectedSpacePoints, m_pfParticleSpacePointMap, selectedPFParticleSpacePointMap );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedClusters   , m_pfParticleClusterMap   , selectedPFParticleClusterMap    );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedVertices   , m_pfParticleVertexMap    , selectedPFParticleVertexMap     );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedTracks     , m_pfParticleTrackMap     , selectedPFParticleTrackMap      );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedShowers    , m_pfParticleShowerMap    , selectedPFParticleShowerMap     );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedSeeds      , m_pfParticleSeedMap      , selectedPFParticleSeedMap       );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedPCAxes     , m_pfParticlePCAxisMap    , selectedPFParticlePCAxisMap     );
    this->GetFilteredAssociationMap( selectedSpacePoints, m_hits, m_spacePointHitMap, selectedSpacePointHitMap ); 
    this->GetFilteredAssociationMap( selectedClusters   , m_hits, m_clusterHitMap   , selectedClusterHitMap    );
    this->GetFilteredAssociationMap( selectedTracks     , m_hits, m_trackHitMap     , selectedTrackHitMap      );
    this->GetFilteredAssociationMap( selectedShowers    , m_hits, m_showerHitMap    , selectedShowerHitMap     );
    this->GetFilteredAssociationMap( selectedSeeds      , m_hits, m_seedHitMap      , selectedSeedHitMap       );
    this->GetFilteredAssociationMap( selectedShowers, selectedPCAxes, m_showerPCAxisMap, selectedShowerPCAxisMap );

    filteredEvent.m_pfParticleSpacePointMap = selectedPFParticleSpacePointMap;
    filteredEvent.m_pfParticleClusterMap    = selectedPFParticleClusterMap;
    filteredEvent.m_pfParticleVertexMap     = selectedPFParticleVertexMap;
    filteredEvent.m_pfParticleTrackMap      = selectedPFParticleTrackMap;
    filteredEvent.m_pfParticleShowerMap     = selectedPFParticleShowerMap;
    filteredEvent.m_pfParticleSeedMap       = selectedPFParticleSeedMap;
    filteredEvent.m_pfParticlePCAxisMap     = selectedPFParticlePCAxisMap;
    filteredEvent.m_spacePointHitMap        = selectedSpacePointHitMap;
    filteredEvent.m_clusterHitMap           = selectedClusterHitMap;
    filteredEvent.m_trackHitMap             = selectedTrackHitMap;
    filteredEvent.m_showerHitMap            = selectedShowerHitMap;
    filteredEvent.m_seedHitMap              = selectedSeedHitMap;
    filteredEvent.m_showerPCAxisMap         = selectedShowerPCAxisMap;
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

inline void LArPandoraEvent::WriteToEvent()
{
    
    this->WriteCollection( m_pfParticles );
    this->WriteCollection( m_spacePoints );
    this->WriteCollection( m_clusters );   
    this->WriteCollection( m_vertices );   
    this->WriteCollection( m_tracks );     
    this->WriteCollection( m_showers );    
    this->WriteCollection( m_seeds );      
    this->WriteCollection( m_pcAxes );     
    
    this->WriteAssociation( m_pfParticleSpacePointMap );
    this->WriteAssociation( m_pfParticleClusterMap );   
    this->WriteAssociation( m_pfParticleVertexMap );    
    this->WriteAssociation( m_pfParticleTrackMap );     
    this->WriteAssociation( m_pfParticleShowerMap );    
    this->WriteAssociation( m_pfParticleSeedMap );      
    this->WriteAssociation( m_pfParticlePCAxisMap );    
    this->WriteAssociation( m_spacePointHitMap );       
    this->WriteAssociation( m_clusterHitMap );          
    this->WriteAssociation( m_trackHitMap );            
    this->WriteAssociation( m_showerHitMap );           
    this->WriteAssociation( m_seedHitMap );             
    this->WriteAssociation( m_showerPCAxisMap );        
      
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

// ---------------------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------------------

inline LArPandoraEvent LArPandoraEvent::Merge( LArPandoraEvent & other )
{
    LArPandoraEvent outputEvent( other );

    this->MergeCollection( other.m_pfParticles, m_pfParticles );
    this->MergeCollection( other.m_spacePoints, m_spacePoints );
    this->MergeCollection( other.m_clusters   , m_clusters );   
    this->MergeCollection( other.m_vertices   , m_vertices );   
    this->MergeCollection( other.m_tracks     , m_tracks );     
    this->MergeCollection( other.m_showers    , m_showers );    
    this->MergeCollection( other.m_seeds      , m_seeds );      
    this->MergeCollection( other.m_pcAxes     , m_pcAxes );     
    
    this->MergeAssociation( other.m_pfParticleSpacePointMap, m_pfParticleSpacePointMap );
    this->MergeAssociation( other.m_pfParticleClusterMap   , m_pfParticleClusterMap );   
    this->MergeAssociation( other.m_pfParticleVertexMap    , m_pfParticleVertexMap );    
    this->MergeAssociation( other.m_pfParticleTrackMap     , m_pfParticleTrackMap );     
    this->MergeAssociation( other.m_pfParticleShowerMap    , m_pfParticleShowerMap );    
    this->MergeAssociation( other.m_pfParticleSeedMap      , m_pfParticleSeedMap );      
    this->MergeAssociation( other.m_pfParticlePCAxisMap    , m_pfParticlePCAxisMap );    
    this->MergeAssociation( other.m_spacePointHitMap       , m_spacePointHitMap );       
    this->MergeAssociation( other.m_clusterHitMap          , m_clusterHitMap );          
    this->MergeAssociation( other.m_trackHitMap            , m_trackHitMap );            
    this->MergeAssociation( other.m_showerHitMap           , m_showerHitMap );           
    this->MergeAssociation( other.m_seedHitMap             , m_seedHitMap );             
    this->MergeAssociation( other.m_showerPCAxisMap        , m_showerPCAxisMap );        

    return outputEvent;
}

// ---------------------------------------------------------------------------------------

template < class T >
inline void LArPandoraEvent::MergeCollection( std::vector< art::Ptr< T > > &  collectionToMerge, 
                                              std::vector< art::Ptr< T > > &  collection )
{
    collectionToMerge.insert( collectionToMerge.end(), collection.begin(), collection.end() );
}

// ---------------------------------------------------------------------------------------

template < class T, class U >
inline void MergeAssociation( std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > &  associationToMerge, 
                              std::map< art::Ptr< T >, std::vector< art::Ptr< U > > > &  association )
{
    associationToMerge.insert( association.begin(), association.end() );
}

// ---------------------------------------------------------------------------------------

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_EVENT_H
