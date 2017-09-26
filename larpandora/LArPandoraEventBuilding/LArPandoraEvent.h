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

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/Hit.h"

#include <memory>
#include <algorithm>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @breif LArPandoraEvent class
 */
class LArPandoraEvent {

private:
    
    // Meta data
    art::Event * m_event;                ///<
    std::string  m_inputProducerLabel;   ///<
    std::string  m_hitProducerLabel;     ///<

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
     *  @breif  Get the collections and associations from m_event with the required labels
     */
    void GetCollections();

    /**
     *  @breif  Gets a given collection from m_event with the label supplied
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

public:

    /**
     *  @breif  Constructor from an art::Event
     */
    LArPandoraEvent( art::Event *  event, 
                     std::string   inputProducerLabel,
                     std::string   hitProducerLabel );

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArPandoraEvent::LArPandoraEvent( art::Event *  event,
                                         std::string   inputProducerLabel,
                                         std::string   hitProducerLabel ) :
    m_event( event ), 
    m_inputProducerLabel( inputProducerLabel ),
    m_hitProducerLabel( hitProducerLabel )
{

    this->GetCollections();
}

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

    // Get the association maps
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
    m_event->getByLabel( inputLabel, outputHandle);   

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
    art::FindManyP< U > assoc( inputHandleT, (*m_event), inputLabel );

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

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_EVENT_H
