/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraEvent.cxx
 *
 *  @brief  A description of all outputs from an instance of pandora with functionality to filter and merge multiple output
 */

#include "larpandora/LArPandoraEventBuilding/LArPandoraEvent.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

LArPandoraEvent::LArPandoraEvent( art::EDProducer *  pProducer,
                                  art::Event *       pEvent,
                                  const Labels &     inputLabels ) :
    m_pProducer( pProducer ),
    m_pEvent( pEvent ), 
    m_labels( inputLabels )
{

    this->GetCollections();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::GetCollections()
{
    // Get each of the collections in turn
    art::Handle< std::vector< recob::PFParticle > > pfParticleHandle;
    this->GetCollection( Labels::PFParticleLabel, pfParticleHandle, m_pfParticles );

    art::Handle< std::vector< recob::SpacePoint > > spacePointHandle;
    this->GetCollection( Labels::SpacePointLabel, spacePointHandle, m_spacePoints ); 

    art::Handle< std::vector< recob::Cluster > > clusterHandle;
    this->GetCollection( Labels::ClusterLabel, clusterHandle, m_clusters ); 

    art::Handle< std::vector< recob::Vertex > > vertexHandle;
    this->GetCollection( Labels::VertexLabel, vertexHandle, m_vertices ); 

    art::Handle< std::vector< recob::Track > > trackHandle;
    this->GetCollection( Labels::TrackLabel, trackHandle, m_tracks ); 

    art::Handle< std::vector< recob::Shower > > showerHandle;
    this->GetCollection( Labels::ShowerLabel, showerHandle, m_showers ); 

    art::Handle< std::vector< recob::Seed > > seedHandle;
    this->GetCollection( Labels::SeedLabel, seedHandle, m_seeds ); 

    art::Handle< std::vector< recob::PCAxis > > pcAxisHandle;
    this->GetCollection( Labels::PCAxisLabel, pcAxisHandle, m_pcAxes ); 

    art::Handle< std::vector< recob::Hit > > hitHandle;
    this->GetCollection( Labels::HitLabel, hitHandle, m_hits ); 

    this->GetAssociationMap( Labels::PFParticleToSpacePointLabel, pfParticleHandle, m_pfParticleSpacePointMap );
    this->GetAssociationMap( Labels::PFParticleToClusterLabel   , pfParticleHandle, m_pfParticleClusterMap    );
    this->GetAssociationMap( Labels::PFParticleToVertexLabel    , pfParticleHandle, m_pfParticleVertexMap     );
    this->GetAssociationMap( Labels::PFParticleToTrackLabel     , pfParticleHandle, m_pfParticleTrackMap      );
    this->GetAssociationMap( Labels::PFParticleToShowerLabel    , pfParticleHandle, m_pfParticleShowerMap     );
    this->GetAssociationMap( Labels::PFParticleToSeedLabel      , pfParticleHandle, m_pfParticleSeedMap       );
    this->GetAssociationMap( Labels::PFParticleToPCAxisLabel    , pfParticleHandle, m_pfParticlePCAxisMap     );

    this->GetAssociationMap( Labels::SpacePointToHitLabel, spacePointHandle, m_spacePointHitMap );
    this->GetAssociationMap( Labels::ClusterToHitLabel   , clusterHandle   , m_clusterHitMap    );
    this->GetAssociationMap( Labels::TrackToHitLabel     , trackHandle     , m_trackHitMap      );
    this->GetAssociationMap( Labels::ShowerToHitLabel    , showerHandle    , m_showerHitMap     );
    this->GetAssociationMap( Labels::SeedToHitLabel      , seedHandle      , m_seedHitMap       );

    this->GetAssociationMap( Labels::ShowerToPCAxisLabel, showerHandle, m_showerPCAxisMap );
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraEvent LArPandoraEvent::FilterByPdgCode( const bool shouldProduceNeutrinos )
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
    LArPandoraEvent filteredEvent( *this, selectedPFParticles );

    return filteredEvent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraEvent LArPandoraEvent::FilterByCRTag( const bool          shouldProduceNeutrinos, 
                                                const std::string & tagProducerLabel )
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
    LArPandoraEvent filteredEvent( *this, selectedPFParticles );

    return filteredEvent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraEvent::LArPandoraEvent ( const LArPandoraEvent &                         event, 
                                   std::vector< art::Ptr< recob::PFParticle > > &  selectedPFParticles ) : 
    m_pProducer( event.m_pProducer ),
    m_pEvent( event.m_pEvent ), 
    m_labels( event.m_labels ),
    m_hits( event.m_hits )
{
    m_pfParticles = selectedPFParticles;

    std::vector< art::Ptr< recob::PFParticle> > selectedParticles;
    std::vector< art::Ptr< recob::SpacePoint> > selectedSpacePoints;
    std::vector< art::Ptr< recob::Cluster>    > selectedClusters;
    std::vector< art::Ptr< recob::Seed>       > selectedSeeds;
    std::vector< art::Ptr< recob::Vertex>     > selectedVertices;
    std::vector< art::Ptr< recob::Track>      > selectedTracks;
    std::vector< art::Ptr< recob::Shower>     > selectedShowers; 
    std::vector< art::Ptr< recob::PCAxis>     > selectedPCAxes;

    for ( art::Ptr< recob::PFParticle > part : selectedPFParticles ) {
        this->CollectAssociated( part, event.m_pfParticleSpacePointMap, selectedSpacePoints );
        this->CollectAssociated( part, event.m_pfParticleClusterMap   , selectedClusters    );
        this->CollectAssociated( part, event.m_pfParticleSeedMap      , selectedSeeds       );
        this->CollectAssociated( part, event.m_pfParticleVertexMap    , selectedVertices    );
        this->CollectAssociated( part, event.m_pfParticleTrackMap     , selectedTracks      );
        this->CollectAssociated( part, event.m_pfParticleShowerMap    , selectedShowers     );
        this->CollectAssociated( part, event.m_pfParticlePCAxisMap    , selectedPCAxes      );
    }

    m_spacePoints = selectedSpacePoints;
    m_clusters    = selectedClusters;
    m_seeds       = selectedSeeds;
    m_vertices    = selectedVertices;
    m_tracks      = selectedTracks;
    m_showers     = selectedShowers;
    m_pcAxes      = selectedPCAxes;

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

    this->GetFilteredAssociationMap( selectedPFParticles, selectedSpacePoints, event.m_pfParticleSpacePointMap, selectedPFParticleSpacePointMap );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedClusters   , event.m_pfParticleClusterMap   , selectedPFParticleClusterMap    );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedVertices   , event.m_pfParticleVertexMap    , selectedPFParticleVertexMap     );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedTracks     , event.m_pfParticleTrackMap     , selectedPFParticleTrackMap      );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedShowers    , event.m_pfParticleShowerMap    , selectedPFParticleShowerMap     );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedSeeds      , event.m_pfParticleSeedMap      , selectedPFParticleSeedMap       );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedPCAxes     , event.m_pfParticlePCAxisMap    , selectedPFParticlePCAxisMap     );
    this->GetFilteredAssociationMap( selectedSpacePoints, event.m_hits, event.m_spacePointHitMap, selectedSpacePointHitMap ); 
    this->GetFilteredAssociationMap( selectedClusters   , event.m_hits, event.m_clusterHitMap   , selectedClusterHitMap    );
    this->GetFilteredAssociationMap( selectedTracks     , event.m_hits, event.m_trackHitMap     , selectedTrackHitMap      );
    this->GetFilteredAssociationMap( selectedShowers    , event.m_hits, event.m_showerHitMap    , selectedShowerHitMap     );
    this->GetFilteredAssociationMap( selectedSeeds      , event.m_hits, event.m_seedHitMap      , selectedSeedHitMap       );
    this->GetFilteredAssociationMap( selectedShowers, selectedPCAxes, event.m_showerPCAxisMap, selectedShowerPCAxisMap );

    m_pfParticleSpacePointMap = selectedPFParticleSpacePointMap;
    m_pfParticleClusterMap    = selectedPFParticleClusterMap;
    m_pfParticleVertexMap     = selectedPFParticleVertexMap;
    m_pfParticleTrackMap      = selectedPFParticleTrackMap;
    m_pfParticleShowerMap     = selectedPFParticleShowerMap;
    m_pfParticleSeedMap       = selectedPFParticleSeedMap;
    m_pfParticlePCAxisMap     = selectedPFParticlePCAxisMap;
    m_spacePointHitMap        = selectedSpacePointHitMap;
    m_clusterHitMap           = selectedClusterHitMap;
    m_trackHitMap             = selectedTrackHitMap;
    m_showerHitMap            = selectedShowerHitMap;
    m_seedHitMap              = selectedSeedHitMap;
    m_showerPCAxisMap         = selectedShowerPCAxisMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::GetPrimaryPFParticles( std::vector< art::Ptr< recob::PFParticle > > & primaryPFParticles )
{
    for ( art::Ptr< recob::PFParticle > part : m_pfParticles ) 
        if ( part->IsPrimary() )
            primaryPFParticles.push_back( part );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::GetFilteredParticlesByPdgCode( const bool                                            shouldProduceNeutrinos, 
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

void LArPandoraEvent::GetFilteredParticlesByCRTag(   const bool                                            shouldProduceNeutrinos, 
                                                     const std::string &                                   tagProducerLabel,
                                                     const std::vector< art::Ptr< recob::PFParticle > > &  inputPFParticles, 
                                                     std::vector< art::Ptr< recob::PFParticle > > &        outputPFParticles )
{

    art::Handle< std::vector< recob::PFParticle > > pfParticleHandle;
    m_pEvent->getByLabel(m_labels.GetLabel( Labels::PFParticleLabel ), pfParticleHandle);

    art::FindManyP< anab::CosmicTag > pfParticleTagAssoc( pfParticleHandle, *m_pEvent, tagProducerLabel );
    
    for ( art::Ptr< recob::PFParticle > part : inputPFParticles ) {
        const std::vector< art::Ptr< anab::CosmicTag > > cosmicTags = pfParticleTagAssoc.at( part.key() );

        if (cosmicTags.size() != 1) 
            throw cet::exception("LArPandora") << " LArPandoraEvent::GetFilteredParticlesByCRTag -- Found " << cosmicTags.size() << " CR tags for a PFParticle (require 1).";

        art::Ptr< anab::CosmicTag > cosmicTag = cosmicTags.front();
        bool isNeutrino = ( cosmicTag->CosmicType() == anab::kNotTagged );

        if ( ( shouldProduceNeutrinos && isNeutrino ) || ( !shouldProduceNeutrinos && !isNeutrino ) ) {
            outputPFParticles.push_back( part );
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::GetIdToPFParticleMap( std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap ) 
{
    for ( art::Ptr< recob::PFParticle > part : m_pfParticles ) 
        if ( !idToPFParticleMap.insert( std::map< size_t, art::Ptr< recob::PFParticle > >::value_type( part->Self(), part ) ).second )
            throw cet::exception("LArPandora") << " LArPandoraEvent::GetIdToPFParticleMap -- Can't insert multiple entries with the same Id";
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::GetDownstreamPFParticles( const std::vector< art::Ptr< recob::PFParticle > > &       inputPFParticles, 
                                                const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                                std::map< size_t, art::Ptr< recob::PFParticle > > &        idToDownstreamPFParticleMap )
{
  for ( art::Ptr< recob::PFParticle > part : inputPFParticles )
    this->GetDownstreamPFParticles( part, idToPFParticleMap, idToDownstreamPFParticleMap );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::GetDownstreamPFParticles( art::Ptr< recob::PFParticle >                              part, 
                                                const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                                std::map< size_t, art::Ptr< recob::PFParticle > > &        idToDownstreamPFParticleMap )
{
    if ( idToDownstreamPFParticleMap.find( part->Self() ) == idToDownstreamPFParticleMap.end() )
        idToDownstreamPFParticleMap.insert( std::map< size_t, art::Ptr< recob::PFParticle > >::value_type( part->Self(), part ) );

    for ( size_t daughterId : part->Daughters() ) {
        std::map< size_t, art::Ptr< recob::PFParticle > >::const_iterator daughterIt = idToPFParticleMap.find( daughterId );

        if ( daughterIt == idToPFParticleMap.end() )
            throw cet::exception("LArPandora") << " LArPandoraEvent::GetDownstreamPFParticles -- Could not find daughter of PFParticle in the supplied map";

        this->GetDownstreamPFParticles( daughterIt->second, idToPFParticleMap, idToDownstreamPFParticleMap );
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::GetParticleVector( const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                         std::vector< art::Ptr< recob::PFParticle > > &             pfParticleVector )
{
    for ( std::map< size_t, art::Ptr< recob::PFParticle > >::const_iterator it=idToPFParticleMap.begin(); it != idToPFParticleMap.end(); ++it ) {
        art::Ptr< recob::PFParticle > part = it->second;
        pfParticleVector.push_back( part );
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::WriteToEvent()
{
    
    this->WriteCollection( m_pfParticles );
    this->WriteCollection( m_spacePoints );
    this->WriteCollection( m_clusters );   
    this->WriteCollection( m_vertices );   
    this->WriteCollection( m_tracks );     
    this->WriteCollection( m_showers );    
    this->WriteCollection( m_seeds );      
    this->WriteCollection( m_pcAxes );     
    
    this->WriteAssociation( m_pfParticleSpacePointMap, m_pfParticles, m_spacePoints );
    this->WriteAssociation( m_pfParticleClusterMap   , m_pfParticles, m_clusters    );   
    this->WriteAssociation( m_pfParticleVertexMap    , m_pfParticles, m_vertices    );    
    this->WriteAssociation( m_pfParticleTrackMap     , m_pfParticles, m_tracks      );     
    this->WriteAssociation( m_pfParticleShowerMap    , m_pfParticles, m_showers     );    
    this->WriteAssociation( m_pfParticleSeedMap      , m_pfParticles, m_seeds       );      
    this->WriteAssociation( m_pfParticlePCAxisMap    , m_pfParticles, m_pcAxes      );     
    this->WriteAssociation( m_spacePointHitMap       , m_spacePoints, m_hits        , m_labels.GetLabel( Labels::HitLabel ) );       
    this->WriteAssociation( m_clusterHitMap          , m_clusters   , m_hits        , m_labels.GetLabel( Labels::HitLabel ) );          
    this->WriteAssociation( m_trackHitMap            , m_tracks     , m_hits        , m_labels.GetLabel( Labels::HitLabel ) );            
    this->WriteAssociation( m_showerHitMap           , m_showers    , m_hits        , m_labels.GetLabel( Labels::HitLabel ) );           
    this->WriteAssociation( m_seedHitMap             , m_seeds      , m_hits        , m_labels.GetLabel( Labels::HitLabel ) );             
    this->WriteAssociation( m_showerPCAxisMap        , m_showers    , m_pcAxes      );        
      
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraEvent LArPandoraEvent::Merge( LArPandoraEvent & other )
{
    LArPandoraEvent outputEvent( other );

    std::vector< art::Ptr< recob::PFParticle > >  adjustedPFParticles;
    std::map< art::Ptr< recob::PFParticle >, art::Ptr< recob::PFParticle > >  adjustedPFParticleMap;

    this->AdjustIds( m_pfParticles, adjustedPFParticles, adjustedPFParticleMap);

    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::SpacePoint > > > adjustedPFParticleSpacePointMap;
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Cluster > > >    adjustedPFParticleClusterMap;
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Vertex > > >     adjustedPFParticleVertexMap;
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Track > > >      adjustedPFParticleTrackMap;
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Shower > > >     adjustedPFParticleShowerMap;
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Seed > > >       adjustedPFParticleSeedMap;
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PCAxis > > >     adjustedPFParticlePCAxisMap;

    this->AdjustAssociation( m_pfParticleSpacePointMap, adjustedPFParticleMap, adjustedPFParticleSpacePointMap );
    this->AdjustAssociation( m_pfParticleClusterMap   , adjustedPFParticleMap, adjustedPFParticleClusterMap   );
    this->AdjustAssociation( m_pfParticleVertexMap    , adjustedPFParticleMap, adjustedPFParticleVertexMap    );
    this->AdjustAssociation( m_pfParticleTrackMap     , adjustedPFParticleMap, adjustedPFParticleTrackMap     );
    this->AdjustAssociation( m_pfParticleShowerMap    , adjustedPFParticleMap, adjustedPFParticleShowerMap    );
    this->AdjustAssociation( m_pfParticleSeedMap      , adjustedPFParticleMap, adjustedPFParticleSeedMap      );
    this->AdjustAssociation( m_pfParticlePCAxisMap    , adjustedPFParticleMap, adjustedPFParticlePCAxisMap    );

    this->MergeCollection( other.m_pfParticles, adjustedPFParticles );
    this->MergeCollection( other.m_spacePoints, m_spacePoints       );
    this->MergeCollection( other.m_clusters   , m_clusters          );   
    this->MergeCollection( other.m_vertices   , m_vertices          );   
    this->MergeCollection( other.m_tracks     , m_tracks            );     
    this->MergeCollection( other.m_showers    , m_showers           );    
    this->MergeCollection( other.m_seeds      , m_seeds             );      
    this->MergeCollection( other.m_pcAxes     , m_pcAxes            );     
    
    this->MergeAssociation( other.m_pfParticleSpacePointMap, adjustedPFParticleSpacePointMap );
    this->MergeAssociation( other.m_pfParticleClusterMap   , adjustedPFParticleClusterMap    );   
    this->MergeAssociation( other.m_pfParticleVertexMap    , adjustedPFParticleVertexMap     );    
    this->MergeAssociation( other.m_pfParticleTrackMap     , adjustedPFParticleTrackMap      );     
    this->MergeAssociation( other.m_pfParticleShowerMap    , adjustedPFParticleShowerMap     );    
    this->MergeAssociation( other.m_pfParticleSeedMap      , adjustedPFParticleSeedMap       );      
    this->MergeAssociation( other.m_pfParticlePCAxisMap    , adjustedPFParticlePCAxisMap     );    
    this->MergeAssociation( other.m_spacePointHitMap       , m_spacePointHitMap              );       
    this->MergeAssociation( other.m_clusterHitMap          , m_clusterHitMap                 );          
    this->MergeAssociation( other.m_trackHitMap            , m_trackHitMap                   );            
    this->MergeAssociation( other.m_showerHitMap           , m_showerHitMap                  );           
    this->MergeAssociation( other.m_seedHitMap             , m_seedHitMap                    );             
    this->MergeAssociation( other.m_showerPCAxisMap        , m_showerPCAxisMap               );        

    return outputEvent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::AdjustIds( const std::vector< art::Ptr< recob::PFParticle > > &                        collection, 
                                 std::vector< art::Ptr< recob::PFParticle > > &                              adjustedCollection,
                                 std::map< art::Ptr< recob::PFParticle >, art::Ptr< recob::PFParticle > > &  adjustedPtrsMap )
{
    // Ensure that we have unique PFParticle IDs

    size_t shift = 100000;
    const lar::PtrMaker< recob::PFParticle > makePtrPFParticle( *m_pEvent, *m_pProducer );

    for ( const art::Ptr< recob::PFParticle > & part : collection ) {

        if ( part->Self() >= shift )
            throw cet::exception("LArPandora") << " LArPandoraEvent::AdjustIds -- PFParticle ID exceeds " << shift << ". Can't merge the collections!";

        size_t adjustedSelf   = part->Self()   + shift;
        size_t adjustedParent = part->Parent() + shift;

        const std::vector< size_t > daughters = part->Daughters();

        std::vector< size_t > adjustedDaughters;
        adjustedDaughters.reserve( daughters.size() );
        std::transform( daughters.begin(), daughters.end(), adjustedDaughters.begin(), [&]( size_t id ) -> size_t { return id + shift; } );

        /* BEGIN_DEBUG */
        for( unsigned int d = 0; d < daughters.size(); d++) {
            std::cout << ( ( adjustedDaughters[d] - daughters[d] == shift ) ? "Adjusting worked" : "Adjusting FAILED!!!!!!!!!!!" ) << std::endl;
        }
        /* END DEBUG */

        recob::PFParticle adjustedPart( part->PdgCode(), adjustedSelf, adjustedParent, adjustedDaughters );
        art::Ptr< recob::PFParticle > pAdjustedPart( makePtrPFParticle( adjustedCollection.size() ) );
        adjustedCollection.emplace_back( pAdjustedPart );

        if ( !adjustedPtrsMap.insert( std::map< art::Ptr< recob::PFParticle >, art::Ptr< recob::PFParticle > >::value_type( part, pAdjustedPart ) ).second )
            throw cet::exception("LArPandora") << " LArPandoraEvent::AdjustIds -- Supplied collection contains repeat PFParticles";
    }

}
//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraEvent::Labels::Labels( std::string pfParticleProducerLabel,
                                 std::string hitProducerLabel )
{
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleLabel, pfParticleProducerLabel) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( SpacePointLabel, pfParticleProducerLabel) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ClusterLabel   , pfParticleProducerLabel) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( VertexLabel    , pfParticleProducerLabel) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( TrackLabel     , pfParticleProducerLabel) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ShowerLabel    , pfParticleProducerLabel) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( SeedLabel      , pfParticleProducerLabel) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PCAxisLabel    , pfParticleProducerLabel) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( HitLabel       , hitProducerLabel       ) );

    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToSpacePointLabel, pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToClusterLabel   , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToVertexLabel    , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToTrackLabel     , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToShowerLabel    , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToSeedLabel      , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToPCAxisLabel    , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( SpacePointToHitLabel       , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ClusterToHitLabel          , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( TrackToHitLabel            , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ShowerToHitLabel           , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( SeedToHitLabel             , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ShowerToPCAxisLabel        , pfParticleProducerLabel ) );
}

//------------------------------------------------------------------------------------------------------------------------------------------
        
LArPandoraEvent::Labels::Labels( std::string pfParticleProducerLabel,
                                 std::string trackProducerLabel,
                                 std::string showerProducerLabel,
                                 std::string hitProducerLabel )
{
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleLabel, pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( SpacePointLabel, pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ClusterLabel   , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( VertexLabel    , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( TrackLabel     , trackProducerLabel      ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ShowerLabel    , showerProducerLabel     ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( SeedLabel      , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PCAxisLabel    , showerProducerLabel     ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( HitLabel       , hitProducerLabel        ) );

    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToSpacePointLabel, pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToClusterLabel   , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToVertexLabel    , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToTrackLabel     , trackProducerLabel      ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToShowerLabel    , showerProducerLabel     ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToSeedLabel      , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToPCAxisLabel    , showerProducerLabel     ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( SpacePointToHitLabel       , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ClusterToHitLabel          , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( TrackToHitLabel            , trackProducerLabel      ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ShowerToHitLabel           , showerProducerLabel     ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( SeedToHitLabel             , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ShowerToPCAxisLabel        , showerProducerLabel     ) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetSpacePointProducerLabel( const std::string & label )
{
    m_labels[SpacePointLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetClusterProducerLabel( const std::string & label )
{
    m_labels[ClusterLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetVertexProducerLabel( const std::string & label )
{
    m_labels[VertexLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetTrackProducerLabel( const std::string & label )
{
    m_labels[TrackLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetShowerProducerLabel( const std::string & label )
{
    m_labels[ShowerLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetSeedProducerLabel( const std::string & label )
{
    m_labels[SeedLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetPCAxisProducerLabel( const std::string & label )
{
    m_labels[PCAxisLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------
        
void LArPandoraEvent::Labels::SetPFParticleToSpacePointProducerLabel( const std::string & label )
{
    m_labels[PFParticleToSpacePointLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetPFParticleToClusterProducerLabel( const std::string & label )
{
    m_labels[PFParticleToClusterLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetPFParticleToVertexProducerLabel( const std::string & label )
{
    m_labels[PFParticleToVertexLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetPFParticleToTrackProducerLabel( const std::string & label )
{
    m_labels[PFParticleToTrackLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetPFParticleToShowerProducerLabel( const std::string & label )
{
    m_labels[PFParticleToShowerLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetPFParticleToSeedProducerLabel( const std::string & label )
{
    m_labels[PFParticleToSeedLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetPFParticleToPCAxisProducerLabel( const std::string & label )
{
    m_labels[PFParticleToPCAxisLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetSpacePointToHitProducerLabel( const std::string & label )
{
    m_labels[SpacePointToHitLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetClusterToHitProducerLabel( const std::string & label )
{
    m_labels[ClusterToHitLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetTrackToHitProducerLabel( const std::string & label )
{
    m_labels[TrackToHitLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetShowerToHitProducerLabel( const std::string & label )
{
    m_labels[ShowerToHitLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetSeedToHitProducerLabel( const std::string & label )
{
    m_labels[SeedToHitLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetShowerToPCAxisProducerLabel( const std::string & label )
{
    m_labels[ShowerToPCAxisLabel] = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string LArPandoraEvent::Labels::GetLabel( const LabelType & type )
{
    return m_labels[ type ];
}

//------------------------------------------------------------------------------------------------------------------------------------------



} // namespace lar_pandora
