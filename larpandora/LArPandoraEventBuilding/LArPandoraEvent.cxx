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
                                         std::string        inputProducerLabel,
                                         std::string        hitProducerLabel ) :
    m_pProducer( pProducer ),
    m_pEvent( pEvent ), 
    m_inputProducerLabel( inputProducerLabel ),
    m_hitProducerLabel( hitProducerLabel )
{

    this->GetCollections();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::GetCollections()
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

    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, m_pfParticleSpacePointMap );
    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, m_pfParticleClusterMap    );
    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, m_pfParticleVertexMap     );
    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, m_pfParticleTrackMap      );
    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, m_pfParticleShowerMap     );
    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, m_pfParticleSeedMap       );
    this->GetAssociationMap( m_inputProducerLabel, pfParticleHandle, m_pfParticlePCAxisMap     );

    this->GetAssociationMap( m_inputProducerLabel, spacePointHandle, m_spacePointHitMap );
    this->GetAssociationMap( m_inputProducerLabel, clusterHandle   , m_clusterHitMap    );
    this->GetAssociationMap( m_inputProducerLabel, trackHandle     , m_trackHitMap      );
    this->GetAssociationMap( m_inputProducerLabel, showerHandle    , m_showerHitMap     );
    this->GetAssociationMap( m_inputProducerLabel, seedHandle      , m_seedHitMap       );

    this->GetAssociationMap( m_inputProducerLabel, showerHandle, m_showerPCAxisMap );
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
    LArPandoraEvent filteredEvent( *this );
    this->FilterByParticleSelection( filteredEvent, selectedPFParticles );

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
    LArPandoraEvent filteredEvent( *this );
    this->FilterByParticleSelection( filteredEvent, selectedPFParticles );

    return filteredEvent;
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

void LArPandoraEvent::GetIdToPFParticleMap( std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap ) 
{
    for ( art::Ptr< recob::PFParticle > part : m_pfParticles ) 
      idToPFParticleMap.insert( std::map< size_t, art::Ptr< recob::PFParticle > >::value_type( part->Self(), part ) );
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

void LArPandoraEvent::GetDownstreamPFParticles( art::Ptr< recob::PFParticle >             part, 
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

void LArPandoraEvent::GetParticleVector( const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                                std::vector< art::Ptr< recob::PFParticle > > &             pfParticleVector )
{
    for ( std::map< size_t, art::Ptr< recob::PFParticle > >::const_iterator it=idToPFParticleMap.begin(); it != idToPFParticleMap.end(); ++it ) {
        art::Ptr< recob::PFParticle > part = it->second;
        pfParticleVector.push_back( part );
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::FilterByParticleSelection( LArPandoraEvent &                               filteredEvent, 
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

LArPandoraEvent LArPandoraEvent::Merge( LArPandoraEvent & other )
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

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraEvent::Labels::Labels( std::string pfParticleProducerLabel,
                                 std::string hitProducerLabel ) :

    m_pfParticleLabel ( pfParticleProducerLabel ),
    m_spacePointLabel ( pfParticleProducerLabel ),
    m_clusterLabel ( pfParticleProducerLabel ),
    m_vertexLabel ( pfParticleProducerLabel ),
    m_trackLabel ( pfParticleProducerLabel ),
    m_showerLabel ( pfParticleProducerLabel ),
    m_seedLabel ( pfParticleProducerLabel ),
    m_pcAxisLabel ( pfParticleProducerLabel ),
    m_hitLabel ( hitProducerLabel )
{}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetSpacePointProducerLabel( const std::string & label )
{
    m_spacePointLabel = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetClusterProducerLabel( const std::string & label )
{
    m_clusterLabel = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetVertexProducerLabel( const std::string & label )
{
    m_vertexLabel = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetTrackProducerLabel( const std::string & label )
{
    m_trackLabel = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetShowerProducerLabel( const std::string & label )
{
    m_showerLabel = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetSeedProducerLabel( const std::string & label )
{
    m_seedLabel = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::Labels::SetPCAxisProducerLabel( const std::string & label )
{
    m_pcAxisLabel = label;
}

//------------------------------------------------------------------------------------------------------------------------------------------


} // namespace lar_pandora
