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
                                  const Labels &     inputLabels,
                                  const bool &       shouldProduceT0s,
                                  const size_t &     shift ) :
    m_pProducer( pProducer ),
    m_pEvent( pEvent ), 
    m_labels( inputLabels ),
    m_shouldProduceT0s( shouldProduceT0s ),
    m_shift( shift )
{

    this->GetCollections();

    // Make the origin of this LArPandoraEvent = 0
    for ( const art::Ptr< recob::PFParticle > & part : m_pfParticles ) 
        if ( !m_pfParticleToOriginIdMap.insert( std::map< art::Ptr< recob::PFParticle >, unsigned int >::value_type( part, 0 ) ).second )
            throw cet::exception("LArPandora") << " LArPandoraEvent::LArPandoraEvent -- Repeated input PFParticles!" << std::endl;

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

    art::Handle< std::vector< recob::PCAxis > > pcAxisHandle;
    this->GetCollection( Labels::PCAxisLabel, pcAxisHandle, m_pcAxes ); 

    art::Handle< std::vector< recob::Hit > > hitHandle;
    this->GetCollection( Labels::HitLabel, hitHandle, m_hits ); 

    if ( m_shouldProduceT0s ) {
        art::Handle< std::vector< anab::T0 > > t0Handle;
        this->GetCollection( Labels::T0Label, t0Handle, m_t0s ); 
    }


    this->GetAssociationMap( Labels::PFParticleToSpacePointLabel, pfParticleHandle, m_pfParticleSpacePointMap );
    this->GetAssociationMap( Labels::PFParticleToClusterLabel   , pfParticleHandle, m_pfParticleClusterMap    );
    this->GetAssociationMap( Labels::PFParticleToVertexLabel    , pfParticleHandle, m_pfParticleVertexMap     );
    this->GetAssociationMap( Labels::PFParticleToTrackLabel     , pfParticleHandle, m_pfParticleTrackMap      );
    this->GetAssociationMap( Labels::PFParticleToShowerLabel    , pfParticleHandle, m_pfParticleShowerMap     );
    this->GetAssociationMap( Labels::PFParticleToPCAxisLabel    , pfParticleHandle, m_pfParticlePCAxisMap     );

    if ( m_shouldProduceT0s ) 
        this->GetAssociationMap( Labels::PFParticleToT0Label        , pfParticleHandle, m_pfParticleT0Map         );

    this->GetAssociationMap( Labels::SpacePointToHitLabel, spacePointHandle, m_spacePointHitMap );
    this->GetAssociationMap( Labels::ClusterToHitLabel   , clusterHandle   , m_clusterHitMap    );
    this->GetAssociationMap( Labels::TrackToHitLabel     , trackHandle     , m_trackHitMap      );
    this->GetAssociationMap( Labels::ShowerToHitLabel    , showerHandle    , m_showerHitMap     );

    this->GetAssociationMap( Labels::ShowerToPCAxisLabel, showerHandle, m_showerPCAxisMap );


    this->GetPFParticleHierarchy();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::GetPFParticleHierarchy()
{
    std::map< size_t, art::Ptr< recob::PFParticle > > idToPFParticleMap;
    this->GetIdToPFParticleMap( idToPFParticleMap );

    for ( const art::Ptr< recob::PFParticle > & part : m_pfParticles ) {
        std::vector< art::Ptr< recob::PFParticle > > daughters;
        if ( !m_pfParticleDaughterMap.insert( std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PFParticle > > >::value_type( part, daughters ) ).second )
            throw cet::exception("LArPandora") << " LArPandoraEvent::GetPFParticleHierarchy -- Repeated PFParticle in heirarchy map!" << std::endl;

        for ( const size_t & daughterId : part->Daughters() ) {
            if ( idToPFParticleMap.find( daughterId ) == idToPFParticleMap.end() )
                throw cet::exception("LArPandora") << " LArPandoraEvent::GetPFParticleHierarchy -- Can't access map entry for daughter of PFParticle supplied." << std::endl;

            art::Ptr< recob::PFParticle > daughter = idToPFParticleMap.at( daughterId );
            if ( std::find( m_pfParticleDaughterMap[ part ].begin(), m_pfParticleDaughterMap[ part ].end(), daughter ) != m_pfParticleDaughterMap[ part ].end() )
                throw cet::exception("LArPandora") << " LArPandoraEvent::GetPFParticleHierarchy -- Can't have the same daughter twice!" << std::endl;

            m_pfParticleDaughterMap[ part ].push_back( daughter );
        }        
    }
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
    std::vector< art::Ptr< recob::PFParticle > > selectedPFParticles;
    this->GetDownstreamPFParticles( filteredPFParticles, selectedPFParticles );

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
    std::vector< art::Ptr< recob::PFParticle > > selectedPFParticles;
    this->GetDownstreamPFParticles( filteredPFParticles, selectedPFParticles );

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
    m_shouldProduceT0s( event.m_shouldProduceT0s ),
    m_shift( event.m_shift ),
    m_hits( event.m_hits )
{


    m_pfParticles = selectedPFParticles;

    this->FillPFParticleToOriginIdMap( event.m_pfParticleToOriginIdMap );

    std::vector< art::Ptr< recob::SpacePoint> > selectedSpacePoints;
    std::vector< art::Ptr< recob::Cluster>    > selectedClusters;
    std::vector< art::Ptr< recob::Vertex>     > selectedVertices;
    std::vector< art::Ptr< recob::Track>      > selectedTracks;
    std::vector< art::Ptr< recob::Shower>     > selectedShowers; 
    std::vector< art::Ptr< recob::PCAxis>     > selectedPCAxes;
    std::vector< art::Ptr< anab::T0 >         > selectedT0s;

    for ( art::Ptr< recob::PFParticle > part : selectedPFParticles ) {
        this->CollectAssociated( part, event.m_pfParticleSpacePointMap, selectedSpacePoints );
        this->CollectAssociated( part, event.m_pfParticleClusterMap   , selectedClusters    );
        this->CollectAssociated( part, event.m_pfParticleVertexMap    , selectedVertices    );
        this->CollectAssociated( part, event.m_pfParticleTrackMap     , selectedTracks      );
        this->CollectAssociated( part, event.m_pfParticleShowerMap    , selectedShowers     );
        this->CollectAssociated( part, event.m_pfParticlePCAxisMap    , selectedPCAxes      );

        if ( m_shouldProduceT0s ) 
            this->CollectAssociated( part, event.m_pfParticleT0Map        , selectedT0s         );
    }

    m_spacePoints = selectedSpacePoints;
    m_clusters    = selectedClusters;
    m_vertices    = selectedVertices;
    m_tracks      = selectedTracks;
    m_showers     = selectedShowers;
    m_pcAxes      = selectedPCAxes;

    if ( m_shouldProduceT0s ) 
        m_t0s         = selectedT0s;

    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::SpacePoint > > >    selectedPFParticleSpacePointMap;
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Cluster > > >       selectedPFParticleClusterMap;  
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Vertex > > >        selectedPFParticleVertexMap;  
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Track > > >         selectedPFParticleTrackMap;  
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Shower > > >        selectedPFParticleShowerMap;
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PCAxis > > >        selectedPFParticlePCAxisMap;
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< anab::T0 > > >             selectedPFParticleT0Map; 

    std::map< art::Ptr< recob::SpacePoint >, std::vector< art::Ptr< recob::Hit > > >           selectedSpacePointHitMap;  
    std::map< art::Ptr< recob::Cluster >   , std::vector< art::Ptr< recob::Hit > > >           selectedClusterHitMap;    
    std::map< art::Ptr< recob::Track >     , std::vector< art::Ptr< recob::Hit > > >           selectedTrackHitMap;     
    std::map< art::Ptr< recob::Shower >    , std::vector< art::Ptr< recob::Hit > > >           selectedShowerHitMap;
    std::map< art::Ptr< recob::Shower >    , std::vector< art::Ptr< recob::PCAxis > > >        selectedShowerPCAxisMap;

    this->GetFilteredAssociationMap( selectedPFParticles, selectedSpacePoints, event.m_pfParticleSpacePointMap, selectedPFParticleSpacePointMap );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedClusters   , event.m_pfParticleClusterMap   , selectedPFParticleClusterMap    );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedVertices   , event.m_pfParticleVertexMap    , selectedPFParticleVertexMap     );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedTracks     , event.m_pfParticleTrackMap     , selectedPFParticleTrackMap      );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedShowers    , event.m_pfParticleShowerMap    , selectedPFParticleShowerMap     );
    this->GetFilteredAssociationMap( selectedPFParticles, selectedPCAxes     , event.m_pfParticlePCAxisMap    , selectedPFParticlePCAxisMap     );

    if ( m_shouldProduceT0s ) 
        this->GetFilteredAssociationMap( selectedPFParticles, selectedT0s        , event.m_pfParticleT0Map        , selectedPFParticleT0Map         );

    this->GetFilteredAssociationMap( selectedSpacePoints, event.m_hits, event.m_spacePointHitMap, selectedSpacePointHitMap ); 
    this->GetFilteredAssociationMap( selectedClusters   , event.m_hits, event.m_clusterHitMap   , selectedClusterHitMap    );
    this->GetFilteredAssociationMap( selectedTracks     , event.m_hits, event.m_trackHitMap     , selectedTrackHitMap      );
    this->GetFilteredAssociationMap( selectedShowers    , event.m_hits, event.m_showerHitMap    , selectedShowerHitMap     );
    this->GetFilteredAssociationMap( selectedShowers, selectedPCAxes, event.m_showerPCAxisMap, selectedShowerPCAxisMap );

    m_pfParticleSpacePointMap = selectedPFParticleSpacePointMap;
    m_pfParticleClusterMap    = selectedPFParticleClusterMap;
    m_pfParticleVertexMap     = selectedPFParticleVertexMap;
    m_pfParticleTrackMap      = selectedPFParticleTrackMap;
    m_pfParticleShowerMap     = selectedPFParticleShowerMap;
    m_pfParticlePCAxisMap     = selectedPFParticlePCAxisMap;

    if ( m_shouldProduceT0s ) 
        m_pfParticleT0Map         = selectedPFParticleT0Map;

    m_spacePointHitMap        = selectedSpacePointHitMap;
    m_clusterHitMap           = selectedClusterHitMap;
    m_trackHitMap             = selectedTrackHitMap;
    m_showerHitMap            = selectedShowerHitMap;
    m_showerPCAxisMap         = selectedShowerPCAxisMap;

    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PFParticle > > >  selectedPFParticleDaughterMap;
    this->GetFilteredHierarchyMap( selectedPFParticles, event.m_pfParticleDaughterMap, selectedPFParticleDaughterMap );

    m_pfParticleDaughterMap = selectedPFParticleDaughterMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::GetFilteredHierarchyMap( const std::vector< art::Ptr< recob::PFParticle > > &                                             filteredParticles,
                                               const std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PFParticle > > > &  unfilteredPFParticleDaughterMap,
                                               std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PFParticle > > > &        outputPFParticleDaughterMap )
{
    for ( std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PFParticle > > >::const_iterator it = unfilteredPFParticleDaughterMap.begin(); it != unfilteredPFParticleDaughterMap.end(); ++it ) {
        if ( std::find( filteredParticles.begin(), filteredParticles.end(), it->first ) == filteredParticles.end() ) continue;

        if ( ! outputPFParticleDaughterMap.insert( std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PFParticle > > >::value_type( it->first, it->second ) ).second )
            throw cet::exception("LArPandora") << " LArPandoraEvent::GetFilteredHierarchyMap -- Can't add multiple map entries for same PFParticle" << std::endl;

    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::FillPFParticleToOriginIdMap( const std::map< art::Ptr< recob::PFParticle >, unsigned int > & existingMap )
{
    for ( const art::Ptr< recob::PFParticle > & part : m_pfParticles ) { 
        if ( existingMap.find(part) == existingMap.end() ) 
            throw cet::exception("LArPandora") << " LArPandoraEvent::FillPFParticleToOriginIdMap -- Can't access map entry for PFParticle supplied." << std::endl;

        if ( !m_pfParticleToOriginIdMap.insert( std::map< art::Ptr< recob::PFParticle >, unsigned int >::value_type( part, existingMap.at(part) ) ).second )
            throw cet::exception("LArPandora") << " LArPandoraEvent::FillPFParticleToOriginIdMap -- Can't add multiple map entries for same PFParticle" << std::endl;
    }
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
            throw cet::exception("LArPandora") << " LArPandoraEvent::GetFilteredParticlesByCRTag -- Found " << cosmicTags.size() << " CR tags for a PFParticle (require 1)." << std::endl;

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
            throw cet::exception("LArPandora") << " LArPandoraEvent::GetIdToPFParticleMap -- Can't insert multiple entries with the same Id" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::GetDownstreamPFParticles( const std::vector< art::Ptr< recob::PFParticle > > &  inputPFParticles,
                                                std::vector< art::Ptr< recob::PFParticle > > &        downstreamPFParticles )
{
  for ( art::Ptr< recob::PFParticle > part : inputPFParticles )
    this->GetDownstreamPFParticles( part, downstreamPFParticles );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::GetDownstreamPFParticles( const art::Ptr< recob::PFParticle > &           part,
                                                std::vector< art::Ptr< recob::PFParticle > > &  downstreamPFParticles )
{
    if ( m_pfParticleDaughterMap.find( part ) == m_pfParticleDaughterMap.end() )
        throw cet::exception("LArPandora") << " LArPandoraEvent::GetDownstreamPFParticles -- Could not find PFParticle in the hierarchy map" << std::endl;

    if ( std::find( downstreamPFParticles.begin(), downstreamPFParticles.end(), part ) == downstreamPFParticles.end() )
        downstreamPFParticles.push_back( part );

    for ( const art::Ptr< recob::PFParticle > & daughter : m_pfParticleDaughterMap.at( part ) )
        this->GetDownstreamPFParticles( daughter, downstreamPFParticles );
}

//------------------------------------------------------------------------------------------------------------------------------------------

/* ATTN now not needed */
void LArPandoraEvent::GetDownstreamPFParticles( const std::vector< art::Ptr< recob::PFParticle > > &       inputPFParticles, 
                                                const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                                std::map< size_t, art::Ptr< recob::PFParticle > > &        idToDownstreamPFParticleMap )
{
  for ( art::Ptr< recob::PFParticle > part : inputPFParticles )
    this->GetDownstreamPFParticles( part, idToPFParticleMap, idToDownstreamPFParticleMap );
}

//------------------------------------------------------------------------------------------------------------------------------------------

/* ATTN now not needed */
void LArPandoraEvent::GetDownstreamPFParticles( art::Ptr< recob::PFParticle >                              part, 
                                                const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                                std::map< size_t, art::Ptr< recob::PFParticle > > &        idToDownstreamPFParticleMap )
{
    if ( idToDownstreamPFParticleMap.find( part->Self() ) == idToDownstreamPFParticleMap.end() )
        idToDownstreamPFParticleMap.insert( std::map< size_t, art::Ptr< recob::PFParticle > >::value_type( part->Self(), part ) );

    for ( size_t daughterId : part->Daughters() ) {
        std::map< size_t, art::Ptr< recob::PFParticle > >::const_iterator daughterIt = idToPFParticleMap.find( daughterId );

        if ( daughterIt == idToPFParticleMap.end() )
            throw cet::exception("LArPandora") << " LArPandoraEvent::GetDownstreamPFParticles -- Could not find daughter of PFParticle in the supplied map" << std::endl;

        this->GetDownstreamPFParticles( daughterIt->second, idToPFParticleMap, idToDownstreamPFParticleMap );
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

/* ATTN now not needed */
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
    this->WriteCollection( m_clusters    );   
    this->WriteCollection( m_vertices    );   
    this->WriteCollection( m_tracks      );     
    this->WriteCollection( m_showers     );    
    this->WriteCollection( m_pcAxes      );     

    if ( m_shouldProduceT0s ) 
        this->WriteCollection( m_t0s         );      
    
    this->WriteAssociation( m_pfParticleSpacePointMap, m_pfParticles, m_spacePoints );
    this->WriteAssociation( m_pfParticleClusterMap   , m_pfParticles, m_clusters    );   
    this->WriteAssociation( m_pfParticleVertexMap    , m_pfParticles, m_vertices    );    
    this->WriteAssociation( m_pfParticleTrackMap     , m_pfParticles, m_tracks      );     
    this->WriteAssociation( m_pfParticleShowerMap    , m_pfParticles, m_showers     );    
    this->WriteAssociation( m_pfParticlePCAxisMap    , m_pfParticles, m_pcAxes      );    

    if ( m_shouldProduceT0s ) 
        this->WriteAssociation( m_pfParticleT0Map        , m_pfParticles, m_t0s         );      

    this->WriteAssociation( m_spacePointHitMap       , m_spacePoints, m_hits        , false );       
    this->WriteAssociation( m_clusterHitMap          , m_clusters   , m_hits        , false );          
    this->WriteAssociation( m_trackHitMap            , m_tracks     , m_hits        , false );            
    this->WriteAssociation( m_showerHitMap           , m_showers    , m_hits        , false );           
    this->WriteAssociation( m_showerPCAxisMap        , m_showers    , m_pcAxes      );        
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraEvent LArPandoraEvent::Merge( LArPandoraEvent & other )
{
    if ( m_shift != other.m_shift )
        throw cet::exception("LArPandora") << " LArPandoraEvent::Merge - Can't merge LArPandoraEvents with differing shift values." << std::endl;

    LArPandoraEvent outputEvent( other );

    this->MergePFParticleToOriginIdMap( outputEvent.m_pfParticleToOriginIdMap, m_pfParticleToOriginIdMap );


    this->MergeCollection( outputEvent.m_pfParticles, m_pfParticles );
    this->MergeCollection( outputEvent.m_spacePoints, m_spacePoints );
    this->MergeCollection( outputEvent.m_clusters   , m_clusters    );   
    this->MergeCollection( outputEvent.m_vertices   , m_vertices    );   
    this->MergeCollection( outputEvent.m_tracks     , m_tracks      );     
    this->MergeCollection( outputEvent.m_showers    , m_showers     );    
    this->MergeCollection( outputEvent.m_pcAxes     , m_pcAxes      );    
    this->MergeCollection( outputEvent.m_hits       , m_hits        );    

    if ( m_shouldProduceT0s ) 
        this->MergeCollection( outputEvent.m_t0s, m_t0s );      
    
    this->MergeAssociation( outputEvent.m_pfParticleSpacePointMap, m_pfParticleSpacePointMap );
    this->MergeAssociation( outputEvent.m_pfParticleClusterMap   , m_pfParticleClusterMap    );   
    this->MergeAssociation( outputEvent.m_pfParticleVertexMap    , m_pfParticleVertexMap     );    
    this->MergeAssociation( outputEvent.m_pfParticleTrackMap     , m_pfParticleTrackMap      );     
    this->MergeAssociation( outputEvent.m_pfParticleShowerMap    , m_pfParticleShowerMap     );    
    this->MergeAssociation( outputEvent.m_pfParticlePCAxisMap    , m_pfParticlePCAxisMap     );   

    if ( m_shouldProduceT0s ) 
        this->MergeAssociation( outputEvent.m_pfParticleT0Map, m_pfParticleT0Map );      

    this->MergeAssociation( outputEvent.m_spacePointHitMap, m_spacePointHitMap );       
    this->MergeAssociation( outputEvent.m_clusterHitMap   , m_clusterHitMap    );          
    this->MergeAssociation( outputEvent.m_trackHitMap     , m_trackHitMap      );            
    this->MergeAssociation( outputEvent.m_showerHitMap    , m_showerHitMap     );           
    this->MergeAssociation( outputEvent.m_showerPCAxisMap , m_showerPCAxisMap  );        

    return outputEvent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEvent::MergePFParticleToOriginIdMap( std::map< art::Ptr< recob::PFParticle >, unsigned int > &         mapToMerge, 
                                                    const std::map< art::Ptr< recob::PFParticle >, unsigned int > &   mapToAdd )
{

    unsigned int maxID = 0;
    if ( mapToMerge.size() != 0 ) {
        maxID = std::max_element( mapToMerge.begin(), mapToMerge.end(), [](const std::pair< art::Ptr< recob::PFParticle >, unsigned int > & p1, const std::pair< art::Ptr< recob::PFParticle >, unsigned int  > & p2) { return p1.second < p2.second; } )->second;
    }

    for ( std::map< art::Ptr< recob::PFParticle >, unsigned int >::const_iterator it=mapToAdd.begin(); it != mapToAdd.end(); ++it ) {
        if ( !mapToMerge.insert( std::map< art::Ptr< recob::PFParticle >, unsigned int >::value_type( it->first, it->second + maxID + 1 ) ).second )
            throw cet::exception("LArPandora") << " LArPandoraEvent::MergePFParticleToOriginIdMap - Can't merge collections containing repeated PFParticles." << std::endl;
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
    m_labels.insert( std::map< LabelType, std::string >::value_type( T0Label        , pfParticleProducerLabel) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PCAxisLabel    , pfParticleProducerLabel) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( HitLabel       , hitProducerLabel       ) );

    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToSpacePointLabel, pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToClusterLabel   , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToVertexLabel    , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToTrackLabel     , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToShowerLabel    , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToT0Label        , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToPCAxisLabel    , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( SpacePointToHitLabel       , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ClusterToHitLabel          , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( TrackToHitLabel            , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ShowerToHitLabel           , pfParticleProducerLabel ) );
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
    m_labels.insert( std::map< LabelType, std::string >::value_type( T0Label        , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PCAxisLabel    , showerProducerLabel     ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( HitLabel       , hitProducerLabel        ) );

    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToSpacePointLabel, pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToClusterLabel   , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToVertexLabel    , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToTrackLabel     , trackProducerLabel      ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToShowerLabel    , showerProducerLabel     ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToT0Label        , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( PFParticleToPCAxisLabel    , showerProducerLabel     ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( SpacePointToHitLabel       , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ClusterToHitLabel          , pfParticleProducerLabel ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( TrackToHitLabel            , trackProducerLabel      ) );
    m_labels.insert( std::map< LabelType, std::string >::value_type( ShowerToHitLabel           , showerProducerLabel     ) );
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

void LArPandoraEvent::Labels::SetT0ProducerLabel( const std::string & label )
{
    m_labels[T0Label] = label;
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

void LArPandoraEvent::Labels::SetPFParticleToT0ProducerLabel( const std::string & label )
{
    m_labels[PFParticleToT0Label] = label;
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
