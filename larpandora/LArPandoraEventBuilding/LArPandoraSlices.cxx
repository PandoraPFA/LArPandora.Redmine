/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraSlices.cxx
 *
 *  @brief  A mapping between two different reconstruction hypotheses on the same input hit collection
 */

#include "larpandora/LArPandoraEventBuilding/LArPandoraSlices.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

LArPandoraSlices::LArPandoraSlices( art::EDProducer *   pProducer,
                                    art::Event *        pEvent, 
                                    const std::string & crRecoProducerLabel,
                                    const std::string & nuRecoProducerLabel,
                                    const std::string & hitProducerLabel ) :
    m_pProducer( pProducer ),
    m_pEvent( pEvent ),
    m_crRecoProducerLabel( crRecoProducerLabel ),
    m_nuRecoProducerLabel( nuRecoProducerLabel ),
    m_hitProducerLabel( hitProducerLabel ),
    m_doesEventContainNeutrino( false )
{
    this->IdentifySlices(); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector< LArPandoraSlices::SliceId > LArPandoraSlices::GetSlices()
{
    std::vector< SliceId > slices( m_crSlicePFParticles.size() );
    std::iota( std::begin(slices), std::end(slices), 0 );
    return slices;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector< art::Ptr< recob::PFParticle > > LArPandoraSlices::GetSliceAsCR( const SliceId & id )
{
    if ( m_crSlicePFParticles.count( id ) == 0 )
        throw cet::exception("LArPandora") << " LArPandoraSlices::GetSliceAsCR -- Slice Id " << id << " is out of bounds.";

    return m_crSlicePFParticles.at( id );
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector< art::Ptr< recob::PFParticle > > LArPandoraSlices::GetSliceAsNu( const SliceId & id )
{
    if ( m_nuSlicePFParticles.count( id ) == 0 )
        throw cet::exception("LArPandora") << " LArPandoraSlices::GetSliceAsNu -- Slice Id " << id << " is out of bounds.";

    return m_nuSlicePFParticles.at( id );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::IdSliceAsNu( const SliceId & id )
{
    if ( m_nuSlicePFParticles.count( id ) == 0 )
        throw cet::exception("LArPandora") << " LArPandoraSlices::IdSliceAsNu -- Can't identify slice " << id << " as the neutrino. Slice Id out of bounds.";

    m_doesEventContainNeutrino = true;
    m_nuSliceId = id;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::WriteTags()
{
    for ( std::map< SliceId, std::vector< art::Ptr< recob::PFParticle > > >::const_iterator it = m_nuSlicePFParticles.begin(); it != m_nuSlicePFParticles.end(); ++it ) {
        SliceId sliceId = it->first;
        std::vector< art::Ptr< recob::PFParticle > > crPFParticles = this->GetSliceAsCR( sliceId );
        std::vector< art::Ptr< recob::PFParticle > > nuPFParticles = this->GetSliceAsNu( sliceId );

        if ( !m_doesEventContainNeutrino ) {
            this->WriteTag( false, crPFParticles );
            this->WriteTag( false, nuPFParticles );
        }
        else  {
            this->WriteTag( (sliceId == m_nuSliceId), crPFParticles );
            this->WriteTag( (sliceId == m_nuSliceId), nuPFParticles );
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::WriteTag( const bool &                                          shouldTagAsNeutrino, 
                                 const std::vector< art::Ptr< recob::PFParticle > > &  pfParticleVector )
{
    for ( const art::Ptr< recob::PFParticle > & part : pfParticleVector ) {

        (void) part; // TEMP
        if ( shouldTagAsNeutrino ) {
            // @todo produce an anab::kNotTagged
        }
        else {
            // @todo produce an anab::kUnknown
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::IdentifySlices()
{
    // Get neutrino PFParticles
    art::Handle< std::vector< recob::PFParticle > > nuPFParticleHandle;
    m_pEvent->getByLabel( m_nuRecoProducerLabel, nuPFParticleHandle); 

    // Get top-level nuPFParticles
    std::vector< art::Ptr< recob::PFParticle > > primaryNuPFParticles;
    this->GetPrimaryPFParticles( nuPFParticleHandle, primaryNuPFParticles);

    // Get mapping between PFParticles and their ID
    std::map< size_t, art::Ptr< recob::PFParticle > > nuIdToPFParticleMap;
    this->GetIdToPFParticleMap( nuPFParticleHandle, nuIdToPFParticleMap );

    // Get top-level nuPFParticles -> downstream nuPFParticles
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PFParticle > > > primaryToDownstreamNuPFParticles;
    this->GetDownstreamPFParticles( primaryNuPFParticles, nuIdToPFParticleMap, primaryToDownstreamNuPFParticles );

    // Get nuPFParticle -> nuClusters map
    std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Cluster > > >  nuPFParticleToClusterMap;
    this->GetPFParticleToClusterMap( nuPFParticleHandle, m_nuRecoProducerLabel, nuPFParticleToClusterMap );
    
    // Get nuCluster -> hits map
    
    // Get top-level nuPFParticle -> hits map

    // Get cosmic PFParticles
    
    // Get top-level crPFParticles

    // Get top-level crPFParticles -> downstream crPFParticles

    // Get crPFParticle -> crClusters map
    
    // Get crCluster -> hits map

    // Get hit -> top-level crPFParticle map

    // Get nuPFParticle -> crPFParticles map

    // Get sliceId to nuPFParticles map
    
    // Get sliceId to crPFParticles map
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::GetPrimaryPFParticles( const art::Handle< std::vector< recob::PFParticle > > &  pfParticleHandle, 
                                              std::vector< art::Ptr< recob::PFParticle > > &           primaryPFParticles )
{
    for( unsigned int i = 0; i != pfParticleHandle->size(); i++ ) {
        art::Ptr< recob::PFParticle > part( pfParticleHandle, i );

        if ( part->IsPrimary() )
            primaryPFParticles.push_back( part );
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::GetIdToPFParticleMap( const art::Handle< std::vector< recob::PFParticle > > &  pfParticleHandle,
                                             std::map< size_t, art::Ptr< recob::PFParticle > > &      idToPFParticleMap )
{
    for( unsigned int i = 0; i != pfParticleHandle->size(); i++ ) {
        art::Ptr< recob::PFParticle > part( pfParticleHandle, i );

        if ( !idToPFParticleMap.insert( std::map< size_t, art::Ptr< recob::PFParticle > >::value_type( part->Self(), part ) ).second )
            throw cet::exception("LArPandora") << " LArPandoraSlices::GetIdToPFParticleMap -- Can't insert multiple entries with the same Id";

    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void LArPandoraSlices::GetDownstreamPFParticles( const std::vector< art::Ptr< recob::PFParticle > > &                                      primaryPFParticles, 
                                                 const std::map< size_t, art::Ptr< recob::PFParticle > > &                                 idToPFParticleMap, 
                                                 std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PFParticle > > > & primaryToDownstreamPFParticles )
{
    
    for ( art::Ptr< recob::PFParticle > part : primaryPFParticles ) {
        
        std::vector< art::Ptr< recob::PFParticle > > downstreamPFParticles;
        this->GetDownstreamPFParticles( part, idToPFParticleMap, downstreamPFParticles );

        if ( !primaryToDownstreamPFParticles.insert( std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::PFParticle > > >::value_type( part, downstreamPFParticles ) ).second )
            throw cet::exception("LArPandora") << " LArPandoraSlices::GetDownstreamPFParticles -- Repeat PFParticles";

    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::GetDownstreamPFParticles( const art::Ptr< recob::PFParticle > &                      part, 
                                                 const std::map< size_t, art::Ptr< recob::PFParticle > > &  idToPFParticleMap, 
                                                 std::vector< art::Ptr< recob::PFParticle > > &             downstreamPFParticles )
{
    if ( std::find( downstreamPFParticles.begin(), downstreamPFParticles.end(), part ) == downstreamPFParticles.end() )
        downstreamPFParticles.push_back( part );

    for ( size_t daughterId : part->Daughters() ) {

        std::map< size_t, art::Ptr< recob::PFParticle > >::const_iterator daughterIt = idToPFParticleMap.find( daughterId );
        if ( daughterIt == idToPFParticleMap.end() )
            throw cet::exception("LArPandora") << " LArPandoraEvent::GetDownstreamPFParticles -- Could not find daughter of PFParticle in the supplied map";
        
        this->GetDownstreamPFParticles( daughterIt->second, idToPFParticleMap, downstreamPFParticles );
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::GetPFParticleToClusterMap( const art::Handle< std::vector< recob::PFParticle > > &                                pfParticleHandle,
                                                  const std::string &                                                                    producerLabel, 
                                                  std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Cluster > > > & pfParticleToClusterMap )
{
    art::FindManyP< recob::Cluster > assoc( pfParticleHandle, (*m_pEvent), producerLabel );

    for ( unsigned int i = 0; i < pfParticleHandle->size(); i++ ) {

        art::Ptr< recob::PFParticle > part( pfParticleHandle, i );

        if ( pfParticleToClusterMap.find( part ) == pfParticleToClusterMap.end() ) {
            std::vector< art::Ptr< recob::Cluster > > emptyVect;
            pfParticleToClusterMap.insert( std::map< art::Ptr< recob::PFParticle >, std::vector< art::Ptr< recob::Cluster > > >::value_type( part, emptyVect ) );
        }

        for ( art::Ptr< recob::Cluster > cluster : assoc.at( part.key() ) ) { 
            pfParticleToClusterMap[ part ].push_back( cluster );        
        }
    } 
}

} // namespace lar_pandora
