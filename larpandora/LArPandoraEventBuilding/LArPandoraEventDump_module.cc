/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraEventDump.cc
 *
 *  @brief  module for lar pandora event dump
 */

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/Hit.h"

namespace lar_pandora
{

class LArPandoraEventDump : public art::EDAnalyzer {
public:
    explicit LArPandoraEventDump(fhicl::ParameterSet const & pset);
    
    LArPandoraEventDump(LArPandoraEventDump const &) = delete;
    LArPandoraEventDump(LArPandoraEventDump &&) = delete;
    LArPandoraEventDump & operator = (LArPandoraEventDump const &) = delete;
    LArPandoraEventDump & operator = (LArPandoraEventDump &&) = delete;

    void analyze(art::Event const & evt) override;

private:
    void PrintParticle( const art::Ptr< recob::PFParticle > &                      part, 
                        const std::map< size_t, art::Ptr< recob::PFParticle > > &  pfParticleIdMap, 
                        const art::FindManyP<recob::SpacePoint> &                  pfPartToSpacePointAssoc,
                        const art::FindManyP<recob::Cluster> &                     pfPartToClusterAssoc,
                        const art::FindManyP<recob::Vertex> &                      pfPartToVertexAssoc,   
                        const art::FindManyP<recob::Track> &                       pfPartToTrackAssoc,   
                        const art::FindManyP<recob::Shower> &                      pfPartToShowerAssoc,    
                        const art::FindManyP<recob::PCAxis> &                      pfPartToPCAxisAssoc,
                        const int &                                                depth ); 

    std::string fPandoraLabel; 
    std::string fTrackLabel; 
    std::string fShowerLabel; 
};

DEFINE_ART_MODULE(LArPandoraEventDump)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows


namespace lar_pandora
{

LArPandoraEventDump::LArPandoraEventDump(fhicl::ParameterSet const & pset ) :
    EDAnalyzer(pset),
    fPandoraLabel(pset.get<std::string>("PandoraLabel")),
    fTrackLabel(pset.get<std::string>("TrackLabel" , fPandoraLabel)),
    fShowerLabel(pset.get<std::string>("ShowerLabel", fPandoraLabel)) {}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEventDump::analyze(art::Event const & evt)
{
  std::cout << std::endl << std::endl; 
  std::cout << std::string(80, '-') << "\r- ";
  std::cout << "Event " << std::endl;
  std::cout << evt.id() << std::endl;
  std::cout << fPandoraLabel << std::endl;
  std::cout << std::string(80, '-') << std::endl;

  // Get the input collections
  art::Handle< std::vector< recob::PFParticle > > pfParticleHandle;
  art::Handle< std::vector< recob::SpacePoint > > spacePointHandle;
  art::Handle< std::vector< recob::Cluster > > clusterHandle;
  art::Handle< std::vector< recob::Vertex > > vertexHandle;
  art::Handle< std::vector< recob::Track > > trackHandle;
  art::Handle< std::vector< recob::Shower > > showerHandle;
  art::Handle< std::vector< recob::PCAxis > > pcAxisHandle;

  evt.getByLabel( fPandoraLabel, pfParticleHandle );
  evt.getByLabel( fPandoraLabel, spacePointHandle );
  evt.getByLabel( fPandoraLabel, clusterHandle );
  evt.getByLabel( fPandoraLabel, vertexHandle );
  evt.getByLabel( fTrackLabel, trackHandle );
  evt.getByLabel( fShowerLabel, showerHandle );
  evt.getByLabel( fShowerLabel, pcAxisHandle );

  // Get the associations
  art::FindManyP<recob::SpacePoint> pfPartToSpacePointAssoc( pfParticleHandle, evt, fPandoraLabel );
  art::FindManyP<recob::Cluster>    pfPartToClusterAssoc(    pfParticleHandle, evt, fPandoraLabel );
  art::FindManyP<recob::Vertex>     pfPartToVertexAssoc(     pfParticleHandle, evt, fPandoraLabel );
  art::FindManyP<recob::Track>      pfPartToTrackAssoc(      pfParticleHandle, evt, fTrackLabel );
  art::FindManyP<recob::Shower>     pfPartToShowerAssoc(     pfParticleHandle, evt, fShowerLabel );
  art::FindManyP<recob::PCAxis>     pfPartToPCAxisAssoc(     pfParticleHandle, evt, fShowerLabel );

  art::FindManyP<recob::Hit> spacePointToHitAssoc( spacePointHandle, evt, fPandoraLabel );
  art::FindManyP<recob::Hit> clusterToHitAssoc(    clusterHandle   , evt, fPandoraLabel );
  art::FindManyP<recob::Hit> trackToHitAssoc(      trackHandle     , evt, fTrackLabel );
  art::FindManyP<recob::Hit> showerToHitAssoc(     showerHandle    , evt, fShowerLabel );

  art::FindManyP<recob::PCAxis> showerToPCAxisAssoc( showerHandle, evt, fShowerLabel );

  // Write out the collection sizes
  std::cout << "N PFParticles : " << pfParticleHandle->size() << std::endl;
  std::cout << "N SpacePoints : " << spacePointHandle->size() << std::endl;
  std::cout << "N Clusters    : " << clusterHandle->size()    << std::endl;
  std::cout << "N Vertices    : " << vertexHandle->size()     << std::endl;
  std::cout << "N Tracks      : " << trackHandle->size()      << std::endl;
  std::cout << "N Showers     : " << showerHandle->size()     << std::endl;
  std::cout << "N PCAxes      : " << pcAxisHandle->size()     << std::endl;
  std::cout << std::string(80, '-') << std::endl;

  // Get the PFParticles ID map
  std::map< size_t, art::Ptr< recob::PFParticle > > pfParticleIdMap;
  for ( unsigned int i = 0; i < pfParticleHandle->size(); ++i ) {
    art::Ptr< recob::PFParticle > part( pfParticleHandle, i );
    pfParticleIdMap[ part->Self() ] = part;
  }

  // Output the PFParticle hierarchy
  for ( auto it = pfParticleIdMap.begin(); it != pfParticleIdMap.end(); ++it ) {

    art::Ptr< recob::PFParticle > part = it->second;

    if ( part->IsPrimary() ) {
        // Display Particle
        this->PrintParticle( part, pfParticleIdMap, pfPartToSpacePointAssoc, pfPartToClusterAssoc, pfPartToVertexAssoc, pfPartToTrackAssoc, pfPartToShowerAssoc, pfPartToPCAxisAssoc, 0 );
    }
  }

  std::cout << std::string(80, '-') << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraEventDump::PrintParticle( const art::Ptr< recob::PFParticle > &                      part, 
                                         const std::map< size_t, art::Ptr< recob::PFParticle > > &  pfParticleIdMap, 
                                         const art::FindManyP<recob::SpacePoint> &                  pfPartToSpacePointAssoc,
                                         const art::FindManyP<recob::Cluster> &                     pfPartToClusterAssoc,
                                         const art::FindManyP<recob::Vertex> &                      pfPartToVertexAssoc,   
                                         const art::FindManyP<recob::Track> &                       pfPartToTrackAssoc,   
                                         const art::FindManyP<recob::Shower> &                      pfPartToShowerAssoc,    
                                         const art::FindManyP<recob::PCAxis> &                      pfPartToPCAxisAssoc,
                                         const int &                                                depth ) 
{
    int w = 16;    
    std::string indent( std::string( depth, ' ') );
    std::string rule( indent + std::string( 2*w, '-') );

    // Output the basic PFParticle information
    std::cout << std::endl << rule << std::endl;
    std::cout << indent << "PFParticle" << std::endl;
    std::cout << rule << std::endl;
    std::cout << indent << std::setw(w) << std::left << "- Key" << part.key() << std::endl;
    std::cout << indent << std::setw(w) << std::left << "- Id"  << part->Self() << std::endl;

    if ( part->IsPrimary() ) 
        std::cout << indent << std::setw(w) << std::left << "- Primary" << std::endl;
    else
        std::cout << indent << std::setw(w) << std::left << "- Parent" << part->Parent() << std::endl;

    std::cout << indent << std::setw(w) << "- PDG" << part->PdgCode() << std::endl;
    std::cout << rule << std::endl;

    // Get the associated objects
    std::vector< art::Ptr< recob::SpacePoint > > spacePoints = pfPartToSpacePointAssoc.at( part.key() );
    std::vector< art::Ptr< recob::Cluster > >    clusters    = pfPartToClusterAssoc.at( part.key() );
    std::vector< art::Ptr< recob::Vertex > >     vertices    = pfPartToVertexAssoc.at( part.key() );
    std::vector< art::Ptr< recob::Track > >      tracks      = pfPartToTrackAssoc.at( part.key() );
    std::vector< art::Ptr< recob::Shower > >     showers     = pfPartToShowerAssoc.at( part.key() );
    std::vector< art::Ptr< recob::PCAxis > >     pcAxes      = pfPartToPCAxisAssoc.at( part.key() );

    // Output the number of associated objects
    std::cout << indent << std::setw(w) << std::left << "- # SpacePoint" << spacePoints.size() << std::endl;
    std::cout << indent << std::setw(w) << std::left << "- # Cluster"    << clusters.size()    << std::endl;
    std::cout << indent << std::setw(w) << std::left << "- # Vertex"     << vertices.size()    << std::endl;
    std::cout << indent << std::setw(w) << std::left << "- # Track"      << tracks.size()      << std::endl;
    std::cout << indent << std::setw(w) << std::left << "- # Shower"     << showers.size()     << std::endl;
    std::cout << indent << std::setw(w) << std::left << "- # PCAxis"     << pcAxes.size()      << std::endl;
    std::cout << rule << std::endl;
    std::cout << indent << std::setw(w) << std::left << "- # Daughters"  << part->NumDaughters() << std::endl;
    std::cout << rule << std::endl;

    for ( const auto & daughterId : part->Daughters() ) {
        art::Ptr< recob::PFParticle > daughter = pfParticleIdMap.at( daughterId );
        this->PrintParticle( daughter, pfParticleIdMap, pfPartToSpacePointAssoc, pfPartToClusterAssoc, pfPartToVertexAssoc, pfPartToTrackAssoc, pfPartToShowerAssoc, pfPartToPCAxisAssoc, depth + 4 );
    }
}

} // namespace lar_pandora

