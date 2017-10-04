////////////////////////////////////////////////////////////////////////
// Class:       LArPandoraEventDump
// Plugin Type: analyzer (art v2_07_03)
// File:        LArPandoraEventDump_module.cc
//
// Generated at Tue Oct  3 12:55:39 2017 by Andrew D. Smith using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
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
#include "lardataobj/RecoBase/Cluster.h"

class LArPandoraEventDump;


class LArPandoraEventDump : public art::EDAnalyzer {
public:
  explicit LArPandoraEventDump(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LArPandoraEventDump(LArPandoraEventDump const &) = delete;
  LArPandoraEventDump(LArPandoraEventDump &&) = delete;
  LArPandoraEventDump & operator = (LArPandoraEventDump const &) = delete;
  LArPandoraEventDump & operator = (LArPandoraEventDump &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;

  void PrintParticle( const art::Ptr< recob::PFParticle > & part, const std::map< size_t, art::Ptr< recob::PFParticle > > & pfParticleIdMap, const art::FindManyP<recob::Cluster> & pfPartToClusterAssoc, const int & depth ) ;

private:

  // Declare member data here.
  std::string fPandoraLabel; 
};


LArPandoraEventDump::LArPandoraEventDump(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
    reconfigure(p);
}

void LArPandoraEventDump::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  
  std::cout << std::string(80, '-') << "\r- ";
  std::cout << "Event " << std::endl;
  std::cout << e.id() << std::endl;
  std::cout << fPandoraLabel << std::endl;
  std::cout << std::string(80, '-') << std::endl;

  // Get the PFParticles
  art::Handle< std::vector< recob::PFParticle > > pfParticleHandle;
  e.getByLabel( fPandoraLabel, pfParticleHandle );

  std::cout << "N PFParticles : " << pfParticleHandle->size() << std::endl;
  std::cout << std::string(80, '-') << std::endl;

  // Get the PFParticles ID map
  std::map< size_t, art::Ptr< recob::PFParticle > > pfParticleIdMap;
  for ( unsigned int i = 0; i < pfParticleHandle->size(); ++i ) {
    art::Ptr< recob::PFParticle > part( pfParticleHandle, i );
    pfParticleIdMap[ part->Self() ] = part;
  }

  // Get the associations
  art::FindManyP<recob::Cluster> pfPartToClusterAssoc( pfParticleHandle, e, fPandoraLabel );

  std::cout << "PFParticle -> Cluster associations size : " << pfPartToClusterAssoc.size() << std::endl;
  for ( unsigned int i = 0; i < pfPartToClusterAssoc.size(); i++ ) {
    auto clusters = pfPartToClusterAssoc.at( i );
    std::cout << "PFParticle index " << i << "  :  " << clusters.size() << " Clusters" << std::endl;
    for( auto cluster : clusters ) {
        std::cout << "    " << cluster->ID() << std::endl;
    }
  }

  // Output the PFParticle hierarchy
  for ( auto it = pfParticleIdMap.begin(); it != pfParticleIdMap.end(); ++it ) {
    art::Ptr< recob::PFParticle > part = it->second;
    if ( part->IsPrimary() ) {
        // Display Particle
        this->PrintParticle( part, pfParticleIdMap, pfPartToClusterAssoc, 0 );
    }
  }

  std::cout << std::string(80, '-') << std::endl;
}

void LArPandoraEventDump::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fPandoraLabel     = p.get<std::string>("PandoraLabel");
}

void LArPandoraEventDump::PrintParticle( const art::Ptr< recob::PFParticle > & part, const std::map< size_t, art::Ptr< recob::PFParticle > > & pfParticleIdMap, const art::FindManyP<recob::Cluster> & pfPartToClusterAssoc, const int & depth ) 
{
        std::cout << std::string( 4*depth, ' ') << "PFParticle ( " << part->Self() << " )" << std::endl;
        std::cout << std::string( 4*depth, ' ') << "- PDG code     : " << part->PdgCode() << std::endl;

        std::vector< art::Ptr< recob::Cluster > > clusters = pfPartToClusterAssoc.at( part.key() );
        std::cout << std::string( 4*depth, ' ') << "- N Clusters   : " << clusters.size() << std::endl;
        for ( auto cluster : clusters ) {
            std::cout << std::string( 4*depth, ' ') << "  - " << cluster->ID() << " ( " << cluster->View() << " )   : " << cluster->NHits() << std::endl;
        }

        if ( !part->IsPrimary() ) {
            std::cout << std::string( 4*depth, ' ') << "- Parent       : " << part->Parent() << std::endl;
        }
        std::cout << std::string( 4*depth, ' ') << "- N Daughters  : " << part->NumDaughters() << std::endl;
        std::cout << std::endl;

        for ( const auto & daughterId : part->Daughters() ) {
            art::Ptr< recob::PFParticle > daughter = pfParticleIdMap.at( daughterId );
            this->PrintParticle( daughter, pfParticleIdMap, pfPartToClusterAssoc, depth + 1 );
        }
}

DEFINE_ART_MODULE(LArPandoraEventDump)
