////////////////////////////////////////////////////////////////////////
// Class:       PandoraSplitting
// Plugin Type: producer (art v2_07_03)
// File:        PandoraSplitting_module.cc
//
// Generated at Thu Sep 21 13:56:07 2017 by Andrew D. Smith using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

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

/*!
 *  \breif   A module that splits collections from a single Pandora instance into multiple collections
 *
 *  This module is used to split the output of the consolidated Pandora reconstruction into separate collections of "neutrinos" and "cosmic-rays".
 *  It can be configured to produce a copy of all Pandora-produced objects associated with a top-level PFParticle that is:
 *
 *  1. ...       identified as a neutrino by Pandora   (set ShouldProduceNeutrinos to true)
 *  2. ... *not* identified as a neutrino by Pandora   (set ShouldProduceCosmics   to true)
 */
class PandoraSplitting;

class PandoraSplitting : public art::EDProducer {
public:
  explicit PandoraSplitting(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PandoraSplitting(PandoraSplitting const &) = delete;
  PandoraSplitting(PandoraSplitting &&) = delete;
  PandoraSplitting & operator = (PandoraSplitting const &) = delete;
  PandoraSplitting & operator = (PandoraSplitting &&) = delete;

  // Required functions.
  virtual void reconfigure(fhicl::ParameterSet const & p);
  void produce(art::Event & e) override;

private:

  /*!
   *  \brief  Filters primary PFParticles from an input list
   *
   *  \param  allPFParticles      input vector of all PFParticles in the event
   *  \param  primaryPFParticles  output vector of all primary PFParticles in the input vector
   */
  void GetPrimaryPFParticles( const art::Handle< std::vector< recob::PFParticle > > & allPFParticles, std::vector< art::Ptr< recob::PFParticle > > & primaryPFParticles );

  /*!
   *  \brief  Determines if the supplied PFParticle should be persisted into the new collection (given ShouldProduceNeutrinos and ShouldProduceCosmics)
   *
   *  \param  part  the PFParticle in question
   *
   *  \return whether the supplied PFParticle should be added to the new collection
   */
  bool ShouldPersist( const art::Ptr< recob::PFParticle > & part );

  /*!
   *  \brief  Filters the PFParticles that should be persisted from an input list
   *
   *  \param  inputPFParticles       input vector of PFParticles
   *  \param  persistentPFParticles  output vector of all PFParticles in the input vector that should be persisted
   */
  void GetPFParticlesToPersist( const std::vector< art::Ptr< recob::PFParticle > > & inputPFParticles, std::vector< art::Ptr< recob::PFParticle > > & persistentPFParticles );

  /*!
   *  \brief  Produce a mapping between PFParticles and their ID
   *
   *  \param  allPFParticles      input vector of all PFParticles in the event
   *  \param  idToPFParticleMap   output mapping between PFParticles and their IDs
   */
  void GetIdToPFParticleMap( const art::Handle< std::vector< recob::PFParticle > > & allPFParticles, std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap );

  /*!
   *  \brief  Collects all PFParticles downstream (children, grandchildren, ...) of a given PFParticle
   *
   *  \param  part                         input PFParticle
   *  \param  idToPFParticleMap            input mapping between PFParticles and their IDs
   *  \param  idToDownstreamPFParticleMap  output mapping between of all downstream PFParticles and their IDs
   */
  void GetDownstreamPFParticles( art::Ptr< recob::PFParticle > part, const std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap, std::map< size_t, art::Ptr< recob::PFParticle > > & idToDownstreamPFParticleMap );

  /*!
   *  \brief  Collects all PFParticles downstream (children, grandchildren, ...) of a given vector of PFParticle
   *
   *  \param  inputPFParticles             input vector of PFParticles
   *  \param  idToPFParticleMap            input mapping between PFParticles and their IDs
   *  \param  idToDownstreamPFParticleMap  output mapping between of all downstream PFParticles and their IDs
   */
  void GetDownstreamPFParticles( const std::vector< art::Ptr< recob::PFParticle > > & inputPFParticles, const std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap, std::map< size_t, art::Ptr< recob::PFParticle > > & idToDownstreamPFParticleMap );

  /*!
   *  \brief  Collects all objects of type U associated to a given object of type T
   *
   *  \param  anObject         an input object of type T with which we want to collect associated objects of type U
   *  \param  associationTtoU  the general input association between objects of type U and T
   *  \param  associatedU      output vector of objects of type U associated with anObject
   */
  template < class T, class U >
  void CollectAssociated( const art::Ptr< T > & anObject, const art::FindManyP< U > & associationTtoU, std::vector< art::Ptr< U > > & associatedU );


  /*!
   *  \brief  Adds the given vector of objects to the event as a new collection
   *
   *  \param  event       the current event
   *  \param  collection  the collection to add
   */
  template < class T >
  void MakeCollection( art::Event & event, const std::vector< art::Ptr< T > > & collection );

  /*!
   *  \brief  Adds associations (type <A, B>) to the event between two collections (A and B) based on the associations between the input collections of the same types
   *
   *  Example use-case: 
   *
   *  Suppose you have an input collection of, say, SpacePoints (A) and PFParticles (B) from a single Pandora instance. 
   *  You also have associations between these collections (inputAssnAtoB)
   *
   *  You have already selected which objects in these collections you want to persist. 
   *  These are in your output collections (collectionA and collectionB). 
   *
   *  Now you need to produce output associations only between the objects that are in collectionA and collectionB.
   *
   *  This function will make an art::Assns< A, B >, and add it to the event.
   *
   *  \param  prod            the producer (usually *this)
   *  \param  event           the current event
   *  \param  collectionA     an collection of type A
   *  \param  collectionB     an collection of type B
   *  \param  inputAssnAtoB   input assocations between the full collections of type A and B
   */
  template < class PRODUCER, class A, class B >
  void MakeAssociation( PRODUCER const & prod, art::Event & event, const std::vector< art::Ptr< A > > & collectionA, const std::vector< art::Ptr< B > > & collectionB, const art::FindManyP< B > & inputAssnAtoB );

  /*!
   *  \brief  Produces a vector of Hits from an input handle
   *
   *  \param  hitHandle  input handle to a vector of hits
   *  \param  hitVect    output vector of art::Ptr to the hits
   */
  void GetHitVector( const art::Handle< std::vector< recob::Hit > > & hitHandle, std::vector< art::Ptr< recob::Hit > > & hitVect );

  // FHicL congifurable parameters
  std::string     fInputProducerLabel;           ///< Label for the Pandora instance that produced the collections we want to split up
  std::string     fHitProducerLabel;             ///< Label for the hit producer that was used as input to the Pandora instance specified
  bool            fShouldProduceNeutrinos;       ///< If we should produce collections related to neutrino top-level PFParticles
  bool            fShouldProduceCosmics;         ///< If we should produce collections related to cosmic (== non-neutrino) top-level PFParticles


  
  // Useful PDG code for readability
  enum Pdg {
    nue   = 12,
    numu  = 14,
    nutau = 16
  };

};

// ---------------------------------------------------------------------------------------

PandoraSplitting::PandoraSplitting(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  reconfigure(p);
  
  // Define which types of collections this the module produces
  produces< std::vector<recob::PFParticle> >();
  produces< std::vector<recob::SpacePoint> >();
  produces< std::vector<recob::Cluster> >();
  produces< std::vector<recob::Seed> >();
  produces< std::vector<recob::Vertex> >();
  produces< std::vector<recob::Track> >(); 
  produces< std::vector<recob::Shower> >();
  produces< std::vector<recob::PCAxis> >();

  produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
  produces< art::Assns<recob::PFParticle, recob::Cluster> >();
  produces< art::Assns<recob::PFParticle, recob::Seed> >();
  produces< art::Assns<recob::PFParticle, recob::Vertex> >();
  produces< art::Assns<recob::PFParticle, recob::Track> >();
  produces< art::Assns<recob::PFParticle, recob::Shower> >();
  produces< art::Assns<recob::PFParticle, recob::PCAxis> >();
  produces< art::Assns<recob::Track, recob::Hit> >();
  produces< art::Assns<recob::Shower, recob::Hit> >();
  produces< art::Assns<recob::Shower, recob::PCAxis> >();
  produces< art::Assns<recob::SpacePoint, recob::Hit> >();
  produces< art::Assns<recob::Cluster, recob::Hit> >();
  produces< art::Assns<recob::Seed, recob::Hit> >();
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::produce(art::Event & e)
{
  /// \todo  Error check for invalid handles etc.

  // ---------------------------------------------------------------------------------------
  // Get all of the input collections related to the Pandora instance specified
  // ---------------------------------------------------------------------------------------

  art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
  art::Handle< std::vector<recob::SpacePoint> > spacePointHandle;
  art::Handle< std::vector<recob::Cluster>    > clusterHandle;
  art::Handle< std::vector<recob::Seed>       > seedHandle;
  art::Handle< std::vector<recob::Vertex>     > vertexHandle;
  art::Handle< std::vector<recob::Track>      > trackHandle;
  art::Handle< std::vector<recob::Shower>     > showerHandle;
  art::Handle< std::vector<recob::PCAxis>     > pcAxisHandle;
  
  art::Handle< std::vector<recob::Hit>        > hitHandle;

  e.getByLabel(fInputProducerLabel, pfParticleHandle);
  e.getByLabel(fInputProducerLabel, spacePointHandle);
  e.getByLabel(fInputProducerLabel, clusterHandle);
  e.getByLabel(fInputProducerLabel, seedHandle);
  e.getByLabel(fInputProducerLabel, vertexHandle);
  e.getByLabel(fInputProducerLabel, trackHandle);
  e.getByLabel(fInputProducerLabel, showerHandle);
  e.getByLabel(fInputProducerLabel, pcAxisHandle);

  e.getByLabel(fHitProducerLabel  , hitHandle);
  std::vector< art::Ptr< recob::Hit > > hitVect;
  this->GetHitVector( hitHandle, hitVect );

  // Get the associations
  art::FindManyP< recob::SpacePoint > assnPFParticleSpacePoint( pfParticleHandle, e, fInputProducerLabel );
  art::FindManyP< recob::Cluster    > assnPFParticleCluster(    pfParticleHandle, e, fInputProducerLabel );
  art::FindManyP< recob::Seed       > assnPFParticleSeed(       pfParticleHandle, e, fInputProducerLabel );
  art::FindManyP< recob::Vertex     > assnPFParticleVertex(     pfParticleHandle, e, fInputProducerLabel );
  art::FindManyP< recob::Track      > assnPFParticleTrack(      pfParticleHandle, e, fInputProducerLabel );
  art::FindManyP< recob::Shower     > assnPFParticleShower(     pfParticleHandle, e, fInputProducerLabel );
  art::FindManyP< recob::PCAxis     > assnPFParticlePCAxis(     pfParticleHandle, e, fInputProducerLabel );

  art::FindManyP< recob::Hit        > assnSpacePointHit( spacePointHandle, e, fInputProducerLabel );
  art::FindManyP< recob::Hit        > assnClusterHit(    clusterHandle   , e, fInputProducerLabel );
  art::FindManyP< recob::Hit        > assnSeedHit(       seedHandle      , e, fInputProducerLabel );
  art::FindManyP< recob::Hit        > assnTrackHit(      trackHandle     , e, fInputProducerLabel );
  art::FindManyP< recob::Hit        > assnShowerHit(     showerHandle    , e, fInputProducerLabel );

  art::FindManyP< recob::PCAxis     > assnShowerPCAxis(  showerHandle    , e, fInputProducerLabel );

  // ---------------------------------------------------------------------------------------
  // Identify the PFParticles to persist
  // ---------------------------------------------------------------------------------------

  // Determine which top-level PFParticles we want to persist
  std::vector< art::Ptr< recob::PFParticle > > primaryPFParticles;
  this->GetPrimaryPFParticles( pfParticleHandle, primaryPFParticles );
  
  std::vector< art::Ptr< recob::PFParticle > > persistentPrimaryPFParticles;
  this->GetPFParticlesToPersist( primaryPFParticles, persistentPrimaryPFParticles );

  // Collect all daughter PFParticles to produce the final list of selected PFParticles to persist
  std::map< size_t, art::Ptr< recob::PFParticle > > idToPFParticleMap;
  this->GetIdToPFParticleMap( pfParticleHandle, idToPFParticleMap );

  std::map< size_t, art::Ptr< recob::PFParticle > > idToSelectedPFParticleMap;
  this->GetDownstreamPFParticles( persistentPrimaryPFParticles, idToPFParticleMap, idToSelectedPFParticleMap);

  // ---------------------------------------------------------------------------------------
  // Find all other objects related to the selected PFParticles
  // ---------------------------------------------------------------------------------------

  std::vector< art::Ptr< recob::PFParticle> > selectedParticles;
  std::vector< art::Ptr< recob::SpacePoint> > selectedSpacePoints;
  std::vector< art::Ptr< recob::Cluster>    > selectedClusters;
  std::vector< art::Ptr< recob::Seed>       > selectedSeeds;
  std::vector< art::Ptr< recob::Vertex>     > selectedVertices;
  std::vector< art::Ptr< recob::Track>      > selectedTracks;
  std::vector< art::Ptr< recob::Shower>     > selectedShowers; 
  std::vector< art::Ptr< recob::PCAxis>     > selectedPCAxes;

  for ( std::map< size_t, art::Ptr< recob::PFParticle > >::const_iterator selectedParticleIt=idToSelectedPFParticleMap.begin(); selectedParticleIt != idToSelectedPFParticleMap.end(); ++selectedParticleIt ) {
      art::Ptr< recob::PFParticle > part = selectedParticleIt->second;
    
      selectedParticles.push_back( part );

      // Collect all other associated objects
      this->CollectAssociated( part, assnPFParticleSpacePoint, selectedSpacePoints );
      this->CollectAssociated( part, assnPFParticleCluster   , selectedClusters    );
      this->CollectAssociated( part, assnPFParticleSeed      , selectedSeeds       );
      this->CollectAssociated( part, assnPFParticleVertex    , selectedVertices    );
      this->CollectAssociated( part, assnPFParticleTrack     , selectedTracks      );
      this->CollectAssociated( part, assnPFParticleShower    , selectedShowers     );
      this->CollectAssociated( part, assnPFParticlePCAxis    , selectedPCAxes      );
  }

  // ---------------------------------------------------------------------------------------
  // Output the selected collections
  // ---------------------------------------------------------------------------------------

  this->MakeCollection( e, selectedParticles   );
  this->MakeCollection( e, selectedSpacePoints );
  this->MakeCollection( e, selectedClusters    );
  this->MakeCollection( e, selectedSeeds       );
  this->MakeCollection( e, selectedVertices    );
  this->MakeCollection( e, selectedTracks      );
  this->MakeCollection( e, selectedShowers     );
  this->MakeCollection( e, selectedPCAxes      );

  this->MakeAssociation( *this, e, selectedParticles,   selectedSpacePoints, assnPFParticleSpacePoint );
  this->MakeAssociation( *this, e, selectedParticles,   selectedClusters   , assnPFParticleCluster    );
  this->MakeAssociation( *this, e, selectedParticles,   selectedSeeds      , assnPFParticleSeed       );
  this->MakeAssociation( *this, e, selectedParticles,   selectedVertices   , assnPFParticleVertex     );
  this->MakeAssociation( *this, e, selectedParticles,   selectedTracks     , assnPFParticleTrack      );
  this->MakeAssociation( *this, e, selectedParticles,   selectedShowers    , assnPFParticleShower     );
  this->MakeAssociation( *this, e, selectedParticles,   selectedPCAxes     , assnPFParticlePCAxis     );

  this->MakeAssociation( *this, e, selectedSpacePoints, hitVect            , assnSpacePointHit        );
  this->MakeAssociation( *this, e, selectedClusters   , hitVect            , assnClusterHit           );
  this->MakeAssociation( *this, e, selectedSeeds      , hitVect            , assnSeedHit              );
  this->MakeAssociation( *this, e, selectedTracks     , hitVect            , assnTrackHit             );
  this->MakeAssociation( *this, e, selectedShowers    , hitVect            , assnShowerHit            );

  this->MakeAssociation( *this, e, selectedShowers    , selectedPCAxes     , assnShowerPCAxis         );

}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::reconfigure(fhicl::ParameterSet const & p)
{
  fInputProducerLabel     = p.get<std::string>("InputProducerLabel");
  fHitProducerLabel       = p.get<std::string>("HitProducerLabel");
  fShouldProduceNeutrinos = p.get<bool>("ShouldProduceNeutrinos", false);
  fShouldProduceCosmics   = p.get<bool>("ShouldProduceCosmics"  , false);
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetPrimaryPFParticles( const art::Handle< std::vector< recob::PFParticle > > & allPFParticles, std::vector< art::Ptr< recob::PFParticle > > & primaryPFParticles )
{
  for( size_t pfPartIdx = 0; pfPartIdx != allPFParticles->size(); pfPartIdx++ ) {
      art::Ptr<recob::PFParticle> part(allPFParticles, pfPartIdx);
      if ( part->IsPrimary() ) 
        primaryPFParticles.push_back( part );
  }
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetPFParticlesToPersist( const std::vector< art::Ptr< recob::PFParticle > > & inputPFParticles, std::vector< art::Ptr< recob::PFParticle > > & persistentPFParticles )
{
  for ( art::Ptr< recob::PFParticle > part : inputPFParticles )
    if ( this->ShouldPersist( part ) )
      persistentPFParticles.push_back( part );
}

// ---------------------------------------------------------------------------------------

bool PandoraSplitting::ShouldPersist( const art::Ptr< recob::PFParticle > & part ) 
{
  unsigned int pdg = std::abs(part->PdgCode());
  bool isNeutrino = ( pdg == nue || pdg == numu || pdg == nutau );

  if (  isNeutrino && fShouldProduceNeutrinos ) return true;
  if ( !isNeutrino && fShouldProduceCosmics   ) return true;

  return false;
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetIdToPFParticleMap( const art::Handle< std::vector< recob::PFParticle > > & allPFParticles, std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap )
{
  for( size_t pfPartIdx = 0; pfPartIdx != allPFParticles->size(); pfPartIdx++ ) {
      art::Ptr<recob::PFParticle> part(allPFParticles, pfPartIdx);
      idToPFParticleMap.insert( std::map< size_t, art::Ptr< recob::PFParticle > >::value_type( part->Self(), part ) );
  }
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetDownstreamPFParticles( art::Ptr< recob::PFParticle > part, const std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap, std::map< size_t, art::Ptr< recob::PFParticle > > & idToDownstreamPFParticleMap )
{
    // Add part to the downstream map
    if ( idToDownstreamPFParticleMap.find( part->Self() ) == idToDownstreamPFParticleMap.end() )
        idToDownstreamPFParticleMap.insert( std::map< size_t, art::Ptr< recob::PFParticle > >::value_type( part->Self(), part ) );

    // And all of its daughters
    for ( size_t daughterId : part->Daughters() ) {

        std::map< size_t, art::Ptr< recob::PFParticle > >::const_iterator daughterIt = idToPFParticleMap.find( daughterId );

        /// \todo handle error properly
        if ( daughterIt == idToPFParticleMap.end() ) std::exit(1);

        this->GetDownstreamPFParticles( daughterIt->second, idToPFParticleMap, idToDownstreamPFParticleMap );
    }
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetDownstreamPFParticles( const std::vector< art::Ptr< recob::PFParticle > > & inputPFParticles, const std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap, std::map< size_t, art::Ptr< recob::PFParticle > > & idToDownstreamPFParticleMap )
{
  for ( art::Ptr< recob::PFParticle > part : inputPFParticles )
    this->GetDownstreamPFParticles( part, idToPFParticleMap, idToDownstreamPFParticleMap );
}

// ---------------------------------------------------------------------------------------

template < class T, class U >
void PandoraSplitting::CollectAssociated( const art::Ptr< T > & anObject, const art::FindManyP< U > & associationTtoU, std::vector< art::Ptr< U > > & associatedU )
{
  std::vector< art::Ptr< U > > associatedObjects = associationTtoU.at( anObject.key() );
  for ( art::Ptr< U > associatedObject : associatedObjects )
      associatedU.push_back( associatedObject );
}

// ---------------------------------------------------------------------------------------

template < class T >
void PandoraSplitting::MakeCollection( art::Event & event, const std::vector< art::Ptr< T > > & collection )
{
  std::unique_ptr< std::vector< T > > output( new std::vector< T > );

  for ( art::Ptr< T > object : collection )
    output->push_back( *object );
  
  event.put( std::move( output ) );
}

// ---------------------------------------------------------------------------------------

template < class PRODUCER, class A, class B >
void PandoraSplitting::MakeAssociation( PRODUCER const & prod, art::Event & event, const std::vector< art::Ptr< A > > & collectionA, const std::vector< art::Ptr< B > > & collectionB, const art::FindManyP< B > & inputAssnAtoB )
{
  
  // Setup the output association
  std::unique_ptr< art::Assns< A, B > > outputAssn( new art::Assns< A, B > );

  for ( art::Ptr< A > objectA : collectionA ) {
    for ( art::Ptr< B > associatedObjectB : inputAssnAtoB.at( objectA.key() ) ) {

        // Check that the associatedObjectB is in collectionB
        typename std::vector< art::Ptr< B > >::const_iterator associatedObjectIter = std::find( collectionB.begin(), collectionB.end(), associatedObjectB );
        if ( associatedObjectIter == collectionB.end() ) continue;
        
        util::CreateAssn( prod, event, associatedObjectB, objectA, *outputAssn );
    }
  }

  event.put( std::move( outputAssn ) );
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetHitVector( const art::Handle< std::vector< recob::Hit > > & hitHandle, std::vector< art::Ptr< recob::Hit > > & hitVect )
{
  for( size_t iHit = 0; iHit != hitHandle->size(); iHit++ ) {
      art::Ptr< recob::Hit > hit( hitHandle, iHit );
      hitVect.push_back( hit );
  }
}

// ---------------------------------------------------------------------------------------

DEFINE_ART_MODULE(PandoraSplitting)
