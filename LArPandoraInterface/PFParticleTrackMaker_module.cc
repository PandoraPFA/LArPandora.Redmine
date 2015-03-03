/**
 *  @file   larpandora/LArPandoraInterface/PFParticleTrackMaker.h
 *
 *  @brief  Producer module for making reco::Tracks from recob::PFParticles
 *
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  PFParticleTrackMaker class
 */
class PFParticleTrackMaker : public art::EDProducer
{
public:    

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    PFParticleTrackMaker(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~PFParticleTrackMaker();

    void beginJob();
    void endJob();
    void produce(art::Event &evt);

private:

    /**
     *  @brief Configure this module
     *  
     *  @param pset  the set of input parameters
     */
    void reconfigure(fhicl::ParameterSet const &pset);

    bool           m_cosmicMode;
    unsigned int   m_minSpacePoints;

    std::string    m_particleLabel;
    std::string    m_spacepointLabel;
};

DEFINE_ART_MODULE(PFParticleTrackMaker)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

// Framework includes
#include "cetlib/exception.h"

// Local includes
#include "LArPandoraCollector.h"
#include "LArPandoraHelper.h"

// LArSoft includes
#include "Utilities/AssociationUtil.h"
#include "Geometry/Geometry.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/Track.h"
#include "RecoBase/Hit.h"

namespace lar_pandora {

PFParticleTrackMaker::PFParticleTrackMaker(fhicl::ParameterSet const &pset) : art::EDProducer()
{
    produces< std::vector<recob::Track> >();    
    produces< art::Assns<recob::Track, recob::PFParticle> >();
    produces< art::Assns<recob::Track, recob::Hit> >();

    this->reconfigure(pset); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleTrackMaker::~PFParticleTrackMaker()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleTrackMaker::reconfigure(fhicl::ParameterSet const &pset)
{
    m_particleLabel   = pset.get<std::string>("PFParticleModuleLabel","pandora");
    m_spacepointLabel = pset.get<std::string>("SpacePointModuleLabel", "pandora");
    m_minSpacePoints  = pset.get<unsigned int>("MinSpacePoints",3); 
    m_cosmicMode      = pset.get<bool>("CosmicMode",true);
}

//------------------------------------------------------------------------------------------------------------------------------------------
  
void PFParticleTrackMaker::beginJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

void PFParticleTrackMaker::endJob()
{
}
   
//------------------------------------------------------------------------------------------------------------------------------------------  

void PFParticleTrackMaker::produce(art::Event &evt)
{
    mf::LogDebug("LArPandora") << " *** PFParticleTrackMaker::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "] *** " << std::endl;

    // Set up ART outputs
    // ==================
    std::unique_ptr< std::vector<recob::Track> > outputTracks( new std::vector<recob::Track> );
    std::unique_ptr< art::Assns<recob::Track, recob::PFParticle> > outputTracksToParticles( new art::Assns<recob::Track, recob::PFParticle> );
    std::unique_ptr< art::Assns<recob::Track, recob::Hit> > outputTracksToHits( new art::Assns<recob::Track, recob::Hit> );


    // Get particles and space points
    // ==============================
    PFParticleMap            particleMap;
    PFParticleVector         particleVector;
    PFParticlesToSpacePoints particlesToSpacePoints;
    SpacePointVector         spacePointVector;
    SpacePointsToHits        spacePointsToHits;
    
    LArPandoraCollector::CollectPFParticles(evt, m_particleLabel, particleVector, particlesToSpacePoints);
    LArPandoraCollector::CollectSpacePoints(evt, m_spacepointLabel, spacePointVector, spacePointsToHits);


    // Build mapping from particle to particle ID for parent/daughter navigation
    // =========================================================================
    for (PFParticleVector::const_iterator iter = particleVector.begin(), iterEnd = particleVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = *iter;
        particleMap[particle->Self()] = particle;
    }


    // Build tracks from space points for primary particles
    // ====================================================
    unsigned int trackCounter(0);

    for (PFParticlesToSpacePoints::const_iterator iter1 = particlesToSpacePoints.begin(), iterEnd1 = particlesToSpacePoints.end(); 
        iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> particle = iter1->first;
        const SpacePointVector &spacepoints = iter1->second;

        if (spacepoints.size() < m_minSpacePoints)
            continue;

        if (!LArPandoraCollector::IsTrack(particle))
            continue;

        if (!LArPandoraCollector::IsFinalState(particleMap, particle))
            continue;

	mf::LogDebug("LArPandora") << "   Building new track [" << trackCounter << "] (spacepoints=" << spacepoints.size() << ")" << std::endl; 

        PFParticleVector particles;
        particles.push_back(particle);

        HitVector hits;
        for (SpacePointVector::const_iterator iter2 = spacepoints.begin(), iterEnd2 = spacepoints.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::SpacePoint> spacepoint = *iter2;
            SpacePointsToHits::const_iterator iter3 = spacePointsToHits.find(spacepoint);
            if (spacePointsToHits.end() == iter3)
                throw cet::exception("LArPandora") << " PFParticleTrackMaker::produce --- Found space point without associated hit ";

            const art::Ptr<recob::Hit> hit = iter3->second;
            hits.push_back(hit);
        }

        try
	{
            recob::Track newTrack(LArPandoraHelper::BuildTrack(trackCounter++, spacepoints, m_cosmicMode));
            outputTracks->push_back(newTrack);

            util::CreateAssn(*this, evt, *(outputTracks.get()), particles, *(outputTracksToParticles.get()));
            util::CreateAssn(*this, evt, *(outputTracks.get()), hits, *(outputTracksToHits.get()));
	}
        catch (cet::exception &e)
	{
	    mf::LogWarning("LArPandora") << " PFParticleTrackMaker::produce --- Warning: Failed to build track " << std::endl;
	}
    }

    mf::LogDebug("LArPandora") << "   Number of new tracks: " << outputTracks->size() << std::endl;
    
    evt.put(std::move(outputTracks));
    evt.put(std::move(outputTracksToParticles));
    evt.put(std::move(outputTracksToHits));
}

} // namespace lar_pandora

