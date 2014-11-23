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

    std::string   m_particleLabel;
};

DEFINE_ART_MODULE(PFParticleTrackMaker)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

// Local includes
#include "LArPandoraCollector.h"
#include "LArPandoraHelper.h"
#include "PFParticleSeed.h"

// LArSoft includes
#include "Utilities/AssociationUtil.h"
#include "Geometry/Geometry.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/Track.h"

namespace lar_pandora {

PFParticleTrackMaker::PFParticleTrackMaker(fhicl::ParameterSet const &pset) : art::EDProducer()
{
    produces< std::vector<recob::Track> >();    
    produces< art::Assns<recob::Track, recob::PFParticle> >();

    this->reconfigure(pset); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleTrackMaker::~PFParticleTrackMaker()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleTrackMaker::reconfigure(fhicl::ParameterSet const &pset)
{
    m_particleLabel = pset.get<std::string>("PFParticleModuleLabel","pandora");
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


    // Get particles and space points
    // ==============================
    PFParticleMap            particleMap;
    PFParticleVector         particleVector;
    PFParticlesToSpacePoints particlesToSpacePoints;
    
    LArPandoraCollector::CollectPFParticles(evt, m_particleLabel, particleVector, particlesToSpacePoints);
    

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

    for (PFParticlesToSpacePoints::const_iterator iter = particlesToSpacePoints.begin(), iterEnd = particlesToSpacePoints.end(); 
        iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = iter->first;
        const SpacePointVector &spacepoints = iter->second;

        if (!LArPandoraCollector::IsTrack(particle))
	    continue;

        if (!LArPandoraCollector::IsFinalState(particleMap, particle))
	    continue;

        if (spacepoints.empty())
	    continue;

        mf::LogDebug("LArPandora") << "   Building new track [" << trackCounter << "] (spacepoints=" << spacepoints.size() << ")" << std::endl; 

        PFParticleVector particles;
        particles.push_back(particle);

        recob::Track newTrack(LArPandoraHelper::BuildTrack(trackCounter++, spacepoints));
        outputTracks->push_back(newTrack);

	util::CreateAssn(*this, evt, *(outputTracks.get()), particles, *(outputTracksToParticles.get()));
    }

    mf::LogDebug("LArPandora") << "   Number of new tracks: " << outputTracks->size() << std::endl;
    
    evt.put(std::move(outputTracks));
    evt.put(std::move(outputTracksToParticles));
}

} // namespace lar_pandora

