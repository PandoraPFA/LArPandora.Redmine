/**
 *  @file   larpandora/LArPandoraEventBuilding/TrackCreation_module.cc
 *
 *  @brief  module for track creation
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Track.h"

#include <memory>

namespace lar_pandora
{

class TrackCreation : public art::EDProducer
{
public:
    explicit TrackCreation(fhicl::ParameterSet const & p);

    TrackCreation(TrackCreation const &) = delete;
    TrackCreation(TrackCreation &&) = delete;
    TrackCreation & operator = (TrackCreation const &) = delete;
    TrackCreation & operator = (TrackCreation &&) = delete;

    void produce(art::Event & e) override;

private:
    /**
     *  @brief Build a recob::Track object
     *
     *  @param id the id code for the track
     *  @param pTrackStateVector the vector of trajectory points for this track
     */
    recob::Track BuildTrack(const int id, const lar_content::LArTrackStateVector *const pTrackStateVector, const IdToHitMap &idToHitMap) const;

    /**
     *  @brief Build a reco::Seed object
     *
     *  @param trackState the trajectory point for this seed
     */
    recob::Seed BuildSeed(const lar_content::LArTrackState &trackState) const;

    /**
     *  @brief Calculate dQ/dL based on an input hit
     *
     *  @param hit the input ART hit
     *  @param trackPosition the local track position
     *  @param trackDirection the local track direction
     *
     *  @return dQ/dL
     */
    double CalculatedQdL(const art::Ptr<recob::Hit> hit, const TVector3 &trackPosition, const TVector3 &trackDirection) const;

    unsigned int    m_minTrajectoryPoints;          ///<
};

DEFINE_ART_MODULE(TrackCreation)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "canvas/Utilities/InputTag.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>

namespace lar_pandora
{

TrackCreation::TrackCreation(fhicl::ParameterSet const & p) :
    m_minTrajectoryPoints(2)
{
    produces< std::vector<recob::Track> >();
    produces< std::vector<recob::Seed> >();
    produces< art::Assns<recob::PFParticle, recob::Track> >();
    produces< art::Assns<recob::PFParticle, recob::Seed> >();
    produces< art::Assns<recob::Track, recob::Hit> >();
    produces< art::Assns<recob::Seed, recob::Hit> >();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackCreation::produce(art::Event & e)
{
    std::unique_ptr< std::vector<recob::Track> > outputTracks( new std::vector<recob::Track> );
    std::unique_ptr< std::vector<recob::Seed> > outputSeeds( new std::vector<recob::Seed> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Track> > outputParticlesToTracks( new art::Assns<recob::PFParticle, recob::Track> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Seed> > outputParticlesToSeeds( new art::Assns<recob::PFParticle, recob::Seed> );
    std::unique_ptr< art::Assns<recob::Track, recob::Hit> > outputTracksToHits( new art::Assns<recob::Track, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Seed, recob::Hit> > outputSeedsToHits( new art::Assns<recob::Seed, recob::Hit> );

    // For all track-like PFParticles
    int trackCounter(0);

    // TODO don't have access to outputParticles

    for ()
    {
        // For now, call static helper functions in larpandoracontent
        
        const lar_content::LArTrackStateVector &trackStateVector = pLArTrackPfo->m_trackStateVector;

        if (trackStateVector.size() < m_minTrajectoryPoints)
        {
            mf::LogDebug("LArPandora") << " LArPandoraOutput::BuildTrack --- Insufficient input trajectory points to build track ";
            continue;
        }

        for (const lar_content::LArTrackState &nextPoint : trackStateVector)
        {
            const art::Ptr<recob::Hit> hit = LArPandoraOutput::GetHit(idToHitMap, nextPoint.GetCaloHit()); // TODO

            HitVector seedHits;
            seedHits.push_back(hit);

            outputSeeds->emplace_back(LArPandoraOutput::BuildSeed(nextPoint));

            util::CreateAssn(*(settings.m_pProducer), evt, *(outputSeeds.get()), seedHits, *(outputSeedsToHits.get()));
            util::CreateAssn(*(settings.m_pProducer), evt, *(outputParticles.get()), *(outputSeeds.get()), *(outputParticlesToSeeds.get()),
                outputSeeds->size() - 1, outputSeeds->size());
        }

        outputTracks->emplace_back(LArPandoraOutput::BuildTrack(trackCounter++, &trackStateVector, idToHitMap));

        util::CreateAssn(*(settings.m_pProducer), evt, *(outputTracks.get()), particleHitsFromSpacePoints, *(outputTracksToHits.get()));
        util::CreateAssn(*(settings.m_pProducer), evt, *(outputParticles.get()), *(outputTracks.get()), *(outputParticlesToTracks.get()), outputTracks->size() - 1, outputTracks->size());
    }

    mf::LogDebug("LArPandora") << "   Number of new tracks: " << outputTracks->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new seeds: " << outputSeeds->size() << std::endl;

    evt.put(std::move(outputTracks));
    evt.put(std::move(outputParticlesToTracks));
    evt.put(std::move(outputTracksToHits));
    evt.put(std::move(outputSeeds));
    evt.put(std::move(outputParticlesToSeeds));
    evt.put(std::move(outputSeedsToHits));
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::Track TrackCreation::BuildTrack(const int id, const lar_content::LArTrackStateVector *const pTrackStateVector, const IdToHitMap &idToHitMap) const
{
    if (pTrackStateVector->empty())
        throw cet::exception("TrackCreation") << " TrackCreation::BuildTrack --- No input trajectory points were provided ";

    // Fill list of track properties
    std::vector<TVector3>               xyz;
    std::vector<TVector3>               pxpypz;
    std::vector< std::vector<double> >  dQdx(3);
    std::vector<double>                 momentum = std::vector<double>(2, util::kBogusD);

    // Loop over trajectory points
    for (const lar_content::LArTrackState &nextPoint : *pTrackStateVector)
    {
        const art::Ptr<recob::Hit> hit = LArPandoraOutput::GetHit(idToHitMap, nextPoint.GetCaloHit());
        const geo::View_t hit_View(hit->View());

        const TVector3 trackPosition(nextPoint.GetPosition().GetX(), nextPoint.GetPosition().GetY(), nextPoint.GetPosition().GetZ());
        const TVector3 trackDirection(nextPoint.GetDirection().GetX(), nextPoint.GetDirection().GetY(), nextPoint.GetDirection().GetZ());
        const double trackdQdx(this->CalculatedQdL(hit, trackPosition, trackDirection));

        const double dQdxU((geo::kU == hit_View) ? trackdQdx : 0.0);
        const double dQdxV((geo::kV == hit_View) ? trackdQdx : 0.0);
        const double dQdxW((geo::kW == hit_View) ? trackdQdx : 0.0);

        xyz.push_back(trackPosition);
        pxpypz.push_back(trackDirection);

        dQdx.at(geo::kU).push_back(dQdxU); dQdx.at(geo::kV).push_back(dQdxV); dQdx.at(geo::kW).push_back(dQdxW);
    }

    // Return a new recob::Track object (of the Bezier variety)
    return recob::Track(xyz, pxpypz, dQdx, momentum, id);
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::Seed TrackCreation::BuildSeed(const lar_content::LArTrackState &trackState) const
{
    double pos[3]     = { trackState.GetPosition().GetX(), trackState.GetPosition().GetY(), trackState.GetPosition().GetZ() };
    double dir[3]     = { trackState.GetDirection().GetX(), trackState.GetDirection().GetY(), trackState.GetDirection().GetZ() };
    double posErr[3]  = { 0.0, 0.0, 0.0 };  // TODO: Fill in errors
    double dirErr[3]  = { 0.0, 0.0, 0.0 };  // TODO: Fill in errors

    return recob::Seed(pos, dir, posErr, dirErr);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double TrackCreation::CalculatedQdL(const art::Ptr<recob::Hit> hit, const TVector3&, const TVector3 &trackDirection) const
{
    // Extract wire pitch and wire direction from geometry
    art::ServiceHandle<geo::Geometry> theGeometry;
    const geo::WireID wireID(hit->WireID());
    const geo::WireGeo &wireGeo = theGeometry->Cryostat(wireID.Cryostat).TPC(wireID.TPC).Plane(wireID.Plane).Wire(wireID.Wire);
    const TVector3 driftDirection(1.0, 0.0, 0.0);
    const TVector3 wireDirection(wireGeo.Direction());
    const TVector3 wireAxis(wireDirection.Cross(driftDirection));
    const double wirePitch(theGeometry->WirePitch(hit->View()));

    // Calculate dQ/dL by projecting track direction onto wire
    const float cosTheta(std::fabs(wireAxis.Dot(trackDirection)));
    const float inverse_dL(cosTheta / wirePitch);
    const float dL((inverse_dL > std::numeric_limits<float>::epsilon()) ? (1.0 / inverse_dL) : std::numeric_limits<float>::max());
    const float dQ(hit->Integral());

    return (dQ/dL);
}

} // namespace lar_pandora
