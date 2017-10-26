/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraTrackCreation_module.cc
 *
 *  @brief  module for lar pandora track creation
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larpandoracontent/LArObjects/LArPfoObjects.h"

#include <memory>

namespace lar_pandora
{

class LArPandoraTrackCreation : public art::EDProducer
{
public:
    explicit LArPandoraTrackCreation(fhicl::ParameterSet const &pset);

    LArPandoraTrackCreation(LArPandoraTrackCreation const &) = delete;
    LArPandoraTrackCreation(LArPandoraTrackCreation &&) = delete;
    LArPandoraTrackCreation & operator = (LArPandoraTrackCreation const &) = delete;
    LArPandoraTrackCreation & operator = (LArPandoraTrackCreation &&) = delete;

    void produce(art::Event &evt) override;

private:
    /**
     *  @brief Build a recob::Track object
     *
     *  @param id the id code for the track
     *  @param trackStateVector the vector of trajectory points for this track
     */
    recob::Track BuildTrack(const int id, const lar_content::LArTrackStateVector &trackStateVector) const;

    std::string     m_pfParticleLabel;              ///< The pf particle label
    unsigned int    m_minTrajectoryPoints;          ///< The minimum number of trajectory points
    unsigned int    m_slidingFitHalfWindow;         ///< The sliding fit half window
};

DEFINE_ART_MODULE(LArPandoraTrackCreation)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "art/Persistency/Common/PtrMaker.h"

#include "canvas/Utilities/InputTag.h"

#include "larcore/Geometry/Geometry.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <iostream>

namespace lar_pandora
{

LArPandoraTrackCreation::LArPandoraTrackCreation(fhicl::ParameterSet const &pset) :
    m_pfParticleLabel(pset.get<std::string>("PFParticleLabel")),
    m_minTrajectoryPoints(pset.get<unsigned int>("MinTrajectoryPoints", 2)),
    m_slidingFitHalfWindow(pset.get<unsigned int>("SlidingFitHalfWindow", 20))
{
    produces< std::vector<recob::Track> >();
    produces< art::Assns<recob::PFParticle, recob::Track> >();
    produces< art::Assns<recob::Track, recob::Hit> >();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraTrackCreation::produce(art::Event &evt)
{
    std::unique_ptr< std::vector<recob::Track> > outputTracks( new std::vector<recob::Track> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Track> > outputParticlesToTracks( new art::Assns<recob::PFParticle, recob::Track> );
    std::unique_ptr< art::Assns<recob::Track, recob::Hit> > outputTracksToHits( new art::Assns<recob::Track, recob::Hit> );

    art::ServiceHandle<geo::Geometry> theGeometry;
    const float wirePitchW((theGeometry->MaxPlanes() > 2) ? theGeometry->WirePitch(geo::kW) : 0.5f * (theGeometry->WirePitch(geo::kU) + theGeometry->WirePitch(geo::kV)));

    int trackCounter(0);
    const art::PtrMaker<recob::Track> makeTrackPtr(evt, *this);

    // Organise inputs
    PFParticleVector pfParticleVector;
    PFParticlesToSpacePoints pfParticlesToSpacePoints;
    LArPandoraHelper::CollectPFParticles(evt, m_pfParticleLabel, pfParticleVector, pfParticlesToSpacePoints);

    VertexVector vertexVector;
    PFParticlesToVertices pfParticlesToVertices;
    LArPandoraHelper::CollectVertices(evt, m_pfParticleLabel, vertexVector, pfParticlesToVertices);

    for (const art::Ptr<recob::PFParticle> pPFParticle : pfParticleVector)
    {
        // Only interested in track-like pfparticles
        if (!LArPandoraHelper::IsTrack(pPFParticle))
            continue;

        // Obtain associated spacepoints
        PFParticlesToSpacePoints::const_iterator particleToSpacePointIter(pfParticlesToSpacePoints.find(pPFParticle));

        if (pfParticlesToSpacePoints.end() == particleToSpacePointIter)
        {
            mf::LogDebug("LArPandoraTrackCreation") << "No spacepoints associated to particle ";
            continue;
        }

        // Obtain associated vertex
        PFParticlesToVertices::const_iterator particleToVertexIter(pfParticlesToVertices.find(pPFParticle));

        if ((pfParticlesToVertices.end() == particleToVertexIter) || (1 != particleToVertexIter->second.size()))
        {
            mf::LogDebug("LArPandoraTrackCreation") << "Unexpected number of vertices for particle ";
            continue;
        }

        // Copy information into expected pandora form
        pandora::CartesianPointVector cartesianPointVector;
        for (const art::Ptr<recob::SpacePoint> spacePoint : particleToSpacePointIter->second)
            cartesianPointVector.emplace_back(pandora::CartesianVector(spacePoint->XYZ()[0], spacePoint->XYZ()[1], spacePoint->XYZ()[2]));

        double vertexXYZ[3] = {0., 0., 0.};
        particleToVertexIter->second.front()->XYZ(vertexXYZ);
        const pandora::CartesianVector vertexPosition(vertexXYZ[0], vertexXYZ[1], vertexXYZ[2]);

        // Call pandora "fast" track fitter
        lar_content::LArTrackStateVector trackStateVector;

        try
        {
            lar_content::LArPfoHelper::GetSlidingFitTrajectory(cartesianPointVector, vertexPosition, m_slidingFitHalfWindow, wirePitchW, trackStateVector);
        }
        catch (const pandora::StatusCodeException &)
        {
            mf::LogDebug("LArPandoraTrackCreation") << "Unable to extract sliding fit trajectory";
            continue;
        }

        if (trackStateVector.size() < m_minTrajectoryPoints)
        {
            mf::LogDebug("LArPandoraTrackCreation") << "Insufficient input trajectory points to build track: " << trackStateVector.size();
            continue;
        }

        HitVector hitsInParticle;
        LArPandoraHelper::GetAssociatedHits(evt, m_pfParticleLabel, particleToSpacePointIter->second, hitsInParticle);

        // Output objects
        outputTracks->emplace_back(LArPandoraTrackCreation::BuildTrack(trackCounter++, trackStateVector));
        art::Ptr<recob::Track> pTrack(makeTrackPtr(outputTracks->size() - 1));

        // Output associations, after output objects are in place
        util::CreateAssn(*this, evt, pTrack, pPFParticle, *(outputParticlesToTracks.get()));
        util::CreateAssn(*this, evt, *(outputTracks.get()), hitsInParticle, *(outputTracksToHits.get()));
    }

    mf::LogDebug("LArPandoraTrackCreation") << "Number of new tracks: " << outputTracks->size() << std::endl;

    evt.put(std::move(outputTracks));
    evt.put(std::move(outputTracksToHits));
    evt.put(std::move(outputParticlesToTracks));
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::Track LArPandoraTrackCreation::BuildTrack(const int id, const lar_content::LArTrackStateVector &trackStateVector) const
{
    if (trackStateVector.empty())
        throw cet::exception("LArPandoraTrackCreation") << "BuildTrack - No input trajectory points provided ";

    std::vector<TVector3> xyz;
    std::vector<TVector3> pxpypz;

    std::vector<double> dummyMomentum(std::vector<double>(2, util::kBogusD));
    std::vector< std::vector<double> > dummyDQDX(std::vector< std::vector<double> >(0));

    for (const lar_content::LArTrackState &trackState : trackStateVector)
    {
        xyz.emplace_back(TVector3(trackState.GetPosition().GetX(), trackState.GetPosition().GetY(), trackState.GetPosition().GetZ()));
        pxpypz.emplace_back(TVector3(trackState.GetDirection().GetX(), trackState.GetDirection().GetY(), trackState.GetDirection().GetZ()));
    }

    return recob::Track(xyz, pxpypz, dummyDQDX, dummyMomentum, id);
}

} // namespace lar_pandora
