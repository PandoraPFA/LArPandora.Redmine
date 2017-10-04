/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraShowerCreation_module.cc
 *
 *  @brief  module for lar pandora shower creation
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PCAxis.h"

#include "larreco/Calorimetry/LinearEnergyAlg.h"

#include "larpandoracontent/LArObjects/LArPfoObjects.h"

#include <memory>

namespace lar_pandora
{

class LArPandoraShowerCreation : public art::EDProducer
{
public:
    explicit LArPandoraShowerCreation(fhicl::ParameterSet const &pset);

    LArPandoraShowerCreation(LArPandoraShowerCreation const &) = delete;
    LArPandoraShowerCreation(LArPandoraShowerCreation &&) = delete;
    LArPandoraShowerCreation & operator = (LArPandoraShowerCreation const &) = delete;
    LArPandoraShowerCreation & operator = (LArPandoraShowerCreation &&) = delete;

    void beginJob();
    void produce(art::Event &evt) override;

private:
    /**
     *  @brief  Build a recob::Shower object
     *
     *  @param  larShowerPCA the lar shower pca parameters extracted from pandora
     *  @param  vertexPosition the shower vertex position
     */
    recob::Shower BuildShower(const lar_content::LArShowerPCA &larShowerPCA, const pandora::CartesianVector &vertexPosition) const;

    /**
     *  @brief  Build a recob::PCAxis object
     *
     *  @param  larShowerPCA the lar shower pca parameters extracted from pandora
     */
    recob::PCAxis BuildPCAxis(const lar_content::LArShowerPCA &larShowerPCA) const;

    const calo::LinearEnergyAlg    *m_pShowerEnergyAlg;         ///<
    std::string                     m_pfParticleLabel;          ///< 
};

DEFINE_ART_MODULE(LArPandoraShowerCreation)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "canvas/Utilities/InputTag.h"

#include "larcore/Geometry/Geometry.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/PtrMaker.h"

#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <iostream>

namespace lar_pandora
{

LArPandoraShowerCreation::LArPandoraShowerCreation(fhicl::ParameterSet const &pset) :
    m_pShowerEnergyAlg(nullptr),
    m_pfParticleLabel(pset.get<std::string>("PFParticleLabel"))
{
    // prepare the optional cluster energy algorithm
    if (pset.has_key("ShowerEnergy") && pset.is_key_to_table("ShowerEnergy"))
    {
        //m_pShowerEnergyAlg = std::make_unique<calo::LinearEnergyAlg>(p.get<fhicl::ParameterSet>("ShowerEnergy"));
    }
    else mf::LogWarning("LArPandora") << "No shower energy calibration set up.";

    produces< std::vector<recob::Shower> >();
    produces< std::vector<recob::PCAxis> >();
    produces< art::Assns<recob::PFParticle, recob::Shower> >();
    produces< art::Assns<recob::PFParticle, recob::PCAxis> >();
    produces< art::Assns<recob::Shower, recob::Hit> >();
    produces< art::Assns<recob::Shower, recob::PCAxis> >();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraShowerCreation::beginJob()
{
    // Print the configuration of the algorithm at the beginning of the job; the algorithm does not need to be set up for this.
    if (m_pShowerEnergyAlg)
    {
        mf::LogInfo log("LArPandoraShowerCreation");
        log << "Energy shower settings: ";
        m_pShowerEnergyAlg->DumpConfiguration(log, "  ", "");
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraShowerCreation::produce(art::Event &evt)
{
    // we set up the algorithm on each new event, in case the services have changed:
    if (m_pShowerEnergyAlg)
    {
        //m_pShowerEnergyAlg->setup(*(lar::providerFrom<detinfo::DetectorPropertiesService>()),
        //    *(lar::providerFrom<detinfo::DetectorClocksService>()), *(lar::providerFrom<geo::Geometry>()));
    }

    std::unique_ptr< std::vector<recob::Shower> > outputShowers( new std::vector<recob::Shower> );
    std::unique_ptr< std::vector<recob::PCAxis> > outputPCAxes( new std::vector<recob::PCAxis> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Shower> > outputParticlesToShowers( new art::Assns<recob::PFParticle, recob::Shower> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::PCAxis> > outputParticlesToPCAxes( new art::Assns<recob::PFParticle, recob::PCAxis> );
    std::unique_ptr< art::Assns<recob::Shower, recob::SpacePoint> > outputShowersToSpacePoints( new art::Assns<recob::Shower, recob::SpacePoint> );
    std::unique_ptr< art::Assns<recob::Shower, recob::PCAxis> > outputShowersToPCAxes( new art::Assns<recob::Shower, recob::PCAxis> );

    const lar::PtrMaker<recob::Shower> makeShowerPtr(evt, *this);
    const lar::PtrMaker<recob::PCAxis> makePCAxisPtr(evt, *this);

    // Organise inputs
    PFParticleVector pfParticleVector;
    PFParticlesToSpacePoints pfParticlesToSpacePoints;
    LArPandoraHelper::CollectPFParticles(evt, m_pfParticleLabel, pfParticleVector, pfParticlesToSpacePoints);

    VertexVector vertexVector;
    PFParticlesToVertices pfParticlesToVertices;
    LArPandoraHelper::CollectVertices(evt, m_pfParticleLabel, vertexVector, pfParticlesToVertices);

    for (const art::Ptr<recob::PFParticle> pPFParticle : pfParticleVector)
    {
        // Only interested in shower-like pfparticles
        if (!LArPandoraHelper::IsShower(pPFParticle))
            continue;

        // Obtain associated spacepoints
        PFParticlesToSpacePoints::const_iterator particleToSpacePointIter(pfParticlesToSpacePoints.find(pPFParticle));

        if (pfParticlesToSpacePoints.end() == particleToSpacePointIter)
        {
            mf::LogDebug("LArPandoraShowerCreation") << "No spacepoints associated to particle ";
            continue;
        }

        // Obtain associated vertex
        PFParticlesToVertices::const_iterator particleToVertexIter(pfParticlesToVertices.find(pPFParticle));

        if ((pfParticlesToVertices.end() == particleToVertexIter) || (1 != particleToVertexIter->second.size()))
        {
            mf::LogDebug("LArPandoraShowerCreation") << "Unexpected number of vertices for particle ";
            continue;
        }

        // Copy information into expected pandora form
        pandora::CartesianPointVector cartesianPointVector;
        for (const art::Ptr<recob::SpacePoint> spacePoint : particleToSpacePointIter->second)
            cartesianPointVector.emplace_back(pandora::CartesianVector(spacePoint->XYZ()[0], spacePoint->XYZ()[1], spacePoint->XYZ()[2]));

        double vertexXYZ[3] = {0., 0., 0.};
        particleToVertexIter->second.front()->XYZ(vertexXYZ);
        const pandora::CartesianVector vertexPosition(vertexXYZ[0], vertexXYZ[1], vertexXYZ[2]);

        // Call pandora "fast" shower fitter - TODO check exception handling
        try
        {
            const lar_content::LArShowerPCA larShowerPCA(lar_content::LArPfoHelper::GetPrincipalComponents(cartesianPointVector, vertexPosition));
            outputShowers->emplace_back(LArPandoraShowerCreation::BuildShower(larShowerPCA, vertexPosition));
            outputPCAxes->emplace_back(LArPandoraShowerCreation::BuildPCAxis(larShowerPCA));
        }
        catch (const pandora::StatusCodeException &)
        {
            mf::LogDebug("LArPandoraShowerCreation") << "Unable to extract shower pca";
            continue;
        }

        // Output objects
        art::Ptr<recob::Shower> pShower(makeShowerPtr(outputShowers->size() - 1));
        art::Ptr<recob::PCAxis> pPCAxis(makePCAxisPtr(outputPCAxes->size() - 1));

        // Output associations, after output objects are in place
        util::CreateAssn(*this, evt, pShower, pPFParticle, *(outputParticlesToShowers.get()));
        util::CreateAssn(*this, evt, pPCAxis, pPFParticle, *(outputParticlesToPCAxes.get()));
        util::CreateAssn(*this, evt, *(outputShowers.get()), particleToSpacePointIter->second, *(outputShowersToSpacePoints.get()));
        util::CreateAssn(*this, evt, pPCAxis, pShower, *(outputShowersToPCAxes.get()));

        // TODO if (m_showerEnergyAlg) ...
    }
    
    mf::LogDebug("LArPandora") << "   Number of new showers: " << outputShowers->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new pcaxes:  " << outputPCAxes->size() << std::endl;

    evt.put(std::move(outputShowers));
    evt.put(std::move(outputPCAxes));
    evt.put(std::move(outputParticlesToShowers));
    evt.put(std::move(outputParticlesToPCAxes));
    evt.put(std::move(outputShowersToSpacePoints));
    evt.put(std::move(outputShowersToPCAxes));
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::Shower LArPandoraShowerCreation::BuildShower(const lar_content::LArShowerPCA &larShowerPCA, const pandora::CartesianVector &vertexPosition) const
{
    const pandora::CartesianVector &showerLength(larShowerPCA.GetAxisLengths());
    const pandora::CartesianVector &showerDirection(larShowerPCA.GetPrimaryAxis());

    const float length(showerLength.GetX());
    const float openingAngle(larShowerPCA.GetPrimaryLength() > 0.f ? std::atan(larShowerPCA.GetSecondaryLength() / larShowerPCA.GetPrimaryLength()) : 0.f);
    const TVector3 direction(showerDirection.GetX(), showerDirection.GetY(), showerDirection.GetZ());
    const TVector3 vertex(vertexPosition.GetX(), vertexPosition.GetY(), vertexPosition.GetZ());

    // TODO
    const TVector3 directionErr;
    const TVector3 vertexErr;
    const std::vector<double> totalEnergyErr;
    const std::vector<double> dEdx;
    const std::vector<double> dEdxErr;
    const std::vector<double> totalEnergy;
    const int bestplane(0);

    return recob::Shower(direction, directionErr, vertex, vertexErr, totalEnergy, totalEnergyErr, dEdx, dEdxErr, bestplane, util::kBogusI, length, openingAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::PCAxis LArPandoraShowerCreation::BuildPCAxis(const lar_content::LArShowerPCA &larShowerPCA) const
{
    const pandora::CartesianVector &showerCentroid(larShowerPCA.GetCentroid());
    const pandora::CartesianVector &showerDirection(larShowerPCA.GetPrimaryAxis());
    const pandora::CartesianVector &showerSecondaryVector(larShowerPCA.GetSecondaryAxis());
    const pandora::CartesianVector &showerTertiaryVector(larShowerPCA.GetTertiaryAxis());
    const pandora::CartesianVector &showerEigenValues(larShowerPCA.GetEigenValues());

    const bool svdOK(true); ///< SVD Decomposition was successful
    const double eigenValues[3] = {showerEigenValues.GetX(), showerEigenValues.GetY(), showerEigenValues.GetZ()}; ///< Eigen values from SVD decomposition
    const double avePosition[3] = {showerCentroid.GetX(), showerCentroid.GetY(), showerCentroid.GetZ()}; ///< Average position of hits fed to PCA

    std::vector< std::vector<double> > eigenVecs = { /// The three principle axes
        { showerDirection.GetX(), showerDirection.GetY(), showerDirection.GetZ() },
        { showerSecondaryVector.GetX(), showerSecondaryVector.GetY(), showerSecondaryVector.GetZ() },
        { showerTertiaryVector.GetX(), showerTertiaryVector.GetY(), showerTertiaryVector.GetZ() }
    };

    // TODO
    const int numHitsUsed(100); ///< Number of hits in the decomposition, not yet ready
    const double aveHitDoca(0.); ///< Average doca of hits used in PCA, not ready yet
    const size_t iD(util::kBogusI); ///< Axis ID, not ready yet

    return recob::PCAxis(svdOK, numHitsUsed, eigenValues, eigenVecs, avePosition, aveHitDoca, iD);
}

} // namespace lar_pandora
