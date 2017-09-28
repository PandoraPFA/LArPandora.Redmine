/**
 *  @file   larpandora/LArPandoraEventBuilding/ShowerCreation_module.cc
 *
 *  @brief  module for shower creation
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PCAxis.h"

#include "larreco/Calorimetry/LinearEnergyAlg.h"

#include <memory>

namespace lar_pandora
{

class ShowerCreation : public art::EDProducer
{
public:
    explicit ShowerCreation(fhicl::ParameterSet const & p);

    ShowerCreation(ShowerCreation const &) = delete;
    ShowerCreation(ShowerCreation &&) = delete;
    ShowerCreation & operator = (ShowerCreation const &) = delete;
    ShowerCreation & operator = (ShowerCreation &&) = delete;

    void beginJob();
    void produce(art::Event & e) override;

private:
    /**
     *  @brief Build a recob::Shower object
     *
     *  @param pLArShowerPfo the object of the shower parameters filled in pandora
     *  @param totalEnergy calibrated energy for each cluster [GeV]
     */
    recob::Shower BuildShower(const lar_content::LArShowerPfo *const pLArShowerPfo, const std::vector<double> &totalEnergy) const;

    /**
     *  @brief Build a recob::PCAxis object
     *
     *  @param pLArShowerPfo the object of the shower parameters filled in pandora
     */
    recob::PCAxis BuildShowerPCA(const lar_content::LArShowerPfo *const pLArShowerPfo) const;

    const calo::LinearEnergyAlg    *m_pShowerEnergyAlg;         ///<
};

DEFINE_ART_MODULE(ShowerCreation)

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

ShowerCreation::ShowerCreation(fhicl::ParameterSet const & p) :
    m_showerEnergyAlg(nullptr)
{
    // prepare the optional cluster energy algorithm
    if (p.has_key("ShowerEnergy") && p.is_key_to_table("ShowerEnergy"))
    {
        m_showerEnergyAlg = std::make_unique<calo::LinearEnergyAlg>(p.get<fhicl::ParameterSet>("ShowerEnergy"));
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

void ShowerCreation::beginJob()
{
    // Print the configuration of the algorithm at the beginning of the job; the algorithm does not need to be set up for this.
    if (m_showerEnergyAlg)
    {
        mf::LogInfo log("ShowerCreation");
        log << "Energy shower settings: ";
        m_showerEnergyAlg->DumpConfiguration(log, "  ", "");
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerCreation::produce(art::Event & e)
{
    // we set up the algorithm on each new event, in case the services have changed:
    if (m_showerEnergyAlg)
    {
        m_showerEnergyAlg->setup(*(lar::providerFrom<detinfo::DetectorPropertiesService>()), *(lar::providerFrom<detinfo::DetectorClocksService>()),
            *(lar::providerFrom<geo::Geometry>()));
    }

    std::unique_ptr< std::vector<recob::Shower> > outputShowers( new std::vector<recob::Shower> );
    std::unique_ptr< std::vector<recob::PCAxis> > outputPCAxes( new std::vector<recob::PCAxis> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Shower> > outputParticlesToShowers( new art::Assns<recob::PFParticle, recob::Shower> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::PCAxis> > outputParticlesToPCAxes( new art::Assns<recob::PFParticle, recob::PCAxis> );
    std::unique_ptr< art::Assns<recob::Shower, recob::Hit> > outputShowersToHits( new art::Assns<recob::Shower, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Shower, recob::PCAxis> > outputShowersToPCAxes( new art::Assns<recob::Shower, recob::PCAxis> );

    // For all shower-like PFParticles
    int showerCounter(0);

    // TODO don't have access to outputParticles, or particleHitsFromSpacePoints

    for ()
    {
        // For now, call static helper functions in larpandoracontent

        outputShowers->emplace_back(LArPandoraOutput::BuildShower(pLArShowerPfo, showerE));
        outputPCAxes->emplace_back(LArPandoraOutput::BuildShowerPCA(pLArShowerPfo));
        outputShowers->back().set_id(outputShowers->size()); // 1-based sequence

        lar::PtrMaker<recob::Shower> makeShowerPtr(evt, *this));
        lar::PtrMaker<recob::PCAxis> makePCAxisPtr(evt, *this);
        lar::PtrMaker<recob::PFParticle> makePfoPtr(evt, *this);
        lar::PtrMaker<recob::Cluster> makeClusterPtr(evt, *this);

        // This bit will take some serious effort to update...
        std::vector<double> showerE;

        if (m_showerEnergyAlg)
        {
            // TODO Need to take implementation from Yun-Tse and Gianluca
        }

        outputParticlesToShowers->addSingle(makePfoPtr(outputParticles->size() - 1), makeShowerPtr(outputShowers->size() - 1));
        outputParticlesToPCAxes->addSingle(makePfoPtr(outputParticles->size() - 1), makePCAxisPtr(outputPCAxes->size() - 1));
        outputShowersToPCAxes->addSingle(makeShowerPtr(outputShowers->size() - 1), makePCAxisPtr(outputPCAxes->size() - 1));

        // Save associations between showers and hits
        util::CreateAssn(*(settings.m_pProducer), evt, *(outputShowers.get()), particleHitsFromSpacePoints, *(outputShowersToHits.get()));
    }
    
    mf::LogDebug("LArPandora") << "   Number of new showers: " << outputShowers->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new pcaxes:  " << outputPCAxes->size() << std::endl;

    evt.put(std::move(outputShowers));
    evt.put(std::move(outputPCAxes));
    evt.put(std::move(outputParticlesToShowers));
    evt.put(std::move(outputParticlesToPCAxes));
    evt.put(std::move(outputShowersToHits));
    evt.put(std::move(outputShowersToPCAxes));
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::Shower ShowerCreation::BuildShower(const lar_content::LArShowerPfo *const pLArShowerPfo, const std::vector<double>& totalEnergy) const
{
    const pandora::CartesianVector &showerLength(pLArShowerPfo->GetShowerLength());
    const pandora::CartesianVector &showerDirection(pLArShowerPfo->GetShowerDirection());
    const pandora::CartesianVector &showerVertex(pLArShowerPfo->GetShowerVertex());

    const float length(showerLength.GetX());
    const float openingAngle(pLArShowerPfo->GetShowerOpeningAngle());
    const TVector3 direction(showerDirection.GetX(), showerDirection.GetY(), showerDirection.GetZ());
    const TVector3 vertex(showerVertex.GetX(), showerVertex.GetY(), showerVertex.GetZ());

    // TODO
    const TVector3 directionErr;
    const TVector3 vertexErr;
    const std::vector<double> totalEnergyErr;
    const std::vector<double> dEdx;
    const std::vector<double> dEdxErr;
    const int bestplane(0);

    return recob::Shower(direction, directionErr, vertex, vertexErr, totalEnergy, totalEnergyErr, dEdx, dEdxErr, bestplane, util::kBogusI, length, openingAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::PCAxis ShowerCreation::BuildShowerPCA(const lar_content::LArShowerPfo *const pLArShowerPfo) const
{
    const pandora::CartesianVector &showerCentroid(pLArShowerPfo->GetShowerCentroid());
    const pandora::CartesianVector &showerDirection(pLArShowerPfo->GetShowerDirection());
    const pandora::CartesianVector &showerSecondaryVector(pLArShowerPfo->GetShowerSecondaryVector());
    const pandora::CartesianVector &showerTertiaryVector(pLArShowerPfo->GetShowerTertiaryVector());
    const pandora::CartesianVector &showerEigenValues(pLArShowerPfo->GetShowerEigenValues());

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
