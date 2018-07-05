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
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/ShowerEnergyAlg.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"
#include "Pandora/PandoraInternal.h"

#include <memory>

namespace lar_pandora
{

class LArPandoraShowerCreation : public art::EDProducer
{
    
public:
    
    typedef std::vector<art::Ptr<recob::Cluster> > ClusterVector;
    
    explicit LArPandoraShowerCreation(fhicl::ParameterSet const &pset);

    LArPandoraShowerCreation(LArPandoraShowerCreation const &) = delete;
    LArPandoraShowerCreation(LArPandoraShowerCreation &&) = delete;
    LArPandoraShowerCreation & operator = (LArPandoraShowerCreation const &) = delete;
    LArPandoraShowerCreation & operator = (LArPandoraShowerCreation &&) = delete;

    void produce(art::Event &evt) override;

private:
    /**
     * @brief 
     * @param larShowerPCA
     * @param vertexPosition
     * @param cartesianPointVector
     * @param clusterVector
     * @param clustersToHits
     * @return 
     */
    recob::Shower BuildShower(const lar_content::LArShowerPCA &larShowerPCA, const pandora::CartesianVector &vertexPosition, 
        const pandora::CartesianPointVector &cartesianPointVector, const ClusterVector &clusterVector, const ClustersToHits &clustersToHits);
    /**
     * @brief 
     * @param plane
     * @param clusterVector
     * @param clustersToHits
     * @param vertexPosition
     * @param showerStartDirection
     * @param dEdx
     * @param totalEnergy
     * @param bestPlane
     */
    void CalculateEnergyVariables(unsigned int plane, const ClusterVector &clusterVector, const ClustersToHits &clustersToHits, 
        const pandora::CartesianVector &vertexPosition, pandora::CartesianVector &showerStartDirection, double &dEdx, double &totalEnergy, int &bestPlane);
    
    /**
     * @brief 
     * @param hitView
     * @param tpc
     * @param cryo
     * @return 
     */
    float GetWireAngle(const geo::View_t hitView, const unsigned int tpc, const unsigned int cryo) const;

    /**
     * @brief 
     * @param wireAngle
     * @param initPosition
     * @return 
     */
    pandora::CartesianVector GetPositionOnPlane(const float wireAngle, const pandora::CartesianVector &initPosition) const;
    
    /**
     * @brief 
     * @param wireAngle
     * @param wirePitch
     * @param direction
     * @return 
     */
    float CalculateDX(const float wireAngle, const float wirePitch, const pandora::CartesianVector &direction) const;
    
    /**
     * @brief 
     * @param startHits
     * @param dx
     * @param plane
     * @return 
     */
    double CalculateDEdx(HitVector startHits, const float dx, unsigned int plane) const;

    /**
     *  @brief  Build a recob::PCAxis object
     *
     *  @param  larShowerPCA the lar shower pca parameters extracted from pandora
     */
    recob::PCAxis BuildPCAxis(const lar_content::LArShowerPCA &larShowerPCA) const;

    std::string     m_pfParticleLabel;              ///< The pf particle label
    bool            m_useAllParticles;              ///< Build a recob::Track for every recob::PFParticle

    shower::ShowerEnergyAlg  const m_showerEnergyAlg;///< Shower energy algorithm 
    calo::CalorimetryAlg const m_caloAlg;           ///< Calorimetry algorithm for dEdx calculation
    float           m_maxDistanceForDEdx;           ///< Maximum distance from vertex to consider for shower dEdx calculation 
    unsigned int    m_minPointsAtStart;             ///< Minimum number of points to consider the new PCA at start
    
};

DEFINE_ART_MODULE(LArPandoraShowerCreation)

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

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include <iostream>

namespace lar_pandora
{
    
    typedef std::vector<art::Ptr<recob::Cluster> > ClusterVector;
    
LArPandoraShowerCreation::LArPandoraShowerCreation(fhicl::ParameterSet const &pset) :
    m_pfParticleLabel(pset.get<std::string>("PFParticleLabel")),
    m_useAllParticles(pset.get<bool>("UseAllParticles", false)),
    m_showerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg")),
    m_caloAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    m_maxDistanceForDEdx(pset.get<float>("MaxDistanceForDEdx", 4.f)),
    m_minPointsAtStart(pset.get<unsigned int>("MinPointsAtStart", 3))
{
    produces< std::vector<recob::Shower> >();
    produces< std::vector<recob::PCAxis> >();
    produces< art::Assns<recob::PFParticle, recob::Shower> >();
    produces< art::Assns<recob::PFParticle, recob::PCAxis> >();
    produces< art::Assns<recob::Shower, recob::Hit> >();
    produces< art::Assns<recob::Shower, recob::PCAxis> >();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraShowerCreation::produce(art::Event &evt)
{
    std::unique_ptr< std::vector<recob::Shower> > outputShowers( new std::vector<recob::Shower> );
    std::unique_ptr< std::vector<recob::PCAxis> > outputPCAxes( new std::vector<recob::PCAxis> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Shower> > outputParticlesToShowers( new art::Assns<recob::PFParticle, recob::Shower> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::PCAxis> > outputParticlesToPCAxes( new art::Assns<recob::PFParticle, recob::PCAxis> );
    std::unique_ptr< art::Assns<recob::Shower, recob::Hit> > outputShowersToHits( new art::Assns<recob::Shower, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Shower, recob::PCAxis> > outputShowersToPCAxes( new art::Assns<recob::Shower, recob::PCAxis> );

    const art::PtrMaker<recob::Shower> makeShowerPtr(evt, *this);
    const art::PtrMaker<recob::PCAxis> makePCAxisPtr(evt, *this);

    // Organise inputs
    PFParticleVector pfParticleVector, extraPfParticleVector;
    PFParticlesToSpacePoints pfParticlesToSpacePoints;
    PFParticlesToClusters pfParticlesToClusters;
    LArPandoraHelper::CollectPFParticles(evt, m_pfParticleLabel, pfParticleVector, pfParticlesToSpacePoints);
    LArPandoraHelper::CollectPFParticles(evt, m_pfParticleLabel, extraPfParticleVector, pfParticlesToClusters);

    ClusterVector clusterVector;
    ClustersToHits clustersToHits;
    LArPandoraHelper::CollectClusters(evt, m_pfParticleLabel, clusterVector, clustersToHits);

    VertexVector vertexVector;
    PFParticlesToVertices pfParticlesToVertices;
    LArPandoraHelper::CollectVertices(evt, m_pfParticleLabel, vertexVector, pfParticlesToVertices);
    
    for (const art::Ptr<recob::PFParticle> pPFParticle : pfParticleVector)
    {
        // Select shower-like pfparticles
        if (!m_useAllParticles && !LArPandoraHelper::IsShower(pPFParticle))
            continue;

        // Obtain associated spacepoints
        PFParticlesToSpacePoints::const_iterator particleToSpacePointIter(pfParticlesToSpacePoints.find(pPFParticle));

        if (pfParticlesToSpacePoints.end() == particleToSpacePointIter)
        {
            mf::LogDebug("LArPandoraShowerCreation") << "No spacepoints associated to particle ";
            continue;
        }

        // Obtain associated clusters
        PFParticlesToClusters::const_iterator particleToClustersIter(pfParticlesToClusters.find(pPFParticle));

        if (pfParticlesToClusters.end() == particleToClustersIter)
        {
            mf::LogDebug("LArPandoraShowerCreation") << "No clusters associated to particle ";
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

        // Call pandora "fast" shower fitter
        try
        {
            // Ensure successful creation of all structures before placing results in output containers
            const lar_content::LArShowerPCA larShowerPCA(lar_content::LArPfoHelper::GetPrincipalComponents(cartesianPointVector, vertexPosition));
            const recob::Shower shower(LArPandoraShowerCreation::BuildShower(larShowerPCA, vertexPosition, cartesianPointVector, 
                                        particleToClustersIter->second, clustersToHits));
            const recob::PCAxis pcAxis(LArPandoraShowerCreation::BuildPCAxis(larShowerPCA));
            outputShowers->emplace_back(shower);
            outputPCAxes->emplace_back(pcAxis);
        }
        catch (const pandora::StatusCodeException &)
        {
            mf::LogDebug("LArPandoraShowerCreation") << "Unable to extract shower pca";
            continue;
        }

        // Output objects
        art::Ptr<recob::Shower> pShower(makeShowerPtr(outputShowers->size() - 1));
        art::Ptr<recob::PCAxis> pPCAxis(makePCAxisPtr(outputPCAxes->size() - 1));

        HitVector hitsInParticle;
        LArPandoraHelper::GetAssociatedHits(evt, m_pfParticleLabel, particleToClustersIter->second, hitsInParticle);

        // Output associations, after output objects are in place
        util::CreateAssn(*this, evt, pShower, pPFParticle, *(outputParticlesToShowers.get()));
        util::CreateAssn(*this, evt, pPCAxis, pPFParticle, *(outputParticlesToPCAxes.get()));
        util::CreateAssn(*this, evt, *(outputShowers.get()), hitsInParticle, *(outputShowersToHits.get()));
        util::CreateAssn(*this, evt, pPCAxis, pShower, *(outputShowersToPCAxes.get()));
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

recob::Shower LArPandoraShowerCreation::BuildShower(const lar_content::LArShowerPCA &larShowerPCA, const pandora::CartesianVector &vertexPosition, 
    const pandora::CartesianPointVector &cartesianPointVector, const ClusterVector &clusterVector, const ClustersToHits &clustersToHits) 
{
    const pandora::CartesianVector &showerLength(larShowerPCA.GetAxisLengths());
    const pandora::CartesianVector &showerDirection(larShowerPCA.GetPrimaryAxis());

    const float length(showerLength.GetX());
    const float openingAngle(larShowerPCA.GetPrimaryLength() > 0.f ? std::atan(larShowerPCA.GetSecondaryLength() / larShowerPCA.GetPrimaryLength()) : 0.f);
    const TVector3 direction(showerDirection.GetX(), showerDirection.GetY(), showerDirection.GetZ());
    const TVector3 vertex(vertexPosition.GetX(), vertexPosition.GetY(), vertexPosition.GetZ());

    // Find 3D points close to shower vertex and calculate PCA direction of the start of shower 
    pandora::CartesianPointVector cartesianPointVectorStart;
    for (const pandora::CartesianVector &pCartesianVector : cartesianPointVector) 
    {
        if ((pCartesianVector-vertexPosition).GetMagnitude() < m_maxDistanceForDEdx)
            cartesianPointVectorStart.emplace_back(pCartesianVector);
    }
 
    pandora::CartesianVector showerStartDirection = larShowerPCA.GetPrimaryAxis();    
    // Run the PCA analysis only on the initial points if a minimum number was found        
    if (cartesianPointVectorStart.size() >= m_minPointsAtStart)
    {
        try
        {
            const lar_content::LArShowerPCA larStartShowerPCA(lar_content::LArPfoHelper::GetPrincipalComponents(cartesianPointVectorStart, vertexPosition));
            showerStartDirection = larStartShowerPCA.GetPrimaryAxis();
        }
        catch (const pandora::StatusCodeException &){}
    }
    
    art::ServiceHandle<geo::Geometry> theGeometry;
    const unsigned int maxPlanes(theGeometry->MaxPlanes());
    std::vector<double> dEdx(maxPlanes), totalEnergy(maxPlanes);
    int bestPlane(-1);
    
    for (unsigned int plane = 0; plane < maxPlanes; ++plane)
        this->CalculateEnergyVariables(plane, clusterVector, clustersToHits, vertexPosition, showerStartDirection, 
                                        dEdx[plane], totalEnergy[plane], bestPlane);
                                        
    // TODO: There are no calculations of errors available yet
    const TVector3 directionErr;
    const TVector3 vertexErr;
    const std::vector<double> totalEnergyErr;
    const std::vector<double> dEdxErr;
    
    return recob::Shower(direction, directionErr, vertex, vertexErr, totalEnergy, totalEnergyErr, dEdx, dEdxErr, bestPlane, util::kBogusI, length, openingAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraShowerCreation::CalculateEnergyVariables(unsigned int plane, const ClusterVector &clusterVector, const ClustersToHits &clustersToHits, 
    const pandora::CartesianVector &vertexPosition, pandora::CartesianVector &showerStartDirection, double &dEdx, double &totalEnergy, int &bestPlane) 
{
    HitVector startHits, totalHits;
    art::ServiceHandle<geo::Geometry> theGeometry;
    auto const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    for (const art::Ptr<recob::Cluster> cluster : clusterVector)
    {
        if (plane != cluster->Plane().Plane)
            continue;

        ClustersToHits::const_iterator clusterToHitsIter(clustersToHits.find(cluster));
        if (clustersToHits.end() == clusterToHitsIter)
            continue;

        const HitVector &hitVector = clusterToHitsIter->second;
        if (hitVector.empty())
            continue;

        // ATTN: Using first hit of cluster, all hits in cluster should have the same view 
        // ATTN: wire angle is repeated here in case a shower crosses different drift volumes 
        HitVector::const_iterator firstHitIter = hitVector.begin();
        const art::Ptr<recob::Hit> firstHit = *firstHitIter;
        const geo::View_t hitView(LArPandoraGeometry::GetGlobalView(firstHit->WireID().Cryostat, firstHit->WireID().TPC, firstHit->View()));
        const float wireAngle(this->GetWireAngle(hitView,firstHit->WireID().TPC, firstHit->WireID().Cryostat));
        
        const pandora::CartesianVector vertexPosOnPlane(this->GetPositionOnPlane(wireAngle,vertexPosition));

        if (firstHit->SignalType() == geo::kCollection)
            bestPlane = plane; 

        for (HitVector::const_iterator hitIter = hitVector.begin(), hitIterEnd = hitVector.end(); hitIter != hitIterEnd; ++hitIter)
        {
            const art::Ptr<recob::Hit> hit = *hitIter;    
            totalHits.push_back(hit);
                
            const double hit_Time(hit->PeakTime());
            double hitXYZ[3];
            const double xpos_cm(theDetector->ConvertTicksToX(hit_Time, hit->WireID().Plane, hit->WireID().TPC, hit->WireID().Cryostat));
            theGeometry->Cryostat(hit->WireID().Cryostat).TPC(hit->WireID().TPC).Plane(hit->WireID().Plane).Wire(hit->WireID().Wire).GetCenter(hitXYZ);

            const pandora::CartesianVector hitPos(xpos_cm,hitXYZ[1],hitXYZ[2]);
            const pandora::CartesianVector hitPosOnPlane(this->GetPositionOnPlane(wireAngle,hitPos));//ATTN: Z here is U/V/W
            // Use 2D  distance (x,u/v/w) from hit to vertex
            if ((hitPosOnPlane-vertexPosOnPlane).GetMagnitude() < m_maxDistanceForDEdx)
                startHits.push_back(hit); 
        }
    }
 
    //totalEnergy = (!totalHits.empty() ? m_showerEnergyAlg.ShowerEnergy(totalHits, plane) : 0.);
    totalEnergy = 0.;
    
    // Calculate dx
    unsigned int vertexCryo(0), vertexTPC(0);
    const double vertexPos[3] = {vertexPosition.GetX(), vertexPosition.GetY(), vertexPosition.GetZ()};
    theGeometry->PositionToTPC(vertexPos, vertexTPC, vertexCryo);
    // this is ugly but necessary - alternative? 
    geo::PlaneGeo const& vertexPlane = theGeometry->Cryostat(vertexCryo).TPC(vertexTPC).Plane(plane);
    const float wireAngleVertex(this->GetWireAngle(vertexPlane.View(),vertexTPC, vertexCryo));
    const float wirePitchVertex(theGeometry->WirePitch(plane,vertexTPC,vertexCryo));

    const float dx(this->CalculateDX(wireAngleVertex, wirePitchVertex, showerStartDirection));
    
    // Calculate dEdx    
    const unsigned int nHits(startHits.size());
    dEdx = ((nHits > m_minPointsAtStart) ? this->CalculateDEdx(startHits, dx, plane) : 0.);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraShowerCreation::GetWireAngle(const geo::View_t view, const unsigned int tpc, const unsigned int cryo) const
{
    art::ServiceHandle<geo::Geometry> theGeometry;
    return (0.5f * M_PI - theGeometry->WireAngleToVertical(view, tpc, cryo));
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector LArPandoraShowerCreation::GetPositionOnPlane(const float wireAngle, const pandora::CartesianVector &initPosition) const
{
        pandora::CartesianVector position(initPosition.GetX(), 0.f, (initPosition.GetZ()*std::cos(wireAngle)-initPosition.GetY()*std::sin(wireAngle)));
        return position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraShowerCreation::CalculateDX(const float wireAngle, const float wirePitch, const pandora::CartesianVector &direction) const
{
    const pandora::CartesianVector unitDirection(direction.GetUnitVector());
    float projection(std::cos(wireAngle)*unitDirection.GetZ() - std::sin(wireAngle)*unitDirection.GetY());
    
    return ( (projection > std::numeric_limits<float>::epsilon()) ? (wirePitch/projection) : wirePitch );
}

//------------------------------------------------------------------------------------------------------------------------------------------
//ATTN: This returns double because the recob::Shower constructor uses doubles 
double LArPandoraShowerCreation::CalculateDEdx(HitVector startHits, const float dx, unsigned int plane) const
{
    float totalCharge(0.f), totalTime(0.f);
    std::vector<float> chargeVector;
    
    for (HitVector::const_iterator hitIter = startHits.begin(), hitIterEnd = startHits.end(); hitIter != hitIterEnd; ++hitIter)
    {
        const art::Ptr<recob::Hit> hit = *hitIter;
            
        chargeVector.push_back(hit->Integral());
        totalCharge += hit->Integral();
        totalTime += hit->PeakTime();
    }
        
    if (totalCharge < std::numeric_limits<float>::epsilon())
        return 0.; 
        
    if (chargeVector.empty())
        return 0.; 
        
    //TODO - allow other options: median, average, truncated mean... 
    float dQdx = TMath::Median(chargeVector.size(), &chargeVector[0])/dx;
    float averageTime(totalTime/(static_cast<float>(startHits.size())));
    
    return m_caloAlg.dEdx_AREA(dQdx, averageTime, plane);

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
