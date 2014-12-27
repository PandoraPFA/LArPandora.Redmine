// Framework includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/cpu_timer.h"
#include "cetlib/exception.h"

// LArSoft includes
#include "Utilities/AssociationUtil.h"
#include "SimulationBase/MCTruth.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/PFParticle.h"

// Local includes 
#include "PFParticleStitcher.h"

// System includes
#include <iostream>
#include <limits>

namespace lar_pandora {

PFParticleStitcher::PFParticleStitcher(fhicl::ParameterSet const &pset) : art::EDProducer()
{
    produces< std::vector<recob::PFParticle> >();    
    produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
    produces< art::Assns<recob::PFParticle, recob::Cluster> >();
    
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleStitcher::~PFParticleStitcher()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::reconfigure(fhicl::ParameterSet const &pset)
{
    m_enableMonitoring = pset.get<bool>("EnableMonitoring", false);
    m_particleLabel = pset.get<std::string>("PFParticleModuleLabel", "pandora");

    m_useXcoordinate = pset.get<bool>("UseXCoordinate", true);
    m_minCosRelativeAngle = pset.get<float>("MaxCosRelativeAngle", 0.966);
    m_maxLongitudinalDisplacementX = pset.get<float>("MaxLongitudinalDisplacementX", 15.f);
    m_maxTransverseDisplacement = pset.get<float>("MaxTransverseDisplacement", 5.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::beginJob()
{  
    if (m_enableMonitoring)
        this->InitializeMonitoring();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::endJob()
{   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::produce(art::Event &evt)
{ 
    // Run/Event Numbers
    // =================
    m_run   = evt.run();
    m_event = evt.id().event();

    mf::LogDebug("LArPandora") << " *** PFParticleStitcher::produce(...)  [Run=" << m_run << ", Event=" << m_event << "] *** " << std::endl;


    // Collect reconstructed PFParticles
    // =================================
    PFParticleVector         inputParticles;
    PFParticleVector         inputParticles2;
    PFParticleVector         parentParticles;

    PFParticlesToClusters    particlesToClusters;
    PFParticlesToSpacePoints particlesToSpacePoints;
    PFParticleSeedMap        particleSeedMap;
    PFParticleVolumeMap      particleVolumeMap;
    PFParticleMap            particleMap;

    SpacePointVector         inputSpacePoints;
    SpacePointsToHits        spacePointsToHits;

    LArPandoraCollector::CollectPFParticles(evt, m_particleLabel, inputParticles, particlesToClusters);
    LArPandoraCollector::CollectPFParticles(evt, m_particleLabel, inputParticles2, particlesToSpacePoints);
    LArPandoraCollector::CollectSpacePoints(evt, m_particleLabel, inputSpacePoints, spacePointsToHits); // (Assume that these were made 
                                                                                                        //  along with the PFParticles)
    
    // Build mapping from particle to particle ID for parent/daughter navigation
    for (PFParticleVector::const_iterator iter = inputParticles.begin(), iterEnd = inputParticles.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = *iter;
        particleMap[particle->Self()] = particle;
    }
   
    // Select final-state particles
    for (PFParticleVector::const_iterator iter = inputParticles.begin(), iterEnd = inputParticles.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = *iter;

        if (LArPandoraCollector::IsFinalState(particleMap, particle))
            parentParticles.push_back(particle);
    }

    // Create particle seeds, and associate them with drift volumes
    for (PFParticlesToSpacePoints::const_iterator iter =  particlesToSpacePoints.begin(), iterEnd = particlesToSpacePoints.end(); 
        iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = iter->first;
        const SpacePointVector &spacepoints = iter->second;

        if(spacepoints.empty())
            continue;

        particleSeedMap.insert(PFParticleSeedMap::value_type(particle, PFParticleSeed(spacepoints)));
        particleVolumeMap.insert(PFParticleVolumeMap::value_type(particle, this->GetVolumeID(spacepoints, spacePointsToHits)));
    }

    // Match Parent PFParticles across drift volumes
    // =============================================
    ParticleAssociationMatrix particleAssociationMatrix;
    PFParticleMergeMap particleMatches, particleMerges;

    this->CreateParticleMatches(parentParticles, particleVolumeMap, particleSeedMap, particleAssociationMatrix);
    this->SelectParticleMatches(particleMap, particleAssociationMatrix, particleMatches);
    this->SelectParticleMerges(parentParticles, particleMatches, particleMerges);

    if (m_enableMonitoring)
        this->WriteParticleMatches(particleMatches, particleSeedMap);

    // Merge PFParticles
    // =================
    this->ProduceArtOutput(evt, particleMap, particleMerges, particlesToClusters, particlesToSpacePoints);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::InitializeMonitoring()
{
    art::ServiceHandle<art::TFileService> tfs;
    m_pRecoTree = tfs->make<TTree>("monitoring", "LAr Reco");
    m_pRecoTree->Branch("run", &m_run, "run/I");
    m_pRecoTree->Branch("event", &m_event, "event/I");
    m_pRecoTree->Branch("particle1", &m_particle1, "particle1/I");
    m_pRecoTree->Branch("particle2", &m_particle2, "particle2/I");
    m_pRecoTree->Branch("cosRelativeAngle", &m_cosRelativeAngle, "cosRelativeAngle/F");
    m_pRecoTree->Branch("transverseDisplacement", &m_transverseDisplacement, "transverseDisplacement/F");
    m_pRecoTree->Branch("longitudinalDisplacement", &m_longitudinalDisplacement, "longitudinalDisplacement/F");
    m_pRecoTree->Branch("longitudinalDisplacementCut", &m_longitudinalDisplacementCut, "longitudinalDisplacementCut/F");
    m_pRecoTree->Branch("deltaX", &m_deltaX, "deltaX/F");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::ProduceArtOutput(art::Event &evt, const PFParticleMap &particleMap, const PFParticleMergeMap &particleMerges, 
    const PFParticlesToClusters &particlesToClusters, const PFParticlesToSpacePoints &particlesToSpacePoints)
{
    mf::LogDebug("LArPandora") << " **** PFParticleStitcher::ProduceArtOutput(...) **** " << std::endl;

    // Set up ART outputs
    // ==================
    std::unique_ptr< std::vector<recob::PFParticle> > outputParticles( new std::vector<recob::PFParticle> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint> > outputParticlesToSpacePoints( new art::Assns<recob::PFParticle, recob::SpacePoint> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster> >    outputParticlesToClusters( new art::Assns<recob::PFParticle, recob::Cluster> );
   

    // Select final-state particles
    // ============================
    PFParticleVector  parentParticles;
    PFParticleVector  daughterParticles;
    
    for (PFParticleMap::const_iterator iter = particleMap.begin(), iterEnd = particleMap.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = iter->second;

        if (LArPandoraCollector::IsFinalState(particleMap, particle))
        {
            parentParticles.push_back(particle);
        }
        else if (LArPandoraCollector::IsNeutrino(particle))
        {
            // TODO: Recover neutrinos
        }
        else
        {
            daughterParticles.push_back(particle);
        }
    }

 
    // Create a new set of particle ID codes
    // =====================================
    typedef std::map<art::Ptr<recob::PFParticle>, size_t> PFParticleIDMap;

    PFParticleIDMap particleIdMap;
    PFParticleIDMap mergeIdMap;
    
    size_t particleCounter(0);
 
    for (PFParticleMergeMap::const_iterator iter1 = particleMerges.begin(), iterEnd1 = particleMerges.end(); iter1 != iterEnd1; ++iter1)
    {
        ++particleCounter;
        const art::Ptr<recob::PFParticle> particle = iter1->first;
        particleIdMap[particle] = particleCounter;

        for (PFParticleList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::PFParticle> mergedParticle = *iter2;

            if (particle->Self() == mergedParticle->Self())
                continue;

            mergeIdMap[mergedParticle] = particleCounter;
        }
    }

    for (PFParticleVector::const_iterator iter1 = daughterParticles.begin(), iterEnd1 = daughterParticles.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> particle = *iter1;
        particleIdMap[particle] = (++particleCounter);
    }


    // Create new particles
    // ====================
    for (PFParticleIDMap::const_iterator iter1 = particleIdMap.begin(), iterEnd1 = particleIdMap.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> oldParticle = iter1->first;

        // Get the old particles
        PFParticleList oldParticleList;
        PFParticleMergeMap::const_iterator iter2A = particleMerges.find(oldParticle);

        if (particleMerges.end() == iter2A)
        {
            oldParticleList.insert(oldParticle);
        }
        else
        {
            oldParticleList.insert(iter2A->second.begin(), iter2A->second.end());
        }

        // Get the new particle codes
        const int newPdgCode(oldParticle->PdgCode());
        const size_t newSelfCode(iter1->second);

        // Get the new parent particle
        size_t newParentCode(recob::PFParticle::kPFParticlePrimary);

        if (!oldParticle->IsPrimary())
        {    
            PFParticleMap::const_iterator iter3A = particleMap.find(oldParticle->Parent());
            if (particleMap.end() == iter3A)
                throw cet::exception("LArPandora") << " PFParticleStitcher::ProduceArtOutput --- No trace of particle in particle maps";

            const art::Ptr<recob::PFParticle> oldParentParticle = iter3A->second;

            if (particleIdMap.find(oldParentParticle) != particleIdMap.end())
            {
                newParentCode = particleIdMap[oldParentParticle];
            }
            else if (mergeIdMap.find(oldParentParticle) != mergeIdMap.end())
            {
                newParentCode = mergeIdMap[oldParentParticle];
            }
        }

        // Get the new daughter particles
        std::vector<size_t> newDaughterCodes;

        const std::vector<size_t> &oldDaughterCodes = oldParticle->Daughters();

        for (std::vector<size_t>::const_iterator iter2B = oldDaughterCodes.begin(), iterEnd2B = oldDaughterCodes.end(); iter2B != iterEnd2B; ++iter2B)
        {
            int oldDaughterCode = static_cast<int>(*iter2B);

            PFParticleMap::const_iterator iter3B = particleMap.find(oldDaughterCode);
            if (particleMap.end() == iter3B)
                throw cet::exception("LArPandora") << " PFParticleStitcher::ProduceArtOutput --- No trace of particle in particle maps";

            const art::Ptr<recob::PFParticle> oldDaughterParticle = iter3B->second;

            if (particleIdMap.find(oldDaughterParticle) != particleIdMap.end())
            {
                newDaughterCodes.push_back(particleIdMap[oldDaughterParticle]);
            }
            else if (mergeIdMap.find(oldDaughterParticle) != mergeIdMap.end())
            {
                newDaughterCodes.push_back(mergeIdMap[oldDaughterParticle]);
            }
        }

        // Build new particle
        recob::PFParticle newParticle(newPdgCode, newSelfCode, newParentCode, newDaughterCodes);
        outputParticles->push_back(newParticle);

        // Build associations between new particle and old clusters, old space points
        ClusterVector clusterVector;
        SpacePointVector spacePointVector;

        for (PFParticleList::const_iterator iter4 = oldParticleList.begin(), iterEnd4 = oldParticleList.end(); iter4 != iterEnd4; ++iter4)
        {
            const art::Ptr<recob::PFParticle> particle = *iter4;
          
            PFParticlesToClusters::const_iterator iter5 = particlesToClusters.find(particle);
            if (particlesToClusters.end() == iter5)
                continue;

            for (ClusterVector::const_iterator iter6 = iter5->second.begin(), iterEnd6 = iter5->second.end(); iter6 != iterEnd6; ++iter6)
            {
                const art::Ptr<recob::Cluster> cluster = *iter6;
                clusterVector.push_back(cluster);
            }

            PFParticlesToSpacePoints::const_iterator iter7 = particlesToSpacePoints.find(particle);
            if (particlesToSpacePoints.end() == iter7)
                continue;

            for (SpacePointVector::const_iterator iter8 = iter7->second.begin(), iterEnd8 = iter7->second.end(); iter8 != iterEnd8; ++iter8)
            {
                const art::Ptr<recob::SpacePoint> spacepoint = *iter8;
                spacePointVector.push_back(spacepoint);
            }
        }

        util::CreateAssn(*this, evt, *(outputParticles.get()), clusterVector, *(outputParticlesToClusters.get()));
        util::CreateAssn(*this, evt, *(outputParticles.get()), spacePointVector, *(outputParticlesToSpacePoints.get()));
    }

    mf::LogDebug("LArPandora") << "   Number of new particles: " << outputParticles->size() << std::endl;
    
    evt.put(std::move(outputParticles));
    evt.put(std::move(outputParticlesToSpacePoints));
    evt.put(std::move(outputParticlesToClusters));
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int PFParticleStitcher::GetVolumeID(const SpacePointVector &spacePoints, const SpacePointsToHits &spacePointsToHits) const
{
    for (SpacePointVector::const_iterator iter1 = spacePoints.begin(), iterEnd1 = spacePoints.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::SpacePoint> spacepoint = *iter1;

        SpacePointsToHits::const_iterator iter2 = spacePointsToHits.find(spacepoint);
        if (spacePointsToHits.end() == iter2)
            throw cet::exception("LArPandora") << " PFParticleStitcher::GetVolumeID --- Found a space point without an associated hit";

        const art::Ptr<recob::Hit> hit = iter2->second;
        const geo::WireID hit_WireID(hit->WireID());

        return this->GetVolumeID(hit_WireID.Cryostat, hit_WireID.TPC);
    }

    throw cet::exception("LArPandora") << " PFParticleStitcher::GetVolumeID --- No volume ID for this collection of hits";
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::WriteParticleMatches(const PFParticleMergeMap &particleMergeMap, const PFParticleSeedMap &particleSeedMap)
{
    for (PFParticleMergeMap::const_iterator iter1 = particleMergeMap.begin(), iterEnd1 = particleMergeMap.end(); iter1 != iterEnd1; ++iter1)    
    {
        const art::Ptr<recob::PFParticle> particle1 = iter1->first;

        for (PFParticleList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::PFParticle> particle2 = *iter2;

            m_particle1 = particle1->Self();
            m_particle2 = particle2->Self();

            PFParticleSeedMap::const_iterator iter3 = particleSeedMap.find(particle1);
            PFParticleSeedMap::const_iterator iter4 = particleSeedMap.find(particle2);

            if (particleSeedMap.end() == iter3 || particleSeedMap.end() == iter4)
                throw cet::exception("LArPandora") << " PFParticleStitcher::WriteParticleMatches --- No trace of particle seed in seed maps";

            const PFParticleSeed seed1(iter3->second);
            const PFParticleSeed seed2(iter4->second);

            // Get closest pair of vertices
            ParticleAssociation::VertexType vertexType1(ParticleAssociation::UNDEFINED);
            ParticleAssociation::VertexType vertexType2(ParticleAssociation::UNDEFINED);

            this->GetClosestVertices(seed1, seed2, vertexType1, vertexType2);

            // Get vertex and direction at closest pair of vertices
            const pandora::CartesianVector vtx1((ParticleAssociation::INNER == vertexType1) ? seed1.GetInnerPosition() :  seed1.GetOuterPosition());
            const pandora::CartesianVector dir1((ParticleAssociation::INNER == vertexType1) ? seed1.GetInnerDirection() : seed1.GetOuterDirection());

            const pandora::CartesianVector vtx2((ParticleAssociation::INNER == vertexType2) ? seed2.GetInnerPosition() :  seed2.GetOuterPosition());
            const pandora::CartesianVector dir2((ParticleAssociation::INNER == vertexType2) ? seed2.GetInnerDirection() : seed2.GetOuterDirection());

            if (m_useXcoordinate)
            {
                this->GetImpactParameters3D(vtx1, dir1, vtx2, m_longitudinalDisplacement, m_transverseDisplacement);
                m_longitudinalDisplacementCut = this->GetLongitudinalDisplacementCut3D(dir1);
            }
            else
            {
                this->GetImpactParameters2D(vtx1, dir1, vtx2, m_longitudinalDisplacement, m_transverseDisplacement);
                m_longitudinalDisplacementCut = this->GetLongitudinalDisplacementCut2D(dir1);
            }
            
            m_cosRelativeAngle = -dir1.GetDotProduct(dir2);
            m_deltaX = this->GetDeltaX(vtx1, dir1, vtx2);

            m_pRecoTree->Fill();
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::CreateParticleMatches(const PFParticleVector &particleVector, const PFParticleVolumeMap &particleVolumeMap,
    const PFParticleSeedMap &particleSeedMap, ParticleAssociationMatrix &particleAssociationMatrix) const
{
    for (PFParticleVector::const_iterator pIter1 = particleVector.begin(), pIterEnd1 = particleVector.end(); pIter1 != pIterEnd1; ++pIter1)
    {
        const art::Ptr<recob::PFParticle> particle1 = *pIter1;

        for (PFParticleVector::const_iterator pIter2 = pIter1, pIterEnd2 = particleVector.end(); pIter2 != pIterEnd2; ++pIter2)
        {
            const art::Ptr<recob::PFParticle> particle2 = *pIter2;

            if (particle1->Self() == particle2->Self())
                continue;

            if (particle1->PdgCode() != particle2->PdgCode())
                continue;

            PFParticleVolumeMap::const_iterator vIter1 = particleVolumeMap.find(particle1);
            PFParticleVolumeMap::const_iterator vIter2 = particleVolumeMap.find(particle2);

            if (particleVolumeMap.end() == vIter1 || particleVolumeMap.end() == vIter2)
                continue;

            const unsigned int volume1(vIter1->second);
            const unsigned int volume2(vIter2->second);

            if (volume1 == volume2)
                continue;

            this->CreateParticleMatches(particle1, particle2, particleSeedMap, particleAssociationMatrix);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::CreateParticleMatches(const art::Ptr<recob::PFParticle> particle1, const art::Ptr<recob::PFParticle> particle2,
    const PFParticleSeedMap &particleSeedMap, ParticleAssociationMatrix &particleAssociationMatrix) const
{
    PFParticleSeedMap::const_iterator iter1 = particleSeedMap.find(particle1);
    PFParticleSeedMap::const_iterator iter2 = particleSeedMap.find(particle2);

    if (particleSeedMap.end() == iter1 || particleSeedMap.end() == iter2)
        return;

    const PFParticleSeed seed1(iter1->second);
    const PFParticleSeed seed2(iter2->second);

    // Get closest pair of vertices
    ParticleAssociation::VertexType vertexType1(ParticleAssociation::UNDEFINED);
    ParticleAssociation::VertexType vertexType2(ParticleAssociation::UNDEFINED);

    try
    {
        this->GetClosestVertices(seed1, seed2, vertexType1, vertexType2);
    }
    catch (pandora::StatusCodeException& )
    {
        return;
    }
    
    // Get vertex and direction at closest pair of vertices
    const pandora::CartesianVector vtx1((ParticleAssociation::INNER == vertexType1) ? seed1.GetInnerPosition() :  seed1.GetOuterPosition());
    const pandora::CartesianVector dir1((ParticleAssociation::INNER == vertexType1) ? seed1.GetInnerDirection() : seed1.GetOuterDirection());

    const pandora::CartesianVector vtx2((ParticleAssociation::INNER == vertexType2) ? seed2.GetInnerPosition() :  seed2.GetOuterPosition());
    const pandora::CartesianVector dir2((ParticleAssociation::INNER == vertexType2) ? seed2.GetInnerDirection() : seed2.GetOuterDirection());

    // Relative Angles
    const float cosRelativeAngle(-dir1.GetDotProduct(dir2));

    if (cosRelativeAngle < m_minCosRelativeAngle)
        return;

    // Impact Parameters
    float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);

    if (m_useXcoordinate)
    {
        this->GetImpactParameters3D(vtx1, dir1, vtx2, rL1, rT1);
        this->GetImpactParameters3D(vtx2, dir2, vtx1, rL2, rT2);
    }
    else
    {
        try
        {
            this->GetImpactParameters2D(vtx1, dir1, vtx2, rL1, rT1);
            this->GetImpactParameters2D(vtx2, dir2, vtx1, rL2, rT2);
        }
        catch (pandora::StatusCodeException& )
        {
            return;
        }
    }

    // Selection Cuts on longitudinal Impact Parameters
    float rCutL1(-std::numeric_limits<float>::max()), rCutL2(-std::numeric_limits<float>::max());

    try
    {
        if (m_useXcoordinate)
        {
            rCutL1 = this->GetLongitudinalDisplacementCut3D(dir1);
            rCutL2 = this->GetLongitudinalDisplacementCut3D(dir2);
        }
        else
        {
            rCutL1 = this->GetLongitudinalDisplacementCut2D(dir1);
            rCutL2 = this->GetLongitudinalDisplacementCut2D(dir2);
        }
    }
    catch (pandora::StatusCodeException& )
    {
    }

    if (rL1 < -1.f || rL1 > rCutL1 || rL2 < -1.f || rL2 > rCutL2 ||
        rT1 > m_maxTransverseDisplacement || rT2 > m_maxTransverseDisplacement)
        return;

    // Store Association 
    const float particleLength1((seed1.GetInnerPosition() - seed1.GetOuterPosition()).GetMagnitudeSquared());
    const float particleLength2((seed2.GetInnerPosition() - seed2.GetOuterPosition()).GetMagnitudeSquared());

    mf::LogDebug("LArPandora")  << "   *** ParticleStitcher::MatchParticles(...) *** " << std::endl
                                << "     id1=" << particle1->Self() << ", id2=" << particle2->Self() << std::endl
                                << "     cosTheta=" << -dir1.GetDotProduct(dir2) << std::endl
                                << "     rL1=" << rL1 << ", rT1=" << rT1 << ", rL2=" << rL2 << ", rT2=" << rT2 << std::endl
                                << "     Length1=" << std::sqrt(particleLength1) << " Length2=" << std::sqrt(particleLength2) << std::endl;

    (void) particleAssociationMatrix[particle1].insert(ParticleAssociationMap::value_type(particle2, 
            ParticleAssociation(vertexType1, vertexType2, particleLength2)));
    (void) particleAssociationMatrix[particle2].insert(ParticleAssociationMap::value_type(particle1, 
            ParticleAssociation(vertexType2, vertexType1, particleLength1)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::SelectParticleMatches(const PFParticleMap &particleMap, const ParticleAssociationMatrix &particleAssociationMatrix,
    PFParticleMergeMap &particleMergeMap) const
{
    // First step: find the best associations A -> X and B -> Y
    // ========================================================
    ParticleAssociationMatrix intermediateAssociationMatrix;

    for (ParticleAssociationMatrix::const_iterator iter1 = particleAssociationMatrix.begin(), iterEnd1 = particleAssociationMatrix.end(); 
        iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> parentParticle(iter1->first);
        const ParticleAssociationMap &particleAssociationMap(iter1->second);

        int bestParticleInner(-1);
        ParticleAssociation bestAssociationInner(ParticleAssociation::UNDEFINED, ParticleAssociation::UNDEFINED, 0.f);

        int bestParticleOuter(-1);
        ParticleAssociation bestAssociationOuter(ParticleAssociation::UNDEFINED, ParticleAssociation::UNDEFINED, 0.f);

        for (ParticleAssociationMap::const_iterator iter2 = particleAssociationMap.begin(), iterEnd2 = particleAssociationMap.end(); 
            iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::PFParticle> daughterParticle(iter2->first);
            const ParticleAssociation &particleAssociation(iter2->second);

            // Inner associations
            if (particleAssociation.GetParent() == ParticleAssociation::INNER)
            {
                if (particleAssociation.GetFigureOfMerit() > bestAssociationInner.GetFigureOfMerit())
                {
                    bestAssociationInner = particleAssociation;
                    bestParticleInner = daughterParticle->Self();
                }
            }

            // Outer associations
            if (particleAssociation.GetParent() == ParticleAssociation::OUTER)
            {
                if (particleAssociation.GetFigureOfMerit() > bestAssociationOuter.GetFigureOfMerit())
                {
                    bestAssociationOuter = particleAssociation;
                    bestParticleOuter = daughterParticle->Self();
                }
            }
        }

        if (bestAssociationInner.GetFigureOfMerit() > std::numeric_limits<float>::epsilon())
        {
            PFParticleMap::const_iterator pIter = particleMap.find(bestParticleInner);
            if (particleMap.end() == pIter)
                throw cet::exception("LArPandora") << " PFParticleStitcher::SelectParticleMatches --- No trace of particle in particle maps";

            const art::Ptr<recob::PFParticle> daughterParticle = pIter->second;
            (void) intermediateAssociationMatrix[parentParticle].insert(ParticleAssociationMap::value_type(daughterParticle, bestAssociationInner));
        }

        if (bestAssociationOuter.GetFigureOfMerit() > std::numeric_limits<float>::epsilon())
        {
            PFParticleMap::const_iterator pIter = particleMap.find(bestParticleOuter);
            if (particleMap.end() == pIter)
                throw cet::exception("LArPandora") << " PFParticleStitcher::SelectParticleMatches --- No trace of particle in particle maps";

            const art::Ptr<recob::PFParticle> daughterParticle = pIter->second;
            (void) intermediateAssociationMatrix[parentParticle].insert(ParticleAssociationMap::value_type(daughterParticle, bestAssociationOuter));
        }
    }

    // Second step: make the merge if A -> X and B -> Y is in fact A -> B and B -> A
    // =============================================================================
    for (ParticleAssociationMatrix::const_iterator iter3 = intermediateAssociationMatrix.begin(), iterEnd3 = intermediateAssociationMatrix.end(); iter3 != iterEnd3; ++iter3)
    {
        const art::Ptr<recob::PFParticle> parentParticle(iter3->first);
        const ParticleAssociationMap &parentAssociationMap(iter3->second);

        for (ParticleAssociationMap::const_iterator iter4 = parentAssociationMap.begin(), iterEnd4 = parentAssociationMap.end(); iter4 != iterEnd4; ++iter4)
        {
            const art::Ptr<recob::PFParticle> daughterParticle(iter4->first);
            const ParticleAssociation &parentToDaughterAssociation(iter4->second);

            ParticleAssociationMatrix::const_iterator iter5 = intermediateAssociationMatrix.find(daughterParticle);

            if (intermediateAssociationMatrix.end() == iter5)
                continue;

            const ParticleAssociationMap &daughterAssociationMap(iter5->second);

            ParticleAssociationMap::const_iterator iter6 = daughterAssociationMap.find(parentParticle);

            if (daughterAssociationMap.end() == iter6)
                continue;

            const ParticleAssociation &daughterToParentAssociation(iter6->second);

            if (parentToDaughterAssociation.GetParent() == daughterToParentAssociation.GetDaughter() &&
                parentToDaughterAssociation.GetDaughter() == daughterToParentAssociation.GetParent())
            {
                particleMergeMap[parentParticle].insert(daughterParticle);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::SelectParticleMerges(const PFParticleVector &particleVector, const PFParticleMergeMap &inputMergeMap,
    PFParticleMergeMap &outputMergeMap) const
{
    PFParticleList vetoList;

    for (PFParticleVector::const_iterator iter1 = particleVector.begin(), iterEnd1 = particleVector.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> seedParticle = *iter1;

        if (vetoList.count(seedParticle))
            continue;

        PFParticleList mergeList;
        this->CollectAssociatedParticles(seedParticle, seedParticle, inputMergeMap, vetoList, mergeList);

        vetoList.insert(seedParticle);
        outputMergeMap[seedParticle].insert(seedParticle);

        for (PFParticleList::const_iterator iter2 = mergeList.begin(), iterEnd2 = mergeList.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::PFParticle> associatedParticle = *iter2;

            if (vetoList.count(associatedParticle))
                throw cet::exception("LArPandora") << " PFParticleStitcher::SelectParticleMerges --- This particle has been counted twice!";

            vetoList.insert(associatedParticle);
            outputMergeMap[seedParticle].insert(associatedParticle);

            mf::LogDebug("LArPandora")  << "   *** ParticleStitcher::SelectMerges(...) *** " << std::endl
                                        << "     Seed=" << seedParticle->Self() << ", Associated=" << associatedParticle->Self() << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::CollectAssociatedParticles(art::Ptr<recob::PFParticle> seedParticle, art::Ptr<recob::PFParticle> currentParticle, 
    const PFParticleMergeMap &particleMergeMap, const PFParticleList &vetoList, PFParticleList &associatedParticleList) const
{
    if (vetoList.count(currentParticle))
        return;

    PFParticleMergeMap::const_iterator iter1 = particleMergeMap.find(currentParticle);

    if (particleMergeMap.end() == iter1)
        return;

    for (PFParticleList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
    {
        const art::Ptr<recob::PFParticle> associatedParticle = *iter2;

        if (associatedParticle->Self() == seedParticle->Self())
            continue;

        if (!associatedParticleList.insert(associatedParticle).second)
            continue;

        this->CollectAssociatedParticles(seedParticle, associatedParticle, particleMergeMap, vetoList, associatedParticleList);
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::GetClosestVertices(const PFParticleSeed& seed1, const PFParticleSeed &seed2, 
    ParticleAssociation::VertexType &vertexType1, ParticleAssociation::VertexType &vertexType2) const
{
    for (unsigned int inner1 = 0; inner1 < 2; ++inner1)
    {
        vertexType1 = ((0 == inner1) ? ParticleAssociation::INNER : ParticleAssociation::OUTER);
        const pandora::CartesianVector vtx1((ParticleAssociation::INNER == vertexType1) ? seed1.GetInnerPosition()  : seed1.GetOuterPosition());
        const pandora::CartesianVector dir1((ParticleAssociation::INNER == vertexType1) ? seed1.GetInnerDirection() : seed1.GetOuterDirection());
        const pandora::CartesianVector end1((ParticleAssociation::INNER == vertexType1) ? seed1.GetOuterPosition()  : seed1.GetInnerPosition());

        for (unsigned int inner2 = 0; inner2 < 2; ++inner2)
        {
            vertexType2 = ((0 == inner2) ? ParticleAssociation::INNER : ParticleAssociation::OUTER);
            const pandora::CartesianVector vtx2((ParticleAssociation::INNER == vertexType2) ? seed2.GetInnerPosition()  : seed2.GetOuterPosition());
            const pandora::CartesianVector dir2((ParticleAssociation::INNER == vertexType2) ? seed2.GetInnerDirection() : seed2.GetOuterDirection());
            const pandora::CartesianVector end2((ParticleAssociation::INNER == vertexType2) ? seed2.GetOuterPosition()  : seed2.GetInnerPosition());

            float vtxT(0.f), vtxL(0.f), endT(0.f), endL(0.f);

            // Project seed1 onto seed2
            if (m_useXcoordinate)
            {
                this->GetImpactParameters3D(vtx1, dir1, vtx2, vtxL, vtxT);
                this->GetImpactParameters3D(vtx1, dir1, end2, endL, endT);
            }
            else
                {
                try
                {
                    this->GetImpactParameters2D(vtx1, dir1, vtx2, vtxL, vtxT);
                    this->GetImpactParameters2D(vtx1, dir1, end2, endL, endT);
                }
                catch (pandora::StatusCodeException& )
                {
                    continue;
                }
            }

            if (vtxL > endL || endL < 0.f)
                continue;
            
            // Project seed2 onto seed1
            if (m_useXcoordinate)
            {
                this->GetImpactParameters3D(vtx2, dir2, vtx1, vtxL, vtxT);
                this->GetImpactParameters3D(vtx2, dir2, end1, endL, endT);
            }
            else
            {
                try
                {
                    this->GetImpactParameters2D(vtx2, dir2, vtx1, vtxL, vtxT);
                    this->GetImpactParameters2D(vtx2, dir2, end1, endL, endT);
                }
                catch (pandora::StatusCodeException& )
                {
                    continue;
                }
            }

            if (vtxL > endL || endL < 0.f)
                continue;

            // These are the closest vertices [return]
            return;
        }
    }

    // Can't find closest vertices [bail out] (ATTN: throw Pandora exception here, as code is based on pandora objects) 
    throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::GetImpactParameters3D(const pandora::CartesianVector &initialPosition, const pandora::CartesianVector &initialDirection, 
    const pandora::CartesianVector &targetPosition, float &longitudinal, float &transverse) const
{
    // sign convention for longitudinal distance:
    // -positive value means initial position is downstream of target position
    transverse = initialDirection.GetCrossProduct(targetPosition-initialPosition).GetMagnitude();
    longitudinal = -initialDirection.GetDotProduct(targetPosition-initialPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleStitcher::GetImpactParameters2D(const pandora::CartesianVector &initialPosition, const pandora::CartesianVector &initialDirection, 
    const pandora::CartesianVector &targetPosition, float &longitudinal, float &transverse) const
{
    if (std::fabs(initialDirection.GetX()) > 1.f - std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    const pandora::CartesianVector initialPosition2D(0.f, initialPosition.GetY(), initialPosition.GetZ());
    const pandora::CartesianVector initialDirection2D(0.f, initialDirection.GetY(), initialDirection.GetZ());
    const pandora::CartesianVector targetPosition2D(0.f, targetPosition.GetY(), targetPosition.GetZ());
   
    this->GetImpactParameters3D(initialPosition2D, initialDirection2D.GetUnitVector(), targetPosition2D, longitudinal, transverse);  
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

float PFParticleStitcher::GetLongitudinalDisplacementCut2D(const pandora::CartesianVector &direction) const
{
    if (std::fabs(direction.GetX()) < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    if (std::fabs(direction.GetX()) > 1.f - std::numeric_limits<float>::epsilon())
        return m_maxLongitudinalDisplacementX;

    const pandora::CartesianVector directionX(direction.GetX(), 0.f, 0.f);
    const pandora::CartesianVector directionYZ(0.f, direction.GetY(), direction.GetZ());

    return (m_maxLongitudinalDisplacementX * directionYZ.GetMagnitude() / directionX.GetMagnitude());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PFParticleStitcher::GetLongitudinalDisplacementCut3D(const pandora::CartesianVector &direction) const
{
    if (std::fabs(direction.GetX()) < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    if (std::fabs(direction.GetX()) > 1.f - std::numeric_limits<float>::epsilon())
        return m_maxLongitudinalDisplacementX;

    const pandora::CartesianVector directionX(direction.GetX(), 0.f, 0.f);
    const pandora::CartesianVector directionYZ(0.f, direction.GetY(), direction.GetZ());

    return (m_maxLongitudinalDisplacementX / directionX.GetMagnitude());    
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PFParticleStitcher::GetDeltaX(const pandora::CartesianVector &initialPosition, const pandora::CartesianVector &initialDirection, 
    const pandora::CartesianVector &targetPosition) const
{
    if (std::fabs(initialDirection.GetX()) > 1.f - std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    const pandora::CartesianVector targetPositionYZ(0.f, targetPosition.GetY(), targetPosition.GetZ());
    const pandora::CartesianVector initialPositionYZ(0.f, initialPosition.GetY(), initialPosition.GetZ());
    const pandora::CartesianVector initialDirectionYZ(0.f, initialDirection.GetY(), initialDirection.GetZ());
    const pandora::CartesianVector initialDirectionX(initialDirection.GetX(), 0.f, 0.f);
    
    const float R(initialDirectionYZ.GetUnitVector().GetDotProduct(initialPositionYZ - targetPositionYZ) / initialDirectionYZ.GetMagnitude());

    return (-1.f * initialDirectionX.GetUnitVector().GetDotProduct(targetPosition - (initialPosition - initialDirection * R)));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleStitcher::ParticleAssociation::ParticleAssociation(const VertexType parent, const VertexType daughter, const float fom) :
    m_parent(parent),
    m_daughter(daughter),
    m_fom(fom)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleStitcher::ParticleAssociation::VertexType PFParticleStitcher::ParticleAssociation::GetParent() const
{
    return m_parent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleStitcher::ParticleAssociation::VertexType PFParticleStitcher::ParticleAssociation::GetDaughter() const
{
    return m_daughter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PFParticleStitcher::ParticleAssociation::GetFigureOfMerit() const
{
    return m_fom;
}

} // namespace lar_pandora
