/**
 *  @file   larpandora/LArPandoraInterface/PFParticleMonitoring_module.cc
 *
 *  @brief  Analysis module for created particles
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

// ROOT includes
#include "TTree.h"

// Local includes
#include "LArPandoraCollector.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  PFParticleMonitoring class
 */
class PFParticleMonitoring : public art::EDAnalyzer
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
     PFParticleMonitoring(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
     virtual ~PFParticleMonitoring();

     void beginJob();
     void endJob();
     void analyze(const art::Event &evt);
     void reconfigure(fhicl::ParameterSet const &pset);

private:

    /**
     *  @brief Performing matching between true and reconstructed particles
     *
     *  @param recoParticlesToHits the mapping from reconstructed particles to hits
     *  @param trueHitsToParticles the mapping from hits to true particles
     *  @param matchedParticles the output matches between reconstructed and true particles
     *  @param matchedHits the output matches between reconstructed particles and hits
     */
     void GetRecoToTrueMatches(const PFParticlesToHits &recoParticlesToHits, const HitsToMCParticles &trueHitsToParticles,
         MCParticlesToPFParticles &matchedParticles, MCParticlesToHits &matchedHits) const;

    /**
     *  @brief Build particle maps for reconstructed particles
     *
     *  @param particleVector the input vector of reconstructed particles
     *  @param particleMap the output mapping between reconstructed particles and particle ID
     */
     void BuildRecoParticleMap(const PFParticleVector &particleVector, PFParticleMap &particleMap) const;

    /**
     *  @brief Build particle maps for true particles
     *
     *  @param particleVector the input vector of true particles
     *  @param particleMap the output mapping between true particle and true track ID
     */
     void BuildTrueParticleMap(const MCParticleVector &particleVector, MCParticleMap &particleMap) const;

    /**
     *  @brief Count the number of reconstructed hits in a given wire plane
     *
     *  @param view the wire plane ID
     *  @param hitVector the input vector of reconstructed hits
     */
     int CountHitsByType(const int view, const HitVector &hitVector) const;

    /**
     *  @brief Calculate length of true particle from its trajectory points
     *
     *  @param trueParticle the input true particle
     */
     double GetLength(const art::Ptr<simb::MCParticle> trueParticle) const;


     TTree       *m_pRecoTree;              ///<

     int          m_run;                    ///<
     int          m_event;                  ///<
     int          m_index;                  ///<

     int          m_nMCParticles;           ///<
     int          m_nNeutrinoPfos;          ///<
     int          m_nPrimaryPfos;           ///<
     int          m_nDaughterPfos;          ///<
     int          m_mcPdg;                  ///<
     int          m_mcNuPdg;                ///<
     int          m_pfoPdg;                 ///<
     int          m_pfoNuPdg;               ///<

     double       m_mcLength;               ///<
     double       m_completeness;           ///<
     double       m_purity;                 ///<

     int          m_nMCHits;                ///<
     int          m_nPfoHits;               ///<
     int          m_nMatchedHits;           ///<

     int          m_nMCHitsU;               ///<
     int          m_nMCHitsV;               ///<
     int          m_nMCHitsW;               ///<

     int          m_nPfoHitsU;              ///<
     int          m_nPfoHitsV;              ///<
     int          m_nPfoHitsW;              ///<

     int          m_nMatchedHitsU;          ///<
     int          m_nMatchedHitsV;          ///<
     int          m_nMatchedHitsW;          ///<

     std::string  m_hitfinderLabel;         ///<
     std::string  m_spacepointLabel;        ///<
     std::string  m_particleLabel;          ///<
     std::string  m_geantModuleLabel;       ///<

     bool         m_useDaughterPFParticles; ///<
     bool         m_useDaughterMCParticles; ///<
};

DEFINE_ART_MODULE(PFParticleMonitoring)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows
/**
 *  @file   LArPandora/PFParticleMonitoring.cxx
 *
 *  @brief  Implementation of the lar pandora analysis producer.
 *
 *  $Log: $
 */

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"

#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/PFParticle.h"

// std includes
#include <iostream>

namespace lar_pandora
{

PFParticleMonitoring::PFParticleMonitoring(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset)
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleMonitoring::~PFParticleMonitoring()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleMonitoring::reconfigure(fhicl::ParameterSet const &pset)
{
    m_particleLabel = pset.get<std::string>("PFParticleModule","pandora");
    m_spacepointLabel = pset.get<std::string>("SpacePointModule","pandora");
    m_hitfinderLabel = pset.get<std::string>("HitFinderModule","gaushit");
    m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");

    m_useDaughterPFParticles = pset.get<bool>("UseDaughterPFParticles",false);
    m_useDaughterMCParticles = pset.get<bool>("UseDaughterMCParticles",true);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleMonitoring::beginJob()
{
    mf::LogDebug("LArPandora") << " *** PFParticleMonitoring::beginJob() *** " << std::endl;

    //
    art::ServiceHandle<art::TFileService> tfs;

    m_pRecoTree = tfs->make<TTree>("pandora", "LAr Reco vs True");
    m_pRecoTree->Branch("run", &m_run,"run/I");
    m_pRecoTree->Branch("event", &m_event,"event/I");
    m_pRecoTree->Branch("index", &m_index,"index/I");
    m_pRecoTree->Branch("nMCParticles", &m_nMCParticles, "nMCParticles/I");
    m_pRecoTree->Branch("nNeutrinoPfos", &m_nNeutrinoPfos, "nNeutrinoPfos/I");
    m_pRecoTree->Branch("nPrimaryPfos", &m_nPrimaryPfos, "nPrimaryPfos/I");
    m_pRecoTree->Branch("nDaughterPfos", &m_nDaughterPfos, "nDaughterPfos/I");
    m_pRecoTree->Branch("mcPdg", &m_mcPdg , "mcPdg/I");
    m_pRecoTree->Branch("mcNuPdg", &m_mcNuPdg, "mcNuPdg/I");
    m_pRecoTree->Branch("pfoPdg", &m_pfoPdg, "pfoPdg/I");
    m_pRecoTree->Branch("pfoNuPdg", &m_pfoNuPdg, "pfoNuPdg/I");
    m_pRecoTree->Branch("mcLength", &m_mcLength, "mcLength/D");
    m_pRecoTree->Branch("completeness", &m_completeness, "completeness/D");
    m_pRecoTree->Branch("purity", &m_purity, "purity/D");
    m_pRecoTree->Branch("nMCHits", &m_nMCHits, "nMCHits/I");
    m_pRecoTree->Branch("nPfoHits", &m_nPfoHits, "nPfoHits/I");
    m_pRecoTree->Branch("nMatchedHits", &m_nMatchedHits, "nMatchedHits/I");
    m_pRecoTree->Branch("nMCHitsU", &m_nMCHitsU, "nMCHitsU/I");
    m_pRecoTree->Branch("nMCHitsV", &m_nMCHitsV, "nMCHitsV/I");
    m_pRecoTree->Branch("nMCHitsW", &m_nMCHitsW, "nMCHitsW/I");
    m_pRecoTree->Branch("nPfoHitsU", &m_nPfoHitsU, "nPfoHitsU/I");
    m_pRecoTree->Branch("nPfoHitsV", &m_nPfoHitsV, "nPfoHitsV/I");
    m_pRecoTree->Branch("nPfoHitsW", &m_nPfoHitsW, "nPfoHitsW/I");
    m_pRecoTree->Branch("nMatchedHitsU", &m_nMatchedHitsU, "nMatchedHitsU/I");
    m_pRecoTree->Branch("nMatchedHitsV", &m_nMatchedHitsV, "nMatchedHitsV/I");
    m_pRecoTree->Branch("nMatchedHitsW", &m_nMatchedHitsW, "nMatchedHitsW/I");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleMonitoring::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleMonitoring::analyze(const art::Event &evt)
{
    std::cout << " *** PFParticleMonitoring::analyze(...) *** " << std::endl;

    m_run = evt.run();
    m_event = evt.id().event();
    m_index = 0;

    m_nMCParticles = 0;
    m_nNeutrinoPfos = 0;
    m_nPrimaryPfos = 0;
    m_nDaughterPfos = 0;

    m_mcPdg = 0;
    m_mcNuPdg = 0;
    m_pfoPdg = 0;
    m_pfoNuPdg = 0;

    m_mcLength = 0.0;

    m_completeness = 0.0;
    m_purity = 0.0;

    m_nMCHits = 0;
    m_nPfoHits = 0;
    m_nMatchedHits = 0;
    m_nMCHitsU = 0;
    m_nMCHitsV = 0;
    m_nMCHitsW = 0;
    m_nPfoHitsU = 0;
    m_nPfoHitsV = 0;
    m_nPfoHitsW = 0;
    m_nMatchedHitsU = 0;
    m_nMatchedHitsV = 0;
    m_nMatchedHitsW = 0;

    std::cout << "  Run: " << m_run << std::endl;
    std::cout << "  Event: " << m_event << std::endl;


    // Match Reco Particles to Hits
    // ============================
    PFParticleVector recoParticleVector;
    PFParticlesToHits recoParticlesToHits;
    HitsToPFParticles recoHitsToParticles;

    LArPandoraCollector::CollectPFParticles(evt, m_particleLabel, recoParticleVector);
    LArPandoraCollector::BuildPFParticleHitMaps(evt, m_particleLabel, m_spacepointLabel, recoParticlesToHits, recoHitsToParticles,
        (m_useDaughterPFParticles ? LArPandoraCollector::kAddDaughters : LArPandoraCollector::kIgnoreDaughters));

    std::cout << "  PFParticles: " << recoParticleVector.size() << std::endl;

    // Match True Particles to Hits
    // ============================
    HitVector hitVector;
    MCTruthToMCParticles truthToParticles;
    MCParticlesToMCTruth particlesToTruth;
    MCParticlesToHits trueParticlesToHits;
    HitsToMCParticles trueHitsToParticles;

    LArPandoraCollector::CollectHits(evt, m_hitfinderLabel, hitVector);
    LArPandoraCollector::CollectMCParticles(evt, m_geantModuleLabel, truthToParticles, particlesToTruth);
    LArPandoraCollector::BuildMCParticleHitMaps(evt, m_geantModuleLabel, hitVector, trueParticlesToHits, trueHitsToParticles,
        (m_useDaughterMCParticles ? LArPandoraCollector::kAddDaughters : LArPandoraCollector::kIgnoreDaughters));

    std::cout << "  MCParticles: " << trueParticlesToHits.size() << std::endl;

    if (trueParticlesToHits.empty())
    {
        m_pRecoTree->Fill();
        return;
    }

    // Match Reco Particles to True Particles
    // ======================================
    MCParticlesToPFParticles matchedParticles;
    MCParticlesToHits matchedHits;
    this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits);


    // Build Reco Particle Maps (for Parent/Daughter Navigation)
    // =========================================================
    PFParticleMap recoParticleMap;
    this->BuildRecoParticleMap(recoParticleVector, recoParticleMap);


    // Loop over True Particles and Write out Reco/True information
    // ============================================================
    m_nMCParticles  = trueParticlesToHits.size();
    m_nNeutrinoPfos = 0;
    m_nPrimaryPfos  = 0;
    m_nDaughterPfos = 0;

    // Count reconstructed particles
    for (PFParticleVector::const_iterator iter = recoParticleVector.begin(), iterEnd = recoParticleVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> recoParticle = *iter;

        if (LArPandoraCollector::IsNeutrino(recoParticle))
        {
            m_nNeutrinoPfos++;
        }
        else if (LArPandoraCollector::IsFinalState(recoParticleMap, recoParticle))
        {
            m_nPrimaryPfos++;
        }
        else
        {
            m_nDaughterPfos++;
        }
    }

    // Compare true and reconstructed particles
    for (MCParticlesToHits::const_iterator iter = trueParticlesToHits.begin(), iterEnd = trueParticlesToHits.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<simb::MCParticle> trueParticle = iter->first;
        const HitVector &trueHitVector = iter->second;

        if (trueHitVector.empty())
            continue;

        m_mcPdg = trueParticle->PdgCode();
        m_mcNuPdg = 0;
        m_pfoPdg = 0;
        m_pfoNuPdg = 0;

        m_mcLength = this->GetLength(trueParticle);

        MCParticlesToMCTruth::const_iterator nuIter = particlesToTruth.find(trueParticle);
        if (particlesToTruth.end() == nuIter)
            throw cet::exception("LArPandora") << " PFParticleMonitoring::analyze --- Found a true particle without any ancestry information ";

        const art::Ptr<simb::MCTruth> trueNeutrino = nuIter->second;

        if (trueNeutrino->NeutrinoSet())
        {
            const simb::MCNeutrino neutrino(trueNeutrino->GetNeutrino());
            m_mcNuPdg = neutrino.Nu().PdgCode();
        }

        m_nMCHits = trueHitVector.size();
        m_nMCHitsU = this->CountHitsByType(geo::kU, trueHitVector);
        m_nMCHitsV = this->CountHitsByType(geo::kV, trueHitVector);
        m_nMCHitsW = this->CountHitsByType(geo::kW, trueHitVector);

        m_nPfoHits = 0;
        m_nPfoHitsU = 0;
        m_nPfoHitsV = 0;
        m_nPfoHitsW = 0;

        m_nMatchedHits = 0;
        m_nMatchedHitsU = 0;
        m_nMatchedHitsV = 0;
        m_nMatchedHitsW = 0;

        MCParticlesToPFParticles::const_iterator pIter1 = matchedParticles.find(trueParticle);
        if (matchedParticles.end() != pIter1)
        {
            const art::Ptr<recob::PFParticle> recoParticle = pIter1->second;
            m_pfoPdg = recoParticle->PdgCode();
            m_pfoNuPdg = LArPandoraCollector::GetParentNeutrino(recoParticleMap, recoParticle);

            PFParticlesToHits::const_iterator pIter2 = recoParticlesToHits.find(recoParticle);
            if (recoParticlesToHits.end() == pIter2)
                throw cet::exception("LArPandora") << " PFParticleMonitoring::analyze --- Found a reco particle without any hits ";

            const HitVector &recoHitVector = pIter2->second;

            MCParticlesToHits::const_iterator pIter3 = matchedHits.find(trueParticle);
            if (matchedHits.end() == pIter3)
                throw cet::exception("LArPandora") << " PFParticleMonitoring::analyze --- Found a matched true particle without matched hits ";

            const HitVector &matchedHitVector = pIter3->second;

            m_nPfoHits = recoHitVector.size();
            m_nPfoHitsU = this->CountHitsByType(geo::kU, recoHitVector);
            m_nPfoHitsV = this->CountHitsByType(geo::kV, recoHitVector);
            m_nPfoHitsW = this->CountHitsByType(geo::kW, recoHitVector);

            m_nMatchedHits = matchedHitVector.size();
            m_nMatchedHitsU = this->CountHitsByType(geo::kU, matchedHitVector);
            m_nMatchedHitsV = this->CountHitsByType(geo::kV, matchedHitVector);
            m_nMatchedHitsW = this->CountHitsByType(geo::kW, matchedHitVector);
        }

        m_purity = ((m_nPfoHits == 0) ? 0.0 : static_cast<double>(m_nMatchedHits) / static_cast<double>(m_nPfoHits));
        m_completeness = ((m_nPfoHits == 0) ? 0.0 : static_cast<double>(m_nMatchedHits) / static_cast<double>(m_nMCHits));

        std::cout << "    MCParticle [" << m_index << "] trueId=" << trueParticle->TrackId()
                  << ", trueNu=" << m_mcNuPdg << ", truePdg=" << m_mcPdg << ", recoNu=" << m_pfoNuPdg << ", recoPdg=" << m_pfoPdg
                  << ", mcHits=" << m_nMCHits << ", pfoHits=" << m_nPfoHits << ", matchedHits=" << m_nMatchedHits << std::endl;

        m_pRecoTree->Fill();
        ++m_index; // Increment index number
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleMonitoring::GetRecoToTrueMatches(const PFParticlesToHits &recoParticlesToHits, const HitsToMCParticles &trueHitsToParticles,
    MCParticlesToPFParticles &matchedParticles, MCParticlesToHits &matchedHits) const
{
    for (PFParticlesToHits::const_iterator iter1 = recoParticlesToHits.begin(), iterEnd1 = recoParticlesToHits.end();
        iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> recoParticle = iter1->first;
        const HitVector &hitVector = iter1->second;

        MCParticlesToHits truthContributionMap;

        for (HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::Hit> hit = *iter2;

            HitsToMCParticles::const_iterator iter3 = trueHitsToParticles.find(hit);
            if (trueHitsToParticles.end() == iter3)
                continue;

            const art::Ptr<simb::MCParticle> trueParticle = iter3->second;

            truthContributionMap[trueParticle].push_back(hit);
        }

        MCParticlesToHits::const_iterator mIter = truthContributionMap.end();

        for (MCParticlesToHits::const_iterator iter4 = truthContributionMap.begin(), iterEnd4 = truthContributionMap.end();
            iter4 != iterEnd4; ++iter4)
        {
            if ((truthContributionMap.end() == mIter) || (iter4->second.size() > mIter->second.size()))
            {
                mIter = iter4;
            }
        }

        if (truthContributionMap.end() != mIter)
        {
            const art::Ptr<simb::MCParticle> trueParticle = mIter->first;

            MCParticlesToHits::const_iterator iter5 = matchedHits.find(trueParticle);

            if ((matchedHits.end() == iter5) || (mIter->second.size() > iter5->second.size()))
            {
                matchedParticles[trueParticle] = recoParticle;
                matchedHits[trueParticle] = mIter->second;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleMonitoring::BuildRecoParticleMap(const PFParticleVector &particleVector, PFParticleMap &particleMap) const
{
    for (PFParticleVector::const_iterator iter = particleVector.begin(), iterEnd = particleVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = *iter;
        particleMap[particle->Self()] = particle;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleMonitoring::BuildTrueParticleMap(const MCParticleVector &particleVector, MCParticleMap &particleMap) const
{
    for (MCParticleVector::const_iterator iter = particleVector.begin(), iterEnd = particleVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<simb::MCParticle> particle = *iter;
        particleMap[particle->TrackId()] = particle;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int PFParticleMonitoring::CountHitsByType(const int view, const HitVector &hitVector) const
{
    int nHits(0);

    for (HitVector::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::Hit> hit = *iter;
        if (hit->View() == view)
            ++nHits;
    }

    return nHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double PFParticleMonitoring::GetLength(const art::Ptr<simb::MCParticle> particle) const
{
    art::ServiceHandle<geo::Geometry> theGeometry;

    bool foundStartPosition(false);
    int startT(0), endT(0);

    const int numTrajectoryPoints(static_cast<int>(particle->NumberTrajectoryPoints()));

    for (int nt = 0; nt < numTrajectoryPoints; ++nt)
    {
        try
        {
            double pos[3] = {particle->Vx(nt), particle->Vy(nt), particle->Vz(nt)};
            unsigned int which_tpc(std::numeric_limits<unsigned int>::max());
            unsigned int which_cstat(std::numeric_limits<unsigned int>::max());
            theGeometry->PositionToTPC(pos, which_tpc, which_cstat);

            endT = nt;
            if (!foundStartPosition)
            {
                startT = endT;
                foundStartPosition = true;
            }
        }
        catch (cet::exception &e){
            continue;
        }
    }

    if (!foundStartPosition)
        return 0.0;

    const double dx(particle->Vx(endT) - particle->Vx(startT));
    const double dy(particle->Vy(endT) - particle->Vy(startT));
    const double dz(particle->Vz(endT) - particle->Vz(startT));

    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} //namespace lar_pandora
