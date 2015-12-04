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

// Local LArPandora includes
#include "LArPandoraInterface/LArPandoraCollector.h"

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

    typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
    typedef std::set< art::Ptr<simb::MCParticle> > MCParticleSet;
    
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
     *  @brief Performing matching between true and reconstructed particles
     *
     *  @param recoParticlesToHits the mapping from reconstructed particles to hits
     *  @param trueHitsToParticles the mapping from hits to true particles
     *  @param matchedParticles the output matches between reconstructed and true particles
     *  @param matchedHits the output matches between reconstructed particles and hits
     *  @param recoVeto the veto list for reconstructed particles
     *  @param trueVeto the veto list for true particles
     */
     void GetRecoToTrueMatches(const PFParticlesToHits &recoParticlesToHits, const HitsToMCParticles &trueHitsToParticles,
         MCParticlesToPFParticles &matchedParticles, MCParticlesToHits &matchedHits, PFParticleSet &recoVeto, MCParticleSet &trueVeto) const;

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
     *  @brief Find the start and end points of the true particle in the active region of detector
     *
     *  @param trueParticle the input true particle
     *  @param startT  the true start point
     *  @param endT  the true end point
     */
     void GetStartAndEndPoints(const art::Ptr<simb::MCParticle> trueParticle, int &startT, int &endT) const;

    /**
     *  @brief Find the length of the true particle trajectory through the active region of the detector
     *
     *  @param trueParticle the input true particle
     *  @param startT  the true start point
     *  @param endT  the true end point
     */
     double GetLength(const art::Ptr<simb::MCParticle> trueParticle, const int startT, const int endT) const;
     

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
     int          m_mcParentPdg;            ///<
     int          m_mcPrimaryPdg;           ///<
     int          m_mcIsPrimary;            ///<
     int          m_mcIsDecay;              ///<
     int          m_pfoPdg;                 ///<
     int          m_pfoNuPdg;               ///<
     int          m_pfoIsPrimary;           ///<

     int          m_pfoTrack;               ///<
     int          m_pfoVertex;              ///<
     double       m_pfoVtxX;                ///<
     double       m_pfoVtxY;                ///<
     double       m_pfoVtxZ;                ///<
     double       m_pfoEndX;                ///<
     double       m_pfoEndY;                ///<
     double       m_pfoEndZ;                ///<
     double       m_pfoDirX;                ///<
     double       m_pfoDirY;                ///<
     double       m_pfoDirZ;                ///<
     double       m_pfoLength;              ///<
     double       m_pfoStraightLength;      ///<

     int          m_mcVertex;               ///<
     double       m_mcVtxX;                 ///<
     double       m_mcVtxY;                 ///<
     double       m_mcVtxZ;                 ///<
     double       m_mcEndX;                 ///<
     double       m_mcEndY;                 ///<
     double       m_mcEndZ;                 ///<
     double       m_mcDirX;                 ///<
     double       m_mcDirY;                 ///<
     double       m_mcDirZ;                 ///<
     double       m_mcEnergy;               ///<
     double       m_mcLength;               ///<
     double       m_mcStraightLength;       ///<

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

     int          m_nAvailableHits;         ///<

     double       m_spacepointsMinX;        ///<
     double       m_spacepointsMaxX;        ///<
     
     std::string  m_hitfinderLabel;         ///<
     std::string  m_spacepointLabel;        ///<
     std::string  m_particleLabel;          ///<
     std::string  m_trackLabel;             ///<
     std::string  m_geantModuleLabel;       ///<

     bool         m_useDaughterPFParticles; ///<
     bool         m_useDaughterMCParticles; ///<
     bool         m_addDaughterPFParticles; ///<
     bool         m_addDaughterMCParticles; ///<
     bool         m_recursiveMatching;      ///<
     bool         m_printDebug;             ///< switch for print statements (TODO: use message service!)
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
#include "RecoBase/Track.h"

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
    m_trackLabel = pset.get<std::string>("TrackModule","pandora");
    m_particleLabel = pset.get<std::string>("PFParticleModule","pandora");
    m_spacepointLabel = pset.get<std::string>("SpacePointModule","pandora");
    m_hitfinderLabel = pset.get<std::string>("HitFinderModule","gaushit");
    m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");

    m_useDaughterPFParticles = pset.get<bool>("UseDaughterPFParticles",false);
    m_useDaughterMCParticles = pset.get<bool>("UseDaughterMCParticles",true);
    m_addDaughterPFParticles = pset.get<bool>("AddDaughterPFParticles",true);
    m_addDaughterMCParticles = pset.get<bool>("AddDaughterMCParticles",true);
    m_recursiveMatching = pset.get<bool>("RecursiveMatching",false);
    m_printDebug = pset.get<bool>("PrintDebug",false);
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
    m_pRecoTree->Branch("mcParentPdg", &m_mcParentPdg, "mcParentPdg/I");
    m_pRecoTree->Branch("mcPrimaryPdg", &m_mcPrimaryPdg, "mcPrimaryPdg/I");
    m_pRecoTree->Branch("mcIsPrimary", &m_mcIsPrimary, "mcIsPrimary/I");
    m_pRecoTree->Branch("mcIsDecay", &m_mcIsDecay, "mcIsDecay/I");
    m_pRecoTree->Branch("pfoPdg", &m_pfoPdg, "pfoPdg/I");
    m_pRecoTree->Branch("pfoNuPdg", &m_pfoNuPdg, "pfoNuPdg/I");
    m_pRecoTree->Branch("pfoIsPrimary", &m_pfoIsPrimary, "pfoIsPrimary/I");
    m_pRecoTree->Branch("pfoTrack", &m_pfoTrack, "pfoTrack/I");
    m_pRecoTree->Branch("pfoVertex", &m_pfoVertex, "pfoVertex/I");
    m_pRecoTree->Branch("pfoVtxX", &m_pfoVtxX, "pfoVtxX/D");
    m_pRecoTree->Branch("pfoVtxY", &m_pfoVtxY, "pfoVtxY/D");
    m_pRecoTree->Branch("pfoVtxZ", &m_pfoVtxZ, "pfoVtxZ/D");
    m_pRecoTree->Branch("pfoEndX", &m_pfoEndX, "pfoEndX/D");
    m_pRecoTree->Branch("pfoEndY", &m_pfoEndY, "pfoEndY/D");
    m_pRecoTree->Branch("pfoEndZ", &m_pfoEndZ, "pfoEndZ/D");
    m_pRecoTree->Branch("pfoDirX", &m_pfoDirX, "pfoDirX/D");
    m_pRecoTree->Branch("pfoDirY", &m_pfoDirY, "pfoDirY/D");
    m_pRecoTree->Branch("pfoDirZ", &m_pfoDirZ, "pfoDirZ/D");
    m_pRecoTree->Branch("pfoLength", &m_pfoLength, "pfoLength/D");
    m_pRecoTree->Branch("pfoStraightLength", &m_pfoStraightLength, "pfoStraightLength/D");
    m_pRecoTree->Branch("mcVertex", &m_mcVertex, "mcVertex/I");
    m_pRecoTree->Branch("mcVtxX", &m_mcVtxX, "mcVtxX/D");
    m_pRecoTree->Branch("mcVtxY", &m_mcVtxY, "mcVtxY/D");
    m_pRecoTree->Branch("mcVtxZ", &m_mcVtxZ, "mcVtxZ/D");
    m_pRecoTree->Branch("mcEndX", &m_mcEndX, "mcEndX/D");
    m_pRecoTree->Branch("mcEndY", &m_mcEndY, "mcEndY/D");
    m_pRecoTree->Branch("mcEndZ", &m_mcEndZ, "mcEndZ/D"); 
    m_pRecoTree->Branch("mcDirX", &m_mcDirX, "mcDirX/D");
    m_pRecoTree->Branch("mcDirY", &m_mcDirY, "mcDirY/D");
    m_pRecoTree->Branch("mcDirZ", &m_mcDirZ, "mcDirZ/D");
    m_pRecoTree->Branch("mcEnergy", &m_mcEnergy, "mcEnergy/D");
    m_pRecoTree->Branch("mcLength", &m_mcLength, "mcLength/D");
    m_pRecoTree->Branch("mcStraightLength", &m_mcStraightLength, "mcStraightLength/D");
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
    m_pRecoTree->Branch("nAvailableHits", &m_nAvailableHits, "nAvailableHits/I");
    m_pRecoTree->Branch("spacepointsMinX", &m_spacepointsMinX, "spacepointsMinX/D");
    m_pRecoTree->Branch("spacepointsMaxX", &m_spacepointsMaxX, "spacepointsMaxX/D");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleMonitoring::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleMonitoring::analyze(const art::Event &evt)
{
    if (m_printDebug)
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
    m_mcParentPdg = 0;
    m_mcPrimaryPdg = 0;
    m_mcIsPrimary = 0;
    m_mcIsDecay = 0;

    m_pfoPdg = 0;
    m_pfoNuPdg = 0;
    m_pfoIsPrimary = 0;
    m_pfoTrack = 0;
    m_pfoVertex = 0;
    m_pfoVtxX = 0.0;
    m_pfoVtxY = 0.0;
    m_pfoVtxZ = 0.0;
    m_pfoEndX = 0.0;
    m_pfoEndY = 0.0;
    m_pfoEndZ = 0.0;
    m_pfoDirX = 0.0;
    m_pfoDirY = 0.0;
    m_pfoDirZ = 0.0;
    m_pfoLength = 0.0;
    m_pfoStraightLength = 0.0;
   
    m_mcVertex = 0;
    m_mcVtxX = 0.0;
    m_mcVtxY = 0.0;
    m_mcVtxZ = 0.0;
    m_mcEndX = 0.0;
    m_mcEndY = 0.0;
    m_mcEndZ = 0.0;
    m_mcDirX = 0.0;
    m_mcDirY = 0.0;
    m_mcDirZ = 0.0;
    m_mcEnergy = 0.0;
    m_mcLength = 0.0;
    m_mcStraightLength = 0.0;

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
    m_nAvailableHits = 0;

    m_spacepointsMinX = 0.0;
    m_spacepointsMaxX = 0.0;
    
    if (m_printDebug)
    {
        std::cout << "  Run: " << m_run << std::endl;
        std::cout << "  Event: " << m_event << std::endl;
    }

    // Collect Hits
    // ============
    HitVector hitVector;
    LArPandoraCollector::CollectHits(evt, m_hitfinderLabel, hitVector);

    if (m_printDebug)
        std::cout << "  Hits: " << hitVector.size() << std::endl;

    // Collect SpacePoints and SpacePoint <-> Hit Associations
    // =======================================================
    SpacePointVector spacePointVector;
    SpacePointsToHits spacePointsToHits;
    HitsToSpacePoints hitsToSpacePoints;
    LArPandoraCollector::CollectSpacePoints(evt, m_spacepointLabel, spacePointVector, spacePointsToHits, hitsToSpacePoints);

    if (m_printDebug)
        std::cout << "  SpacePoints: " << spacePointVector.size() << std::endl;

    // Collect Tracks and PFParticle <-> Track Associations
    // ====================================================
    TrackVector recoTrackVector;
    PFParticlesToTracks recoParticlesToTracks;
    LArPandoraCollector::CollectTracks(evt, m_trackLabel, recoTrackVector, recoParticlesToTracks);
    
    if (m_printDebug)
        std::cout << "  Tracks: " << recoTrackVector.size() << std::endl;

    // Collect Vertices and PFParticle <-> Vertex Associations
    // =======================================================
    VertexVector recoVertexVector;
    PFParticlesToVertices recoParticlesToVertices;
    LArPandoraCollector::CollectVertices(evt, m_particleLabel, recoVertexVector, recoParticlesToVertices);

    if (m_printDebug)
        std::cout << "  Vertices: " << recoVertexVector.size() << std::endl;

    // Collect PFParticles and match Reco Particles to Hits
    // ====================================================
    PFParticleVector recoParticleVector;
    PFParticlesToHits recoParticlesToHits;
    HitsToPFParticles recoHitsToParticles;

    LArPandoraCollector::CollectPFParticles(evt, m_particleLabel, recoParticleVector);
    LArPandoraCollector::BuildPFParticleHitMaps(evt, m_particleLabel, m_spacepointLabel, recoParticlesToHits, recoHitsToParticles,
        (m_useDaughterPFParticles ? (m_addDaughterPFParticles ? LArPandoraCollector::kAddDaughters : LArPandoraCollector::kUseDaughters) : LArPandoraCollector::kIgnoreDaughters)); 

    if (m_printDebug)
        std::cout << "  PFParticles: " << recoParticleVector.size() << std::endl;

    // Collect MCParticles and match True Particles to Hits
    // ====================================================
    MCParticleVector trueParticleVector;
    MCTruthToMCParticles truthToParticles;
    MCParticlesToMCTruth particlesToTruth;
    MCParticlesToHits trueParticlesToHits;
    HitsToMCParticles trueHitsToParticles;

    if (!evt.isRealData())
    {
        LArPandoraCollector::CollectMCParticles(evt, m_geantModuleLabel, trueParticleVector);
        LArPandoraCollector::CollectMCParticles(evt, m_geantModuleLabel, truthToParticles, particlesToTruth);
        LArPandoraCollector::BuildMCParticleHitMaps(evt, m_geantModuleLabel, hitVector, trueParticlesToHits, trueHitsToParticles,
            (m_useDaughterMCParticles ? (m_addDaughterMCParticles ? LArPandoraCollector::kAddDaughters : LArPandoraCollector::kUseDaughters) : LArPandoraCollector::kIgnoreDaughters));
    }

    if (m_printDebug)
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


    // Build True Particle Maps (for Parent/Daughter Navigation)
    // =========================================================
    MCParticleMap trueParticleMap;
    this->BuildTrueParticleMap(trueParticleVector, trueParticleMap);


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
        m_mcParentPdg = 0;
        m_mcPrimaryPdg = 0;
        m_mcIsPrimary = 0;
        m_mcIsDecay = 0;
 
        m_pfoPdg = 0;
        m_pfoNuPdg = 0;
        m_pfoIsPrimary = 0;
        m_pfoTrack = 0;
        m_pfoVertex = 0;
        m_pfoVtxX = 0.0;
        m_pfoVtxY = 0.0;
        m_pfoVtxZ = 0.0;
        m_pfoEndX = 0.0;
        m_pfoEndY = 0.0;
        m_pfoEndZ = 0.0;
        m_pfoDirX = 0.0;
        m_pfoDirY = 0.0;
        m_pfoDirZ = 0.0;
        m_pfoLength = 0.0;
        m_pfoStraightLength = 0.0;

        m_mcVertex = 0;
        m_mcVtxX = 0.0;
        m_mcVtxY = 0.0;
        m_mcVtxZ = 0.0;
        m_mcEndX = 0.0;
        m_mcEndY = 0.0;
        m_mcEndZ = 0.0;
        m_mcDirX = 0.0;
        m_mcDirY = 0.0;
        m_mcDirZ = 0.0;
        m_mcEnergy = 0.0;
        m_mcLength = 0.0;
        m_mcStraightLength = 0.0;

        m_completeness = 0.0;
        m_purity = 0.0;

        m_nMCHits = 0;
        m_nMCHitsU = 0;
        m_nMCHitsV = 0;
        m_nMCHitsW = 0;

        m_nPfoHits = 0;
        m_nPfoHitsU = 0;
        m_nPfoHitsV = 0;
        m_nPfoHitsW = 0;

        m_nMatchedHits = 0;
        m_nMatchedHitsU = 0;
        m_nMatchedHitsV = 0;
        m_nMatchedHitsW = 0;
        
        m_nAvailableHits = 0;
 
        m_spacepointsMinX = 0.0;
        m_spacepointsMaxX = 0.0;

        // Set true properties
        try
	{
            int startT(-1);
            int endT(-1);
            this->GetStartAndEndPoints(trueParticle, startT, endT);

	    // vertex and end positions
            m_mcVertex = 1;
            m_mcVtxX = trueParticle->Vx(startT);
            m_mcVtxY = trueParticle->Vy(startT);
            m_mcVtxZ = trueParticle->Vz(startT);
            m_mcEndX = trueParticle->Vx(endT);
            m_mcEndY = trueParticle->Vy(endT); 
            m_mcEndZ = trueParticle->Vz(endT);
            
            const double dx(m_mcEndX - m_mcVtxX);
            const double dy(m_mcEndY - m_mcVtxY);
            const double dz(m_mcEndZ - m_mcVtxZ);

            m_mcStraightLength = std::sqrt(dx * dx + dy *dy + dz * dz);
            m_mcLength = this->GetLength(trueParticle, startT, endT);
            
            // energy and momentum
            const double Ptot(trueParticle->P(startT));

            if (Ptot > 0.0)
	    {
	        m_mcDirX = trueParticle->Px(startT) / Ptot;
                m_mcDirY = trueParticle->Py(startT) / Ptot;
                m_mcDirZ = trueParticle->Pz(startT) / Ptot;
                m_mcEnergy = trueParticle->E(startT);
	    }
	}
        catch (cet::exception &e){
	}

        // Get the true parent neutrino
        MCParticlesToMCTruth::const_iterator nuIter = particlesToTruth.find(trueParticle);
        if (particlesToTruth.end() == nuIter)
            throw cet::exception("LArPandora") << " PFParticleMonitoring::analyze --- Found a true particle without any ancestry information ";

        const art::Ptr<simb::MCTruth> trueNeutrino = nuIter->second;

        if (trueNeutrino->NeutrinoSet())
        {
            const simb::MCNeutrino neutrino(trueNeutrino->GetNeutrino());
            m_mcNuPdg = neutrino.Nu().PdgCode();
        }

	// Get the true 'parent' and 'primary' particles
        try
	{
	    const art::Ptr<simb::MCParticle> parentParticle(LArPandoraCollector::GetParentMCParticle(trueParticleMap, trueParticle));
            const art::Ptr<simb::MCParticle> primaryParticle(LArPandoraCollector::GetPrimaryMCParticle(trueParticleMap, trueParticle));
            m_mcParentPdg = ((parentParticle != trueParticle) ? parentParticle->PdgCode() : 0);
            m_mcPrimaryPdg = primaryParticle->PdgCode();
            m_mcIsPrimary = (primaryParticle == trueParticle);
            m_mcIsDecay = ("Decay" == trueParticle->Process()); 
	}
        catch (cet::exception &e){
	}

        // Find min and max X positions of space points
        bool foundSpacePoints(false);

        for (HitVector::const_iterator hIter1 = trueHitVector.begin(), hIterEnd1 = trueHitVector.end(); hIter1 != hIterEnd1; ++hIter1)
	{
	    const art::Ptr<recob::Hit> hit = *hIter1;
            
	    HitsToSpacePoints::const_iterator hIter2 = hitsToSpacePoints.find(hit);
            if (hitsToSpacePoints.end() == hIter2)
	        continue;

            const art::Ptr<recob::SpacePoint> spacepoint = hIter2->second;
            const double X(spacepoint->XYZ()[0]);

            if (!foundSpacePoints)
	    {
                m_spacepointsMinX = X;
                m_spacepointsMaxX = X;
                foundSpacePoints = true;
	    }
            else
	    {
	        m_spacepointsMinX = std::min(m_spacepointsMinX, X);
                m_spacepointsMaxX = std::max(m_spacepointsMaxX, X);
	    }
	}
  
        // Count number of available hits
        for (HitVector::const_iterator hIter = trueHitVector.begin(), hIterEnd = trueHitVector.end(); hIter != hIterEnd; ++hIter)
	{            
	    if (recoHitsToParticles.find(*hIter) == recoHitsToParticles.end())
                ++m_nAvailableHits;
	}

        // Match true and reconstructed hits
        m_nMCHits = trueHitVector.size();
        m_nMCHitsU = this->CountHitsByType(geo::kU, trueHitVector);
        m_nMCHitsV = this->CountHitsByType(geo::kV, trueHitVector);
        m_nMCHitsW = this->CountHitsByType(geo::kW, trueHitVector);

        MCParticlesToPFParticles::const_iterator pIter1 = matchedParticles.find(trueParticle);
        if (matchedParticles.end() != pIter1)
        {
            const art::Ptr<recob::PFParticle> recoParticle = pIter1->second;
            m_pfoPdg = recoParticle->PdgCode();
            m_pfoNuPdg = LArPandoraCollector::GetParentNeutrino(recoParticleMap, recoParticle);
            m_pfoIsPrimary = LArPandoraCollector::IsFinalState(recoParticleMap, recoParticle);

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

	    PFParticlesToVertices::const_iterator pIter4 = recoParticlesToVertices.find(recoParticle);
            if (recoParticlesToVertices.end() != pIter4)
	    {
                const VertexVector &vertexVector = pIter4->second;
                if (!vertexVector.empty())
                {
                    if (vertexVector.size() !=1 && m_printDebug)
                        std::cout << " Warning: Found particle with more than one associated vertex " << std::endl;

                    const art::Ptr<recob::Vertex> recoVertex = *(vertexVector.begin());
                    double xyz[3] = {0.0, 0.0, 0.0} ;
                    recoVertex->XYZ(xyz);

                    m_pfoVertex = 1;
                    m_pfoVtxX = xyz[0];
                    m_pfoVtxY = xyz[1];
                    m_pfoVtxZ = xyz[2];
                }
	    }

	    PFParticlesToTracks::const_iterator pIter5 = recoParticlesToTracks.find(recoParticle);
            if (recoParticlesToTracks.end() != pIter5)
	    {
                const TrackVector &trackVector = pIter5->second;
                if (!trackVector.empty())
                {
                    if (trackVector.size() !=1 && m_printDebug)
                        std::cout << " Warning: Found particle with more than one associated track " << std::endl;
	        
                    const art::Ptr<recob::Track> recoTrack = *(trackVector.begin());
                    const TVector3 &vtxPosition = recoTrack->Vertex();
                    const TVector3 &endPosition = recoTrack->End();
                    const TVector3 &vtxDirection = recoTrack->VertexDirection();

                    m_pfoTrack = 1;
                    m_pfoVtxX = vtxPosition.x();
                    m_pfoVtxY = vtxPosition.y();
                    m_pfoVtxZ = vtxPosition.z();
                    m_pfoEndX = endPosition.x();
                    m_pfoEndY = endPosition.y();
                    m_pfoEndZ = endPosition.z();
                    m_pfoDirX = vtxDirection.x();
                    m_pfoDirY = vtxDirection.y();
                    m_pfoDirZ = vtxDirection.z();
		    m_pfoStraightLength = (endPosition - vtxPosition).Mag();
                    m_pfoLength = recoTrack->Length();
		}
	    }
        }

        m_purity = ((m_nPfoHits == 0) ? 0.0 : static_cast<double>(m_nMatchedHits) / static_cast<double>(m_nPfoHits));
        m_completeness = ((m_nPfoHits == 0) ? 0.0 : static_cast<double>(m_nMatchedHits) / static_cast<double>(m_nMCHits));

        if (m_printDebug)
          std::cout << "    MCParticle [" << m_index << "] trueId=" << trueParticle->TrackId()
                    << ", trueNu=" << m_mcNuPdg << ", truePdg=" << m_mcPdg << ", recoNu=" << m_pfoNuPdg << ", recoPdg=" << m_pfoPdg
                    << ", mcHits=" << m_nMCHits << ", pfoHits=" << m_nPfoHits << ", matchedHits=" << m_nMatchedHits 
                    << ", availableHits=" << m_nAvailableHits << std::endl;

        m_pRecoTree->Fill();
        ++m_index; // Increment index number
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleMonitoring::GetRecoToTrueMatches(const PFParticlesToHits &recoParticlesToHits, const HitsToMCParticles &trueHitsToParticles,
    MCParticlesToPFParticles &matchedParticles, MCParticlesToHits &matchedHits) const
{
    PFParticleSet recoVeto; MCParticleSet trueVeto;

    this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, recoVeto, trueVeto);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleMonitoring::GetRecoToTrueMatches(const PFParticlesToHits &recoParticlesToHits, const HitsToMCParticles &trueHitsToParticles,
    MCParticlesToPFParticles &matchedParticles, MCParticlesToHits &matchedHits, PFParticleSet &vetoReco, MCParticleSet &vetoTrue) const
{
    bool foundMatches(false);

    for (PFParticlesToHits::const_iterator iter1 = recoParticlesToHits.begin(), iterEnd1 = recoParticlesToHits.end();
        iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> recoParticle = iter1->first;
        if (vetoReco.count(recoParticle) > 0)
	    continue;

        const HitVector &hitVector = iter1->second;

        MCParticlesToHits truthContributionMap;

        for (HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::Hit> hit = *iter2;

            HitsToMCParticles::const_iterator iter3 = trueHitsToParticles.find(hit);
            if (trueHitsToParticles.end() == iter3)
                continue;

            const art::Ptr<simb::MCParticle> trueParticle = iter3->second;
            if (vetoTrue.count(trueParticle) > 0)
	        continue;
	    
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
                foundMatches = true;
            }
        }
    }

    if (!foundMatches)
        return;

    for (MCParticlesToPFParticles::const_iterator pIter = matchedParticles.begin(), pIterEnd = matchedParticles.end(); 
        pIter != pIterEnd; ++pIter)
    {
        vetoTrue.insert(pIter->first);
        vetoReco.insert(pIter->second);
    }
 
    if (m_recursiveMatching)
        this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, vetoReco, vetoTrue);
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

void PFParticleMonitoring::GetStartAndEndPoints(const art::Ptr<simb::MCParticle> particle, int &startT, int &endT) const
{
    art::ServiceHandle<geo::Geometry> theGeometry;

    bool foundStartPosition(false);

    const int numTrajectoryPoints(static_cast<int>(particle->NumberTrajectoryPoints()));

    for (int nt = 0; nt < numTrajectoryPoints; ++nt)
    {
        try
        {
            double pos[3] = {particle->Vx(nt), particle->Vy(nt), particle->Vz(nt)};
            unsigned int which_tpc(std::numeric_limits<unsigned int>::max());
            unsigned int which_cstat(std::numeric_limits<unsigned int>::max());
            theGeometry->PositionToTPC(pos, which_tpc, which_cstat);

            // TODO: Apply fiducial cut due to readout window

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
        throw cet::exception("LArPandora");
}

//------------------------------------------------------------------------------------------------------------------------------------------

double PFParticleMonitoring::GetLength(const art::Ptr<simb::MCParticle> particle, const int startT, const int endT) const
{
    if (endT <= startT)
        return 0.0;

    double length(0.0);

    for (int nt = startT; nt < endT; ++nt)
    {
        const double dx(particle->Vx(nt+1) - particle->Vx(nt));
        const double dy(particle->Vy(nt+1) - particle->Vy(nt));
        const double dz(particle->Vz(nt+1) - particle->Vz(nt));
        length += sqrt(dx * dx + dy * dy + dz * dz);
    }

    return length;
}

} //namespace lar_pandora
