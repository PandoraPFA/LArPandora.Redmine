/**
 *  @file   larpandora/LArPandoraAnalysis/PFParticleAnalysis_module.cc
 *
 *  @brief  Analysis module for created particles
 */

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "TTree.h"
#include "TVector3.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  PFParticleAnalysis class
 */
class PFParticleAnalysis : public art::EDAnalyzer
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
     PFParticleAnalysis(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
     virtual ~PFParticleAnalysis();

     void beginJob();
     void endJob();
     void analyze(const art::Event &evt);
     void reconfigure(fhicl::ParameterSet const &pset);

private:

    /**
     *  @brief Build particle maps for reconstructed particles
     *
     *  @param particleVector the input vector of reconstructed particles
     *  @param particleMap the output mapping between reconstructed particles and particle ID
     */
     void BuildParticleMap(const PFParticleVector &particleVector, PFParticleMap &particleMap) const;

     TTree       *m_pRecoTree;             ///< 

     int          m_run;                   ///< 
     int          m_event;                 ///< 
     int          m_index;                 ///<

     int          m_self;                  ///<
     int          m_pdgcode;               ///<
     int          m_primary;               ///<
     int          m_parent;                ///<
     int          m_daughters;             ///<
     int          m_generation;            ///<
     int          m_neutrino;              ///<
     int          m_finalstate;            ///<
     int          m_vertex;                ///<
     int          m_track;                 ///<
     int          m_trackid;               ///<
     int          m_shower;                ///<
     int          m_showerid;              ///<

     int          m_seeds;                 ///<
     int          m_clusters;              ///<
     int          m_spacepoints;           ///<
     int          m_hits;                  ///<
     int          m_trajectorypoints;      ///<
     int          m_trackhits;             ///<
     int          m_showerhits;            ///<
     int          m_seedhits;              ///<

     double       m_seedvtxx;               ///< 
     double       m_seedvtxy;               ///<
     double       m_seedvtxz;               ///<
     double       m_seedvtxdirx;            ///< 
     double       m_seedvtxdiry;            ///<
     double       m_seedvtxdirz;            ///<

     double       m_pfovtxx;               ///< 
     double       m_pfovtxy;               ///<
     double       m_pfovtxz;               ///<

     double       m_trkvtxx;               ///< 
     double       m_trkvtxy;               ///<
     double       m_trkvtxz;               ///<
     double       m_trkvtxdirx;            ///< 
     double       m_trkvtxdiry;            ///<
     double       m_trkvtxdirz;            ///<
     double       m_trkendx;               ///< 
     double       m_trkendy;               ///<
     double       m_trkendz;               ///<
     double       m_trkenddirx;            ///< 
     double       m_trkenddiry;            ///<
     double       m_trkenddirz;            ///<
     double       m_trklength;             ///<
     double       m_trkstraightlength;     ///<
      
     double       m_shwvtxx;               ///< 
     double       m_shwvtxy;               ///<
     double       m_shwvtxz;               ///<
     double       m_shwvtxdirx;            ///< 
     double       m_shwvtxdiry;            ///<
     double       m_shwvtxdirz;            ///<
     double       m_shwlength;             ///<
     double       m_shwopenangle;          ///<
     double       m_shwbestplane;          ///<

     double       m_t0;                    ///<

     std::string  m_particleLabel;         ///<
     bool         m_printDebug;            ///< switch for print statements (TODO: use message service!)
};

DEFINE_ART_MODULE(PFParticleAnalysis)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "lardataobj/AnalysisBase/T0.h"

#include <iostream>

namespace lar_pandora
{

PFParticleAnalysis::PFParticleAnalysis(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset)
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleAnalysis::~PFParticleAnalysis()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleAnalysis::reconfigure(fhicl::ParameterSet const &pset)
{
    m_particleLabel = pset.get<std::string>("PFParticleModule","pandora");
    m_printDebug = pset.get<bool>("PrintDebug",false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleAnalysis::beginJob()
{
    mf::LogDebug("LArPandora") << " *** PFParticleAnalysis::beginJob() *** " << std::endl; 

    // 
    art::ServiceHandle<art::TFileService> tfs;

    m_pRecoTree = tfs->make<TTree>("pandora", "LAr PFParticles");
    m_pRecoTree->Branch("run", &m_run, "run/I");
    m_pRecoTree->Branch("event", &m_event, "event/I");
    m_pRecoTree->Branch("index", &m_index, "index/I");
    m_pRecoTree->Branch("self", &m_self, "self/I");
    m_pRecoTree->Branch("pdgcode", &m_pdgcode, "pdgcode/I");
    m_pRecoTree->Branch("primary", &m_primary, "primary/I");
    m_pRecoTree->Branch("parent", &m_parent, "parent/I");  
    m_pRecoTree->Branch("daughters", &m_daughters, "daughters/I"); 
    m_pRecoTree->Branch("generation", &m_generation, "generation/I"); 
    m_pRecoTree->Branch("neutrino", &m_neutrino, "neutrino/I");  
    m_pRecoTree->Branch("finalstate", &m_finalstate, "finalstate/I"); 
    m_pRecoTree->Branch("vertex", &m_vertex, "vertex/I"); 
    m_pRecoTree->Branch("track", &m_track, "track/I"); 
    m_pRecoTree->Branch("trackid", &m_trackid, "trackid/I"); 
    m_pRecoTree->Branch("shower", &m_shower, "shower/I"); 
    m_pRecoTree->Branch("showerid", &m_showerid, "showerid/I"); 
    m_pRecoTree->Branch("seeds", &m_seeds, "seeds/I");
    m_pRecoTree->Branch("clusters", &m_clusters, "clusters/I");
    m_pRecoTree->Branch("spacepoints", &m_spacepoints, "spacepoints/I"); 
    m_pRecoTree->Branch("hits", &m_hits, "hits/I"); 
    m_pRecoTree->Branch("trackhits", &m_trackhits, "trackhits/I");
    m_pRecoTree->Branch("trajectorypoints", &m_trajectorypoints, "trajectorypoints/I");
    m_pRecoTree->Branch("seedhits", &m_seedhits, "seedhits/I");   
    m_pRecoTree->Branch("showerhits", &m_showerhits, "showerhits/I");  
    m_pRecoTree->Branch("seedvtxx", &m_seedvtxx, "seedvtxx/D");
    m_pRecoTree->Branch("seedvtxy", &m_seedvtxy, "seedvtxy/D");
    m_pRecoTree->Branch("seedvtxz", &m_seedvtxz, "seedvtxz/D");
    m_pRecoTree->Branch("seedvtxdirx", &m_seedvtxdirx, "seedvtxdirx/D");
    m_pRecoTree->Branch("seedvtxdiry", &m_seedvtxdiry, "seedvtxdiry/D");
    m_pRecoTree->Branch("seedvtxdirz", &m_seedvtxdirz, "seedvtxdirz/D");
    m_pRecoTree->Branch("pfovtxx", &m_pfovtxx, "pfovtxx/D");
    m_pRecoTree->Branch("pfovtxy", &m_pfovtxy, "pfovtxy/D");
    m_pRecoTree->Branch("pfovtxz", &m_pfovtxz, "pfovtxz/D");
    m_pRecoTree->Branch("trkvtxx", &m_trkvtxx, "trkvtxx/D");
    m_pRecoTree->Branch("trkvtxy", &m_trkvtxy, "trkvtxy/D");
    m_pRecoTree->Branch("trkvtxz", &m_trkvtxz, "trkvtxz/D");
    m_pRecoTree->Branch("trkvtxdirx", &m_trkvtxdirx, "trkvtxdirx/D");
    m_pRecoTree->Branch("trkvtxdiry", &m_trkvtxdiry, "trkvtxdiry/D");
    m_pRecoTree->Branch("trkvtxdirz", &m_trkvtxdirz, "trkvtxdirz/D");
    m_pRecoTree->Branch("trkendx", &m_trkendx, "trkendx/D");
    m_pRecoTree->Branch("trkendy", &m_trkendy, "trkendy/D");
    m_pRecoTree->Branch("trkendz", &m_trkendz, "trkendz/D");
    m_pRecoTree->Branch("trkenddirx", &m_trkenddirx, "trkenddirx/D");
    m_pRecoTree->Branch("trkenddiry", &m_trkenddiry, "trkenddiry/D");
    m_pRecoTree->Branch("trkenddirz", &m_trkenddirz, "trkenddirz/D");
    m_pRecoTree->Branch("trklength", &m_trklength, "trklength/D");
    m_pRecoTree->Branch("trkstraightlength", &m_trkstraightlength, "trkstraightlength/D");   
    m_pRecoTree->Branch("shwvtxx", &m_shwvtxx, "shwvtxx/D");
    m_pRecoTree->Branch("shwvtxy", &m_shwvtxy, "shwvtxy/D");
    m_pRecoTree->Branch("shwvtxz", &m_shwvtxz, "shwvtxz/D");
    m_pRecoTree->Branch("shwvtxdirx", &m_shwvtxdirx, "shwvtxdirx/D");
    m_pRecoTree->Branch("shwvtxdiry", &m_shwvtxdiry, "shwvtxdiry/D");
    m_pRecoTree->Branch("shwvtxdirz", &m_shwvtxdirz, "shwvtxdirz/D");    
    m_pRecoTree->Branch("shwlength", &m_shwlength, "shwlength/D");
    m_pRecoTree->Branch("shwopenangle", &m_shwopenangle, "shwopenangle/D");
    m_pRecoTree->Branch("shwbestplane", &m_shwbestplane, "shwbestplane/D");
    m_pRecoTree->Branch("t0", &m_t0, "t0/D");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleAnalysis::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleAnalysis::analyze(const art::Event &evt)
{
    if (m_printDebug)
        std::cout << " *** PFParticleAnalysis::analyze(...) *** " << std::endl;

    m_run = evt.run();
    m_event = evt.id().event();
    m_index = 0;

    m_self = 0;
    m_pdgcode = 0;
    m_primary = 0;
    m_parent = 0;
    m_daughters = 0;
    m_generation = 0;
    m_neutrino = 0;
    m_finalstate = 0;
    m_vertex = 0;
    m_track = 0;
    m_trackid = -999;
    m_shower = 0;
    m_showerid = -999;

    m_seeds = 0;
    m_clusters = 0;
    m_spacepoints = 0;
    m_hits = 0;
    m_trajectorypoints = 0;
    m_trackhits = 0;
    m_seedhits = 0;
    m_showerhits = 0;
      
    m_seedvtxx = 0.0;
    m_seedvtxy = 0.0;
    m_seedvtxz = 0.0;
    m_seedvtxdirx = 0.0;
    m_seedvtxdiry = 0.0;
    m_seedvtxdirz = 0.0;

    m_pfovtxx = 0.0;
    m_pfovtxy = 0.0;
    m_pfovtxz = 0.0;

    m_trkvtxx = 0.0;
    m_trkvtxy = 0.0;
    m_trkvtxz = 0.0;
    m_trkvtxdirx = 0.0;
    m_trkvtxdiry = 0.0;
    m_trkvtxdirz = 0.0;
    m_trkendx = 0.0;
    m_trkendy = 0.0;
    m_trkendz = 0.0;
    m_trkenddirx = 0.0;
    m_trkenddiry = 0.0;
    m_trkenddirz = 0.0;
    m_trklength = 0.0;
    m_trkstraightlength = 0.0;
     
    m_shwvtxx = 0.0;
    m_shwvtxy = 0.0;
    m_shwvtxz = 0.0;
    m_shwvtxdirx = 0.0;
    m_shwvtxdiry = 0.0;
    m_shwvtxdirz = 0.0;
    m_shwlength = 0.0;
    m_shwopenangle = 0.0;
    m_shwbestplane = 0.0;

    m_t0 = 0;

    if (m_printDebug)
    {
        std::cout << "  Run: " << m_run << std::endl;
        std::cout << "  Event: " << m_event << std::endl; 
    }

    // Get the reconstructed PFParticles
    // =================================
    PFParticleVector particleVector;
    PFParticleVector particles1, particles2;
    PFParticlesToClusters particlesToClusters;
    PFParticlesToSpacePoints particlesToSpacePoints;
    PFParticlesToHits particlesToHits;
    HitsToPFParticles hitsToParticles;

    LArPandoraHelper::CollectPFParticles(evt, m_particleLabel, particleVector);
    LArPandoraHelper::CollectPFParticles(evt, m_particleLabel, particles1, particlesToClusters);
    LArPandoraHelper::CollectPFParticles(evt, m_particleLabel, particles2, particlesToSpacePoints);
    LArPandoraHelper::BuildPFParticleHitMaps(evt, m_particleLabel, particlesToHits, hitsToParticles);

    if (m_printDebug)
        std::cout << "  PFParticles: " << particleVector.size() << std::endl;

    if (particleVector.empty())
    {
        m_pRecoTree->Fill();
        return;
    }

    // Get the reconstructed vertices
    // ==============================
    VertexVector vertexVector;
    PFParticlesToVertices particlesToVertices;
    LArPandoraHelper::CollectVertices(evt, m_particleLabel, vertexVector, particlesToVertices);

    // Get the reconstructed seeds
    // ===========================
    SeedVector seedVector, seedVector2;
    PFParticlesToSeeds particlesToSeeds;
    SeedsToHits seedsToHits;
    LArPandoraHelper::CollectSeeds(evt, m_particleLabel, seedVector, particlesToSeeds);
    LArPandoraHelper::CollectSeeds(evt, m_particleLabel, seedVector2, seedsToHits);

    // Get the reconstructed tracks
    // ============================
    TrackVector trackVector, trackVector2;
    PFParticlesToTracks particlesToTracks;
    TracksToHits tracksToHits;
    LArPandoraHelper::CollectTracks(evt, m_particleLabel, trackVector, particlesToTracks);
    LArPandoraHelper::CollectTracks(evt, m_particleLabel, trackVector2, tracksToHits);    

    // Get the reconstructed showers
    // ============================
    ShowerVector showerVector, showerVector2;
    PFParticlesToShowers particlesToShowers;
    ShowersToHits showersToHits;
    LArPandoraHelper::CollectShowers(evt, m_particleLabel, showerVector, particlesToShowers);
    LArPandoraHelper::CollectShowers(evt, m_particleLabel, showerVector2, showersToHits);

    // Get the reconstructed T0 objects
    // ================================
    T0Vector t0Vector;
    TracksToT0s tracksToT0s;
    LArPandoraHelper::CollectT0s(evt, m_particleLabel, t0Vector, tracksToT0s);

    // Build an indexed map of the PFParticles
    // =======================================
    PFParticleMap particleMap;
    this->BuildParticleMap(particleVector, particleMap);


    // Write PFParticle properties to ROOT file
    // ========================================
    for (unsigned int n = 0; n < particleVector.size(); ++n)
    {
        const art::Ptr<recob::PFParticle> particle = particleVector.at(n);

        m_index = n;
        m_self = particle->Self();
        m_pdgcode = particle->PdgCode();
        m_primary = particle->IsPrimary();
        m_parent = (particle->IsPrimary() ? -1 : particle->Parent());
        m_daughters = particle->NumDaughters();
        m_generation = LArPandoraHelper::GetGeneration(particleMap, particle);
        m_neutrino = LArPandoraHelper::GetParentNeutrino(particleMap, particle);
        m_finalstate = LArPandoraHelper::IsFinalState(particleMap, particle);
        m_vertex = 0;
        m_track = 0;
        m_trackid = -999;
        m_shower = 0;
        m_showerid = -999;

        m_seeds = 0;
        m_clusters = 0;
        m_spacepoints = 0;
        m_hits = 0;
        m_trajectorypoints = 0;
        m_trackhits = 0;
        m_seedhits = 0;
        m_showerhits = 0;
              
        m_seedvtxx = 0.0;
        m_seedvtxy = 0.0;
        m_seedvtxz = 0.0;
        m_seedvtxdirx = 0.0;
        m_seedvtxdiry = 0.0;
        m_seedvtxdirz = 0.0;

        m_pfovtxx = 0.0;
        m_pfovtxy = 0.0;
        m_pfovtxz = 0.0;

        m_trkvtxx = 0.0;
        m_trkvtxy = 0.0;
        m_trkvtxz = 0.0;
        m_trkvtxdirx = 0.0;
        m_trkvtxdiry = 0.0;
        m_trkvtxdirz = 0.0;
        m_trkendx = 0.0;
        m_trkendy = 0.0;
        m_trkendz = 0.0;
        m_trkenddirx = 0.0;
        m_trkenddiry = 0.0;
        m_trkenddirz = 0.0;
        m_trklength = 0.0;
        m_trkstraightlength = 0.0;
	      
        m_shwvtxx = 0.0;
        m_shwvtxy = 0.0;
        m_shwvtxz = 0.0;
        m_shwvtxdirx = 0.0;
        m_shwvtxdiry = 0.0;
        m_shwvtxdirz = 0.0;
        m_shwlength = 0.0;
        m_shwopenangle = 0.0;
        m_shwbestplane = 0.0;

        m_t0 = 0.0;

        // Particles <-> Clusters
        PFParticlesToClusters::const_iterator cIter = particlesToClusters.find(particle);
        if (particlesToClusters.end() != cIter)
            m_clusters = cIter->second.size();

        // Particles <-> SpacePoints
        PFParticlesToSpacePoints::const_iterator pIter = particlesToSpacePoints.find(particle);
        if (particlesToSpacePoints.end() != pIter)
            m_spacepoints = pIter->second.size();

        // Particles <-> Hits
        PFParticlesToHits::const_iterator hIter = particlesToHits.find(particle);
        if (particlesToHits.end() != hIter)
            m_hits = hIter->second.size();

        // Particles <-> Vertices
        PFParticlesToVertices::const_iterator vIter = particlesToVertices.find(particle);
        if (particlesToVertices.end() != vIter)
        {
            const VertexVector &vertexVector = vIter->second;
            if (!vertexVector.empty())
            {
                if (vertexVector.size() !=1 && m_printDebug)
                    std::cout << " Warning: Found particle with more than one associated vertex " << std::endl;

                const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
                double xyz[3] = {0.0, 0.0, 0.0} ;
                vertex->XYZ(xyz);

                m_vertex  = 1;
                m_pfovtxx = xyz[0];
                m_pfovtxy = xyz[1];
                m_pfovtxz = xyz[2];
            }
        }

        // Particles <-> Seeds <-> Hits
        PFParticlesToSeeds::const_iterator seedIter = particlesToSeeds.find(particle);
        if (particlesToSeeds.end() != seedIter)
        {
            const SeedVector &seedVector = seedIter->second;
            m_seeds = seedVector.size();

            if (!seedVector.empty())
            {
                const art::Ptr<recob::Seed> firstSeed = *(seedVector.begin());
                double vxvyvz[3] = {0.0, 0.0, 0.0} ;
                double pxpypz[3] = {0.0, 0.0, 0.0} ;
                double err[3] = {0.0, 0.0, 0.0} ;
                firstSeed->GetPoint(vxvyvz, err);
                firstSeed->GetDirection(pxpypz, err);

                m_seedvtxx = vxvyvz[0];
                m_seedvtxy = vxvyvz[1];
                m_seedvtxz = vxvyvz[2];
                m_seedvtxdirx = pxpypz[0];
                m_seedvtxdiry = pxpypz[1];
                m_seedvtxdirz = pxpypz[2];

                for (SeedVector::const_iterator seedIter1 = seedVector.begin(), seedIterEnd1 = seedVector.end(); seedIter1 != seedIterEnd1; ++seedIter1)
                {
                    SeedsToHits::const_iterator seedIter2 = seedsToHits.find(*seedIter1);
                    if (seedsToHits.end() != seedIter2)
                        ++m_seedhits;
                }
            }
        }

        // Particles <-> Tracks <-> Hits, T0s
        PFParticlesToTracks::const_iterator trkIter = particlesToTracks.find(particle);
        if (particlesToTracks.end() != trkIter)
        {
            const TrackVector &trackVector = trkIter->second;
            if (!trackVector.empty())
            {
                if (trackVector.size() !=1 && m_printDebug)
                    std::cout << " Warning: Found particle with more than one associated track " << std::endl;
 
                const art::Ptr<recob::Track> track = *(trackVector.begin());
                const TVector3 &trackVtxPosition = track->Vertex();
                const TVector3 &trackVtxDirection = track->VertexDirection();
                const TVector3 &trackEndPosition = track->End();
                const TVector3 &trackEndDirection = track->EndDirection();
		
                m_track = 1;
                m_trackid = track->ID();
                m_trajectorypoints = track->NumberTrajectoryPoints();
                m_trkvtxx = trackVtxPosition.x();
                m_trkvtxy = trackVtxPosition.y();
                m_trkvtxz = trackVtxPosition.z();
                m_trkvtxdirx = trackVtxDirection.x();
                m_trkvtxdiry = trackVtxDirection.y();
                m_trkvtxdirz = trackVtxDirection.z();
		m_trkendx = trackEndPosition.x();
                m_trkendy = trackEndPosition.y();
                m_trkendz = trackEndPosition.z();
                m_trkenddirx = trackEndDirection.x();
                m_trkenddiry = trackEndDirection.y();
                m_trkenddirz = trackEndDirection.z();
                m_trklength = track->Length();
                m_trkstraightlength = (trackEndPosition - trackVtxPosition).Mag();

	        TracksToT0s::const_iterator t0Iter = tracksToT0s.find(track);
                if (tracksToT0s.end() != t0Iter)
                {
	            const art::Ptr<anab::T0> t0 = t0Iter->second;
                    m_t0 = t0->Time();
		}

                TracksToHits::const_iterator trkIter2 = tracksToHits.find(track);
                if (tracksToHits.end() != trkIter2)
		    m_trackhits = trkIter2->second.size();
            }
        }

        // Particles <-> Showers <-> Hits
        PFParticlesToShowers::const_iterator shwIter = particlesToShowers.find(particle);
        if (particlesToShowers.end() != shwIter)
        {
            const ShowerVector &showerVector = shwIter->second;
            if (!showerVector.empty())
            {
                if (showerVector.size() !=1 && m_printDebug)
                    std::cout << " Warning: Found particle with more than one associated shower " << std::endl;
 
                const art::Ptr<recob::Shower> shower = *(showerVector.begin());
                const TVector3 &showerVtxPosition = shower->ShowerStart();
                const TVector3 &showerVtxDirection = shower->Direction();

                m_shower = 1;
                m_showerid = shower->ID();
 
                m_shwvtxx = showerVtxPosition.x();
                m_shwvtxy = showerVtxPosition.y();
                m_shwvtxz = showerVtxPosition.z();
                m_shwvtxdirx = showerVtxDirection.x();
                m_shwvtxdiry = showerVtxDirection.y();
                m_shwvtxdirz = showerVtxDirection.z();

                m_shwlength = shower->Length();
                m_shwopenangle = shower->OpenAngle();
                m_shwbestplane = shower->best_plane();

                ShowersToHits::const_iterator shwIter2 = showersToHits.find(shower);
                if (showersToHits.end() != shwIter2)
		    m_showerhits = shwIter2->second.size();
            }
        }

        if (m_printDebug)
            std::cout << "    PFParticle [" << n << "] Primary=" << m_primary << " FinalState=" << m_finalstate 
                      << " Pdg=" << m_pdgcode << " NuPdg=" << m_neutrino
                      << " (Self=" << m_self << ", Parent=" << m_parent << ")"
                      << " (Vertex=" << m_vertex << ", Seeds=" << m_seeds << ", Track=" << m_track << ", Shower=" << m_shower
                      << ", Clusters=" << m_clusters << ", SpacePoints=" << m_spacepoints << ", Hits=" << m_hits << ") " << std::endl;

        m_pRecoTree->Fill();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleAnalysis::BuildParticleMap(const PFParticleVector &particleVector, PFParticleMap &particleMap) const
{
    for (PFParticleVector::const_iterator iter = particleVector.begin(), iterEnd = particleVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = *iter;
        particleMap[particle->Self()] = particle;
    }
}

} //namespace lar_pandora
