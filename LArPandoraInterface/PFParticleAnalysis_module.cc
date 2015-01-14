/**
 *  @file   larpandora/LArPandoraInterface/PFParticleAnalysis_module.cc
 *
 *  @brief  Analysis module for created particles
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
     int          m_neutrino;              ///<
     int          m_finalstate;            ///<
     int          m_vertex;                ///<

     int          m_clusters;              ///<
     int          m_spacepoints;           ///<

     double       m_vtxx;                  ///< 
     double       m_vtxy;                  ///<
     double       m_vtxz;                  ///<
     double       m_px;                    ///< 
     double       m_py;                    ///<
     double       m_pz;                    ///<
     double       m_ptot;                  ///<

     std::string  m_particleLabel;         ///<
};

DEFINE_ART_MODULE(PFParticleAnalysis)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows
/**
 *  @file   LArPandora/PFParticleAnalysis.cxx
 * 
 *  @brief  Implementation of PFParticle analysis module
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

// std includes
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
    m_pRecoTree->Branch("neutrino", &m_neutrino, "neutrino/I");  
    m_pRecoTree->Branch("finalstate", &m_finalstate, "finalstate/I"); 
    m_pRecoTree->Branch("vertex", &m_vertex, "vertex/I"); 
    m_pRecoTree->Branch("clusters", &m_clusters, "clusters/I");
    m_pRecoTree->Branch("spacepoints", &m_spacepoints, "spacepoints/I"); 
    m_pRecoTree->Branch("vtxx", &m_vtxx, "vtxx/D");
    m_pRecoTree->Branch("vtxy", &m_vtxy, "vtxy/D");
    m_pRecoTree->Branch("vtxz", &m_vtxz, "vtxz/D");
    m_pRecoTree->Branch("px", &m_px, "px/D");
    m_pRecoTree->Branch("py", &m_py, "py/D");
    m_pRecoTree->Branch("pz", &m_pz, "pz/D");
    m_pRecoTree->Branch("ptot", &m_ptot, "ptot/D");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleAnalysis::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleAnalysis::analyze(const art::Event &evt)
{
    std::cout << " *** PFParticleAnalysis::analyze(...) *** " << std::endl;

    m_run = evt.run();
    m_event = evt.id().event();
    m_index = 0;

    m_self = 0;
    m_pdgcode = 0;
    m_primary = 0;
    m_parent = 0;
    m_neutrino = 0;
    m_finalstate = 0;
    m_vertex = 0;

    m_clusters = 0;
    m_spacepoints = 0;
      
    m_vtxx = 0.0;
    m_vtxy = 0.0;
    m_vtxz = 0.0;

    m_px = 0.0;
    m_py = 0.0;
    m_pz = 0.0;
    m_ptot = 0.0;

    std::cout << "  Run: " << m_run << std::endl;
    std::cout << "  Event: " << m_event << std::endl; 

    // Get the reconstructed PFParticles
    // =================================
    PFParticleVector particleVector;
    PFParticleVector particles1, particles2;
    PFParticlesToClusters particlesToClusters;
    PFParticlesToSpacePoints particlesToSpacePoints;

    LArPandoraCollector::CollectPFParticles(evt, m_particleLabel, particleVector);
    LArPandoraCollector::CollectPFParticles(evt, m_particleLabel, particles1, particlesToClusters);
    LArPandoraCollector::CollectPFParticles(evt, m_particleLabel, particles2, particlesToSpacePoints);
    
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
    LArPandoraCollector::CollectVertices(evt, m_particleLabel, vertexVector, particlesToVertices);

    // Get the reconstructed seeds
    // ===========================
    SeedVector seedVector;
    PFParticlesToSeeds particlesToSeeds;
    LArPandoraCollector::CollectSeeds(evt, m_particleLabel, seedVector, particlesToSeeds);

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
        m_neutrino = LArPandoraCollector::GetParentNeutrino(particleMap, particle);
        m_finalstate = LArPandoraCollector::IsFinalState(particleMap, particle);
        m_vertex = 0;

        m_clusters = 0;
        m_spacepoints = 0;
      
        m_vtxx = 0.0;
        m_vtxy = 0.0;
        m_vtxz = 0.0;
     
        m_px = 0.0;
        m_py = 0.0;
        m_pz = 0.0;
        m_ptot = 0.0;

        PFParticlesToClusters::const_iterator iter1 = particlesToClusters.find(particle);
        if (particlesToClusters.end() != iter1)
            m_clusters = iter1->second.size();

        PFParticlesToSpacePoints::const_iterator iter2 = particlesToSpacePoints.find(particle);
        if (particlesToSpacePoints.end() != iter2)
            m_spacepoints = iter2->second.size();

        PFParticlesToVertices::const_iterator iter3 = particlesToVertices.find(particle);
        if (particlesToVertices.end() != iter3)
        {
            const VertexVector &vertexVector = iter3->second;
            if (!vertexVector.empty())
            {
                if (vertexVector.size() !=1 )
                    std::cout << " Warning: Found particle with more than one associated vertex " << std::endl;

                const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
                double xyz[3] = {0.0, 0.0, 0.0} ;
                vertex->XYZ(xyz);

                m_vertex = 1;
                m_vtxx = xyz[0];
                m_vtxy = xyz[1];
                m_vtxz = xyz[2];
            }
        }

        PFParticlesToSeeds::const_iterator iter4 = particlesToSeeds.find(particle);
        if (particlesToSeeds.end() != iter4)
        {
            const SeedVector &seedVector = iter4->second;
            if (!seedVector.empty())
            {
                if (seedVector.size() !=1 )
                  std::cout << " Warning: Found particle with more than one associated seed " << std::endl;

                const art::Ptr<recob::Seed> seed = *(seedVector.begin());
                double pxpypz[3] = {0.0, 0.0, 0.0} ;
                double err[3] = {0.0, 0.0, 0.0} ;
                seed->GetDirection(pxpypz, err);

                m_px = pxpypz[0];
                m_py = pxpypz[1];
                m_pz = pxpypz[2];
                m_ptot = std::sqrt(m_px * m_px + m_py * m_py + m_pz * m_pz);
            }
        }

        std::cout << "    PFParticle [" << n << "] Primary=" << m_primary << " FinalState=" << m_finalstate 
                  << " Pdg=" << m_pdgcode << " NuPdg=" << m_neutrino
                  << " (Self=" << m_self << ", Parent=" << m_parent << ")"
                  << " (Clusters=" << m_clusters << ", SpacePoints=" << m_spacepoints << ") " << std::endl;

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
