/**
 *  @file   larpandora/LArPandoraInterface/PFParticleHitDumper_module.cc
 *
 *  @brief  Analysis module for created particles
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

// ROOT includes
#include "TTree.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  PFParticleHitDumper class
 */
class PFParticleHitDumper : public art::EDAnalyzer
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
     PFParticleHitDumper(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
     virtual ~PFParticleHitDumper();

     void beginJob();
     void endJob();
     void analyze(const art::Event &evt);
     void reconfigure(fhicl::ParameterSet const &pset);

private:
     TTree       *m_pReco3D;         ///< 
     TTree       *m_pReco2D;         ///<

     int          m_run;             ///< 
     int          m_event;           ///< 
     int          m_particle;        ///<
     int          m_primary;         ///<
     
     int          m_cstat;           ///<
     int          m_tpc;             ///<
     int          m_plane;           ///<
     int          m_wire;            ///<

     double       m_w;               ///<
     double       m_x;               ///<
     double       m_y;               ///<
     double       m_z;               ///<
     double       m_q;               ///<

     std::string  m_hitfinderLabel;  ///<
     std::string  m_spacepointLabel; ///< 
     std::string  m_particleLabel;   ///<
};

DEFINE_ART_MODULE(PFParticleHitDumper)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows
/**
 *  @file   LArPandora/PFParticleHitDumper.cxx
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

// Local includes
#include "LArPandoraCollector.h"

// std includes
#include <iostream>

namespace lar_pandora
{

PFParticleHitDumper::PFParticleHitDumper(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset)
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleHitDumper::~PFParticleHitDumper()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHitDumper::reconfigure(fhicl::ParameterSet const &pset)
{
    m_particleLabel   = pset.get<std::string>("PFParticleModule","pandora");
    m_spacepointLabel = pset.get<std::string>("SpacePointModule","pandora");
    m_hitfinderLabel  = pset.get<std::string>("HitFinderModule","gaushit");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHitDumper::beginJob()
{
    mf::LogDebug("LArPandora") << " *** PFParticleHitDumper::beginJob() *** " << std::endl; 

    // 
    art::ServiceHandle<art::TFileService> tfs;

    m_pReco3D = tfs->make<TTree>("pandora3D", "LAr Reco 3D");
    m_pReco3D->Branch("run", &m_run,"run/I");
    m_pReco3D->Branch("event", &m_event,"event/I");
    m_pReco3D->Branch("particle", &m_particle, "particle/I");
    m_pReco3D->Branch("primary", &m_primary, "primary/I");
    m_pReco3D->Branch("cstat", &m_cstat, "cstat/I");
    m_pReco3D->Branch("tpc", &m_tpc, "tpc/I");
    m_pReco3D->Branch("plane", &m_plane, "plane/I");
    m_pReco3D->Branch("x", &m_x, "x/D");
    m_pReco3D->Branch("y", &m_y, "y/D");
    m_pReco3D->Branch("z", &m_z, "z/D");

    m_pReco2D = tfs->make<TTree>("pandora2D", "LAr Reco 2D");
    m_pReco2D->Branch("run", &m_run,"run/I");
    m_pReco2D->Branch("event", &m_event,"event/I");
    m_pReco2D->Branch("particle", &m_particle, "particle/I");
    m_pReco2D->Branch("cstat", &m_cstat, "cstat/I");
    m_pReco2D->Branch("tpc", &m_tpc, "tpc/I");
    m_pReco2D->Branch("plane", &m_plane, "plane/I");
    m_pReco2D->Branch("wire", &m_wire, "wire/I");
    m_pReco2D->Branch("x", &m_x, "x/D");
    m_pReco2D->Branch("w", &m_w, "w/D");
    m_pReco2D->Branch("q", &m_q, "q/D");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHitDumper::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHitDumper::analyze(const art::Event &evt)
{
    std::cout << " *** PFParticleHitDumper::analyze(...) *** " << std::endl;

    m_run = evt.run();
    m_event = evt.id().event();

    m_particle = -1;
    m_primary = 0;
     
    m_cstat = 0;
    m_tpc = 0;
    m_plane = 0;
    m_wire = 0;

    m_x = 0.0;
    m_y = 0.0;
    m_z = 0.0;
    m_w = 0.0;
    m_q = 0.0;

    std::cout << "  Run: " << m_run << std::endl;
    std::cout << "  Event: " << m_event << std::endl; 


    // Get particles, space points and hits
    // ====================================
    PFParticleVector         particleVector;
    SpacePointVector         spacePointVector;
    HitVector                hitVector;

    SpacePointsToHits        spacePointsToHits;
    PFParticlesToSpacePoints particlesToSpacePoints;
    HitsToPFParticles        hitsToParticles;
    PFParticlesToHits        particlesToHits;

    LArPandoraCollector::CollectHits(evt, m_hitfinderLabel, hitVector);
    LArPandoraCollector::CollectSpacePoints(evt, m_spacepointLabel, spacePointVector, spacePointsToHits);
    LArPandoraCollector::CollectPFParticles(evt, m_particleLabel, particleVector, particlesToSpacePoints);
    LArPandoraCollector::BuildPFParticleHitMaps(particleVector, particlesToSpacePoints, spacePointsToHits, particlesToHits, hitsToParticles);

    std::cout << "  PFParticles (with SpacePoints): " << particleVector.size() << std::endl; 


    // Get geometry and detector properties
    // ====================================
    art::ServiceHandle<geo::Geometry>            theGeometry;
    art::ServiceHandle<util::DetectorProperties> theDetector;


    // Loop over Hits (Fill 2D Tree)
    // =============================
    if (hitVector.empty())
    {
        m_pReco2D->Fill();
    }
    
    for (unsigned int i = 0; i<hitVector.size(); ++i)
    {
        const art::Ptr<recob::Hit> hit = hitVector.at(i);

        m_particle = -1;

        HitsToPFParticles::const_iterator pIter = hitsToParticles.find(hit);
        if (hitsToParticles.end() != pIter)
	{
	    const art::Ptr<recob::PFParticle> particle = pIter->second;
	    m_particle = particle->Self();
        }
                
        const geo::WireID &wireID(hit->WireID());
        m_cstat = wireID.Cryostat;
        m_tpc   = wireID.TPC;
        m_plane = wireID.Plane;
        m_wire  = wireID.Wire; 

        m_q = hit->Charge();
        m_x = theDetector->ConvertTicksToX(hit->PeakTime(), wireID.Plane, wireID.TPC, wireID.Cryostat);

        // define UVW as closest distance from (0,0) to wire
        double xyzStart[3];
        theGeometry->TPC(wireID.TPC).Plane(wireID.Plane).Wire(wireID.Wire).GetStart(xyzStart);
        const double ay(xyzStart[1]);
        const double az(xyzStart[2]);

        double xyzEnd[3];
        theGeometry->TPC(wireID.TPC).Plane(wireID.Plane).Wire(wireID.Wire).GetEnd(xyzEnd);
        const double by(xyzEnd[1]);
        const double bz(xyzEnd[2]);

        const double ny(by - ay);
        const double nz(bz - az);
        const double N2(ny * ny + nz * nz);

        const double ry(ay - (ay * ny + az * nz) * ny / N2);
        const double rz(az - (ay * ny + az * nz) * nz / N2);
        const double sign((rz >0.0) ? +1.0 : -1.0);

        m_w = sign * std::sqrt(ry * ry + rz * rz);

        m_pReco2D->Fill();
    }


    // Loop over PFParticles (Fill 3D Tree)
    // ====================================
    if (particleVector.empty())
    {
        m_pReco3D->Fill();
    }

    // Store associations between particle and particle ID
    PFParticleMap theParticleMap;

    for (unsigned int i = 0; i < particleVector.size(); ++i )
    {
        const art::Ptr<recob::PFParticle> particle = particleVector.at(i);
        theParticleMap[particle->Self()] = particle;
    }

    // Loop over particles
    for (PFParticlesToSpacePoints::const_iterator iter1 = particlesToSpacePoints.begin(), iterEnd1 = particlesToSpacePoints.end();
        iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> particle = iter1->first;
        const SpacePointVector &spacepoints = iter1->second;

        m_particle = particle->Self();
        m_primary = 0;

        if (particle->IsPrimary())
	{
	    m_primary = 1;
	}
        else
	{
	    const size_t parentID(particle->Parent());
	    PFParticleMap::const_iterator pIter = theParticleMap.find(parentID);
            if (theParticleMap.end() == pIter)
                throw cet::exception("LArPandora") << " PFParticleHitDumper::analyze --- Found particle with ID code";

	    const art::Ptr<recob::PFParticle> particleParent = pIter->second;
            if (LArPandoraCollector::IsNeutrino(particleParent))
                m_primary = 1;
	}

        std::cout << "    PFPARTICLE [" << m_particle << "] [Primary=" << m_primary << "] (" << spacepoints.size() << " Space Points)" << std::endl;

        for (unsigned int j=0; j<spacepoints.size(); ++j)
	{
	    const art::Ptr<recob::SpacePoint> spacepoint = spacepoints.at(j);

	    m_x = spacepoint->XYZ()[0];
            m_y = spacepoint->XYZ()[1];
            m_z = spacepoint->XYZ()[2];

	    SpacePointsToHits::const_iterator iter2 = spacePointsToHits.find(spacepoint);
            if (spacePointsToHits.end() == iter2)
                throw cet::exception("LArPandora") << " PFParticleHitDumper::analyze --- Found space point without associated hit";

	    const art::Ptr<recob::Hit> hit = iter2->second;
            const geo::WireID &wireID(hit->WireID());

            m_cstat = wireID.Cryostat;
            m_tpc   = wireID.TPC;
            m_plane = wireID.Plane;

            m_pReco3D->Fill();
	}
    }
}

} //namespace lar_pandora
