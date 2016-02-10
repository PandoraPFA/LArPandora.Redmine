/**
 *  @file   larpandora/LArPandoraAnalysis/PFParticleHitDumper_module.cc
 *
 *  @brief  Analysis module for created particles
 *
 */

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "TTree.h"

#include "LArPandoraInterface/LArPandoraHelper.h"

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

    /**
     *  @brief Store 3D track hits
     *
     *  @param particlesToTracks mapping between 3D track hits and PFParticles
     */
     void FillRecoTracks(const PFParticlesToTracks &particlesToTracks);

    /**
     *  @brief Store 3D hits
     *
     *  @param particleVector the input vector of PFParticles
     *  @param particlesToSpacePoints mapping between 3D hits and PFParticles
     *  @param spacePointsToHits mapping between 3D hits and 2D hits
     */
    void FillReco3D(const PFParticleVector &particleVector, const PFParticlesToSpacePoints &particlesToSpacePoints,
        const SpacePointsToHits &spacePointsToHits);
    
    /** 
     *  @brief Store 2D hits
     *
     *  @param hitVector the input vector of 2D hits
     *  @param hitsToParticles mapping between 2D hits and PFParticles
     */
     void FillReco2D(const HitVector &hitVector, const HitsToPFParticles &hitsToParticles);
    
    /**
     *  @brief Store raw data
     *
     *  @param wireVector the input vector of reconstructed wires
     */
     void FillRecoWires(const WireVector &wireVector);

    /**
     *  @brief Conversion from wire ID to U/V/W coordinate
     *
     *  @param wireID the input wire ID
     */
     double GetUVW(const geo::WireID &wireID) const;
  
    /**
     *  @brief Convert from (Y,Z) to U coordinate
     *
     *  @param cstat the cryostat
     *  @param tpc the tpc
     *  @param y the y coordinate 
     *  @param z the z coordinate
     */
     double YZtoU(const unsigned int cstat, const unsigned int tpc, const double y, const double z) const;

    /**
     *  @brief Convert from (Y,Z) to V coordinate
     *
     *  @param cstat the crystat
     *  @param tpc the tpc
     *  @param y the y coordinate
     *  @param z the z coordinate
     */
     double YZtoV(const unsigned int cstat, const unsigned int tpc, const double y, const double z) const;

     TTree       *m_pRecoTracks;     ///<
     TTree       *m_pReco3D;         ///< 
     TTree       *m_pReco2D;         ///<
     TTree       *m_pRecoWire;       ///<

     int          m_run;             ///< 
     int          m_event;           ///< 
     int          m_particle;        ///<
     int          m_primary;         ///<
     int          m_pdgcode;         ///<
     
     int          m_cstat;           ///<
     int          m_tpc;             ///<
     int          m_plane;           ///<
     int          m_wire;            ///<

     double       m_u;               ///<
     double       m_v;               ///<
     double       m_w;               ///<
     double       m_x;               ///<
     double       m_y;               ///<
     double       m_z;               ///<
     double       m_q;               ///<

     std::string  m_calwireLabel;    ///<
     std::string  m_hitfinderLabel;  ///<
     std::string  m_spacepointLabel; ///< 
     std::string  m_particleLabel;   ///<
     std::string  m_trackLabel;      ///<

     bool         m_storeWires;      ///<
     bool         m_printDebug;      ///< switch for print statements (TODO: use message service!)
};

DEFINE_ART_MODULE(PFParticleHitDumper)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows
/**
 *  @file   larpandora/LArPandoraAnalysis/PFParticleHitDumper.cxx
 * 
 *  @brief  Implementation of the lar pandora analysis producer.
 * 
 *  $Log: $
 */

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Track.h"
#include "RecoBase/Wire.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

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
    m_storeWires      = pset.get<bool>("StoreWires", false);
    m_trackLabel      = pset.get<std::string>("TrackModule", "pandora");
    m_particleLabel   = pset.get<std::string>("PFParticleModule", "pandora");
    m_spacepointLabel = pset.get<std::string>("SpacePointModule", "pandora");
    m_hitfinderLabel  = pset.get<std::string>("HitFinderModule", "gaushit");
    m_calwireLabel    = pset.get<std::string>("CalWireModule", "caldata");
    m_printDebug      = pset.get<bool>("PrintDebug",false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHitDumper::beginJob()
{
    mf::LogDebug("LArPandora") << " *** PFParticleHitDumper::beginJob() *** " << std::endl; 

    // 
    art::ServiceHandle<art::TFileService> tfs;

    m_pRecoTracks = tfs->make<TTree>("pandoraTracks", "LAr Reco Tracks");
    m_pRecoTracks->Branch("run", &m_run,"run/I");
    m_pRecoTracks->Branch("event", &m_event,"event/I");
    m_pRecoTracks->Branch("particle", &m_particle, "particle/I");
    m_pRecoTracks->Branch("x", &m_x, "x/D");
    m_pRecoTracks->Branch("y", &m_y, "y/D");
    m_pRecoTracks->Branch("z", &m_z, "z/D");

    m_pReco3D = tfs->make<TTree>("pandora3D", "LAr Reco 3D");
    m_pReco3D->Branch("run", &m_run,"run/I");
    m_pReco3D->Branch("event", &m_event,"event/I");
    m_pReco3D->Branch("particle", &m_particle, "particle/I");
    m_pReco3D->Branch("primary", &m_primary, "primary/I");
    m_pReco3D->Branch("pdgcode", &m_pdgcode, "pdgcode/I");
    m_pReco3D->Branch("cstat", &m_cstat, "cstat/I");
    m_pReco3D->Branch("tpc", &m_tpc, "tpc/I");
    m_pReco3D->Branch("plane", &m_plane, "plane/I");
    m_pReco3D->Branch("x", &m_x, "x/D");
    m_pReco3D->Branch("y", &m_y, "y/D");
    m_pReco3D->Branch("u", &m_u, "u/D");
    m_pReco3D->Branch("v", &m_v, "v/D");
    m_pReco3D->Branch("z", &m_z, "z/D");

    m_pReco2D = tfs->make<TTree>("pandora2D", "LAr Reco 2D");
    m_pReco2D->Branch("run", &m_run,"run/I");
    m_pReco2D->Branch("event", &m_event,"event/I");
    m_pReco2D->Branch("particle", &m_particle, "particle/I");
    m_pReco2D->Branch("pdgcode", &m_pdgcode, "pdgcode/I");
    m_pReco2D->Branch("cstat", &m_cstat, "cstat/I");
    m_pReco2D->Branch("tpc", &m_tpc, "tpc/I");
    m_pReco2D->Branch("plane", &m_plane, "plane/I");
    m_pReco2D->Branch("wire", &m_wire, "wire/I");
    m_pReco2D->Branch("x", &m_x, "x/D");
    m_pReco2D->Branch("w", &m_w, "w/D");
    m_pReco2D->Branch("q", &m_q, "q/D");

    m_pRecoWire = tfs->make<TTree>("rawdata", "LAr Reco Wires");
    m_pRecoWire->Branch("run", &m_run,"run/I");
    m_pRecoWire->Branch("event", &m_event,"event/I");
    m_pRecoWire->Branch("cstat", &m_cstat, "cstat/I");
    m_pRecoWire->Branch("tpc", &m_tpc, "tpc/I");
    m_pRecoWire->Branch("plane", &m_plane, "plane/I");
    m_pRecoWire->Branch("wire", &m_wire, "wire/I");
    m_pRecoWire->Branch("x", &m_x, "x/D");
    m_pRecoWire->Branch("w", &m_w, "w/D");
    m_pRecoWire->Branch("q", &m_q, "q/D");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHitDumper::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHitDumper::analyze(const art::Event &evt)
{
    if (m_printDebug)
        std::cout << " *** PFParticleHitDumper::analyze(...) *** " << std::endl;

    m_run = evt.run();
    m_event = evt.id().event();

    m_particle = -1;
    m_primary = 0;
    m_pdgcode = 0;
     
    m_cstat = 0;
    m_tpc = 0;
    m_plane = 0;
    m_wire = 0;

    m_x = 0.0;
    m_y = 0.0;
    m_u = 0.0;
    m_v = 0.0;
    m_z = 0.0;
    m_w = 0.0;
    m_q = 0.0;

    if (m_printDebug)
    {
        std::cout << "  Run: " << m_run << std::endl;
        std::cout << "  Event: " << m_event << std::endl; 
    }

    // Need DetectorProperties service to convert from ticks to X
    art::ServiceHandle<util::DetectorProperties> theDetector;

    // Need geometry service to convert channel to wire ID
    art::ServiceHandle<geo::Geometry> theGeometry;

    // Get particles, tracks, space points, hits (and wires)
    // ====================================================
    TrackVector              trackVector;
    PFParticleVector         particleVector;
    SpacePointVector         spacePointVector;
    HitVector                hitVector;
    WireVector               wireVector;

    PFParticlesToTracks      particlesToTracks;
    PFParticlesToSpacePoints particlesToSpacePoints;
    PFParticlesToHits        particlesToHits;
    HitsToPFParticles        hitsToParticles;
    SpacePointsToHits        spacePointsToHits;

    LArPandoraHelper::CollectHits(evt, m_hitfinderLabel, hitVector);
    LArPandoraHelper::CollectSpacePoints(evt, m_spacepointLabel, spacePointVector, spacePointsToHits);
    LArPandoraHelper::CollectTracks(evt, m_trackLabel, trackVector, particlesToTracks);
    LArPandoraHelper::CollectPFParticles(evt, m_particleLabel, particleVector, particlesToSpacePoints);
    LArPandoraHelper::BuildPFParticleHitMaps(evt, m_particleLabel, m_spacepointLabel, particlesToHits, hitsToParticles);

    if (m_storeWires)
        LArPandoraHelper::CollectWires(evt, m_calwireLabel, wireVector);

    if (m_printDebug)
        std::cout << "  PFParticles: " << particleVector.size() << std::endl; 
   
    // Loop over Tracks (Fill 3D Track Tree)
    // =====================================
    if (m_printDebug)
        std::cout << "   PFParticleHitDumper::FillRecoTracks(...) " << std::endl;
    this->FillRecoTracks(particlesToTracks);

    // Loop over PFParticles (Fill 3D Reco Tree)
    // =========================================
    if (m_printDebug)
        std::cout << "   PFParticleHitDumper::FillReco3D(...) " << std::endl;
    this->FillReco3D(particleVector, particlesToSpacePoints, spacePointsToHits);

    // Loop over Hits (Fill 2D Reco Tree)
    // ==================================
    if (m_printDebug)
        std::cout << "   PFParticleHitDumper::FillReco2D(...) " << std::endl;
    this->FillReco2D(hitVector, hitsToParticles);

    // Loop over Wires (Fill Reco Wire Tree)
    // =====================================
    if (m_printDebug)
        std::cout << "   PFParticleHitDumper::FillRecoWires(...) " << std::endl;
    this->FillRecoWires(wireVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHitDumper::FillRecoTracks(const PFParticlesToTracks &particlesToTracks)
{ 
    // Initialise variables
    m_particle = -1;
    m_x = 0.0;
    m_y = 0.0;
    m_z = 0.0;

    // Create dummy entry if there are no particles
    if (particlesToTracks.empty())
    {
        m_pRecoTracks->Fill();
    }

    // Loop over tracks
    for (PFParticlesToTracks::const_iterator iter = particlesToTracks.begin(), iterEnd = particlesToTracks.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = iter->first;
        const TrackVector &trackVector = iter->second;
  
        m_particle = particle->Self();

        if (!trackVector.empty())
        {
            if (trackVector.size() != 1 && m_printDebug)
                std::cout << " Warning: Found particle with more than one associated track " << std::endl;
	        
            const art::Ptr<recob::Track> track = *(trackVector.begin());

            if (m_printDebug)
                std::cout << "    PFPARTICLE [" << m_particle << "] (" << track->NumberTrajectoryPoints() << " Trajectory Points)" << std::endl;

	    for (unsigned int p = 0; p < track->NumberTrajectoryPoints(); ++p)
            {
	        const TVector3 position(track->LocationAtPoint(p));
                m_x = position.x();
                m_y = position.y();
                m_z = position.z();

                m_pRecoTracks->Fill();
	    }
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHitDumper::FillReco3D(const PFParticleVector &particleVector, const PFParticlesToSpacePoints &particlesToSpacePoints,
    const SpacePointsToHits &spacePointsToHits)
{ 
    // Initialise variables
    m_particle = -1;
    m_primary = 0;
    m_pdgcode = 0;
    m_cstat = 0;
    m_tpc = 0;
    m_plane = 0;
    m_x = 0.0;
    m_u = 0.0;
    m_v = 0.0;
    m_y = 0.0;
    m_z = 0.0;

    // Create dummy entry if there are no particles
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
        m_pdgcode = particle->PdgCode();
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
            if (LArPandoraHelper::IsNeutrino(particleParent))
                m_primary = 1;
        }

        if (m_printDebug)
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

            m_u = this->YZtoU(m_cstat, m_tpc, m_y, m_z);
            m_v = this->YZtoV(m_cstat, m_tpc, m_y, m_z);

            m_pReco3D->Fill();
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHitDumper::FillReco2D(const HitVector &hitVector, const HitsToPFParticles &hitsToParticles)
{ 
    // Initialise variables
    m_particle = -1;
    m_pdgcode = 0;
    m_cstat = 0;
    m_tpc = 0;
    m_plane = 0;
    m_wire = 0;
    m_x = 0.0;
    m_w = 0.0;
    m_q = 0.0;

    // Create dummy entry if there are no 2D hits 
    if (hitVector.empty())
    {
        m_pReco2D->Fill();
    }

    // Need DetectorProperties service to convert from ticks to X
    art::ServiceHandle<util::DetectorProperties> theDetector;

    // Loop over 2D hits
    for (unsigned int i = 0; i<hitVector.size(); ++i)
    {
        const art::Ptr<recob::Hit> hit = hitVector.at(i);

        m_particle = -1;
        m_pdgcode = 0;

        HitsToPFParticles::const_iterator pIter = hitsToParticles.find(hit);
        if (hitsToParticles.end() != pIter)
        {
            const art::Ptr<recob::PFParticle> particle = pIter->second;
            m_particle = particle->Self();
            m_pdgcode = particle->PdgCode();
        }
                
        const geo::WireID &wireID(hit->WireID());
        m_cstat = wireID.Cryostat;
        m_tpc   = wireID.TPC;
        m_plane = wireID.Plane;
        m_wire  = wireID.Wire; 

        m_q = hit->Integral();
        m_x = theDetector->ConvertTicksToX(hit->PeakTime(), wireID.Plane, wireID.TPC, wireID.Cryostat);
        m_w = this->GetUVW(wireID);
     
        m_pReco2D->Fill();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleHitDumper::FillRecoWires(const WireVector &wireVector)
{

    // Create dummy entry if there are no wires
    if (wireVector.empty())
    {
        m_pRecoWire->Fill();
    }

    // Need geometry service to convert channel to wire ID
    art::ServiceHandle<geo::Geometry> theGeometry;

    // Need DetectorProperties service to convert from ticks to X
    art::ServiceHandle<util::DetectorProperties> theDetector;

    // Loop over wires
    int signalCounter(0);

    for (unsigned int i = 0; i<wireVector.size(); ++i)
    {
        const art::Ptr<recob::Wire> wire = wireVector.at(i);
       
        const std::vector<float> &signals(wire->Signal());
        const std::vector<geo::WireID> wireIds = theGeometry->ChannelToWire(wire->Channel());

        if ((signalCounter++) < 10 && m_printDebug)
          std::cout << "    numWires=" << wireVector.size() << " numSignals=" << signals.size() << std::endl;

        double time(0.0);

        m_q = 0.0;

        for (std::vector<float>::const_iterator tIter = signals.begin(), tIterEnd = signals.end(); tIter != tIterEnd; ++tIter)
        {
            time += 1.0;
            m_q = *tIter;

            if (m_q < 2.0) // seems to remove most noise
                continue;

            for (std::vector<geo::WireID>::const_iterator wIter = wireIds.begin(), wIterEnd = wireIds.end(); wIter != wIterEnd; ++wIter)
            {
                const geo::WireID &wireID = *wIter; 
                m_cstat = wireID.Cryostat;
                m_tpc   = wireID.TPC;
                m_plane = wireID.Plane;
                m_wire  = wireID.Wire; 

                m_x = theDetector->ConvertTicksToX(time, wireID.Plane, wireID.TPC, wireID.Cryostat);
                m_w = this->GetUVW(wireID);

                m_pRecoWire->Fill();
            }
        }
    } 
}

//------------------------------------------------------------------------------------------------------------------------------------------

double PFParticleHitDumper::GetUVW(const geo::WireID &wireID) const
{
    // define UVW as closest distance from (0,0) to wire axis
    art::ServiceHandle<geo::Geometry> theGeometry;
     
    double xyzStart[3];
    theGeometry->Cryostat(wireID.Cryostat).TPC(wireID.TPC).Plane(wireID.Plane).Wire(wireID.Wire).GetStart(xyzStart);
    const double ay(xyzStart[1]);
    const double az(xyzStart[2]);

    double xyzEnd[3];
    theGeometry->Cryostat(wireID.Cryostat).TPC(wireID.TPC).Plane(wireID.Plane).Wire(wireID.Wire).GetEnd(xyzEnd);
    const double by(xyzEnd[1]);
    const double bz(xyzEnd[2]);

    const double ny(by - ay);
    const double nz(bz - az);
    const double N2(ny * ny + nz * nz);

    const double ry(ay - (ay * ny + az * nz) * ny / N2);
    const double rz(az - (ay * ny + az * nz) * nz / N2);
    const double sign((rz > 0.0) ? +1.0 : -1.0);

    return sign * std::sqrt(ry * ry + rz * rz);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double PFParticleHitDumper::YZtoU(const unsigned int cstat, const unsigned int tpc, const double y, const double z) const
{
    // TODO: Check that this stills works in DUNE
    art::ServiceHandle<geo::Geometry> theGeometry;
    const double m_theta(theGeometry->WireAngleToVertical(geo::kU, tpc, cstat));
    return z * std::sin(m_theta) - y * std::cos(m_theta);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double PFParticleHitDumper::YZtoV(const unsigned int cstat, const unsigned int tpc, const double y, const double z) const
{
    // TODO; Check that this still works in DUNE
    art::ServiceHandle<geo::Geometry> theGeometry;
    const double m_theta(theGeometry->WireAngleToVertical(geo::kV, tpc, cstat));
    return z * std::sin(m_theta) - y * std::cos(m_theta);
}

//------------------------------------------------------------------------------------------------------------------------------------------

} //namespace lar_pandora
