/**
 *  @file   LArPandora/LArPandora_module.cc
 * 
 *  @brief  lar pandora module.
 * 
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

// LArSoft includes
#include "MCCheater/BackTracker.h"
#include "RecoBase/Hit.h"

// Pandora includes
#include "Api/PandoraApi.h"

// ROOT includes
#include "TTree.h"

// std includes
#include <string>

namespace pandora {class Pandora;}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  LArPandora class
 */
class LArPandora : public art::EDProducer
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    LArPandora(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~LArPandora();

    void beginJob();
    void endJob();
    void produce(art::Event &evt);
    void reconfigure(fhicl::ParameterSet const &pset);

private:
    typedef std::vector< art::Ptr<simb::MCParticle> > ParticleVector;
    typedef std::vector< art::Ptr<recob::Hit> >  HitVector;
    typedef std::map< art::Ptr<recob::Hit>, std::vector<cheat::TrackIDE> > HitToParticleMap;
    typedef std::map< art::Ptr<simb::MCTruth>, std::vector<int> > TruthToParticleMap;
    typedef std::map< int, art::Ptr<simb::MCParticle> > ParticleMap;
    typedef std::map< int, art::Ptr<recob::Hit> > HitMap;

    /**
     *  @brief  Event Preparation 
     * 
     *  @param  evt  the ART event 
     */
    void PrepareEvent(const art::Event &evt);  

    /**
     *  @brief Register the Pandora algorithms, helper functions and geometry
     */
    void InitializePandora() const;

    /**
     *  @brief  Extract the ART hits and the ART hit-particle relationships
     * 
     *  @param  evt  the ART event 
     *  @param  hits  the ART hits for this event
     *  @param  hitToParticleMap  mapping from each ART hit to its underlying G4 track ID
     */
    void CollectArtHits(const art::Event &evt, HitVector &hits, HitToParticleMap &hitToParticleMap) const;

    /**
     *  @brief Extract the ART MC particles
     * 
     *  @param evt  the ART event 
     *  @param particleMap  mapping from each G4 track ID to its corresponding ART MC particle
     *  @param truthMap  mapping from each G4 track ID to its corresponding ART MC primary particle
     */
    void CollectArtParticles(const art::Event &evt, ParticleMap &particleMap, TruthToParticleMap &truthToParticleMap) const;

    /**
     *  @brief Create the Pandora hits from the ART hits
     * 
     *  @param hits  the ART hits for this event
     *  @param hitMap  mapping from Pandora hit ID to ART hit 
     */
    void CreatePandoraHits(const HitVector &hits, HitMap &hitMap) const;

    /**
     *  @brief Create the Pandora MC particles from the MC particles
     * 
     *  @param particleMap  mapping from each G4 track ID to its corresponding ART MC particle 
     *  @param truthMap  mapping from each G4 track ID to its corresponding ART MC primary particle
     */
    void CreatePandoraParticles(const ParticleMap &particleMap, const TruthToParticleMap &truthToParticleMap) const;

    /**
     *  @brief Create the Pandora hit-particle links
     * 
     *  @param hitMap  mapping from Pandora hit ID to ART hit  
     *  @param hitToParticleMap  mapping from each ART hit to its underlying G4 track ID
     */
    void CreatePandoraLinks(const HitMap &hitMap, const HitToParticleMap &hitToParticleMap) const;

    /**
     *  @brief Process the event using Pandora
     */
    void RunPandora() const;

    /**
     *  @brief Convert the Pandora PFOs into ART clusters and write into ART event
     * 
     *  @param evt  the ART event  
     *  @param hitMap  mapping from Pandora hit ID to ART hit   
     */
    void ProduceArtClusters(art::Event &evt, const HitMap &hitMap) const;

    /**
     *  @brief Reset Pandora
     */
    void ResetPandora() const; 

    /**
     *  @brief Initialize the internal monitoring
     */
    void InitializeMonitoring();

    /**
     *  @brief Loop over MC trajectory points and identify start and end points within detector
     *
     *  @param particle  the true particle
     *  @param startT  the first trajectory point in the detector
     *  @param endT  the last trajectory point in the detector
     */
    void GetStartAndEndPoints(const art::Ptr<simb::MCParticle> &particle, int &startT, int& endT) const;



    bool               m_enableProduction;      ///< 
    bool               m_enableMCParticles;     ///< 
    bool               m_enableMonitoring;      ///<
    std::string        m_configFile;            ///< 
    std::string        m_geantModuleLabel;      ///< 
    std::string        m_hitfinderModuleLabel;  ///< 
    pandora::Pandora  *m_pPandora;              ///< 

    TTree             *m_pRecoTree;             ///< 
    int                m_run;                   ///< 
    int                m_event;                 ///< 
    int                m_hits;                  ///<
    float              m_time;                  ///< 

};

DEFINE_ART_MODULE(LArPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

// Framework includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

// LArSoft includes 
#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationBase/MCTruth.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"

// Pandora includes
#include "Api/PandoraApi.h"
#include "Objects/ParticleFlowObject.h"

// Pandora LArContent includes
#include "LArContent.h"

// Local includes
#include "MicroBooNEPseudoLayerCalculator.h"
#include "MicroBooNETransformationCalculator.h"

// System includes
#include <iostream>

namespace lar_pandora {

LArPandora::LArPandora(fhicl::ParameterSet const &pset)
{
    this->reconfigure(pset);
    m_pPandora = new pandora::Pandora();

    produces< std::vector<recob::Cluster> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPandora::~LArPandora()
{
    delete m_pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::reconfigure(fhicl::ParameterSet const &pset)
{
    m_enableProduction = pset.get<bool>("EnableProduction",true);
    m_enableMCParticles = pset.get<bool>("EnableMCParticles",false);
    m_enableMonitoring = pset.get<bool>("EnableMonitoring",false);

    m_configFile = pset.get<std::string>("ConfigFile");
    m_geantModuleLabel = pset.get<std::string>("GeantModuleLabel","largeant");
    m_hitfinderModuleLabel = pset.get<std::string>("HitFinderModuleLabel","ffthit");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::beginJob()
{
    this->InitializePandora();  

    if (m_enableMonitoring)
        this->InitializeMonitoring();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::produce(art::Event &evt)
{ 
    mf::LogInfo("LArPandora") << " *** LArPandora::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "] *** " << std::endl;

    cet::cpu_timer theClock;

    HitVector theArtHits;
    HitToParticleMap theHitToParticleMap;
    TruthToParticleMap theTruthToParticleMap;
    ParticleMap theParticleMap;
    HitMap thePandoraHits;

    this->PrepareEvent(evt);
    this->CollectArtHits(evt, theArtHits, theHitToParticleMap);
    this->CreatePandoraHits(theArtHits, thePandoraHits);

    if (m_enableMCParticles && !evt.isRealData())
    {
        this->CollectArtParticles(evt, theParticleMap, theTruthToParticleMap);
        this->CreatePandoraParticles(theParticleMap, theTruthToParticleMap);
        this->CreatePandoraLinks(thePandoraHits, theHitToParticleMap);
    }

    if (m_enableMonitoring)
        theClock.start();

    this->RunPandora();

    if (m_enableMonitoring)
        theClock.stop();

    if (m_enableProduction)
        this->ProduceArtClusters(evt, thePandoraHits);

    this->ResetPandora();

    if (m_enableMonitoring)
    {
        m_run   = evt.run();
        m_event = evt.id().event();
        m_time  = theClock.accumulated_real_time();
        m_hits  = static_cast<int>(theArtHits.size());
        m_pRecoTree->Fill();
    }
   
    mf::LogDebug("LArPandora") << " *** LArPandora::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "]  Done! *** " << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::InitializePandora() const
{ 
    mf::LogDebug("LArPandora") << " *** LArPandora::InitializePandora(...) *** " << std::endl;

    // Find the Pandora settings file (must be within 'FW_SEARCH_PATH')
    cet::search_path sp("FW_SEARCH_PATH");
    std::string configFileName("");

    mf::LogDebug("LArPandora") << "   Load Pandora settings: " << m_configFile << std::endl;
    mf::LogDebug("LArPandora") << "   Search path: " << sp.to_string() << std::endl;

    if (!sp.find_file(m_configFile, configFileName))
    {
        mf::LogError("LArPandora") << "   Failed to find: " << m_configFile << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    }
    else
    {
        mf::LogDebug("LArPandora") << "   Found it: " <<  configFileName << std::endl;
    }
    
    // Identify the Geometry and load the calculators
    art::ServiceHandle<geo::Geometry> theGeometry;

    if (theGeometry->DetId() == geo::kMicroBooNE)
    {
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerCalculator(*m_pPandora, new MicroBooNEPseudoLayerCalculator));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationCalculator(*m_pPandora, new MicroBooNETransformationCalculator));
    }
    else
    {
        mf::LogError("LArPandora") << " Geometry helpers not yet available for detector: " << theGeometry->GetDetectorName() << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);   
    }

    // Register the algorithms and read the settings
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*m_pPandora));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterHelperFunctions(*m_pPandora));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterResetFunctions(*m_pPandora));

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPandora, configFileName));

    mf::LogDebug("LArPandora") << " *** LArPandora::InitializePandora(...)  Done! *** " << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::InitializeMonitoring()
{
    art::ServiceHandle<art::TFileService> tfs;
    m_pRecoTree = tfs->make<TTree>("monitoring", "LAr Reco");
    m_pRecoTree->Branch("run", &m_run, "run/I");
    m_pRecoTree->Branch("event", &m_event, "event/I");
    m_pRecoTree->Branch("hits", &m_hits, "hits/I");
    m_pRecoTree->Branch("time", &m_time, "time/F");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::PrepareEvent(const art::Event &evt)
{
    m_run   = evt.run();
    m_event = evt.id().event();
    m_hits  = 0;
    m_time  = 0.f;
 
    if (m_enableMCParticles && !evt.isRealData())
    {
        art::ServiceHandle<cheat::BackTracker> theBackTracker; 

	// Bail out if there is no back-tracking information
        if( theBackTracker->GetSetOfTrackIDs().size() == 0 )
	{
	    mf::LogError("LArPandora") << "   Failed to load back-tracking data " << std::endl;
	    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);  
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::CollectArtHits(const art::Event &evt, HitVector &hitVector, HitToParticleMap &hitToParticleMap) const
{
    mf::LogDebug("LArPandora") << " *** LArPandora::CollectArtHits(...) *** " << std::endl;

    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(m_hitfinderModuleLabel, hitHandle);

    art::ServiceHandle<cheat::BackTracker> theBackTracker; 

    for (unsigned int iHit = 0, iHitEnd = hitHandle->size(); iHit < iHitEnd; ++iHit)
    {
        art::Ptr<recob::Hit> hit(hitHandle, iHit);
        hitVector.push_back(hit);

        if (m_enableMCParticles && !evt.isRealData())
        {
            const std::vector<cheat::TrackIDE> trackCollection(theBackTracker->HitToTrackID(hit));

            for (unsigned int iTrack = 0, iTrackEnd = trackCollection.size(); iTrack < iTrackEnd; ++iTrack)
            {
                cheat::TrackIDE trackIDE = trackCollection.at(iTrack);
                hitToParticleMap[hit].push_back(trackIDE);
	    }
	}
    }

    mf::LogDebug("LArPandora") << "   Number of ART hits: " << hitVector.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::CollectArtParticles(const art::Event &evt, ParticleMap &particleMap, TruthToParticleMap &truthToParticleMap) const
{
    mf::LogDebug("LArPandora") << " *** LArPandora::CollectArtParticles(...) *** " << std::endl;

    art::Handle< std::vector<simb::MCParticle> > mcParticleHandle;
    evt.getByLabel(m_geantModuleLabel, mcParticleHandle);

    art::ServiceHandle<cheat::BackTracker> theBackTracker; 

    for (unsigned int i = 0, iEnd = mcParticleHandle->size(); i < iEnd; ++i)
    {
        art::Ptr<simb::MCParticle> particle(mcParticleHandle, i);
        particleMap[particle->TrackId()] = particle;

        art::Ptr<simb::MCTruth> truth(theBackTracker->TrackIDToMCTruth(particle->TrackId()));
        truthToParticleMap[truth].push_back(particle->TrackId());
    }

    mf::LogDebug("LArPandora") << "   Number of ART particles: " << particleMap.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::CreatePandoraHits(const HitVector &hitVector, HitMap &hitMap) const
{
    mf::LogDebug("LArPandora") << " *** LArPandora::CreatePandoraHits(...) *** " << std::endl;

    static const double dx_cm(0.5);          // cm
    static const double int_cm(84.0);        // cm
    static const double rad_cm(14.0);        // cm
    static const double dEdX_max(25.0);      // MeV/cm
    static const double dEdX_mip(2.0);       // MeV/cm (for now)
    static const double mips_to_gev(3.5e-4); // from 100 single-electrons

    art::ServiceHandle<util::DetectorProperties> theDetector;
    static const double us_per_tdc(1.0e-3 * theDetector->SamplingRate()); // ns->us
    static const double tdc_offset(theDetector->TriggerOffset());

    art::ServiceHandle<geo::Geometry> theGeometry;
    art::ServiceHandle<util::LArProperties> theLiquidArgon;

    // Calculate offsets in coordinate system
    // (until we find the wireID->Upos and wireID->Vpos methods in LArSoft!)
    double y0(0.f); double z0(0.f);
    theGeometry->IntersectionPoint(0, 0, geo::kU,geo::kV, 0, 0, y0, z0);

    const double wire_pitch_cm(lar::LArGeometryHelper::GetLArPseudoLayerCalculator()->GetZPitch());
    const double u0(lar::LArGeometryHelper::GetLArTransformationCalculator()->YZtoU(y0, z0));
    const double v0(lar::LArGeometryHelper::GetLArTransformationCalculator()->YZtoV(y0, z0));
    const double w0(wire_pitch_cm);
    const double x0(0.f); // if necessary, make the hits line up with the MC particles

    // Loop over hits
    int hitCounter(0);

    for (HitVector::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::Hit> hit = *iter;

        const int hit_View(hit->View());
        const double hit_Time(hit->PeakTime());
        const double hit_Charge(hit->Charge(true));
        const double hit_TimeStart(hit->StartTime());
        const double hit_TimeEnd(hit->EndTime());

        const geo::WireID hit_WireID(hit->WireID());

        const double wpos_cm(hit_WireID.Wire * wire_pitch_cm);
        const double xpos_cm(theDetector->ConvertTicksToX(hit_Time, hit_WireID.Plane, hit_WireID.TPC, hit_WireID.Cryostat));
        const double dxpos_cm(theDetector->ConvertTicksToX(hit_TimeEnd, hit_WireID.Plane, hit_WireID.TPC, hit_WireID.Cryostat) -
                              theDetector->ConvertTicksToX(hit_TimeStart, hit_WireID.Plane, hit_WireID.TPC, hit_WireID.Cryostat));

        const double t_us((hit_Time - tdc_offset) * us_per_tdc);
        const double dQdX(hit_Charge / wire_pitch_cm); // ADC/cm
        const double dQdX_e(dQdX / (theDetector->ElectronsToADC() * exp(-t_us / theLiquidArgon->ElectronLifetime()))); // e/cm

        double dEdX(theLiquidArgon->BirksCorrection(dQdX_e));

        if ((dEdX < 0) || (dEdX > dEdX_max))
            dEdX = dEdX_max;

        const double mips(dEdX / dEdX_mip); // TODO: Check if calibration procedure is correct

        hitMap[++hitCounter] = hit;

        // Create Pandora CaloHit
        PandoraApi::CaloHit::Parameters caloHitParameters;
        caloHitParameters.m_expectedDirection = pandora::CartesianVector(0., 0., 1.);
        caloHitParameters.m_cellNormalVector = pandora::CartesianVector(0., 0., 1.);
        caloHitParameters.m_cellSizeU = dx_cm;
        caloHitParameters.m_cellSizeV = dxpos_cm;
        caloHitParameters.m_cellThickness = wire_pitch_cm;
        caloHitParameters.m_time = 0.;
        caloHitParameters.m_nCellRadiationLengths = dx_cm / rad_cm;
        caloHitParameters.m_nCellInteractionLengths = dx_cm / int_cm;
        caloHitParameters.m_isDigital = false;
        caloHitParameters.m_detectorRegion = pandora::ENDCAP;
        caloHitParameters.m_layer = 0;
        caloHitParameters.m_isInOuterSamplingLayer = false;
        caloHitParameters.m_inputEnergy = hit_Charge;
        caloHitParameters.m_mipEquivalentEnergy = mips;
        caloHitParameters.m_electromagneticEnergy = mips * mips_to_gev; 
        caloHitParameters.m_hadronicEnergy = mips * mips_to_gev; 
        caloHitParameters.m_pParentAddress = (void*)((intptr_t)hitCounter); 

        if (hit_View == geo::kW)
        {
            caloHitParameters.m_hitType = pandora::TPC_VIEW_W;
            caloHitParameters.m_positionVector = pandora::CartesianVector(xpos_cm + x0, 0., wpos_cm + w0);
        }
        else if(hit_View == geo::kU)
        {
            caloHitParameters.m_hitType = pandora::TPC_VIEW_U;
            caloHitParameters.m_positionVector = pandora::CartesianVector(xpos_cm + x0, 0., wpos_cm + u0);
        }
        else if(hit_View == geo::kV)
        {
            caloHitParameters.m_hitType = pandora::TPC_VIEW_V;
            caloHitParameters.m_positionVector = pandora::CartesianVector(xpos_cm + x0, 0., wpos_cm + v0);
        }
        else
        {
            mf::LogError("LArPandora") << " --- WARNING: UNKNOWN VIEW !!!  (View=" << hit_View << ")" << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }

        if (std::isnan(mips))
        {
            mf::LogError("LArPandora") << " --- WARNING: UNPHYSICAL PULSEHEIGHT !!! (MIPs=" << mips << ")" << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters)); 
    }

    mf::LogDebug("LArPandora") << "   Number of Pandora hits: " << hitCounter << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::CreatePandoraParticles(const ParticleMap& particleMap, const TruthToParticleMap &truthToParticleMap) const
{
    mf::LogDebug("LArPandora") << " *** LArPandora::CreatePandoraParticles(...) *** " << std::endl;
    static const int id_offset = 100000000;

    // Loop over incident particles
    int neutrinoCounter(0);

    for (TruthToParticleMap::const_iterator iter = truthToParticleMap.begin(), iterEnd = truthToParticleMap.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<simb::MCTruth> truth = iter->first;

        if (truth->NeutrinoSet())
        {
            const simb::MCNeutrino neutrino(truth->GetNeutrino());
            ++neutrinoCounter;

            if (neutrinoCounter > id_offset)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            const int neutrinoID(neutrinoCounter + 4 * id_offset);

            // Create Pandora 3D MC Particle
            PandoraApi::MCParticle::Parameters mcParticleParameters;
            mcParticleParameters.m_energy = neutrino.Nu().E();
            mcParticleParameters.m_momentum = pandora::CartesianVector(neutrino.Nu().Px(), neutrino.Nu().Py(), neutrino.Nu().Pz());
            mcParticleParameters.m_vertex = pandora::CartesianVector(neutrino.Nu().Vx(), neutrino.Nu().Vy(), neutrino.Nu().Vz());
            mcParticleParameters.m_endpoint = pandora::CartesianVector(neutrino.Nu().Vx(), neutrino.Nu().Vy(), neutrino.Nu().Vz());
            mcParticleParameters.m_particleId = neutrino.Nu().PdgCode();
            mcParticleParameters.m_mcParticleType = pandora::MC_STANDARD;
            mcParticleParameters.m_pParentAddress = (void*)((intptr_t)neutrinoID);
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*m_pPandora, mcParticleParameters));

            // Loop over associated particles
            const std::vector<int> particleCollection = iter->second;

            for (unsigned int j = 0; j < particleCollection.size(); ++j)
            {
                const int trackID = particleCollection.at(j);

                ParticleMap::const_iterator iterJ = particleMap.find(trackID);

                if( iterJ == particleMap.end() )
                    throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

                const art::Ptr<simb::MCParticle> particle = iterJ->second;

                // Mother/Daughter Links
                if (particle->Mother() == 0) 
                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*m_pPandora,
                        (void*)((intptr_t)neutrinoID), (void*)((intptr_t)trackID)));
            }
        }
    }

    mf::LogDebug("LArPandora") << "   Number of Pandora neutrinos: " << neutrinoCounter << std::endl;


    // Loop over G4 particles
    int particleCounter(0);
    
    for (ParticleMap::const_iterator iterI = particleMap.begin(), iterEndI = particleMap.end(); iterI != iterEndI; ++iterI )
    {
        const art::Ptr<simb::MCParticle> particle = iterI->second;

        if (particle->TrackId() != iterI->first)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        if (particle->TrackId() > id_offset)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        ++particleCounter;

        // Find Start and End Points
        int startT(0), endT(0);
        this->GetStartAndEndPoints(particle, startT, endT);

        const float vtxX(particle->Vx(startT));
        const float vtxY(particle->Vy(startT));
        const float vtxZ(particle->Vz(startT));

        const float endX(particle->Vx(endT));
        const float endY(particle->Vy(endT));
        const float endZ(particle->Vz(endT));

        const float pX(particle->Px(startT));
        const float pY(particle->Py(startT));
        const float pZ(particle->Pz(startT));
        const float E(particle->E(startT));

        // Create 3D Pandora MC Particle
        PandoraApi::MCParticle::Parameters mcParticleParameters;
        mcParticleParameters.m_energy = E;
        mcParticleParameters.m_particleId = particle->PdgCode();
        mcParticleParameters.m_momentum = pandora::CartesianVector(pX, pY, pZ);
        mcParticleParameters.m_vertex = pandora::CartesianVector(vtxX, vtxY, vtxZ);
        mcParticleParameters.m_endpoint = pandora::CartesianVector(endX, endY, endZ);
        mcParticleParameters.m_mcParticleType = pandora::MC_STANDARD;
        mcParticleParameters.m_pParentAddress = (void*)((intptr_t)particle->TrackId());
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*m_pPandora, mcParticleParameters));

        // Create Mother/Daughter Links between 3D MC Particles
        const int id_mother(particle->Mother());
        ParticleMap::const_iterator iterJ = particleMap.find(id_mother);

        if (iterJ != particleMap.end())
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*m_pPandora,
                (void*)((intptr_t)id_mother), (void*)((intptr_t)particle->TrackId())));

   
        // Create 2D Pandora MC Particles for Event Display
        const float dx(endX - vtxX);
        const float dy(endY - vtxY);
        const float dz(endZ - vtxZ);
        const float dw(lar::LArGeometryHelper::GetLArPseudoLayerCalculator()->GetZPitch());

        if (dx * dx + dy * dy + dz * dz < 0.5 * dw * dw)
            continue;

        // Create U projection
        mcParticleParameters.m_momentum = pandora::CartesianVector(pX, 0.f, 
            lar::LArGeometryHelper::GetLArTransformationCalculator()->PYPZtoPU(pY, pZ));
        mcParticleParameters.m_vertex = pandora::CartesianVector(vtxX, 0.f, 
            lar::LArGeometryHelper::GetLArTransformationCalculator()->YZtoU(vtxY, vtxZ));
        mcParticleParameters.m_endpoint = pandora::CartesianVector(endX,  0.f, 
            lar::LArGeometryHelper::GetLArTransformationCalculator()->YZtoU(endY, endZ));
        mcParticleParameters.m_mcParticleType = pandora::MC_VIEW_U;
        mcParticleParameters.m_pParentAddress = (void*)((intptr_t)(particle->TrackId() + 1 * id_offset)); 
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*m_pPandora, mcParticleParameters));

        // Create V projection
        mcParticleParameters.m_momentum = pandora::CartesianVector(pX, 0.f, 
            lar::LArGeometryHelper::GetLArTransformationCalculator()->PYPZtoPV(pY, pZ));
        mcParticleParameters.m_vertex = pandora::CartesianVector(vtxX, 0.f, 
            lar::LArGeometryHelper::GetLArTransformationCalculator()->YZtoV(vtxY, vtxZ));
        mcParticleParameters.m_endpoint = pandora::CartesianVector(endX,  0.f, 
            lar::LArGeometryHelper::GetLArTransformationCalculator()->YZtoV(endY, endZ));
        mcParticleParameters.m_mcParticleType = pandora::MC_VIEW_V;
        mcParticleParameters.m_pParentAddress = (void*)((intptr_t)(particle->TrackId() + 2 * id_offset)); 
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*m_pPandora, mcParticleParameters));

        // Create W projection
        mcParticleParameters.m_momentum = pandora::CartesianVector(pX, 0.f, pZ);
        mcParticleParameters.m_vertex = pandora::CartesianVector(vtxX, 0.f, vtxZ);
        mcParticleParameters.m_endpoint = pandora::CartesianVector(endX,  0.f, endZ);
        mcParticleParameters.m_mcParticleType = pandora::MC_VIEW_W;
        mcParticleParameters.m_pParentAddress = (void*)((intptr_t)(particle->TrackId() + 3 * id_offset)); 
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*m_pPandora, mcParticleParameters));
    }

    mf::LogDebug("LArPandora") << "   Number of Pandora particles: " << particleCounter << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::CreatePandoraLinks(const HitMap &hitMap, const HitToParticleMap &hitToParticleMap) const
{
    mf::LogDebug("LArPandora") << " *** LArPandora::CreatePandoraLinks(...) *** " << std::endl;

    for (HitMap::const_iterator iterI = hitMap.begin(), iterEndI = hitMap.end(); iterI != iterEndI ; ++iterI)
    {
        const int hitID(iterI->first);
        const art::Ptr<recob::Hit> hit(iterI->second);

        HitToParticleMap::const_iterator iterJ = hitToParticleMap.find(hit);

        if (iterJ == hitToParticleMap.end())
            continue;

        std::vector<cheat::TrackIDE> trackCollection = iterJ->second;

        if (trackCollection.size() == 0)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        for (unsigned int k = 0; k < trackCollection.size(); ++k)
        {
            const cheat::TrackIDE trackIDE(trackCollection.at(k));
            const int trackID(std::abs(trackIDE.trackID)); // TODO: Find out why std::abs is needed
            const float energyFrac(trackIDE.energyFrac);

            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(*m_pPandora,
                (void*)((intptr_t)hitID), (void*)((intptr_t)trackID), energyFrac));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::RunPandora() const
{
    mf::LogDebug("LArPandora") << " *** LArPandora::RunPandora() *** " << std::endl;

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPandora));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::ResetPandora() const
{
    mf::LogDebug("LArPandora") << " *** LArPandora::ResetPandora() *** " << std::endl;

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPandora));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::ProduceArtClusters(art::Event &evt, const HitMap &hitMap) const
{
    mf::LogDebug("LArPandora") << " *** LArPandora::ProduceArtClusters() *** " << std::endl;

    const pandora::PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pPandora, pPfoList));

    std::unique_ptr< std::vector<recob::Cluster> > artClusterVector( new std::vector<recob::Cluster> );
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > artClusterAssociations( new art::Assns<recob::Cluster, recob::Hit> );

    for (pandora::PfoList::iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
    {
        const pandora::ParticleFlowObject *const pPfo = *iter;
        const pandora::ClusterAddressList clusterAddressList(pPfo->GetClusterAddressList()); 

        HitVector artClusterHits;

        for (pandora::ClusterAddressList::const_iterator cIter = clusterAddressList.begin(), cIterEnd = clusterAddressList.end();
            cIter != cIterEnd; ++cIter)
        {
            const pandora::CaloHitAddressList &caloHitAddressList(*cIter);

            for (pandora::CaloHitAddressList::const_iterator hIter = caloHitAddressList.begin(), hIterEnd = caloHitAddressList.end();
                hIter != hIterEnd; ++hIter)
            {
                const void *pHitAddress(*hIter);
                const intptr_t hitID_temp((intptr_t)(pHitAddress)); // TODO
                const int hitID((int)(hitID_temp));

                HitMap::const_iterator artIter = hitMap.find(hitID);

                if (artIter == hitMap.end())
                    throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

                art::Ptr<recob::Hit> hit = artIter->second;
                artClusterHits.push_back(hit);
            }
        }

        if (!artClusterHits.empty())
        {
            recob::Cluster artCluster;
            artClusterVector->push_back(artCluster);
            util::CreateAssn(*this, evt, *(artClusterVector.get()), artClusterHits, *(artClusterAssociations.get()));
        }
    }

    evt.put(std::move(artClusterVector));
    evt.put(std::move(artClusterAssociations));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::GetStartAndEndPoints(const art::Ptr<simb::MCParticle> &particle, int& startT, int& endT) const
{
    art::ServiceHandle<geo::Geometry> theGeometry;

    bool foundStartPosition(false);

    const int numTrajectoryPoints(static_cast<int>(particle->NumberTrajectoryPoints()));

    for (int nt = 0; nt < numTrajectoryPoints; ++nt)
    {
        try{
            double pos[3] = {particle->Vx(nt), particle->Vy(nt), particle->Vz(nt)};
            unsigned int tpc   = 0;
            unsigned int cstat = 0;
            theGeometry->PositionToTPC(pos, tpc, cstat);

            endT = nt;
            if (!foundStartPosition)
            {
                startT = endT;
                foundStartPosition = true;
            }   
        }
        catch(cet::exception &e){
            continue;
        }
    }
}

} // namespace lar_pandora
