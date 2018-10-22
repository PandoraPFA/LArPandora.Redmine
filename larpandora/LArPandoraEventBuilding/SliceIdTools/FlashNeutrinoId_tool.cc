/**
 *  @file   larpandora/LArPandoraEventBuilding/SliceIdTools/FlashNeutrinoId_tool.cc
 *
 *  @brief  implementation of the flash based neutrino id tool
 */

#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

#include "ubana/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "ubana/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraSliceIdHelper.h"
#include "larpandora/LArPandoraEventBuilding/SliceIdBaseTool.h"
#include "larpandora/LArPandoraEventBuilding/Slice.h"

#include "Objects/CartesianVector.h"

#include "TFile.h"
#include "TTree.h"

#include <numeric>

namespace lar_pandora
{

/**
 *  @brief  Neutrino ID tool that selects the most likely neutrino slice using PMT information
 */
class FlashNeutrinoId : SliceIdBaseTool
{
public:
    /**
     *  @brief  Default constructor
     *
     *  @param  pset FHiCL parameter set
     */
    FlashNeutrinoId(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Classify slices as neutrino or cosmic
     *
     *  @param  slices the input vector of slices to classify
     *  @param  evt the art event
     */
    void ClassifySlices(SliceVector &slices, const art::Event &evt) override;

private:

    /**
     *  @brief  A description of the reason the tool couldn't find a neutrino candidate
     */
    class FailureMode
    {
    public:
        /**
         *  @brief  Default constructor
         *
         *  @param  reason the reason for the failure
         */
        FailureMode(const std::string &reason);

        /**
         *  @brief  Default destructor - explains the failure
         */
        ~FailureMode();

    private:
        std::string m_reason;  ///< The reason for the failure
    };
    
    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Class to hold information about the event for monitoring
     */
    class OutputEvent
    {
    public:
        /**
         *  @brief  Reset the variables to default dummy values
         *
         *  @param  event the art event
         */
        void Reset(const art::Event &event);

        int  m_run;                   ///< The run number
        int  m_subRun;                ///< The subRun number
        int  m_event;                 ///< The event number
        int  m_nFlashes;              ///< The number of flashes
        int  m_nFlashesInBeamWindow;  ///< The number of flashes in the beam window
        bool m_hasBeamFlash;          ///< If a beam flash was found
        int  m_nSlices;               ///< The number of slices
        int  m_nSlicesAfterPrecuts;   ///< The number of slices remaining after the preselection cuts
        bool m_foundATargetSlice;     ///< If a slice was identified as the target (neutrino)
    };
    
    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  A candidate for the beam flash
     */
    class FlashCandidate
    {
    public:
        /**
         *  @brief  Default constructor
         */
        FlashCandidate();

        /**
         *  @brief  Parametrized constructor
         *
         *  @param  event the art event
         *  @param  flash the flash
         */
        FlashCandidate(const art::Event &event, const recob::OpFlash &flash);

        /**
         *  @brief  Check if the time of the flash is in the beam window, and save the information for later
         *
         *  @param  beamWindowStart the starting time of the beam window
         *  @param  beamWindowWend the end time of the beam window
         */
        bool IsInBeamWindow(const float beamWindowStart, const float beamWindowEnd);

        /**
         *  @brief  Check if the flash passes the minimum PE threshold
         *
         *  @param  minBeamFlashPE the minimum number of photo electrons to pass
         */
        bool PassesPEThreshold(const float minBeamFlashPE) const;
    
        /**
         *  @breif  Convert to a flashana::Flash_t
         *
         *  @param  opDetVector the ordered vector of optical detector IDs
         *
         *  @return the flashana::Flash_t
         */
        flashana::Flash_t ConvertFlashFormat(const std::vector<unsigned int> &opDetVector) const;

        // Features of the flash are used when writing to file is enabled
        int                 m_run;                  ///< The run number
        int                 m_subRun;               ///< The subRun number
        int                 m_event;                ///< The event number
        float               m_time;                 ///< Time of the flash
        std::vector<float>  m_peSpectrum;           ///< The number of PEs on each PMT
        float               m_totalPE;              ///< The total number of photoelectrons over all PMTs in the flash
        float               m_centerY;              ///< The PE weighted center Y position of the flash
        float               m_centerZ;              ///< The PE weighted center Z position of the flash
        float               m_widthY;               ///< The PE weighted width of the flash in Y
        float               m_widthZ;               ///< The PE weighted width of the flash in Z
        bool                m_inBeamWindow;         ///< If the flash is in time with the beam window
        bool                m_isBrightestInWindow;  ///< If the flash is the brightest in the event
        bool                m_isBeamFlash;          ///< If the flash has been selected as the beam flash
    };
    
    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  A candidate for the target slice
     */
    class SliceCandidate
    {
    public:
        /**
         *  @brief  Data to describe an amount of charge deposited in a given 3D position
         */
        class Deposition
        {
        public:
            /**
             *  @brief  Default constructor
             *
             *  @param  x the x-component of the charge position
             *  @param  y the z-component of the charge position
             *  @param  z the z-component of the charge position
             *  @param  charge the charge deposited
             *  @param  nPhotons the estimated numer of photons produced
             */
            Deposition(const float x, const float y, const float z, const float charge, const float nPhotons);

            float m_x;         ///< The x-component of the charge position
            float m_y;         ///< The z-component of the charge position
            float m_z;         ///< The z-component of the charge position
            float m_charge;    ///< The charge deposited
            float m_nPhotons;  ///< The estimated numer of photons produced
        };
    
        typedef std::vector<Deposition> DepositionVector;

        // ---------------------------------------------------------------------------------------------------------------------------------

        /**
         *  @brief  Default constructor
         */
        SliceCandidate();

        /**
         *  @brief  Parametrized constructor
         *
         *  @param  event the art event
         *  @param  slice the slice
         *  @param  pfParticleMap the input mapping from PFParticle ID to PFParticle
         *  @param  pfParticleToSpacePointMap the input mapping from PFParticles to SpacePoints
         *  @param  spacePointToHitMap the input mapping from SpacePoints to Hits
         */
        SliceCandidate(const art::Event &event, const Slice &slice, const PFParticleMap &pfParticleMap,
            const PFParticlesToSpacePoints &pfParticleToSpacePointMap, const SpacePointsToHits &spacePointToHitMap,
            const float chargeToNPhotonsTrack, const float chargeToNPhotonsShower);

        /**
         *  @brief  Parametrized constructor for slices that aren't considered by the flash neutrino ID module - monitoring only
         *
         *  @param  event the art event
         *  @param  slice the slice
         */
        SliceCandidate(const art::Event &event, const Slice &slice);
        
        /**
         *  @breif  Determine if a given slice is compatible with the beam flash by applying pre-selection cuts
         *
         *  @param  beamFlash the beam flash
         *  @param  maxDeltaY the maximum difference in Y between the beam flash center and the weighted charge center
         *  @param  maxDeltaZ the maximum difference in Z between the beam flash center and the weighted charge center
         *  @param  maxDeltaYSigma as for maxDeltaY, but measured in units of the flash width in Y
         *  @param  maxDeltaZSigma as for maxDeltaZ, but measured in units of the flash width in Z
         *  @param  minChargeToLightRatio the minimum ratio between the total charge and the total PE
         *  @param  maxChargeToLightRatio the maximum ratio between the total charge and the total PE
         *
         *  @return if the slice is compatible with the beamFlash
         */
        bool IsCompatibleWithBeamFlash(const FlashCandidate &beamFlash, const float maxDeltaY, const float maxDeltaZ,
            const float maxDeltaYSigma, const float maxDeltaZSigma, const float minChargeToLightRatio, const float maxChargeToLightRatio);

        /** 
         *  @brief  Get the flash matching score between this slice and the beam flash
         *
         *  @param  beamFlash the beam flash
         *  @param  flashMatchManager the flash matching manager
         *  @param  opDetVector the ordered vector of optical detector IDs
         *
         *  @return the flash matching score
         */
        float GetFlashMatchScore(const FlashCandidate &beamFlash, flashana::FlashMatchManager &flashMatchManager,
            const std::vector<unsigned int> &opDetVector);

    private:
        /**
         *  @breif  Get the 3D spacepoints (with charge) associated with the PFParticles in the slice that are produced from hits in the W view
         *
         *  @param  pfParticleMap the input mapping from PFParticle ID to PFParticle
         *  @param  pfParticleToSpacePointMap the input mapping from PFParticles to SpacePoints
         *  @param  spacePointToHitMap the input mapping from SpacePoints to Hits
         *  @param  slice the input slice
         *
         *  @return the output depositionVector
         */
        DepositionVector GetDepositionVector(const PFParticleMap &pfParticleMap, const PFParticlesToSpacePoints &pfParticleToSpacePointMap,
            const SpacePointsToHits &spacePointToHitMap, const Slice &slice) const;
    
        /**
         *  @breif  Collect all downstream particles of those in the input vector
         *
         *  @param  pfParticleMap the mapping from PFParticle ID to PFParticle
         *  @param  parentPFParticles the input vector of PFParticles
         *  @param  downstreamPFParticle the output vector of PFParticles including those downstream of the input
         */
        void CollectDownstreamPFParticles(const PFParticleMap &pfParticleMap, const PFParticleVector &parentPFParticles, 
            PFParticleVector &downstreamPFParticles) const;
    
        /**
         *  @breif  Collect all downstream particles of a given particle
         *
         *  @param  pfParticleMap the mapping from PFParticle ID to PFParticle
         *  @param  particle the input PFParticle
         *  @param  downstreamPFParticle the output vector of PFParticles including those downstream of the input
         */
        void CollectDownstreamPFParticles(const PFParticleMap &pfParticleMap, const art::Ptr<recob::PFParticle> &particle,
            PFParticleVector &downstreamPFParticles) const;
    
        /**
         *  @brief  Convert from deposited charge to number of photons for a given particle
         *
         *  @param  charge the input charge
         *  @param  particle the input particle
         *
         *  @return the number of photons
         */
        float GetNPhotons(const float charge, const art::Ptr<recob::PFParticle> &particle) const;
    
        /**
         *  @brief  Get the centroid of the input charge cluster, weighted by charge
         *
         *  @param  depositionVector the input charge cluster
         *
         *  @return the charge weighted centroid
         */
        pandora::CartesianVector GetChargeWeightedCenter(const DepositionVector &depositionVector) const;
    
        /**
         *  @brief  Get the total charge from an input charge cluster
         *
         *  @param  depositionVector the input charge cluster
         *
         *  @return the total charge
         */
        float GetTotalCharge(const DepositionVector &depositionVector) const;

        /**
         *  @brief  Get the minimum X-position of all deposition points given
         *
         *  @param  depositionVector the input charge cluster
         *
         *  @return the minimum X-position
         */
        float GetMinimumXPosition(const DepositionVector &depositionVector) const;

        /**
         *  @brief  Convert a charge deposition into a light cluster by applying the chargeToPhotonFactor to every point
         *
         *  @param  depositionVector the input charge cluster
         *
         *  @return the output light cluster
         */
        flashana::QCluster_t GetLightCluster(const DepositionVector &depositionVector) const;
    
    public:
        // Features of the slice are used when writing to file is enabled
        int                  m_run;                      ///< The run number
        int                  m_subRun;                   ///< The subRun number
        int                  m_event;                    ///< The event number
        bool                 m_hasDeposition;            ///< If the slice has any charge deposited on the collection plane which produced a spacepoint
        float                m_totalCharge;              ///< The total charge deposited on the collection plane by hits that produced spacepoints
        float                m_centerX;                  ///< The charge weighted center of the slice in X
        float                m_centerY;                  ///< The charge weighted center of the slice in Y
        float                m_centerZ;                  ///< The charge weighted center of the slice in Z
        float                m_minX;                     ///< The minimum X-coordinate of all spacepoints in the slice
        float                m_deltaY;                   ///< The distance of the slice centroid from the flash centroid in Y
        float                m_deltaZ;                   ///< The distance of the slice centroid from the flash centroid in Z
        float                m_deltaYSigma;              ///< deltaY but in units of the flash width in Y
        float                m_deltaZSigma;              ///< deltaZ but in units of the flash width in Z
        float                m_chargeToLightRatio;       ///< The ratio between the total charge and the total PE of the beam flash
        bool                 m_passesPrecuts;            ///< If the slice passes the preselection cuts
        float                m_flashMatchScore;          ///< The flash matching score between the slice and the beam flash
        float                m_flashMatchX;              ///< The etimated X coordinate of the flashmatching
        float                m_totalPEHypothesis;        ///< The total PE of the hypothesized flash for this slice
        std::vector<float>   m_peHypothesisSpectrum;     ///< The PE of the hypothesized flash of this slice 
        bool                 m_isTaggedAsTarget;         ///< If the slice has been tagged as the target (neutrino)
        bool                 m_isConsideredByFlashId;    ///< If the slice was considered by the flash ID tool - this will be false if there wasn't a beam flash found in the event
        float                m_topologicalNeutrinoScore; ///< The topological-information-only neutrino ID score from Pandora
        bool                 m_hasBestTopologicalScore;  ///< If this slice has the highest topological score in the event
        
        float                m_chargeToNPhotonsTrack;    ///< The conversion factor between charge and number of photons for tracks
        float                m_chargeToNPhotonsShower;   ///< The conversion factor between charge and number of photons for showers
        flashana::QCluster_t m_lightCluster;             ///< The hypothesised light produced - used by flashmatching
    };
    
    // -------------------------------------------------------------------------------------------------------------------------------------

    typedef std::vector<recob::OpFlash> FlashVector;
    typedef std::vector<SliceCandidate> SliceCandidateVector;
    typedef std::vector<FlashCandidate> FlashCandidateVector;

    /**
     *  @brief  Get the ordered vector of optical detector IDs, use the remapping provided in FHiCL if required
     *
     *  @param  pset FHiCL parameter set
     */
    void GetOrderedOpDetVector(fhicl::ParameterSet const &pset);
    
    /**
     *  @brief  Get the candidate flashes in the event
     *
     *  @param  event the art event
     *  @param  flashCandidates the output vector of flash candidates
     */
    void GetFlashCandidates(const art::Event &event, FlashCandidateVector &flashCandidates);

    /**
     *  @breif  Try to find the brightest flash with sufficent photoelectons that is in time with the beam
     *
     *  @param  flashCandidates the input vector of slice candidates
     *
     *  @return the beam flash
     */
    FlashCandidate& GetBeamFlash(FlashCandidateVector &flashCandidates);
    
    /**
     *  @brief  Get the candidate slices in the event
     *
     *  @param  event the art event
     *  @param  slices the input vector of slices
     *  @param  sliceCandidates the output vector of slice candidates
     */
    void GetSliceCandidates(const art::Event &event, SliceVector &slices, SliceCandidateVector &sliceCandidates);

    /**
     *  @brief  Get the index of the slice which should be tagged as a neutrino
     *
     *  @param  beamFlash the beam flash
     *  @param  sliceCandidates the neutrino slice candidates
     */
    unsigned int GetBestSliceIndex(const FlashCandidate &beamFlash, SliceCandidateVector &sliceCandidates);

    /**
     *  @brief  Fill the event tree
     */
    void FillEventTree();
    
    /**
     *  @brief  Fill the flash tree
     *
     *  @param  flashCandidate the candidate flashes
     */
    void FillFlashTree(const FlashCandidateVector &flashCandidates);
    
    /**
     *  @brief  Fill the slice tree
     *
     *  @param  evt the art event
     *  @param  slices the input vector of slices
     *  @param  sliceCandidates the candidate slices
     */
    void FillSliceTree(const art::Event &evt, const SliceVector &slices, const SliceCandidateVector &sliceCandidates);

    /**
     *  @brief  Identify the slice which has the largest topological neutrino ID score, and flag it
     *
     *  @param  sliceCandidates the candidate slices
     */
    void IdentifySliceWithBestTopologicalScore(SliceCandidateVector &sliceCandidates) const;

    // Producer labels
    std::string  m_flashLabel;    ///< The label of the flash producer
    std::string  m_pandoraLabel;  ///< The label of the allOutcomes pandora producer

    // Cuts for selecting the beam flash
    float        m_beamWindowStart;  ///< The start time of the beam window
    float        m_beamWindowEnd;    ///< The end time of the beam window
    float        m_minBeamFlashPE;   ///< The minimum number of photoelectrons required to consider a flash as the beam flash

    // Pre-selection cuts to determine if a slice is compatible with the beam flash
    float        m_maxDeltaY;              ///< The maximum difference in Y between the beam flash center and the weighted charge center
    float        m_maxDeltaZ;              ///< The maximum difference in Z between the beam flash center and the weighted charge center
    float        m_maxDeltaYSigma;         ///< As for maxDeltaY, but measured in units of the flash width in Y
    float        m_maxDeltaZSigma;         ///< As for maxDeltaZ, but measured in units of the flash width in Z
    float        m_minChargeToLightRatio;  ///< The minimum ratio between the total charge and the total PE
    float        m_maxChargeToLightRatio;  ///< The maximum ratio between the total charge and the total PE

    // Variables required for flash matching
    float                          m_chargeToNPhotonsTrack;   ///< The conversion factor between charge and number of photons for tracks
    float                          m_chargeToNPhotonsShower;  ///< The conversion factor between charge and number of photons for showers
    flashana::FlashMatchManager    m_flashMatchManager;       ///< The flash match manager
    std::vector<unsigned int>      m_opDetVector;             ///< The ordered vector of optical detector IDs

    // Debugging / testing
    bool                                    m_shouldWriteToFile;   ///< If we should write interesting information to a root file
    bool                                    m_hasMCNeutrino;       ///< If there is an MC neutrino we can use to get truth information
    int                                     m_nuInteractionType;   ///< The interaction type code from MCTruth
    float                                   m_nuEnergy;            ///< The true neutrino energy
    float                                   m_nuVertexX;           ///< The true neutrino vertex X position
    float                                   m_nuVertexY;           ///< The true neutrino vertex Y position
    float                                   m_nuVertexZ;           ///< The true neutrino vertex Z position
    float                                   m_nuTime;              ///< The time of the true neutrino interaction
    std::string                             m_truthLabel;          ///< The MCTruth producer label
    std::string                             m_mcParticleLabel;     ///< The MCParticle producer label
    std::string                             m_hitLabel;            ///< The Hit producer label
    std::string                             m_backtrackLabel;      ///< The Hit -> MCParticle producer label
    OutputEvent                             m_outputEvent;         ///< The output event whose address is used by the output branch
    FlashCandidate                          m_outputFlash;         ///< The output flash whose address is used by the output branch
    SliceCandidate                          m_outputSlice;         ///< The output slice whose address is used by the output branch
    LArPandoraSliceIdHelper::SliceMetadata  m_outputSliceMetadata; ///< The output slice metadata whose address is used by the output branch
    TTree                                  *m_pEventTree;          ///< The event tree
    TTree                                  *m_pFlashTree;          ///< The flash tree
    TTree                                  *m_pSliceTree;          ///< The slice tree
};

DEFINE_ART_CLASS_TOOL(FlashNeutrinoId)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_pandora
{
    
FlashNeutrinoId::FlashNeutrinoId(fhicl::ParameterSet const &pset) :
    m_flashLabel(pset.get<std::string>("FlashLabel")),
    m_pandoraLabel(pset.get<std::string>("PandoraAllOutcomesLabel")),
    m_beamWindowStart(pset.get<float>("BeamWindowStartTime")),
    m_beamWindowEnd(pset.get<float>("BeamWindowEndTime")),
    m_minBeamFlashPE(pset.get<float>("BeamFlashPEThreshold")),
    m_maxDeltaY(pset.get<float>("MaxDeltaY")),
    m_maxDeltaZ(pset.get<float>("MaxDeltaZ")),
    m_maxDeltaYSigma(pset.get<float>("MaxDeltaYSigma")),
    m_maxDeltaZSigma(pset.get<float>("MaxDeltaZSigma")),
    m_minChargeToLightRatio(pset.get<float>("MinChargeToLightRatio")),
    m_maxChargeToLightRatio(pset.get<float>("MaxChargeToLightRatio")),
    m_chargeToNPhotonsTrack(pset.get<float>("ChargeToNPhotonsTrack")),
    m_chargeToNPhotonsShower(pset.get<float>("ChargeToNPhotonsShower")),
    m_shouldWriteToFile(pset.get<bool>("ShouldWriteToFile", false)),
    m_hasMCNeutrino(m_shouldWriteToFile ? pset.get<bool>("HasMCNeutrino") : false),
    m_truthLabel(m_hasMCNeutrino ? pset.get<std::string>("MCTruthLabel") : ""),
    m_mcParticleLabel(m_hasMCNeutrino ? pset.get<std::string>("MCParticleLabel") : ""),
    m_hitLabel(m_hasMCNeutrino ? pset.get<std::string>("HitLabel") : ""),
    m_backtrackLabel(m_hasMCNeutrino ? pset.get<std::string>("BacktrackerLabel") : ""),
    m_pEventTree(nullptr),
    m_pFlashTree(nullptr),
    m_pSliceTree(nullptr)
{
    m_flashMatchManager.Configure(pset.get<flashana::Config_t>("FlashMatchConfig")); 
    this->GetOrderedOpDetVector(pset);

    if (!m_shouldWriteToFile)
        return;

    // Set up the output branches
    art::ServiceHandle<art::TFileService> fileService;

    TTree *pMetadataTree = fileService->make<TTree>("metadata","");
    pMetadataTree->Branch("beamWindowStart"       , &m_beamWindowStart       , "beamWindowStart/F");
    pMetadataTree->Branch("beamWindowEnd"         , &m_beamWindowEnd         , "beamWindowEnd/F");
    pMetadataTree->Branch("minBeamFlashPE"        , &m_minBeamFlashPE        , "minBeamFlashPE/F");
    pMetadataTree->Branch("maxDeltaY"             , &m_maxDeltaY             , "maxDeltaY/F");
    pMetadataTree->Branch("maxDeltaZ"             , &m_maxDeltaZ             , "maxDeltaZ/F");
    pMetadataTree->Branch("maxDeltaYSigma"        , &m_maxDeltaYSigma        , "maxDeltaYSigma/F");
    pMetadataTree->Branch("maxDeltaZSigma"        , &m_maxDeltaZSigma        , "maxDeltaZSigma/F");
    pMetadataTree->Branch("minChargeToLightRatio" , &m_minChargeToLightRatio , "minChargeToLightRatio/F");
    pMetadataTree->Branch("maxChargeToLightRatio" , &m_maxChargeToLightRatio , "maxChargeToLightRatio/F");
    pMetadataTree->Branch("chargeToNPhotonsTrack" , &m_chargeToNPhotonsTrack , "chargeToNPhotonsTrack/F");
    pMetadataTree->Branch("chargeToNPhotonsShower", &m_chargeToNPhotonsShower, "chargeToNPhotonsShower/F");
    pMetadataTree->Fill();

    m_pEventTree = fileService->make<TTree>("events","");
    m_pEventTree->Branch("run"                 , &m_outputEvent.m_run                 , "run/I");
    m_pEventTree->Branch("subRun"              , &m_outputEvent.m_subRun              , "subRun/I");
    m_pEventTree->Branch("event"               , &m_outputEvent.m_event               , "event/I");
    m_pEventTree->Branch("nFlashes"            , &m_outputEvent.m_nFlashes            , "nFlashes/I");
    m_pEventTree->Branch("nFlashesInBeamWindow", &m_outputEvent.m_nFlashesInBeamWindow, "nFlashesInBeamWindow/I");
    m_pEventTree->Branch("hasBeamFlash"        , &m_outputEvent.m_hasBeamFlash        , "hasBeamFlash/O");
    m_pEventTree->Branch("nSlices"             , &m_outputEvent.m_nSlices             , "nSlices/I");
    m_pEventTree->Branch("nSlicesAfterPrecuts" , &m_outputEvent.m_nSlicesAfterPrecuts , "nSlicesAfterPrecuts/I");
    m_pEventTree->Branch("foundATargetSlice"   , &m_outputEvent.m_foundATargetSlice   , "foundATarget/O");

    m_pFlashTree = fileService->make<TTree>("flashes","");
    m_pFlashTree->Branch("run"                , &m_outputFlash.m_run                , "run/I");
    m_pFlashTree->Branch("subRun"             , &m_outputFlash.m_subRun             , "subRun/I");
    m_pFlashTree->Branch("event"              , &m_outputFlash.m_event              , "event/I");
    m_pFlashTree->Branch("time"               , &m_outputFlash.m_time               , "time/F");
    m_pFlashTree->Branch("centerY"            , &m_outputFlash.m_centerY            , "centerY/F");
    m_pFlashTree->Branch("centerZ"            , &m_outputFlash.m_centerZ            , "centerZ/F");
    m_pFlashTree->Branch("widthY"             , &m_outputFlash.m_widthY             , "widthY/F");
    m_pFlashTree->Branch("widthZ"             , &m_outputFlash.m_widthZ             , "widthZ/F");
    m_pFlashTree->Branch("totalPE"            , &m_outputFlash.m_totalPE            , "totalPE/F");
    m_pFlashTree->Branch("peSpectrum"         , "std::vector< float >"              , &m_outputFlash.m_peSpectrum);
    m_pFlashTree->Branch("inBeamWindow"       , &m_outputFlash.m_inBeamWindow       , "inBeamWindow/O");
    m_pFlashTree->Branch("isBrightestInWindow", &m_outputFlash.m_isBrightestInWindow, "isBrightestInWindow/O");
    m_pFlashTree->Branch("isBeamFlash"        , &m_outputFlash.m_isBeamFlash        , "isBeamFlash/O");

    m_pSliceTree = fileService->make<TTree>("slices","");
    m_pSliceTree->Branch("run"                    , &m_outputSlice.m_run                     , "run/I");
    m_pSliceTree->Branch("subRun"                 , &m_outputSlice.m_subRun                  , "subRun/I");
    m_pSliceTree->Branch("event"                  , &m_outputSlice.m_event                   , "event/I");
    m_pSliceTree->Branch("hasDeposition"          , &m_outputSlice.m_hasDeposition           , "hasDeposition/O");
    m_pSliceTree->Branch("totalCharge"            , &m_outputSlice.m_totalCharge             ,  "totalCharge/F");
    m_pSliceTree->Branch("centerX"                , &m_outputSlice.m_centerX                 , "centerX/F");
    m_pSliceTree->Branch("centerY"                , &m_outputSlice.m_centerY                 , "centerY/F");
    m_pSliceTree->Branch("centerZ"                , &m_outputSlice.m_centerZ                 , "centerZ/F");
    m_pSliceTree->Branch("minX"                   , &m_outputSlice.m_minX                    , "minX/F");
    m_pSliceTree->Branch("deltaY"                 , &m_outputSlice.m_deltaY                  , "deltaY/F");
    m_pSliceTree->Branch("deltaZ"                 , &m_outputSlice.m_deltaZ                  , "deltaZ/F");
    m_pSliceTree->Branch("deltaYSigma"            , &m_outputSlice.m_deltaYSigma             , "deltaYSigma/F");
    m_pSliceTree->Branch("deltaZSigma"            , &m_outputSlice.m_deltaZSigma             , "deltaZSigma/F");
    m_pSliceTree->Branch("chargeToLightRatio"     , &m_outputSlice.m_chargeToLightRatio      , "chargeToLightRatio/F");
    m_pSliceTree->Branch("passesPreCuts"          , &m_outputSlice.m_passesPrecuts           , "passesPrecuts/O");
    m_pSliceTree->Branch("flashMatchScore"        , &m_outputSlice.m_flashMatchScore         , "flashMatchScore/F");
    m_pSliceTree->Branch("flashMatchX"            , &m_outputSlice.m_flashMatchX             , "flashMatchX/F");
    m_pSliceTree->Branch("totalPEHypothesis"      , &m_outputSlice.m_totalPEHypothesis       , "totalPEHypothesis/F");
    m_pSliceTree->Branch("peHypothesisSpectrum"   , "std::vector< float >"                   , &m_outputSlice.m_peHypothesisSpectrum);
    m_pSliceTree->Branch("isTaggedAsTarget"       , &m_outputSlice.m_isTaggedAsTarget        , "isTaggedAsTarget/O");
    m_pSliceTree->Branch("isConsideredByFlashId"  , &m_outputSlice.m_isConsideredByFlashId   , "isConsideredByFlashId/O");
    m_pSliceTree->Branch("topologicalScore"       , &m_outputSlice.m_topologicalNeutrinoScore, "topologicalScore/F");
    m_pSliceTree->Branch("hasBestTopologicalScore", &m_outputSlice.m_hasBestTopologicalScore , "hasBestTopologicalScore/O");

    if (m_hasMCNeutrino)
    {
        // Truth MC information about the slice
        m_pSliceTree->Branch("purity"           , &m_outputSliceMetadata.m_purity        , "purity/F");
        m_pSliceTree->Branch("completeness"     , &m_outputSliceMetadata.m_completeness  , "completeness/F");
        m_pSliceTree->Branch("isMostComplete"   , &m_outputSliceMetadata.m_isMostComplete, "isMostComplete/O");
        m_pSliceTree->Branch("nHits"            , &m_outputSliceMetadata.m_nHits         , "nHits/I");
        m_pSliceTree->Branch("nuInteractionType", &m_nuInteractionType                   , "nuInteractionType/F");
        m_pSliceTree->Branch("nuEnergy"         , &m_nuEnergy                            , "nuEnergy/F");
        m_pSliceTree->Branch("nuInteractionTime", &m_nuTime                              , "nuInteractionTime/F");
        m_pSliceTree->Branch("nuVertexX"        , &m_nuVertexX                           , "nuVertexX/F");
        m_pSliceTree->Branch("nuVertexY"        , &m_nuVertexY                           , "nuVertexY/F");
        m_pSliceTree->Branch("nuVertexZ"        , &m_nuVertexZ                           , "nuVertexZ/F");
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::GetOrderedOpDetVector(fhicl::ParameterSet const &pset)
{
    art::ServiceHandle<geo::Geometry> geometry;
    const auto nOpDets(geometry->NOpDets());

    // Get the OpDets in their default order
    std::vector<unsigned int> opDetVector;
    for (unsigned int iChannel = 0; iChannel < nOpDets; ++iChannel)
        opDetVector.push_back(geometry->OpDetFromOpChannel(iChannel));

    // Get the remapped OpDets if required
    if (pset.get<bool>("ShouldRemapPMTs"))
    {
        const auto pmtMapping(pset.get<std::vector<unsigned int> >("OrderedPMTList"));
        
        // Ensure there are the correct number of OpDets
        if (pmtMapping.size() != nOpDets)
            throw cet::exception("FlashNeutrinoId") << "The input PMT remapping vector has the wrong size. Expected " << nOpDets << " elements." << std::endl;

        for (const auto &opDet : pmtMapping)
        {
            // Each OpDet in the default list must be listed once and only once
            if (std::count(opDetVector.begin(), opDetVector.end(), opDet) != 1)
                throw cet::exception("FlashNeutrinoId") << "Unknown or repeated PMT ID: " << opDet << std::endl;

            m_opDetVector.push_back(opDet);
        }
    }
    else
    {
        m_opDetVector = opDetVector;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::ClassifySlices(SliceVector &slices, const art::Event &evt) 
{
    // Reset the output addresses in case we are writing monitoring details to an outpu file
    m_outputEvent.Reset(evt);

    FlashCandidateVector flashCandidates;
    SliceCandidateVector sliceCandidates;

    try
    {
        // Find the flash, if any, in time with the beam with the largest number of photoelectrons that is sufficiently bright
        this->GetFlashCandidates(evt, flashCandidates);
        const auto beamFlash(this->GetBeamFlash(flashCandidates));
        
        // Find the slice - if any that matches best with the beamFlash 
        this->GetSliceCandidates(evt, slices, sliceCandidates);
        const auto bestSliceIndex(this->GetBestSliceIndex(beamFlash, sliceCandidates));
    
        // Tag the choesn slice as a neutrino
        slices.at(bestSliceIndex).TagAsTarget();
    }
    catch (const FailureMode &)
    {
    }

    if (!m_shouldWriteToFile)
        return;

    this->FillFlashTree(flashCandidates);
    this->FillSliceTree(evt, slices, sliceCandidates);
    this->FillEventTree();
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void FlashNeutrinoId::GetFlashCandidates(const art::Event &event, FlashCandidateVector &flashCandidates)
{
    // Collect all flashes from the event
    art::InputTag flashTag(m_flashLabel); 
    const auto flashes(*event.getValidHandle<FlashVector>(flashTag));
   
    for (const auto &flash : flashes)
        flashCandidates.emplace_back(event, flash);

    m_outputEvent.m_nFlashes = flashCandidates.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::FlashCandidate& FlashNeutrinoId::GetBeamFlash(FlashCandidateVector &flashCandidates)
{
    bool foundFlashInBeamWindow(false);
    unsigned int brightestFlashIndex(std::numeric_limits<unsigned int>::max());
    float maxTotalPE(-std::numeric_limits<float>::max());
    m_outputEvent.m_nFlashesInBeamWindow = 0;

    // Find the brightest flash in the beam window
    for (unsigned int flashIndex = 0; flashIndex < flashCandidates.size(); ++flashIndex)
    {
        // ATTN non const reference is required since monitoring variables are stored in the slice candidate
        auto &flashCandidate(flashCandidates.at(flashIndex));

        if (!flashCandidate.IsInBeamWindow(m_beamWindowStart, m_beamWindowEnd))
            continue;
    
        m_outputEvent.m_nFlashesInBeamWindow++;
       
        const auto totalPE(flashCandidate.m_totalPE);
        if (totalPE < maxTotalPE)
            continue;
        
        foundFlashInBeamWindow = true;
        maxTotalPE = totalPE;
        brightestFlashIndex = flashIndex;
    }

    if (!foundFlashInBeamWindow)
        throw FailureMode("There were no flashes in the beam window");

    // Ensure it is sufficiently bright
    auto &brightestFlash(flashCandidates.at(brightestFlashIndex));
    brightestFlash.m_isBrightestInWindow = true;

    if (!brightestFlash.PassesPEThreshold(m_minBeamFlashPE))
        throw FailureMode("No flashes in the beam window passed the PE threshold");

    // Save the monitoring information
    brightestFlash.m_isBeamFlash = true;
    m_outputEvent.m_hasBeamFlash = true;

    return brightestFlash;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void FlashNeutrinoId::GetSliceCandidates(const art::Event &event, SliceVector &slices, SliceCandidateVector &sliceCandidates)
{
    m_outputEvent.m_nSlices = slices.size();

    if (slices.empty())
        throw FailureMode("No slices to choose from");

    // Collect the PFParticles and their associations to SpacePoints and Hits
    PFParticleVector pfParticles;
    SpacePointVector spacePoints;
    SpacePointsToHits spacePointToHitMap;
    PFParticleMap pfParticleMap;

    PFParticlesToSpacePoints pfParticleToSpacePointMap;
    LArPandoraHelper::CollectPFParticles(event, m_pandoraLabel, pfParticles, pfParticleToSpacePointMap);
    LArPandoraHelper::CollectSpacePoints(event, m_pandoraLabel, spacePoints, spacePointToHitMap);
    LArPandoraHelper::BuildPFParticleMap(pfParticles, pfParticleMap);
    
    for (const auto &slice : slices)
        sliceCandidates.emplace_back(event, slice, pfParticleMap, pfParticleToSpacePointMap, spacePointToHitMap, m_chargeToNPhotonsTrack, m_chargeToNPhotonsShower);
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int FlashNeutrinoId::GetBestSliceIndex(const FlashCandidate &beamFlash, SliceCandidateVector &sliceCandidates)
{
    bool foundViableSlice(false);
    unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());
    float maxScore(-std::numeric_limits<float>::max());
    m_outputEvent.m_nSlicesAfterPrecuts = 0;

    for (unsigned int sliceIndex = 0; sliceIndex < sliceCandidates.size(); ++sliceIndex)
    {
        auto &sliceCandidate(sliceCandidates.at(sliceIndex));

        // Apply the pre-selection cuts to ensure that the slice is compatible with the beam flash
        if (!sliceCandidate.IsCompatibleWithBeamFlash(beamFlash, m_maxDeltaY, m_maxDeltaZ, m_maxDeltaYSigma, m_maxDeltaZSigma,
            m_minChargeToLightRatio, m_maxChargeToLightRatio))
            continue;
        
        m_outputEvent.m_nSlicesAfterPrecuts++;

        // ATTN if there is only one slice that passes the pre-selection cuts, then the score won't be used
        const auto &score(sliceCandidate.GetFlashMatchScore(beamFlash, m_flashMatchManager, m_opDetVector));
        if (score < maxScore)
            continue;

        foundViableSlice = true;
        bestSliceIndex = sliceIndex;
        maxScore = score;
    }

    if (!foundViableSlice)
        throw FailureMode("None of the slices passed the pre-selection cuts");

    m_outputEvent.m_foundATargetSlice = true;
    sliceCandidates.at(bestSliceIndex).m_isTaggedAsTarget = true;

    return bestSliceIndex;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void FlashNeutrinoId::FillEventTree()
{
    if (!m_pEventTree)
        throw cet::exception("FlashNeutrinoId") << "Trying to fill the event tree which hasn't been configured" << std::endl;

    m_pEventTree->Fill();
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void FlashNeutrinoId::FillFlashTree(const FlashCandidateVector &flashCandidates)
{
    if (!m_pFlashTree)
        throw cet::exception("FlashNeutrinoId") << "Trying to fill the flash tree which hasn't been configured" << std::endl;

    for (const auto &flashCandidate : flashCandidates)
    {
        m_outputFlash = flashCandidate;
        m_pFlashTree->Fill();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void FlashNeutrinoId::FillSliceTree(const art::Event &evt, const SliceVector &slices, const SliceCandidateVector &sliceCandidates)
{
    if (!m_pSliceTree)
        throw cet::exception("FlashNeutrinoId") << "Trying to fill the slice tree which hasn't been configured" << std::endl;

    // We won't ever have any slice candidates if there wasn't a beam flash
    SliceCandidateVector allSliceCandidates(sliceCandidates);
    if (!m_outputEvent.m_hasBeamFlash)
    {
        if (!allSliceCandidates.empty())
            throw cet::exception("FlashNeutrinoId") << "There were slice candidates made even though there wasn't a beam flash!" << std::endl;

        // ATTN this code is only required for monitoring to compare with the topological score
        for (const auto &slice : slices)
            allSliceCandidates.emplace_back(evt, slice);
    }

    if (slices.size() != allSliceCandidates.size())
        throw cet::exception("FlashNeutrinoId") << "The number of slice candidates doesn't match the number of slices" << std::endl;

    this->IdentifySliceWithBestTopologicalScore(allSliceCandidates);

    // If available, get the information from the MC neutrino
    LArPandoraSliceIdHelper::SliceMetadataVector sliceMetadata;
    if (m_hasMCNeutrino)
    {
        simb::MCNeutrino mcNeutrino;
        LArPandoraSliceIdHelper::GetSliceMetadata(slices, evt, m_truthLabel, m_mcParticleLabel, m_hitLabel, m_backtrackLabel,
            m_pandoraLabel, sliceMetadata, mcNeutrino);

        m_nuInteractionType = mcNeutrino.InteractionType();
        const auto nuMCParticle(mcNeutrino.Nu());

        m_nuEnergy = nuMCParticle.E();
        m_nuVertexX = nuMCParticle.Vx();
        m_nuVertexY = nuMCParticle.Vy();
        m_nuVertexZ = nuMCParticle.Vz();
        m_nuTime = nuMCParticle.T();
    }

    if (slices.size() != sliceMetadata.size())
        throw cet::exception("FlashNeutrinoId") << "The number of slice metadata doesn't match the number of slices" << std::endl;

    // Output the info for each slice
    for (unsigned int sliceIndex = 0; sliceIndex < slices.size(); ++sliceIndex)
    {
        m_outputSlice = allSliceCandidates.at(sliceIndex);

        if (m_hasMCNeutrino)
            m_outputSliceMetadata = sliceMetadata.at(sliceIndex);

        m_pSliceTree->Fill();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::IdentifySliceWithBestTopologicalScore(SliceCandidateVector &sliceCandidates) const
{
    float bestTopologicalScore(-std::numeric_limits<float>::max());
    unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());

    for (unsigned int sliceIndex = 0; sliceIndex < sliceCandidates.size(); ++sliceIndex)
    {
        const float topologicalScore(sliceCandidates.at(sliceIndex).m_topologicalNeutrinoScore);
        if (topologicalScore < bestTopologicalScore)
            continue;

        bestTopologicalScore = topologicalScore;
        bestSliceIndex = sliceIndex;
    }

    if (bestSliceIndex > sliceCandidates.size())
        throw cet::exception("FlashNeutrinoId") << "Couldn't find slice the best topological score" << std::endl;

    sliceCandidates.at(bestSliceIndex).m_hasBestTopologicalScore = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::FailureMode::FailureMode(const std::string &reason) :
    m_reason(reason)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::FailureMode::~FailureMode()
{
    std::cout << "Flash neutrino ID - failed to find a viable neutrino slice." << std::endl;
    std::cout << m_reason << std::endl << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::Deposition::Deposition(const float x, const float y, const float z, const float charge, const float nPhotons) :
    m_x(x),
    m_y(y),
    m_z(z),
    m_charge(charge),
    m_nPhotons(nPhotons)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::OutputEvent::Reset(const art::Event &event)
{
    m_run = event.run();
    m_subRun = event.subRun();
    m_event = event.event();
    m_nFlashes = -std::numeric_limits<int>::max();
    m_nFlashesInBeamWindow = -std::numeric_limits<int>::max();
    m_hasBeamFlash = false;
    m_nSlices = -std::numeric_limits<int>::max();
    m_nSlicesAfterPrecuts = -std::numeric_limits<int>::max();
    m_foundATargetSlice = false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::FlashCandidate::FlashCandidate() :
    m_run(-std::numeric_limits<int>::max()),
    m_subRun(-std::numeric_limits<int>::max()),
    m_event(-std::numeric_limits<int>::max()),
    m_time(-std::numeric_limits<float>::max()),
    m_totalPE(-std::numeric_limits<float>::max()),
    m_centerY(-std::numeric_limits<float>::max()),
    m_centerZ(-std::numeric_limits<float>::max()),
    m_widthY(-std::numeric_limits<float>::max()),
    m_widthZ(-std::numeric_limits<float>::max()),
    m_inBeamWindow(false),
    m_isBrightestInWindow(false),
    m_isBeamFlash(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::FlashCandidate::FlashCandidate(const art::Event &event, const recob::OpFlash &flash) :
    m_run(event.run()),
    m_subRun(event.subRun()),
    m_event(event.event()),
    m_time(flash.Time()),
    m_peSpectrum(flash.PEs().begin(), flash.PEs().end()),
    m_totalPE(flash.TotalPE()),
    m_centerY(flash.YCenter()),
    m_centerZ(flash.ZCenter()),
    m_widthY(flash.YWidth()),
    m_widthZ(flash.ZWidth()),
    m_inBeamWindow(false),
    m_isBrightestInWindow(false),
    m_isBeamFlash(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool FlashNeutrinoId::FlashCandidate::IsInBeamWindow(const float beamWindowStart, const float beamWindowEnd)
{
    m_inBeamWindow = (m_time > beamWindowStart && m_time < beamWindowEnd);
    return m_inBeamWindow;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool FlashNeutrinoId::FlashCandidate::PassesPEThreshold(const float minBeamFlashPE) const
{
    return (m_totalPE > minBeamFlashPE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

flashana::Flash_t FlashNeutrinoId::FlashCandidate::ConvertFlashFormat(const std::vector<unsigned int> &opDetVector) const
{
    // Ensure the input flash is valid
    const auto nOpDets(opDetVector.size());
    if (m_peSpectrum.size() != nOpDets)
        throw cet::exception("FlashNeutrinoId") << "Number of channels in beam flash doesn't match the number of OpDets!" << std::endl;

    // Set the flash properties
    flashana::Flash_t flash;
    flash.x = 0;
    flash.x_err = 0;
    flash.y = m_centerY;
    flash.y_err = m_widthY;
    flash.z = m_centerZ;
    flash.z_err = m_widthZ;
    flash.time = m_time;
    flash.pe_v.resize(nOpDets);
    flash.pe_err_v.resize(nOpDets);

    // Fill the flash with the PE spectrum
    for (unsigned int i = 0; i < nOpDets; ++i)
    {
        const auto opDet(opDetVector.at(i));
        if (opDet < 0 || opDet >= nOpDets)
            throw cet::exception("FlashNeutrinoId") << "OpDet ID, " << opDet << ", is out of range: 0 - " << (nOpDets-1) << std::endl;

        const auto PE(m_peSpectrum.at(i));
        flash.pe_v.at(opDet) = PE;
        flash.pe_err_v.at(opDet) = std::sqrt(PE);
    }

    return flash;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::SliceCandidate() :
    m_run(-std::numeric_limits<int>::max()),
    m_subRun(-std::numeric_limits<int>::max()),
    m_event(-std::numeric_limits<int>::max()),
    m_hasDeposition(false),
    m_totalCharge(-std::numeric_limits<float>::max()),
    m_centerX(-std::numeric_limits<float>::max()),
    m_centerY(-std::numeric_limits<float>::max()),
    m_centerZ(-std::numeric_limits<float>::max()),
    m_minX(-std::numeric_limits<float>::max()),
    m_deltaY(-std::numeric_limits<float>::max()),
    m_deltaZ(-std::numeric_limits<float>::max()),
    m_deltaYSigma(-std::numeric_limits<float>::max()),
    m_deltaZSigma(-std::numeric_limits<float>::max()),
    m_chargeToLightRatio(-std::numeric_limits<float>::max()),
    m_passesPrecuts(false),
    m_flashMatchScore(-std::numeric_limits<float>::max()),
    m_totalPEHypothesis(-std::numeric_limits<float>::max()),
    m_isTaggedAsTarget(false),
    m_isConsideredByFlashId(false),
    m_topologicalNeutrinoScore(-std::numeric_limits<float>::max()),
    m_hasBestTopologicalScore(false),
    m_chargeToNPhotonsTrack(-std::numeric_limits<float>::max()),
    m_chargeToNPhotonsShower(-std::numeric_limits<float>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::SliceCandidate(const art::Event &event, const Slice &slice) :
    m_run(event.run()),
    m_subRun(event.subRun()),
    m_event(event.event()),
    m_hasDeposition(false),
    m_totalCharge(-std::numeric_limits<float>::max()),
    m_centerX(-std::numeric_limits<float>::max()),
    m_centerY(-std::numeric_limits<float>::max()),
    m_centerZ(-std::numeric_limits<float>::max()),
    m_minX(-std::numeric_limits<float>::max()),
    m_deltaY(-std::numeric_limits<float>::max()),
    m_deltaZ(-std::numeric_limits<float>::max()),
    m_deltaYSigma(-std::numeric_limits<float>::max()),
    m_deltaZSigma(-std::numeric_limits<float>::max()),
    m_chargeToLightRatio(-std::numeric_limits<float>::max()),
    m_passesPrecuts(false),
    m_flashMatchScore(-std::numeric_limits<float>::max()),
    m_totalPEHypothesis(-std::numeric_limits<float>::max()),
    m_isTaggedAsTarget(false),
    m_isConsideredByFlashId(false),
    m_topologicalNeutrinoScore(slice.GetTopologicalScore()),
    m_hasBestTopologicalScore(false),
    m_chargeToNPhotonsTrack(-std::numeric_limits<float>::max()),
    m_chargeToNPhotonsShower(-std::numeric_limits<float>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::SliceCandidate(const art::Event &event, const Slice &slice, const PFParticleMap &pfParticleMap,
    const PFParticlesToSpacePoints &pfParticleToSpacePointMap, const SpacePointsToHits &spacePointToHitMap,
    const float chargeToNPhotonsTrack, const float chargeToNPhotonsShower) :
    m_run(event.run()),
    m_subRun(event.subRun()),
    m_event(event.event()),
    m_hasDeposition(false),
    m_totalCharge(-std::numeric_limits<float>::max()),
    m_centerX(-std::numeric_limits<float>::max()),
    m_centerY(-std::numeric_limits<float>::max()),
    m_centerZ(-std::numeric_limits<float>::max()),
    m_minX(-std::numeric_limits<float>::max()),
    m_deltaY(-std::numeric_limits<float>::max()),
    m_deltaZ(-std::numeric_limits<float>::max()),
    m_deltaYSigma(-std::numeric_limits<float>::max()),
    m_deltaZSigma(-std::numeric_limits<float>::max()),
    m_chargeToLightRatio(-std::numeric_limits<float>::max()),
    m_passesPrecuts(false),
    m_flashMatchScore(-std::numeric_limits<float>::max()),
    m_totalPEHypothesis(-std::numeric_limits<float>::max()),
    m_isTaggedAsTarget(false),
    m_isConsideredByFlashId(true),
    m_topologicalNeutrinoScore(slice.GetTopologicalScore()),
    m_hasBestTopologicalScore(false),
    m_chargeToNPhotonsTrack(chargeToNPhotonsTrack),
    m_chargeToNPhotonsShower(chargeToNPhotonsShower)
{
    const auto chargeDeposition(this->GetDepositionVector(pfParticleMap, pfParticleToSpacePointMap, spacePointToHitMap, slice));
    m_lightCluster = this->GetLightCluster(chargeDeposition);
    
    m_totalCharge = this->GetTotalCharge(chargeDeposition);
    m_hasDeposition = (m_totalCharge > std::numeric_limits<float>::epsilon());

    if (!m_hasDeposition)
        return;

    const auto chargeCenter(this->GetChargeWeightedCenter(chargeDeposition));
    m_centerX = chargeCenter.GetX();
    m_centerY = chargeCenter.GetY();
    m_centerZ = chargeCenter.GetZ();

    m_minX = this->GetMinimumXPosition(chargeDeposition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::SliceCandidate::DepositionVector FlashNeutrinoId::SliceCandidate::GetDepositionVector(const PFParticleMap &pfParticleMap,
    const PFParticlesToSpacePoints &pfParticleToSpacePointMap, const SpacePointsToHits &spacePointToHitMap, const Slice &slice) const
{
    // Collect all PFParticles in the slice, including those downstream of the primaries
    // ATTN here we only use the neutrino hypothesis, in theory this should work with either (or indeed both with some thought)
    PFParticleVector allParticlesInSlice;
    this->CollectDownstreamPFParticles(pfParticleMap, slice.GetTargetHypothesis(), allParticlesInSlice);

    DepositionVector depositionVector;
    for (const auto &particle : allParticlesInSlice)
    {
        // Get the associated spacepoints
        const auto &partToSpacePointIter(pfParticleToSpacePointMap.find(particle));
        if (partToSpacePointIter == pfParticleToSpacePointMap.end())
            continue;

        for (const auto &spacePoint : partToSpacePointIter->second)
        {
            // Get the associated hit
            const auto &spacePointToHitIter(spacePointToHitMap.find(spacePoint));
            if (spacePointToHitIter == spacePointToHitMap.end())
                continue;

            // Only use hits from the collection plane
            const auto &hit(spacePointToHitIter->second);
            if (hit->View() != geo::kZ)
                continue;
            
            // Add the charged point to the vector
            const auto &position(spacePoint->XYZ());
            const auto charge(hit->Integral());

            depositionVector.emplace_back(position[0], position[1], position[2], charge, this->GetNPhotons(charge, particle));
        }
    }

    return depositionVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::SliceCandidate::CollectDownstreamPFParticles(const PFParticleMap &pfParticleMap, const PFParticleVector &parentPFParticles,
    PFParticleVector &downstreamPFParticles) const
{
    for (const auto &particle : parentPFParticles)
        this->CollectDownstreamPFParticles(pfParticleMap, particle, downstreamPFParticles);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::SliceCandidate::CollectDownstreamPFParticles(const PFParticleMap &pfParticleMap, const art::Ptr<recob::PFParticle> &particle,
    PFParticleVector &downstreamPFParticles) const
{
    if (std::find(downstreamPFParticles.begin(), downstreamPFParticles.end(), particle) == downstreamPFParticles.end())
        downstreamPFParticles.push_back(particle);

    for (const auto &daughterId : particle->Daughters())
    {
        const auto iter(pfParticleMap.find(daughterId));
        if (iter == pfParticleMap.end())
            throw cet::exception("FlashNeutrinoId") << "Scrambled PFParticle IDs" << std::endl;

        this->CollectDownstreamPFParticles(pfParticleMap, iter->second, downstreamPFParticles);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float FlashNeutrinoId::SliceCandidate::GetNPhotons(const float charge, const art::Ptr<recob::PFParticle> &particle) const
{
    return charge * (LArPandoraHelper::IsTrack(particle) ? m_chargeToNPhotonsTrack : m_chargeToNPhotonsShower);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector FlashNeutrinoId::SliceCandidate::GetChargeWeightedCenter(const DepositionVector &depositionVector) const
{
    pandora::CartesianVector center(0.f, 0.f, 0.f);
    float totalCharge(0.f);

    for (const auto &chargePoint : depositionVector)
    {
        center += pandora::CartesianVector(chargePoint.m_x, chargePoint.m_y, chargePoint.m_z) * chargePoint.m_charge;
        totalCharge += chargePoint.m_charge;
    }

    if (totalCharge <= std::numeric_limits<float>::epsilon())
        throw cet::exception("FlashNeutrinoId") << "Can't find charge weighted center of slice with zero total charge" << std::endl;

    center *= (1.f / totalCharge);

    return center;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float FlashNeutrinoId::SliceCandidate::GetTotalCharge(const DepositionVector &depositionVector) const
{
    float totalCharge(0.f);

    for (const auto &chargePoint : depositionVector)
        totalCharge += chargePoint.m_charge;

    return totalCharge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float FlashNeutrinoId::SliceCandidate::GetMinimumXPosition(const DepositionVector &depositionVector) const
{
    float minX(std::numeric_limits<float>::max());

    for (const auto &chargePoint : depositionVector)
        minX = std::min(chargePoint.m_x, minX);

    return minX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

flashana::QCluster_t FlashNeutrinoId::SliceCandidate::GetLightCluster(const DepositionVector &depositionVector) const
{
    flashana::QCluster_t lightCluster;

    for (const auto &point : depositionVector)
        lightCluster.emplace_back(point.m_x, point.m_y, point.m_z, point.m_nPhotons);

    return lightCluster;
}        

//------------------------------------------------------------------------------------------------------------------------------------------

bool FlashNeutrinoId::SliceCandidate::IsCompatibleWithBeamFlash(const FlashCandidate &beamFlash, const float maxDeltaY,
    const float maxDeltaZ, const float maxDeltaYSigma, const float maxDeltaZSigma, const float minChargeToLightRatio,
    const float maxChargeToLightRatio)
{
    // Check the flash is usable
    if (beamFlash.m_totalPE <= std::numeric_limits<float>::epsilon())
        return false;
    
    if (beamFlash.m_widthY <= std::numeric_limits<float>::epsilon())
        return false;
    
    if (beamFlash.m_widthZ <= std::numeric_limits<float>::epsilon())
        return false;

    if (m_totalCharge <= std::numeric_limits<float>::epsilon())
        return false;
    
    // Calculate the pre-selection variables
    m_deltaY = std::abs(m_centerY - beamFlash.m_centerY);
    m_deltaZ = std::abs(m_centerZ - beamFlash.m_centerZ);
    m_deltaYSigma = m_deltaY / beamFlash.m_widthY;
    m_deltaZSigma = m_deltaZ / beamFlash.m_widthZ;
    m_chargeToLightRatio = m_totalCharge / beamFlash.m_totalPE;  // TODO ATTN check if this should be total PE or max PE. Code differs from technote
    
    // Check if the slice passes the pre-selection cuts
    m_passesPrecuts = (m_deltaY < maxDeltaY                           &&
                       m_deltaZ < maxDeltaZ                           &&
                       m_deltaYSigma < maxDeltaYSigma                 &&
                       m_deltaZSigma < maxDeltaZSigma                 &&
                       m_chargeToLightRatio > minChargeToLightRatio   &&
                       m_chargeToLightRatio < maxChargeToLightRatio   );

    return m_passesPrecuts;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float FlashNeutrinoId::SliceCandidate::GetFlashMatchScore(const FlashCandidate &beamFlash, flashana::FlashMatchManager &flashMatchManager,
    const std::vector<unsigned int> &opDetVector)
{
    flashMatchManager.Reset();

    // Convert the flash and the charge cluster into the required format for flash matching
    auto flash(beamFlash.ConvertFlashFormat(opDetVector));

    // Perform the match
    flashMatchManager.Emplace(std::move(flash));
    flashMatchManager.Emplace(std::move(m_lightCluster));
    const auto matches(flashMatchManager.Match());

    // Unable to match
    if (matches.empty())
        return -1.f;

    if (matches.size() != 1)
        throw cet::exception("FlashNeutrinoId") << "Flash matching returned multiple matches!" << std::endl;
  
    // Fill the slice candidate with the details of the matching
    const auto match(matches.front());

    m_flashMatchScore = match.score;
    m_flashMatchX = match.tpc_point.x;
    m_totalPEHypothesis = std::accumulate(match.hypothesis.begin(), match.hypothesis.end(), 0.f);

    // Fill the slice with the hypothesized PE spectrum
    if (!m_peHypothesisSpectrum.empty())
        throw cet::exception("FlashNeutrinoId") << "Hypothesized PE spectrum already set for this flash" << std::endl;

    const unsigned int nOpDets(opDetVector.size());
    if (match.hypothesis.size() != nOpDets)
        throw cet::exception("FlashNeutrinoId") << "Hypothesized PE spectrum has the wrong size" << std::endl;

    for (unsigned int i = 0; i < nOpDets; ++i)
        m_peHypothesisSpectrum.push_back(static_cast<float>(match.hypothesis.at(i)));

    return m_flashMatchScore;
}

} // namespace lar_pandora
