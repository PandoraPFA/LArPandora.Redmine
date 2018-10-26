/**
 *  @file   larpandora/LArPandoraEventBuilding/SliceIdTools/FlashNeutrinoId_tool.h
 *
 *  @brief  header of the flash based neutrino id tool
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
    int                                     m_nuCCNC;                ///< Charged current or neutral current?
    float                                   m_nuEnergy;            ///< The true neutrino energy
    float                                   m_leptonEnergy;        ///< The true energy of the lepton coming from the CC interaction
    float                                   m_nuVertexX;           ///< The true neutrino vertex X position
    float                                   m_nuVertexY;           ///< The true neutrino vertex Y position
    float                                   m_nuVertexZ;           ///< The true neutrino vertex Z position
    float                                   m_nuTime;              ///< The time of the true neutrino interaction
    int                                     m_nuPdgCode;           ///< The true neutrino pdg code
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
    if (m_hasMCNeutrino)
    {
        // Truth MC information about the neutrino
        m_pEventTree->Branch("nuInteractionType", &m_nuInteractionType                   , "nuInteractionType/I");
        m_pEventTree->Branch("nuCCNC"           , &m_nuCCNC                              , "nuCCNC/I");
        m_pEventTree->Branch("nuEnergy"         , &m_nuEnergy                            , "nuEnergy/F");
        m_pEventTree->Branch("leptonEnergy"     , &m_leptonEnergy                        , "leptonEnergy/F");
        m_pEventTree->Branch("nuInteractionTime", &m_nuTime                              , "nuInteractionTime/F");
        m_pEventTree->Branch("nuPdgCode"        , &m_nuPdgCode                           , "nuPdgCode/I");
        m_pEventTree->Branch("nuVertexX"        , &m_nuVertexX                           , "nuVertexX/F");
        m_pEventTree->Branch("nuVertexY"        , &m_nuVertexY                           , "nuVertexY/F");
        m_pEventTree->Branch("nuVertexZ"        , &m_nuVertexZ                           , "nuVertexZ/F");
    }

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
        m_pSliceTree->Branch("nuInteractionType", &m_nuInteractionType                   , "nuInteractionType/I");
        m_pSliceTree->Branch("nuCCNC"           , &m_nuCCNC                              , "nuCCNC/I");
        m_pSliceTree->Branch("nuEnergy"         , &m_nuEnergy                            , "nuEnergy/F");
        m_pSliceTree->Branch("leptonEnergy"     , &m_leptonEnergy                        , "leptonEnergy/F");
        m_pSliceTree->Branch("nuInteractionTime", &m_nuTime                              , "nuInteractionTime/F");
        m_pSliceTree->Branch("nuPdgCode"        , &m_nuPdgCode                           , "nuPdgCode/I");
        m_pSliceTree->Branch("nuVertexX"        , &m_nuVertexX                           , "nuVertexX/F");
        m_pSliceTree->Branch("nuVertexY"        , &m_nuVertexY                           , "nuVertexY/F");
        m_pSliceTree->Branch("nuVertexZ"        , &m_nuVertexZ                           , "nuVertexZ/F");
    }
}