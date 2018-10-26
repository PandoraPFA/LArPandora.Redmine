/**
 *  @file   larpandora/LArPandoraEventBuilding/SliceIdTools/FlashNeutrinoId_tool.cc
 *
 *  @brief  implementation of the flash based neutrino id tool
 */

#include "larpandora/LArPandoraEventBuilding/SliceIdTools/FlashNeutrinoId_tool.h"

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
        m_nuCCNC = mcNeutrino.CCNC();
        const auto nuMCParticle(mcNeutrino.Nu());
        const auto leptonMCParticle(mcNeutrino.Lepton());

        m_nuEnergy = nuMCParticle.E();
        m_leptonEnergy = leptonMCParticle.E();
        m_nuVertexX = nuMCParticle.Vx();
        m_nuVertexY = nuMCParticle.Vy();
        m_nuVertexZ = nuMCParticle.Vz();
        m_nuTime = nuMCParticle.T();
        m_nuPdgCode = nuMCParticle.PdgCode();
    
        if (slices.size() != sliceMetadata.size())
            throw cet::exception("FlashNeutrinoId") << "The number of slice metadata doesn't match the number of slices" << std::endl;
    }

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
    if (sliceCandidates.empty())
        return;

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
    m_deltaY = (m_centerY - beamFlash.m_centerY);
    m_deltaZ = (m_centerZ - beamFlash.m_centerZ);
    m_deltaYSigma = m_deltaY / beamFlash.m_widthY;
    m_deltaZSigma = m_deltaZ / beamFlash.m_widthZ;
    m_chargeToLightRatio = m_totalCharge / beamFlash.m_totalPE;  // TODO ATTN check if this should be total PE or max PE. Code differs from technote
    
    // Check if the slice passes the pre-selection cuts
    m_passesPrecuts = (std::abs(m_deltaY) < maxDeltaY                           &&
                       std::abs(m_deltaZ) < maxDeltaZ                           &&
                       std::abs(m_deltaYSigma) < maxDeltaYSigma                 &&
                       std::abs(m_deltaZSigma) < maxDeltaZSigma                 &&
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
