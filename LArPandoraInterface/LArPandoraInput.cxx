/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraInput.cxx
 *
 *  @brief  Helper functions for providing inputs to pandora
 */

namespace lar_pandora
{

void LArPandoraInput::CreatePandoraHits2D(const HitVector &hitVector, HitMap &hitMap)
{
    mf::LogDebug("LArPandora") << " *** LArPandora::CreatePandoraHits2D(...) *** " << std::endl;

    // Set up ART services
    art::ServiceHandle<geo::Geometry> theGeometry;
    art::ServiceHandle<util::DetectorProperties> theDetector;

    // Loop over ART hits
    int hitCounter(0);

    PandoraAddressList pandoraAddressList;

    for (HitVector::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::Hit> hit = *iter;
        const geo::WireID hit_WireID(hit->WireID());
        const geo::View_t hit_View(hit->View());

        PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.end();

        try
        {
            const unsigned int volumeID(this->GetPandoraVolumeID(hit_WireID.Cryostat, hit_WireID.TPC));
            pIter = m_pandoraInstanceMap.find(volumeID);
        }
        catch (pandora::StatusCodeException&)
        {
            continue;
        }

        if (m_pandoraInstanceMap.end() == pIter)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        const pandora::Pandora *const pPandora = pIter->second;

        const double hit_Time(hit->PeakTime());
        const double hit_Charge(hit->Integral());
        const double hit_TimeStart(hit->PeakTimeMinusRMS());
        const double hit_TimeEnd(hit->PeakTimePlusRMS());

        double xyz[3];
        theGeometry->Cryostat(hit_WireID.Cryostat).TPC(hit_WireID.TPC).Plane(hit_WireID.Plane).Wire(hit_WireID.Wire).GetCenter(xyz);
        const double y0_cm(xyz[1]);
        const double z0_cm(xyz[2]);

        const double wire_pitch_cm(theGeometry->WirePitch(hit_View)); // cm

        const double xpos_cm(theDetector->ConvertTicksToX(hit_Time, hit_WireID.Plane, hit_WireID.TPC, hit_WireID.Cryostat));
        const double dxpos_cm(std::fabs(theDetector->ConvertTicksToX(hit_TimeEnd, hit_WireID.Plane, hit_WireID.TPC, hit_WireID.Cryostat) -
	    theDetector->ConvertTicksToX(hit_TimeStart, hit_WireID.Plane, hit_WireID.TPC, hit_WireID.Cryostat)));

        const double mips(this->GetMips(hit_Charge, hit_View));

        // Create Pandora CaloHit
        PandoraApi::CaloHit::Parameters caloHitParameters;
        caloHitParameters.m_expectedDirection = pandora::CartesianVector(0., 0., 1.);
        caloHitParameters.m_cellNormalVector = pandora::CartesianVector(0., 0., 1.);
        caloHitParameters.m_cellSize0 = m_dx_cm;
        caloHitParameters.m_cellSize1 = (m_useHitWidths ? dxpos_cm : m_dx_cm);
        caloHitParameters.m_cellThickness = wire_pitch_cm;
        caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
        caloHitParameters.m_time = 0.;
        caloHitParameters.m_nCellRadiationLengths = m_dx_cm / m_rad_cm;
        caloHitParameters.m_nCellInteractionLengths = m_dx_cm / m_int_cm;
        caloHitParameters.m_isDigital = false;
        caloHitParameters.m_hitRegion = pandora::SINGLE_REGION;
        caloHitParameters.m_layer = 0;
        caloHitParameters.m_isInOuterSamplingLayer = false;
        caloHitParameters.m_inputEnergy = hit_Charge;
        caloHitParameters.m_mipEquivalentEnergy = mips;
        caloHitParameters.m_electromagneticEnergy = mips * m_mips_to_gev;
        caloHitParameters.m_hadronicEnergy = mips * m_mips_to_gev;
        caloHitParameters.m_pParentAddress = (void*)((intptr_t)(++hitCounter));

        if (hit_View == geo::kW)
        {
            caloHitParameters.m_hitType = pandora::TPC_VIEW_W;
            const double wpos_cm(z0_cm);
            caloHitParameters.m_positionVector = pandora::CartesianVector(xpos_cm, 0., wpos_cm);
        }
        else if(hit_View == geo::kU)
        {
            caloHitParameters.m_hitType = pandora::TPC_VIEW_U;
            const double upos_cm(lar_content::LArGeometryHelper::GetLArTransformationPlugin(*pPandora)->YZtoU(y0_cm, z0_cm));
            caloHitParameters.m_positionVector = pandora::CartesianVector(xpos_cm, 0., upos_cm);
        }
        else if(hit_View == geo::kV)
        {
            caloHitParameters.m_hitType = pandora::TPC_VIEW_V;
            const double vpos_cm(lar_content::LArGeometryHelper::GetLArTransformationPlugin(*pPandora)->YZtoV(y0_cm, z0_cm));
            caloHitParameters.m_positionVector = pandora::CartesianVector(xpos_cm, 0., vpos_cm);
        }
        else
        {
            mf::LogError("LArPandora") << " --- WARNING: UNKNOWN VIEW !!!  (View=" << hit_View << ")" << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        } 

        // Check for unphysical pulse heights
        if (std::isnan(mips))
        {
            mf::LogError("LArPandora") << " --- WARNING: UNPHYSICAL PULSEHEIGHT !!! (MIPs=" << mips << ")" << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }

        // Store the hit address
        if (hitCounter >= m_uidOffset)
        {
            mf::LogError("LArPandora") << " --- WARNING: TOO MANY HITS !!! (hitCounter=" << hitCounter << ")" << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }

        hitMap[hitCounter] = hit;

        // Create the Pandora hit
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPandora, caloHitParameters));
        pandoraAddressList[pPandora].push_back(hitCounter);
    }

    // Print out results
    for (PandoraInstanceMap::const_iterator pIter1 = m_pandoraInstanceMap.begin(), pIterEnd1 = m_pandoraInstanceMap.end();
        pIter1 != pIterEnd1; ++pIter1)
    {
        const unsigned int volID = pIter1->first;
        const pandora::Pandora *const pPandora = pIter1->second;

        PandoraAddressList::const_iterator pIter2 = pandoraAddressList.find(pPandora);
        const unsigned int numHits((pandoraAddressList.end() == pIter2) ? 0 : pIter2->second.size());

        mf::LogDebug("LArPandora") << "   Number of Pandora 2D Hits [" << volID << "] " << numHits << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraInput::CreatePandoraHits3D(const SpacePointVector &spacePointVector, const SpacePointsToHits &spacePointsToHits,
    SpacePointMap &spacePointMap)
{
    mf::LogDebug("LArPandora") << " *** LArPandora::CreatePandoraHits3D(...) *** " << std::endl;

    // Set up ART services
    art::ServiceHandle<geo::Geometry> theGeometry;
    art::ServiceHandle<util::DetectorProperties> theDetector;
    art::ServiceHandle<util::LArProperties> theLiquidArgon;

    // Loop over ART SpacePoints
    int spacePointCounter(m_uidOffset);

    PandoraAddressList pandoraAddressList;

    for (SpacePointVector::const_iterator iter1 = spacePointVector.begin(), iterEnd1 = spacePointVector.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::SpacePoint> spacepoint = *iter1;

        const double xpos_cm(spacepoint->XYZ()[0]);
        const double ypos_cm(spacepoint->XYZ()[1]);
        const double zpos_cm(spacepoint->XYZ()[2]);

        SpacePointsToHits::const_iterator iter2 = spacePointsToHits.find(spacepoint);
        if (spacePointsToHits.end() == iter2)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        const art::Ptr<recob::Hit> hit = iter2->second;
 
        const geo::View_t hit_View(hit->View());
        const double hit_Charge(hit->Integral());
        
        const double wire_pitch_cm(theGeometry->WirePitch(hit_View));
        const double mips(this->GetMips(hit_Charge, hit_View));

        // Create Pandora CaloHit
        PandoraApi::CaloHit::Parameters caloHitParameters;
        caloHitParameters.m_hitType = pandora::TPC_3D;
        caloHitParameters.m_positionVector = pandora::CartesianVector(xpos_cm, ypos_cm, zpos_cm);
        caloHitParameters.m_cellNormalVector = pandora::CartesianVector(0., 0., 1.);
        caloHitParameters.m_expectedDirection = pandora::CartesianVector(0., 0., 1.);
        caloHitParameters.m_cellSize0 = m_dx_cm;
        caloHitParameters.m_cellSize1 = m_dx_cm;
        caloHitParameters.m_cellThickness = wire_pitch_cm;
        caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
        caloHitParameters.m_time = 0.;
        caloHitParameters.m_nCellRadiationLengths = m_dx_cm / m_rad_cm;
        caloHitParameters.m_nCellInteractionLengths = m_dx_cm / m_int_cm;
        caloHitParameters.m_isDigital = false;
        caloHitParameters.m_hitRegion = pandora::SINGLE_REGION;
        caloHitParameters.m_layer = 0;
        caloHitParameters.m_isInOuterSamplingLayer = false;
        caloHitParameters.m_inputEnergy = hit_Charge;
        caloHitParameters.m_mipEquivalentEnergy = mips;
        caloHitParameters.m_electromagneticEnergy = mips * m_mips_to_gev;
        caloHitParameters.m_hadronicEnergy = mips * m_mips_to_gev;
        caloHitParameters.m_pParentAddress = (void*)((intptr_t)(++spacePointCounter));

        // Check for unphysical pulse heights
        if (std::isnan(mips))
        {
            mf::LogError("LArPandora") << " --- WARNING: UNPHYSICAL PULSEHEIGHT !!! (MIPs=" << mips << ")" << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }

        // Store the hit address
        spacePointMap[spacePointCounter] = spacepoint;

        // Create the Pandora hit
        for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end();
            pIter != pIterEnd; ++pIter)
        {
            const pandora::Pandora *const pPandora = pIter->second; 

            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pPandora, caloHitParameters));
            pandoraAddressList[pPandora].push_back(spacePointCounter);
        }
    } 
 
    // Print out results
    for (PandoraInstanceMap::const_iterator pIter1 = m_pandoraInstanceMap.begin(), pIterEnd1 = m_pandoraInstanceMap.end();
        pIter1 != pIterEnd1; ++pIter1)
    {
        const unsigned int volID = pIter1->first;
        const pandora::Pandora *const pPandora = pIter1->second;

        PandoraAddressList::const_iterator pIter2 = pandoraAddressList.find(pPandora);
        const unsigned int numHits((pandoraAddressList.end() == pIter2) ? 0 : pIter2->second.size());

        mf::LogDebug("LArPandora") << "   Number of Pandora 3D Hits [" << volID << "] " << numHits << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraInput::CreatePandoraParticles(const MCTruthToMCParticles &truthToParticleMap, const MCParticlesToMCTruth &particleToTruthMap)
{
    mf::LogDebug("LArPandora") << " *** LArPandora::CreatePandoraParticles(...) *** " << std::endl;

    // Make indexed list of MC particles
    MCParticleMap particleMap;

    for (MCParticlesToMCTruth::const_iterator iter = particleToTruthMap.begin(), iterEnd = particleToTruthMap.end(); 
        iter != iterEnd; ++iter)
    {
        const art::Ptr<simb::MCParticle> particle = iter->first;
        particleMap[particle->TrackId()] = particle;
    }

    // Loop over MC truth objects
    int neutrinoCounter(0);

    lar_content::LArMCParticleFactory mcParticleFactory;

    for (MCTruthToMCParticles::const_iterator iter1 = truthToParticleMap.begin(), iterEnd1 = truthToParticleMap.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<simb::MCTruth> truth = iter1->first;

        if (truth->NeutrinoSet())
        {
            const simb::MCNeutrino neutrino(truth->GetNeutrino());
            ++neutrinoCounter;

            if (neutrinoCounter >= m_uidOffset)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            const int neutrinoID(neutrinoCounter + 4 * m_uidOffset);

            // Create Pandora 3D MC Particle
            lar_content::LArMCParticleParameters mcParticleParameters;
            mcParticleParameters.m_nuanceCode = neutrino.InteractionType();
            mcParticleParameters.m_energy = neutrino.Nu().E();
            mcParticleParameters.m_momentum = pandora::CartesianVector(neutrino.Nu().Px(), neutrino.Nu().Py(), neutrino.Nu().Pz());
            mcParticleParameters.m_vertex = pandora::CartesianVector(neutrino.Nu().Vx(), neutrino.Nu().Vy(), neutrino.Nu().Vz());
            mcParticleParameters.m_endpoint = pandora::CartesianVector(neutrino.Nu().Vx(), neutrino.Nu().Vy(), neutrino.Nu().Vz());
            mcParticleParameters.m_particleId = neutrino.Nu().PdgCode();
            mcParticleParameters.m_mcParticleType = pandora::MC_3D;
            mcParticleParameters.m_pParentAddress = (void*)((intptr_t)neutrinoID);

            for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end();
                pIter != pIterEnd; ++pIter)
            {
                const pandora::Pandora *const pPandora = pIter->second;
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPandora, mcParticleParameters, mcParticleFactory));
            }

            // Loop over associated particles
            const MCParticleVector &particleVector = iter1->second;
            
            for (MCParticleVector::const_iterator iter2 = particleVector.begin(), iterEnd2 = particleVector.end(); iter2 != iterEnd2; ++iter2)
            {
                const art::Ptr<simb::MCParticle> particle = *iter2;
                const int trackID(particle->TrackId());

                // Mother/Daughter Links
                if (particle->Mother() == 0)
                {
                    for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end();
                        pIter != pIterEnd; ++pIter)
                    {
                        const pandora::Pandora *const pPandora = pIter->second;
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*pPandora,
                            (void*)((intptr_t)neutrinoID), (void*)((intptr_t)trackID)));
                    }
                }
            }
        }
    }

    mf::LogDebug("LArPandora") << "   Number of Pandora neutrinos: " << neutrinoCounter << std::endl;


    // Loop over G4 particles
    int particleCounter(0);

    for (MCParticleMap::const_iterator iterI = particleMap.begin(), iterEndI = particleMap.end(); iterI != iterEndI; ++iterI)
    {
        const art::Ptr<simb::MCParticle> particle = iterI->second;

        if (particle->TrackId() != iterI->first)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        if (particle->TrackId() >= m_uidOffset)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        ++particleCounter;

        // Loop over drift volumes
        for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end();
            pIter != pIterEnd; ++pIter)
        {
            const unsigned int pVolume  = pIter->first;
            const pandora::Pandora *const pPandora = pIter->second;

            // Find start and end trajectory points
            int firstT(-1), lastT(-1);
            this->GetTrueStartAndEndPoints(pVolume, particle, firstT, lastT);

            if (firstT < 0 && lastT < 0)
            {
                firstT = 0; lastT = 0;
            }

            // Lookup position and kinematics at start and end points
            const float vtxX(particle->Vx(firstT));
            const float vtxY(particle->Vy(firstT));
            const float vtxZ(particle->Vz(firstT));

            const float endX(particle->Vx(lastT));
            const float endY(particle->Vy(lastT));
            const float endZ(particle->Vz(lastT));

            const float pX(particle->Px(firstT));
            const float pY(particle->Py(firstT));
            const float pZ(particle->Pz(firstT));
            const float E(particle->E(firstT));

            // Create 3D Pandora MC Particle
            lar_content::LArMCParticleParameters mcParticleParameters;
            mcParticleParameters.m_nuanceCode = 0;
            mcParticleParameters.m_energy = E;
            mcParticleParameters.m_particleId = particle->PdgCode();
            mcParticleParameters.m_momentum = pandora::CartesianVector(pX, pY, pZ);
            mcParticleParameters.m_vertex = pandora::CartesianVector(vtxX, vtxY, vtxZ);
            mcParticleParameters.m_endpoint = pandora::CartesianVector(endX, endY, endZ);
            mcParticleParameters.m_mcParticleType = pandora::MC_3D;
            mcParticleParameters.m_pParentAddress = (void*)((intptr_t)particle->TrackId());
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPandora, mcParticleParameters, mcParticleFactory));

            // Create Mother/Daughter Links between 3D MC Particles
            const int id_mother(particle->Mother());
            MCParticleMap::const_iterator iterJ = particleMap.find(id_mother);

            if (iterJ != particleMap.end())
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*pPandora,
                    (void*)((intptr_t)id_mother), (void*)((intptr_t)particle->TrackId())));  
        }
    }

    mf::LogDebug("LArPandora") << "   Number of Pandora particles: " << particleCounter << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraInput::CreatePandoraParticles2D(const MCParticleVector &particleVector)
{
    lar_content::LArMCParticleFactory mcParticleFactory;

    for (MCParticleVector::const_iterator iter = particleVector.begin(), iterEnd = particleVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<simb::MCParticle> particle = *iter;

        if (particle->TrackId() >= m_uidOffset)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        // Loop over drift volumes
        for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end();
            pIter != pIterEnd; ++pIter)
        {
            const unsigned int pVolume  = pIter->first;
            const pandora::Pandora *const pPandora = pIter->second;

            // Find start and end trajectory points
            int firstT(-1), lastT(-1);
            bool foundStartAndEndPoints(false);
            this->GetTrueStartAndEndPoints(pVolume, particle, firstT, lastT);

            if (firstT >= 0 && lastT >= 0)
            {
                foundStartAndEndPoints = true;
            }
            else
            {
                firstT = 0; lastT = 0;
            }

            // Lookup position and kinematics at start and end points
            const float vtxX(particle->Vx(firstT));
            const float vtxY(particle->Vy(firstT));
            const float vtxZ(particle->Vz(firstT));

            const float endX(particle->Vx(lastT));
            const float endY(particle->Vy(lastT));
            const float endZ(particle->Vz(lastT));

            const float pX(particle->Px(firstT));
            const float pY(particle->Py(firstT));
            const float pZ(particle->Pz(firstT));
            const float E(particle->E(firstT));

            // Create 2D Pandora MC Particles for Event Display
            if (!foundStartAndEndPoints)
                continue;

            const float dx(endX - vtxX);
            const float dy(endY - vtxY);
            const float dz(endZ - vtxZ);
            const float dw(lar_content::LArGeometryHelper::GetWireZPitch(*pPandora));

            if (dx * dx + dy * dy + dz * dz < 0.5 * dw * dw)
                continue;

            // Add in T0 to 2D projections
            const float vtxX0(this->GetTrueX0(particle, firstT));
            const float endX0(this->GetTrueX0(particle, lastT));

            // Create 2D Pandora MC Particles for each view
            lar_content::LArMCParticleParameters mcParticleParameters;
            mcParticleParameters.m_nuanceCode = 0;
            mcParticleParameters.m_energy = E;
            mcParticleParameters.m_particleId = particle->PdgCode();

            // Create U projection
            mcParticleParameters.m_momentum = pandora::CartesianVector(pX, 0.f,
                lar_content::LArGeometryHelper::GetLArTransformationPlugin(*pPandora)->PYPZtoPU(pY, pZ));
            mcParticleParameters.m_vertex = pandora::CartesianVector(vtxX + vtxX0, 0.f,
                lar_content::LArGeometryHelper::GetLArTransformationPlugin(*pPandora)->YZtoU(vtxY, vtxZ));
            mcParticleParameters.m_endpoint = pandora::CartesianVector(endX + endX0,  0.f,
                lar_content::LArGeometryHelper::GetLArTransformationPlugin(*pPandora)->YZtoU(endY, endZ));
            mcParticleParameters.m_mcParticleType = pandora::MC_VIEW_U;
            mcParticleParameters.m_pParentAddress = (void*)((intptr_t)(particle->TrackId() + 1 * m_uidOffset));
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPandora, mcParticleParameters, mcParticleFactory));

            // Create V projection
            mcParticleParameters.m_momentum = pandora::CartesianVector(pX, 0.f,
                lar_content::LArGeometryHelper::GetLArTransformationPlugin(*pPandora)->PYPZtoPV(pY, pZ));
            mcParticleParameters.m_vertex = pandora::CartesianVector(vtxX + vtxX0, 0.f,
                lar_content::LArGeometryHelper::GetLArTransformationPlugin(*pPandora)->YZtoV(vtxY, vtxZ));
            mcParticleParameters.m_endpoint = pandora::CartesianVector(endX + endX0,  0.f,
                lar_content::LArGeometryHelper::GetLArTransformationPlugin(*pPandora)->YZtoV(endY, endZ));
            mcParticleParameters.m_mcParticleType = pandora::MC_VIEW_V;
            mcParticleParameters.m_pParentAddress = (void*)((intptr_t)(particle->TrackId() + 2 * m_uidOffset));
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPandora, mcParticleParameters, mcParticleFactory));

            // Create W projection
            mcParticleParameters.m_momentum = pandora::CartesianVector(pX, 0.f, pZ);
            mcParticleParameters.m_vertex = pandora::CartesianVector(vtxX + vtxX0, 0.f, vtxZ);
            mcParticleParameters.m_endpoint = pandora::CartesianVector(endX + endX0,  0.f, endZ);
            mcParticleParameters.m_mcParticleType = pandora::MC_VIEW_W;
            mcParticleParameters.m_pParentAddress = (void*)((intptr_t)(particle->TrackId() + 3 * m_uidOffset));
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPandora, mcParticleParameters, mcParticleFactory));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraInput::CreatePandoraLinks2D(const HitMap &hitMap, const HitsToTrackIDEs &hitToParticleMap)
{
    mf::LogDebug("LArPandora") << " *** LArPandora::CreatePandoraLinks(...) *** " << std::endl;

    for (HitMap::const_iterator iterI = hitMap.begin(), iterEndI = hitMap.end(); iterI != iterEndI ; ++iterI)
    {        
        // Get the ART hit
        const int hitID(iterI->first);
        const art::Ptr<recob::Hit> hit(iterI->second);
        const geo::WireID hit_WireID(hit->WireID());

        // Get the Pandora instance and Pandora hit
        PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.end();

        try
        {
            const unsigned int volumeID(this->GetPandoraVolumeID(hit_WireID.Cryostat, hit_WireID.TPC));
            pIter = m_pandoraInstanceMap.find(volumeID);
        }
        catch (pandora::StatusCodeException&)
        {
            continue;
        }

        if (m_pandoraInstanceMap.end() == pIter)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        const pandora::Pandora *const pPandora = pIter->second;

        // Get list of associated MC particles
        HitsToTrackIDEs::const_iterator iterJ = hitToParticleMap.find(hit);

        if (hitToParticleMap.end() == iterJ)
            continue;

        const TrackIDEVector &trackCollection = iterJ->second;

        if (trackCollection.size() == 0)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        // Create links between hits and MC particles
        for (unsigned int k = 0; k < trackCollection.size(); ++k)
        {
            const sim::TrackIDE trackIDE(trackCollection.at(k));
            const int trackID(std::abs(trackIDE.trackID)); // TODO: Find out why std::abs is needed
            const float energyFrac(trackIDE.energyFrac);

            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(*pPandora,
                (void*)((intptr_t)hitID), (void*)((intptr_t)trackID), energyFrac));
        }
    }
}

} // namespace lar_pandora
