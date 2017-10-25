/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraSlices.cxx
 *
 *  @brief  A mapping between two different reconstruction hypotheses on the same input hit collection
 */

#include "larpandora/LArPandoraEventBuilding/LArPandoraSlices.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

LArPandoraSlices::LArPandoraSlices(art::EDProducer *pProducer, art::Event *pEvent, const std::string & crRecoProducerLabel, const std::string & nuRecoProducerLabel, const std::string & hitProducerLabel) :
    m_pProducer(pProducer),
    m_pEvent(pEvent),
    m_crRecoProducerLabel(crRecoProducerLabel),
    m_nuRecoProducerLabel(nuRecoProducerLabel),
    m_hitProducerLabel(hitProducerLabel),
    m_doesEventContainNeutrino(false)
{
    this->IdentifySlices(); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector< LArPandoraSlices::SliceId > LArPandoraSlices::GetSlices()
{
    std::vector< SliceId > slices(m_crSlicePFParticles.size());
    std::iota(std::begin(slices), std::end(slices), 0);
    return slices;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector LArPandoraSlices::GetSliceAsCR(const SliceId & id)
{
    if (m_crSlicePFParticles.count(id) == 0)
        throw cet::exception("LArPandora") << " LArPandoraSlices::GetSliceAsCR -- Slice Id " << id << " is out of bounds.";

    return m_crSlicePFParticles.at(id);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector LArPandoraSlices::GetSliceAsNu(const SliceId & id)
{
    if (m_nuSlicePFParticles.count(id) == 0)
        throw cet::exception("LArPandora") << " LArPandoraSlices::GetSliceAsNu -- Slice Id " << id << " is out of bounds.";

    return m_nuSlicePFParticles.at(id);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::IdSliceAsNu(const SliceId & id)
{
    if (m_nuSlicePFParticles.count(id) == 0)
        throw cet::exception("LArPandora") << " LArPandoraSlices::IdSliceAsNu -- Can't identify slice " << id << " as the neutrino. Slice Id out of bounds.";

    m_doesEventContainNeutrino = true;
    m_nuSliceId = id;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::WriteTags()
{

    if (m_nuSlicePFParticles.size() != m_crSlicePFParticles.size())
        throw cet::exception("LArPandora") << " LArPandoraSlices::WriteTags -- Malformed slices.";

    std::unique_ptr< std::vector< anab::CosmicTag > >                    outputTags(new std::vector< anab::CosmicTag >);
    std::unique_ptr< art::Assns< recob::PFParticle, anab::CosmicTag > >  outputAssn(new art::Assns< recob::PFParticle, anab::CosmicTag >);

    for (std::map< SliceId, PFParticleVector >::const_iterator it = m_nuSlicePFParticles.begin(); it != m_nuSlicePFParticles.end(); ++it) 
    {
        SliceId sliceId = it->first;
        PFParticleVector crPFParticles = this->GetSliceAsCR(sliceId);
        PFParticleVector nuPFParticles = this->GetSliceAsNu(sliceId);

        if (!m_doesEventContainNeutrino)
        {
            this->WriteTag(false, crPFParticles, outputTags, outputAssn);
            this->WriteTag(false, nuPFParticles, outputTags, outputAssn);
        }
        else 
        {
            this->WriteTag((sliceId == m_nuSliceId), crPFParticles, outputTags, outputAssn);
            this->WriteTag((sliceId == m_nuSliceId), nuPFParticles, outputTags, outputAssn);
        }
    }
    
    m_pEvent->put(std::move(outputTags));
    m_pEvent->put(std::move(outputAssn));

}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::WriteTag(const bool & shouldTagAsNeutrino, const PFParticleVector & pfParticleVector, std::unique_ptr< std::vector< anab::CosmicTag > > & outputTags, std::unique_ptr< art::Assns< recob::PFParticle, anab::CosmicTag > > &  outputAssn) 
{
    const lar::PtrMaker< anab::CosmicTag > makeCRTagPtr(*m_pEvent, *m_pProducer);

    for (const art::Ptr< recob::PFParticle > & part : pfParticleVector)
    {
        std::vector< float > endPoints(3, std::numeric_limits< float >::max());
        anab::CosmicTag tag(endPoints, endPoints, std::numeric_limits< float >::max(), (shouldTagAsNeutrino ? anab::kNotTagged : anab::kUnknown));

        outputTags->emplace_back(tag);
        art::Ptr< anab::CosmicTag > pTag(makeCRTagPtr(outputTags->size() - 1));
        util::CreateAssn(*m_pProducer, *m_pEvent, pTag, part, *outputAssn);  
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::IdentifySlices()
{
    // Collect inputs
    PFParticleVector           nuParticles;
    SpacePointVector           nuSpacePoints;
    PFParticlesToSpacePoints   nuParticlesToSpacePoints;
    SpacePointsToHits          nuSpacePointsToHits;
    LArPandoraHelper::CollectPFParticles(*m_pEvent, m_nuRecoProducerLabel, nuParticles, nuParticlesToSpacePoints);
    LArPandoraHelper::CollectSpacePoints(*m_pEvent, m_nuRecoProducerLabel, nuSpacePoints, nuSpacePointsToHits);

    PFParticleVector           crParticles;
    SpacePointVector           crSpacePoints;
    PFParticlesToSpacePoints   crParticlesToSpacePoints;
    SpacePointsToHits          crSpacePointsToHits;
    LArPandoraHelper::CollectPFParticles(*m_pEvent, m_crRecoProducerLabel, crParticles, crParticlesToSpacePoints);
    LArPandoraHelper::CollectSpacePoints(*m_pEvent, m_crRecoProducerLabel, crSpacePoints, crSpacePointsToHits);
    
    // Get Final state PFParticle -> Hit maps
    PFParticlesToHits  nuParticlesToHits;
    HitsToPFParticles  nuHitsToParticles;
    LArPandoraHelper::BuildPFParticleHitMaps(nuParticles, nuParticlesToSpacePoints, nuSpacePointsToHits, nuParticlesToHits, nuHitsToParticles, LArPandoraHelper::kAddDaughters);
    
    PFParticlesToHits  crParticlesToHits;
    HitsToPFParticles  crHitsToParticles;
    LArPandoraHelper::BuildPFParticleHitMaps(crParticles, crParticlesToSpacePoints, crSpacePointsToHits, crParticlesToHits, crHitsToParticles, LArPandoraHelper::kAddDaughters);
    
    // Get top-level and final state PFParticles
    PFParticleVector  nuTopLevelParticles;
    LArPandoraHelper::SelectNeutrinoPFParticles(nuParticles, nuTopLevelParticles);
    
    PFParticleVector  nuFinalStateParticles;
    LArPandoraHelper::SelectFinalStatePFParticles(nuParticles, nuFinalStateParticles);
    
    PFParticleVector  crFinalStateParticles;
    LArPandoraHelper::SelectFinalStatePFParticles(crParticles, crFinalStateParticles);

    // Get PFParticle maps
    PFParticleMap     nuPFParticleIdMap;
    PFParticleMap     crPFParticleIdMap;
    this->GetPFParticleIdMap(nuParticles, nuPFParticleIdMap);
    this->GetPFParticleIdMap(crParticles, crPFParticleIdMap);

    // Make a new slice for each top-level neutrino particle
    std::map< art::Ptr< recob::PFParticle>, SliceId > nuFinalStateParticlesToSlice;
    this->MakeSlicePerNeutrino(nuPFParticleIdMap, nuTopLevelParticles, nuFinalStateParticles, nuFinalStateParticlesToSlice);

    // Add CR PFParticles to slices if they share any hits, and make new slice otherwise
    this->AddCRParticlesToSlices(crPFParticleIdMap, crFinalStateParticles, crParticlesToHits, nuHitsToParticles, nuFinalStateParticlesToSlice);

}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::GetPFParticleIdMap(const PFParticleVector & inputParticles, PFParticleMap & outputMap)
{
    for (const art::Ptr< recob::PFParticle > & part : inputParticles) 
        outputMap.insert(PFParticleMap::value_type(part->Self(), part));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::MakeSlicePerNeutrino(const PFParticleMap & nuPFParticleIdMap, const PFParticleVector & nuTopLevelParticles, const PFParticleVector & nuFinalStateParticles, std::map< art::Ptr< recob::PFParticle>, SliceId > &  nuFinalStateParticlesToSlice)
{
    for (SliceId id = 0; id < nuTopLevelParticles.size(); id++) 
    {
        const art::Ptr< recob::PFParticle > topLevelParticle = nuTopLevelParticles[id];

        PFParticleVector  nuParticlesInSlice;
        PFParticleVector  crParticlesInSlice;
        this->CollectDaughters(nuPFParticleIdMap, topLevelParticle, nuParticlesInSlice);

        m_nuSlicePFParticles.insert(std::map< SliceId, PFParticleVector >::value_type(id, nuParticlesInSlice)); 
        m_crSlicePFParticles.insert(std::map< SliceId, PFParticleVector >::value_type(id, crParticlesInSlice)); 
        
        for (const art::Ptr< recob::PFParticle > & part : nuFinalStateParticles)
        {

            if (LArPandoraHelper::GetParentPFParticle(nuPFParticleIdMap, part) != topLevelParticle) continue;

            // Now have a final state daughter of the topLevelParticle
            nuFinalStateParticlesToSlice.insert(std::map< art::Ptr< recob::PFParticle>, SliceId >::value_type(part, id));
        }
    } 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::CollectDaughters(const PFParticleMap & pfParticleMap, const art::Ptr< recob::PFParticle > & part, PFParticleVector & daughterParticles)
{
    if (std::find(daughterParticles.begin(), daughterParticles.end(), part) == daughterParticles.end()) 
        daughterParticles.push_back(part);
    
    for (const size_t & daughterId : part->Daughters())
    {
        if (pfParticleMap.find(daughterId) == pfParticleMap.end())
            throw cet::exception("LArPandora") << " LArPandoraSlices::CollectDaughters - Can't find any daughter PFParticle in ID map." << std::endl;
            
        this->CollectDaughters(pfParticleMap, pfParticleMap.at(daughterId), daughterParticles);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSlices::AddCRParticlesToSlices(const PFParticleMap & crPFParticleIdMap, const PFParticleVector & crFinalStateParticles, const PFParticlesToHits & crParticlesToHits, const HitsToPFParticles & nuHitsToParticles, const std::map< art::Ptr< recob::PFParticle>, SliceId > & nuFinalStateParticlesToSlice)
{
    for (const art::Ptr< recob::PFParticle> & part : crFinalStateParticles)
    {

        if (crParticlesToHits.find(part) == crParticlesToHits.end()) 
            throw cet::exception("LArPandora") << " LArPandoraSlices::AddCRParticlesToSlices - Can't find any hits for supplied PFParticle." << std::endl;

        PFParticleVector daughterParticles;
        this->CollectDaughters(crPFParticleIdMap, part, daughterParticles);

        // Assume this final state does not belong to an existing slice
        SliceId id = m_nuSlicePFParticles.size();

        // Use its hits to find a matching exisiting slice
        for (const art::Ptr< recob::Hit > & hit : crParticlesToHits.at(part))
        {

            if (nuHitsToParticles.find(hit) == nuHitsToParticles.end()) continue;

            if (nuFinalStateParticlesToSlice.find(nuHitsToParticles.at(hit)) == nuFinalStateParticlesToSlice.end())
                throw cet::exception("LArPandora") << " LArPandoraSlices::AddCRParticlesToSlices - Can't find slice associated with supplied final state PFParticle." << std::endl;

            // Found the slice 
            id = nuFinalStateParticlesToSlice.at(nuHitsToParticles.at(hit));
            break;
        }

        // New slice
        if (id == m_nuSlicePFParticles.size())
        {
            PFParticleVector  nuParticlesInSlice;
            m_nuSlicePFParticles.insert(std::map< SliceId, PFParticleVector >::value_type(id, nuParticlesInSlice)); 
            m_crSlicePFParticles.insert(std::map< SliceId, PFParticleVector >::value_type(id, daughterParticles)); 
        }
        // Existing slice
        else
        {
            if (m_nuSlicePFParticles.find(id) == m_nuSlicePFParticles.end() || m_crSlicePFParticles.find(id) == m_crSlicePFParticles.end())
                throw cet::exception("LArPandora") << " LArPandoraSlices::AddCRParticlesToSlices - Invalid slice found : " << id << "." << std::endl;

            m_crSlicePFParticles.at(id).insert(m_crSlicePFParticles.at(id).end(), daughterParticles.begin(), daughterParticles.end());
        }
    }
}

} // namespace lar_pandora
