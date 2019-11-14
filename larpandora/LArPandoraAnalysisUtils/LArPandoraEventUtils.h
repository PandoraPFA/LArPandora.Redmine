/**
 *
 * @file larpandora/LArPandoraAnalysisUtils/LArPandoraEventUtils.h
 *
 * @brief Utility containing helpful functions for end users to access products from events
 *
 * @author leigh.howard.whitehead@cern.ch
*/

#ifndef LAR_PANDORA_EVENT_UTILS_H
#define LAR_PANDORA_EVENT_UTILS_H

#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraUtilsBase.h"

#include <string>
#include <vector>

namespace lar_pandora
{

/**
 *
 * @brief LArPandoraEventUtils class
 *
*/
class LArPandoraEventUtils:LArPandoraUtilsBase
{

public:

    static const std::vector<art::Ptr<recob::PFParticle>> GetPFParticles(art::Event const &evt, const std::string &label);

    static const std::vector<art::Ptr<recob::Track>> GetTracks(art::Event const &evt, const std::string &label);
    
    static const std::vector<art::Ptr<recob::Shower>> GetShowers(art::Event const &evt, const std::string &label);

    static const std::vector<art::Ptr<recob::Vertex>> GetVertices(art::Event const &evt, const std::string &label);

    static const std::vector<art::Ptr<recob::SpacePoint>> GetSpacePoints(art::Event const &evt, const std::string &label);

    static const std::vector<art::Ptr<recob::Slice>> GetSlices(art::Event const &evt, const std::string &label);

    static const std::vector<art::Ptr<recob::PFParticle>> GetClearCosmics(art::Event const &evt, const std::string &label);

    static const art::Ptr<recob::PFParticle> GetNeutrino(art::Event const &, const std::string &label);

    static const bool HasNeutrino(art::Event const &evt, const std::string &label);
    

private:

};

} // namespace lar_pandora


#endif // LAR_PANDORA_EVENT_UTILS_H

