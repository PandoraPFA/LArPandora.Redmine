/**
 *
 * @file larpandora/LArPandoraAnalysisUtils/LArPandoraUtilsBase.h
 *
 * @brief Base class containing functionality to extract products from the event
 *
 * @author leigh.howard.whitehead@cern.ch
*/

#ifndef LAR_PANDORA_UTILS_BASE_H
#define LAR_PANDORA_UTILS_BASE_H

#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

// Access the type defs defined in the helper
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <string>
#include <vector>

namespace lar_pandora
{

/**
 *
 * @brief LArPandoraUtilsBase class containing some template functions
 *
*/
class LArPandoraUtilsBase
{

protected:

    template <typename T> static void GetProductVector(art::Event const &evt, const std::string &label, std::vector<art::Ptr<T>> &productVector);
    template <typename T, typename U> static void GetAssocProductVector(const art::Ptr<U> &part, art::Event const &evt, const std::string &label, const std::string &assocLabel, std::vector<art::Ptr<T>> &productVector);

};

    // Implementation of the template function to get the products from the event
    template <typename T> void LArPandoraUtilsBase::GetProductVector(art::Event const &evt, const std::string &label, std::vector<art::Ptr<T>> &productVector)
    {

        art::Handle<std::vector<T>> theseProds;
        evt.getByLabel(label,theseProds);

        if (!theseProds.isValid())
        {
            mf::LogError("LArPandora") << " Failed to find product with label " << label << " ... returning empty vector" << std::endl;
            productVector = std::vector<art::Ptr<T>>();
            return;
        }

        // We need to convert these to art pointers
        std::vector<art::Ptr<T>> pointerVersion;
        for(unsigned int i = 0; theseProds->size(); ++i){
          art::Ptr<T> ptr(theseProds,i);
          pointerVersion.push_back(ptr);
        }

        productVector = pointerVersion;
    }

    // Implementation of the template function to get the associated products from the event
    template <typename T, typename U> void LArPandoraUtilsBase::GetAssocProductVector(const art::Ptr<U> &prod, art::Event const &evt, const std::string &label, const std::string &assocLabel, std::vector<art::Ptr<T>> &productVector)
    {

        art::Handle<std::vector<U>> products;
        evt.getByLabel(label,products);

        if (!products.isValid())
        {
            mf::LogError("LArPandora") << " Failed to find product with label " << label << " ... returning empty vector" << std::endl;
            productVector = std::vector<art::Ptr<T>>();
            return;
        }

        const art::FindManyP<T> findParticleAssocs(products,evt,assocLabel);

        const std::vector<art::Ptr<T>> theseAssocs = findParticleAssocs.at(prod.key());

        productVector = theseAssocs;
    }


} // namespace lar_pandora


#endif // LAR_PANDORA_PFPARTICLE_UTILS_H

