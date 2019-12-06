/**
 *
 * @file larpandora/LArPandoraAnalysisUtils/LArPandoraUtilsBase.h
 *
 * @brief Base class containing functionality to extract products from the event
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
//    template <typename T> static void GetProductVector(const art::Event &evt, const std::string &label, std::vector<art::Ptr<T>> &productVector);
//    template <typename T, typename U> static void GetAssocProductVector(const art::Ptr<U> &part, const art::Event &evt, const std::string &label, const std::string &assocLabel, std::vector<art::Ptr<T>> &productVector);

    template <typename T> static std::vector<art::Ptr<T>> GetProductVector(const art::Event &evt, const std::string &label);
    template <typename T, typename U> static std::vector<art::Ptr<T>> GetAssocProductVector(const art::Ptr<U> &part, const art::Event &evt, const std::string &label, const std::string &assocLabel);
    template <typename T, typename U> static art::Ptr<T> GetAssocProduct(const art::Ptr<U> &part, const art::Event &evt, const std::string &label, const std::string &assocLabel); 

};

// Implementation of the template function to get the products from the event
template <typename T> std::vector<art::Ptr<T>> LArPandoraUtilsBase::GetProductVector(const art::Event &evt, const std::string &label)
{
    art::Handle<std::vector<T>> theseProds;
    bool success = evt.getByLabel(label,theseProds);

    if (!success)
    {
        mf::LogError("LArPandora") << " Failed to find product with label " << label << " ... returning empty vector" << std::endl;
        return std::vector<art::Ptr<T>>();
    }

    // We need to convert these to art pointers
    std::vector<art::Ptr<T>> productVector;
    for(unsigned int productIndex = 0; productIndex < theseProds->size(); ++productIndex){
      art::Ptr<T> ptr(theseProds,productIndex);
      productVector.push_back(ptr);
    }
    return productVector;
}

// Implementation of the template function to get the associated products from the event
template <typename T, typename U> std::vector<art::Ptr<T>> LArPandoraUtilsBase::GetAssocProductVector(const art::Ptr<U> &pProd, const art::Event &evt, const std::string &label, const std::string &assocLabel)
{
    art::Handle<std::vector<U>> products;
    bool success = evt.getByLabel(label,products);

    if (!success)
    {
        mf::LogError("LArPandora") << " Failed to find product with label " << label << " ... returning empty vector" << std::endl;
        return std::vector<art::Ptr<T>>();
    }

    const art::FindManyP<T> findParticleAssocs(products,evt,assocLabel);

    return findParticleAssocs.at(pProd.key());
}

// Implementation of the template function to get the associated product from the event
template <typename T, typename U> art::Ptr<T> LArPandoraUtilsBase::GetAssocProduct(const art::Ptr<U> &pProd, const art::Event &evt, const std::string &label, const std::string &assocLabel)
{   
    std::vector<art::Ptr<T>> associatedProducts = LArPandoraUtilsBase::GetAssocProductVector<T>(pProd,evt,label,assocLabel); 
    if (associatedProducts.empty())
    {
        throw cet::exception("LArPandora") << "LArPandoraUtilsBase::GetShower --- No associated object found";
    }
    return associatedProducts.at(0);
}


} // namespace lar_pandora

#endif // LAR_PANDORA_PFPARTICLE_UTILS_H

