/**
 *
 * @file dune/AnaUtils/DUNEAnaUtilsBase.h
 *
 * @brief Base class containing functionality to extract products from the event
*/

#ifndef DUNE_ANA_UTILS_BASE_H
#define DUNE_ANA_UTILS_BASE_H

#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <string>
#include <vector>

namespace dune_ana
{
/**
 *
 * @brief DUNEAnaUtilsBase class containing some template functions
 *
*/
class DUNEAnaUtilsBase
{
protected:
    template <typename T> static std::vector<art::Ptr<T>> GetProductVector(const art::Event &evt, const std::string &label);
    template <typename T, typename U> static std::vector<art::Ptr<T>> GetAssocProductVector(const art::Ptr<U> &part, const art::Event &evt, const std::string &label, const std::string &assocLabel);
    template <typename T, typename U> static art::Ptr<T> GetAssocProduct(const art::Ptr<U> &part, const art::Event &evt, const std::string &label, const std::string &assocLabel); 
};

// Implementation of the template function to get the products from the event
template <typename T> std::vector<art::Ptr<T>> DUNEAnaUtilsBase::GetProductVector(const art::Event &evt, const std::string &label)
{
    auto theseProds = evt.getHandle<std::vector<T>>(label);
    bool success = theseProds.isValid();

    if (!success)
    {
        mf::LogError("DUNEAna") << " Failed to find product with label " << label << " ... returning empty vector" << std::endl;
        return std::vector<art::Ptr<T>>();
    }

    // We need to convert these to art pointers
    std::vector<art::Ptr<T>> productVector;
    art::fill_ptr_vector(productVector,theseProds);

    return productVector;
}

// Implementation of the template function to get the associated products from the event
template <typename T, typename U> std::vector<art::Ptr<T>> DUNEAnaUtilsBase::GetAssocProductVector(const art::Ptr<U> &pProd, const art::Event &evt, const std::string &label, const std::string &assocLabel)
{
    auto products = evt.getHandle<std::vector<U>>(label);
    bool success = products.isValid();

    if (!success)
    {
        mf::LogError("DUNEAna") << " Failed to find product with label " << label << " ... returning empty vector" << std::endl;
        return std::vector<art::Ptr<T>>();
    }

    const art::FindManyP<T> findParticleAssocs(products,evt,assocLabel);

    return findParticleAssocs.at(pProd.key());
}

// Implementation of the template function to get the associated product from the event
template <typename T, typename U> art::Ptr<T> DUNEAnaUtilsBase::GetAssocProduct(const art::Ptr<U> &pProd, const art::Event &evt, const std::string &label, const std::string &assocLabel)
{   
    std::vector<art::Ptr<T>> associatedProducts = DUNEAnaUtilsBase::GetAssocProductVector<T>(pProd,evt,label,assocLabel); 
    if (associatedProducts.empty())
    {
        throw cet::exception("DUNEAna") << "DUNEAnaUtilsBase::GetAssocProduct --- No associated object found";
    }
    return associatedProducts.at(0);
}


} // namespace dune_ana

#endif // DUNE_ANA_PFPARTICLE_UTILS_H

