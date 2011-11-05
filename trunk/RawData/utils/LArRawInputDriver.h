////////////////////////////////////////////////////////////////////////
/// \file  LArRawInputDriver.h
/// \brief Source to convert raw binary files to root files
///
/// \version $Id: T962ConvertBinaryToROOT.h,v 1.7 2010/01/14 19:20:33 brebel Exp $
/// \author  brebel@fnal.gov, soderber@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/Core/PrincipalMaker.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Persistency/Provenance/SubRunID.h"

#include <fstream>
#include <string>
#include <vector>

///Conversion of binary data to root files
namespace lris {
  class LArRawInputDriver;
}

class lris::LArRawInputDriver {
  /// Class to fill the constraints on a template argument to the class,
  /// FileReaderSource
 public:
  // Required constructor
  LArRawInputDriver(fhicl::ParameterSet const &pset,
                    art::ProductRegistryHelper &helper,
                    art::PrincipalMaker const &pm);

  // Required by FileReaderSource:
  void closeCurrentFile();
  void readFile(std::string const &name,
                art::FileBlock* &fb);
  bool readNext(art::RunPrincipal* const &inR,
                art::SubRunPrincipal* const &inSR,
                art::RunPrincipal* &outR,
                art::SubRunPrincipal* &outSR,
                art::EventPrincipal* &outE);

 private:
  // --- data members:
  typedef  std::vector<std::string>  stringvec_t;

  art::PrincipalMaker            principalMaker_;
  std::string                    currentDir_;
  stringvec_t                    inputfiles_;
  stringvec_t::const_iterator    nextfile_;
  stringvec_t::const_iterator    filesdone_;
  art::SubRunID                  currentSubRunID_; 
};  // LArRawInputDriver
