////////////////////////////////////////////////////////////////////////
//
/// \file Muon_elifetime_module.cc
/// \version $Id: NeutrinoAna.h,v 1.3 2010/06/17 12:06:38 antonm Exp $
///
/// \author joshua.spitz@yale.edu
//
//  This class is designed to find electron lifetime from long tracks.
//  The electron lifetime is found using the ADC values and drift times 
//  of hits that have been associated with HoughClusters. 
////////////////////////////////////////////////////////////////////////
#include "Monitoring/Muon_elifetime.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace mon{

  DEFINE_ART_MODULE(Muon_elifetime);

} //end namespace
