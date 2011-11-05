////////////////////////////////////////////////////////////////////////
/// \file  AngleCluster_module.cc
/// \brief Create a Cluster with an angle 
///
/// \version $Id: AngleCluster.cxx,v 0.1 19/07/2011 12:45:16 PM  andrzejs $
/// \author andrzej.szelc@yale.edu
/// \author based on code from brossi and tstrauss and in the future from ShowerFinder by Roxanne Guenette  
////////////////////////////////////////////////////////////////////////
// This class solves the following problem:
//
// Create a cluster saving its slope angle

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

#include "ClusterFinder/ShowerAngleCluster.h"

namespace cluster {

  DEFINE_ART_MODULE(ShowerAngleCluster);

}
