/**
 * @file    RootGraphicsEnablingService.h
 * @brief   Service linked to RootGraphicsEnabler library
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    March 3, 2016
 */


#ifndef LARPANDORA_LARPANDORA_ROOTGRAPHICS_ROOTGRAPHICSENABLINGSERVICE_H
#define LARPANDORA_LARPANDORA_ROOTGRAPHICS_ROOTGRAPHICSENABLINGSERVICE_H

// framework libraries
#include "art/Framework/Services/Registry/ServiceMacros.h"


namespace fhicl { class ParameterSet; }
namespace art { class ActivityRegistry; }


namespace util {
  
  /**
   * @brief Trojan service to inject initialization code
   *
   * Run this service to ensure that a graphical (interactive) ROOT session
   * can be initialized (e.g., for TEve-based visualization.
   * Just add it to the services in your FHiCL configuration file with e.g.,
   *     
   *     services.RootGraphicsEnablingService: {}
   *     
   *
   * This service does very much nothing.
   * But it is linked to the RootGraphicsEnabler library that does magics on loading.
   * 
   * Configuration parameters
   * =========================
   *
   * No configuration.
   *
   */
  class RootGraphicsEnablingService {
      public:
    RootGraphicsEnablingService
      (fhicl::ParameterSet const&, art::ActivityRegistry&) {}
  }; // class RootGraphicsEnablingService
  
} // namespace util


DECLARE_ART_SERVICE(util::RootGraphicsEnablingService, LEGACY)


#endif // LARPANDORA_LARPANDORA_ROOTGRAPHICS_ROOTGRAPHICSENABLINGSERVICE_H

