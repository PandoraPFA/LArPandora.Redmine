////////////////////////////////////////////////////////////////////////
// $Id: RecoDrawingOption.h,v 1.15 2010/08/30 21:33:24 spitz7 Exp $
//
// Display parameters for the raw data
//
// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef RECODRAWINGOPTIONS_H
#define RECODRAWINGOPTIONS_H
#ifndef __CINT__
#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace evd {
  class RecoDrawingOptions 
  {
  public:
    RecoDrawingOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~RecoDrawingOptions();
    
    void reconfigure(fhicl::ParameterSet const& pset);

    int fDrawHits;
    int fDrawClusters;
    int fDrawProngs;
    int fDrawTracks;
    int fDrawShowers;
    int fDrawVertices;
    int fDrawEvents;
    int fDraw2DEndPoints;

    std::vector<std::string> fWireLabels;    ///< module labels that produced wires
    std::vector<std::string> fHitLabels;     ///< module labels that produced hits
    std::vector<std::string> fClusterLabels; ///< module labels that produced clusters
    std::vector<std::string> fProngLabels;   ///< module labels that produced prongs
    std::vector<std::string> fTrackLabels;   ///< module labels that produced tracks
    std::vector<std::string> fShowerLabels;  ///< module labels that produced showers
    std::vector<std::string> fVertexLabels;  ///< module labels that produced vertices
    std::vector<std::string> fEventLabels;   ///< module labels that produced events

  };
}//namespace
#endif // __CINT__
#endif

