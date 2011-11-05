////////////////////////////////////////////////////////////////////////
//
// ScanFilter class:
// Tells the downstream rreconstruction to not process events that
// do not pass certain hand scan criteria
//
// joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef SCANFILTER_H
#define SCANFILTER_H

#include "art/Framework/Core/EDFilter.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TH2D.h"
#include "TTree.h"
class TTree;

///T962 filter for hand scan results
namespace filt {

  class ScanFilter : public art::EDFilter  {
    
  public:
    
    explicit ScanFilter(fhicl::ParameterSet const& ); 
    virtual ~ScanFilter();
    void beginJob();
    bool filter(art::Event& evt);
   
  private: 
    TTree* ftree;
    std::string fScanModuleLabel;
    int fNeutrino_req, fMinShowers_req, fMinTracks_req, fMaxShowers_req, fMaxTracks_req;
    float fFidVolume_cut;
    float fm_vertexx;         
    float fm_vertexy;
    float fm_vertexz;
    int fm_run,fm_event,fm_tracks,fm_showers,fm_neutrino,fm_maybeneutrino;
    int fCreateTTree;
  protected: 
    
  }; // class ScanFilter

}

#endif // SCANFILTER_H
