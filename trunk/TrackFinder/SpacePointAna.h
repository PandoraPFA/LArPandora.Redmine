#ifndef SPACEPOINTANA_H
#define SPACEPOINTANA_H

//
// Name: SpacePointAna.h
//
// Purpose: Header file for module SpacePointAna.
//
// Created: 2-Aug-2011  H. Greenlee
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "TH1F.h"

namespace trkf {

  class SpacePointAna : public art::EDAnalyzer
  {
  public:
 
    // Constructors, destructor

    explicit SpacePointAna(fhicl::ParameterSet const& pset);
    virtual ~SpacePointAna();

    // Book histograms.

    void bookHistograms(bool mc);

    // Overrides.

    void analyze(const art::Event& evt);

  private:

    // Fcl Attributes.

    std::string fHitModuleLabel;
    std::string fClusterModuleLabel;
    std::string fG4ModuleLabel;    // For SimChannel.
    bool fUseClusterHits;
    double fMaxDT;
    double fMaxS;

    // Histograms.

    bool fBooked;   // Have histograms been booked yet?
    TH1F* fHDTUE;   // U-drift electrons time difference.
    TH1F* fHDTVE;   // V-drift electrons time difference.
    TH1F* fHDTWE;   // W-drift electrons time difference.
    TH1F* fHDTUV;   // U-V time difference.
    TH1F* fHDTVW;   // V-W time difference.
    TH1F* fHDTWU;   // W-U time difference.
    TH1F* fHS;      // Spatial separation.
    TH1F* fHx;      // X position.
    TH1F* fHy;      // Y position.
    TH1F* fHz;      // Z position.
    TH1F* fHMCdx;   // X residual (reco vs. mc truth).
    TH1F* fHMCdy;   // Y residual (reco vs. mc truth).
    TH1F* fHMCdz;   // Z residual (reco vs. mc truth).

    // Statistics.

    int fNumEvent;
  };
}

#endif // SPACEPOINTANA_H
