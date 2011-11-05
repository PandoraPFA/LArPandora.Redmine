//
// Name: FFTTest.css
//
// Purpose: Implementation file for module FFTTest.
//
// Created: 29-Aug-2011  H. Greenlee
//

#include <iostream>
#include <cmath>
#include "CalData/test/FFTTest.h"
#include "Utilities/LArFFT.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/search_path.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

// Local functions.

namespace {

  // Fill vector from histogram (ignore underflow/overflow bins).

  void hist_to_vector(const TH1D* h, std::vector<double>& v)
  {
    assert(h != 0);
    int n = h->GetNbinsX();
    v.resize(n);
    for(int i=0; i<n; ++i)
      v[i] = h->GetBinContent(i+1);
  }

  // Fill histogram from vector (set underflow/overflow bins to zero).

  void vector_to_hist(const std::vector<double>& v, TH1D* h)
  {
    assert(h != 0);
    int nvec = v.size();
    int nbins = h->GetNbinsX();
    int nfill = std::min(nvec, nbins);
    h->SetBinContent(0, 0.);
    for(int i=0; i<nfill; ++i)
      h->SetBinContent(i+1, v[i]);
    for(int i=nfill+1; i<=nbins+1; ++i)
      h->SetBinContent(i, 0.);
  }

  // Fill vector with initial delta-function at bin d.

  void fill_delta(std::vector<double>& v, int d)
  {
    int n = v.size();
    assert(d >= 0 && d < n);
    for(int i=0; i<n; ++i)
      v[i] = 0.;
    v[d] = 1.;
  }
}

namespace caldata {

  FFTTest::FFTTest(const fhicl::ParameterSet& pset)
  {
    // Get file service.

    art::ServiceHandle<art::TFileService> tfs;

    // Get FFT service.

    art::ServiceHandle<util::LArFFT> fFFT;
    fNTicks = fFFT->FFTSize();
    std::cout << "Number of ticks = " << fNTicks << std::endl;

    // Get simulation (convolution) response functions.

    fSimFile = pset.get<std::string>("simwire_file");
    std::cout << "SimWire file = " << fSimFile << std::endl;

    TFile fsim(fSimFile.c_str());

    TH1D* hSimElect = dynamic_cast<TH1D*>(fsim.Get("daq/ElectronicsResponse"));
    hist_to_vector(hSimElect, fSimElect);
    fSimElect.resize(fNTicks, 0.);
    fSimElectF.resize(fNTicks/2+1);
    fFFT->DoFFT(fSimElect, fSimElectF);

    TH1D* hSimColField = dynamic_cast<TH1D*>(fsim.Get("daq/CollectionFieldResponse"));
    hist_to_vector(hSimColField, fSimColField);
    fSimColField.resize(fNTicks, 0.);
    fSimColFieldF.resize(fNTicks/2+1);
    fFFT->DoFFT(fSimColField, fSimColFieldF);

    TH1D* hSimIndField = dynamic_cast<TH1D*>(fsim.Get("daq/InductionFieldResponse"));
    hist_to_vector(hSimIndField, fSimIndField);
    fSimIndField.resize(fNTicks, 0.);
    fSimIndFieldF.resize(fNTicks/2+1);
    fFFT->DoFFT(fSimIndField, fSimIndFieldF);

    TH1D* hSimColConv = dynamic_cast<TH1D*>(fsim.Get("daq/ConvolutedCollection"));
    hist_to_vector(hSimColConv, fSimColConv);
    fSimColConv.resize(fNTicks, 0.);
    fSimColConvF.resize(fNTicks/2+1);
    fFFT->DoFFT(fSimColConv, fSimColConvF);

    TH1D* hSimIndConv = dynamic_cast<TH1D*>(fsim.Get("daq/ConvolutedInduction"));
    hist_to_vector(hSimIndConv, fSimIndConv);
    fSimIndConv.resize(fNTicks, 0.);
    fSimIndConvF.resize(fNTicks/2+1);
    fFFT->DoFFT(fSimIndConv, fSimIndConvF);

    // Get reco (deconvolution) response function.

    fhicl::ParameterSet calwire_pset = pset.get<fhicl::ParameterSet>("calwire");
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(calwire_pset.get<std::string>("ResponseFile"), fCalFile);
    std::cout << "CalWire file = " << fCalFile << std::endl;

    TFile fcal(fCalFile.c_str());

    TH2D* respRe = dynamic_cast<TH2D*>(fcal.Get("sim/RespRe"));
    TH2D* respIm = dynamic_cast<TH2D*>(fcal.Get("sim/RespIm"));
    int nx = respRe->GetNbinsX();
    int ny = respRe->GetNbinsY();
    assert(nx == respIm->GetNbinsX());
    assert(ny == respIm->GetNbinsY());
    assert(nx == 2);   // 1=induction, 2=collection.

    fColDeconvF.resize(ny);
    fIndDeconvF.resize(ny);

    for(int i=0; i<ny; ++i) {
      double ac = respRe->GetBinContent(2, i+1);
      double bc = respIm->GetBinContent(2, i+1);
      TComplex zc(ac, bc);
      fColDeconvF[i] = zc;

      double ai = respRe->GetBinContent(1, i+1);
      double bi = respIm->GetBinContent(1, i+1);
      TComplex zi(ai, bi);
      fIndDeconvF[i] = zi;
    }

    // Calculate response of delta function to collection field + electronics.

    art::TFileDirectory dirc = tfs->mkdir("Collection", "Collection");
    int nhist = std::min(200, fNTicks);

    std::vector<double> twork(fNTicks, 0.);
    fill_delta(twork, nhist/2);
    TH1D* hinc = dirc.make<TH1D>("input", "Collection Input", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(twork, hinc);
    hinc->Write();

    fFFT->Convolute(twork, fSimElectF);
    TH1D* helectc = dirc.make<TH1D>("elect", "Collection Electronics", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(twork, helectc);
    helectc->Write();

    fill_delta(twork, nhist/2);
    fFFT->Convolute(twork, fSimColFieldF);
    TH1D* hfieldc = dirc.make<TH1D>("field", "Collection Field", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(twork, hfieldc);
    hfieldc->Write();

    fFFT->Convolute(twork, fSimElectF);
    TH1D* hbothc = dirc.make<TH1D>("both", "Collection Field+Electronics", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(twork, hbothc);
    hbothc->Write();

    // Now do the same thing using the convoluted response function.

    fill_delta(twork, nhist/2);
    fFFT->Convolute(twork, fSimColConvF);
    TH1D* hconvc = dirc.make<TH1D>("conv", "Collection Convoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(twork, hconvc);
    hconvc->Write();

    // Deconvolution.

    fFFT->Convolute(twork, fColDeconvF);
    TH1D* hdeconvc = dirc.make<TH1D>("deconv", "Collection Deconvoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(twork, hdeconvc);
    hdeconvc->Write();

    // Calculate response of delta function to induction field + electronics.

    art::TFileDirectory diri = tfs->mkdir("Induction", "Induction");

    fill_delta(twork, nhist/2);
    TH1D* hini = diri.make<TH1D>("input", "Induction Input", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(twork, hini);
    hini->Write();

    fFFT->Convolute(twork, fSimElectF);
    TH1D* helecti = diri.make<TH1D>("elect", "Induction Electronics", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(twork, helecti);
    helecti->Write();

    fill_delta(twork, nhist/2);
    fFFT->Convolute(twork, fSimIndFieldF);
    TH1D* hfieldi = diri.make<TH1D>("field", "Induction Field", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(twork, hfieldi);
    hfieldi->Write();

    fFFT->Convolute(twork, fSimElectF);
    TH1D* hbothi = diri.make<TH1D>("both", "Induction Field+Electronics", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(twork, hbothi);
    hbothi->Write();

    // Now do the same thing using the convoluted response function.

    fill_delta(twork, nhist/2);
    fFFT->Convolute(twork, fSimIndConvF);
    TH1D* hconvi = diri.make<TH1D>("conv", "Induction Convoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(twork, hconvi);
    hconvi->Write();

    // Deconvolution.

    fFFT->Convolute(twork, fIndDeconvF);
    TH1D* hdeconvi = diri.make<TH1D>("deconv", "Induction Deconvoluted", nhist+1, -0.5, nhist+0.5);
    vector_to_hist(twork, hdeconvi);
    hdeconvi->Write();
  }

  FFTTest::~FFTTest()
  {}

  void FFTTest::analyze(const art::Event& evt)
 {}
}
