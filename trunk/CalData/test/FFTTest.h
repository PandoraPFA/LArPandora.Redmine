#ifndef FFTTEST_H
#define FFTTEST_H

//
// Name:  FFTTest.h
//
// Purpose: FFTTest module.  Test convolution/deconvolution.
//
// Created:  29-Aug-2011  H. Greenlee

#include "art/Framework/Core/EDAnalyzer.h"
#include <vector>
#include "TComplex.h"

namespace caldata
{
  class FFTTest : public art::EDAnalyzer
  {
  public:

    // Constructor, destructor.

    explicit FFTTest(fhicl::ParameterSet const& pset);
    virtual ~FFTTest();

    // Overrides.

    void analyze(const art::Event& evt);

  private:

    // Attributes.

    std::string fSimFile;   // SimWire response file name.
    std::string fCalFile;   // CalWire response file name.
    int fNTicks;                // Number of ticks.

    // Time domain response functions.

    std::vector<double> fSimElect;    // response function for the electronics
    std::vector<double> fSimColField; // response function for the field @ collection plane
    std::vector<double> fSimIndField; // response function for the field @ induction plane
    std::vector<double> fSimColConv;  // Collection plane convoluted response function.
    std::vector<double> fSimIndConv;  // Induction plane convoluted response function.

    // Frequency domain response functions.

    std::vector<TComplex> fSimElectF;    // response function for the electronics
    std::vector<TComplex> fSimColFieldF; // response function for the field @ collection plane
    std::vector<TComplex> fSimIndFieldF; // response function for the field @ induction plane
    std::vector<TComplex> fSimColConvF;  // Collection plane convoluted response function.
    std::vector<TComplex> fSimIndConvF;  // Induction plane convoluted response function.
    std::vector<TComplex> fColDeconvF;   // Collection plane deconvolution.
    std::vector<TComplex> fIndDeconvF;   // Collection plane deconvolution.
  };
}

#endif
