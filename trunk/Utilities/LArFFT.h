////////////////////////////////////////////////////////////////////////
/// \file LArFFT.h
///
/// Utility FFT functions
///
/// \version $Id: SingleGen.h,v 1.2 2010/02/15 19:10:40 brebel Exp $
/// \author  Brian Page
////////////////////////////////////////////////////////////////////////
#ifndef LARFFT_H
#define LARFFT_H

#include "TComplex.h"
#include "TFFTRealComplex.h"
#include "TFFTComplexReal.h"
#include "TF1.h"
#include "TH1D.h"
#include <vector>
#include <string>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

///General LArSoft Utilities
namespace util{
    class LArFFT {
    public:
      LArFFT(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~LArFFT();
      
      void         DoFFT(std::vector<double> & input,
			 std::vector<TComplex> & output);            
      
      void         DoInvFFT(std::vector<TComplex> & input,
                            std::vector<double> & output);
      
      void        Deconvolute(std::vector<double> & input,
			      std::vector<double> & respFunc);
      
      void        Deconvolute(std::vector<double> & input,
			      std::vector<TComplex> & kern);
      
      void        Convolute(std::vector<double> & input,
			    std::vector<double> & respFunc);
      
      void        Convolute(std::vector<double> & input,
			    std::vector<TComplex> & kern);
      
      void        Correlate(std::vector<double> & input,
			    std::vector<double> & respFunc);
      
      void        Correlate(std::vector<double> & input,
			    std::vector<TComplex> & kern);
      
      void        AlignedSum(std::vector<double> & input,
			     std::vector<double> &output,
			     bool add = true);
      
      void        ShiftData(std::vector<TComplex> & input,
			    double shift);
  
      void        ShiftData(std::vector<double> & input,
                            double shift);
      
      double      PeakCorrelation(std::vector<double> &shape1,
				  std::vector<double> &shape2);       
      
      const int   FFTSize() const { return fSize; }
	
	private:
      
      int                    fSize;       //size of transform
      int                    fFreqSize;   //size of frequency space
      std::string            fOption;     //FFTW setting
      int                    fFitBins;    //Bins used for peak fit
      TF1                   *fPeakFit;    //Gaussian peak function
      TH1D                  *fConvHist;   //Fit data histogram
      std::vector<TComplex>  fCompTemp;   //temporary complex data
      std::vector<TComplex>  fKern;       //transformed response function
       
      TFFTRealComplex       *fFFT;        ///< object to do FFT
      TFFTComplexReal       *fInverseFFT; ///< object to do Inverse FF   
      
    }; // class LArFFT
} //namespace utils
#endif // LARFFT_H
