////////////////////////////////////////////////////////////////////////
// $Id: LArFFT.cxx,v 1.11 2010/07/02 20:33:09 bpage Exp $
//
// \file LArFFT class
//
//  This class simplifies implementation of Fourier transforms. 
//  Because all data inputs and outputs are purely real,  the 
//  transforms implemented in this way get a substantial performance
//  increase ~2x.
//
// \author pagebri3@msu.edu
//
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include "Utilities/LArFFT.h"

//-----------------------------------------------
util::LArFFT::LArFFT(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg) :
  fSize    (pset.get< int         >("FFTSize")),
  fOption  (pset.get< std::string >("FFTOption")),
  fFitBins (pset.get< int         >("FitBins"))
{
  fFreqSize = fSize/2+1;  

  //allocate and setup Transform objects  
  fFFT        = new TFFTRealComplex(fSize, false);  
  fInverseFFT = new TFFTComplexReal(fSize, false);  

  int dummy[1] = {0}; 
  ///< appears to be dummy argument from root page  
  fFFT->Init(fOption.c_str(),-1,dummy);  
  fInverseFFT->Init(fOption.c_str(),1,dummy);  

  fPeakFit = new TF1("fPeakFit","gaus"); //allocate function used for peak fitting  
  fConvHist = new TH1D("fConvHist","Convolution Peak Data",fFitBins,0,fFitBins);  //allocate histogram for peak fitting
  //allocate other data vectors  
  fCompTemp.resize(fFreqSize);  
  fKern.resize(fFreqSize);}

//------------------------------------------------
util::LArFFT::~LArFFT() 
{
  delete fFFT;
  delete fInverseFFT;
  delete fPeakFit;
}

//-------------------------------------------------
//For the sake of efficiency, as all transforms should
//be of the same size, all functions expect vectors 
//to be of the desired size.
//Note that because the transforms are real to real
//there is a redundancy in the information in the
//complex part in the positive and negative
//components of the FFT, thus the size of the
//frequency vectors are input_fFreqSize 
//--see the FFTW3 or Root docmentation for details

////"Forward" Fourier Transform
//--------------------------------------------------------
void util::LArFFT::DoFFT(std::vector<double> & input,		   
			  std::vector<TComplex> & output)
{  
  double real = 0.;    //real value holder  
  double imaginary = 0.;  //imaginary value hold   
  fFFT->SetPoints(&input[0]);    
  
  fFFT->Transform();    
  
  for(int i = 0; i < fFreqSize; ++i){    
    fFFT->GetPointComplex(i, real, imaginary);    
    output[i]=TComplex(real, imaginary);  
  }  
  
  return;
}

//Inverse Fourier Transform
//-------------------------------------------------
void util::LArFFT::DoInvFFT(std::vector<TComplex> & input,		      
			     std::vector<double> & output)
{  
  for(int i = 0; i < fFreqSize; ++i)    
    fInverseFFT->SetPointComplex(i, input[i]);   

  fInverseFFT->Transform();  
  double factor = 1.0/(double) fSize;  
  
  for(int i = 0; i < fSize; ++i)    
    output[i]=factor*fInverseFFT->GetPointReal(i,false);  

  return;
}

//Deconvolution scheme taking all time-domain
//information
//--------------------------------------------------
void util::LArFFT::Deconvolute(std::vector<double> & input,			 
				std::vector<double> & respFunction)
{  
  DoFFT(respFunction, fKern);  
  DoFFT(input, fCompTemp);  

  for(int i = 0; i < fFreqSize; i++) 
    fCompTemp[i]/=fKern[i];  

  DoInvFFT(fCompTemp, input);  

  return;
}

//Deconvolution scheme using an already transformed
//response function
//saves cpu time if same response function is used
//for many consecutive transforms
//--------------------------------------------------
void util::LArFFT::Deconvolute(std::vector<double> & input,			 
				std::vector<TComplex> & kern)
{    
  DoFFT(input, fCompTemp);  

  for(int i = 0; i < fFreqSize; i++) 
    fCompTemp[i]/=kern[i];  

  DoInvFFT(fCompTemp, input);  
  
  return;
}

//Convolution scheme taking all time-domain
//information
//--------------------------------------------------
void util::LArFFT::Convolute(std::vector<double> & shape1,
			      std::vector<double> & shape2)
{  
  DoFFT(shape1, fKern);  
  DoFFT(shape2, fCompTemp);  

  for(int i = 0; i < fFreqSize; i++) 
    fCompTemp[i]*=fKern[i];  

  DoInvFFT(fCompTemp, shape1);  

  return;
}

//Convolution scheme using an already transformed
//response function
//saves cpu time if same response function is used
//for many consecutive transforms
//--------------------------------------------------
void util::LArFFT::Convolute(std::vector<double> & input,
			      std::vector<TComplex> & kern)
{  
  DoFFT(input, fCompTemp);  
  
  for(int i = 0; i < fFreqSize; i++) 
    fCompTemp[i]*=kern[i];  
  
  DoInvFFT(fCompTemp, input);    
  
  return;
}

//Correlation taking all time domain data
//--------------------------------------------------
void util::LArFFT::Correlate(std::vector<double> & shape1,
			      std::vector<double> & shape2)
{  
  DoFFT(shape1, fKern);  
  DoFFT(shape2, fCompTemp);  
  
  for(int i = 0; i < fFreqSize; i++)     
    fCompTemp[i]*=TComplex::Conjugate(fKern[i]);  
  
  DoInvFFT(fCompTemp, shape1);  

  return;
}

//Convolution scheme using an already transformed
//response function
//saves cpu time if same response function is used
//for many consecutive transforms
//--------------------------------------------------
void util::LArFFT::Correlate(std::vector<double> & input,	
			      std::vector<TComplex> & kern)
{  
  DoFFT(input, fCompTemp);  

  for(int i = 0; i < fFreqSize; i++)     
    fCompTemp[i]*=TComplex::Conjugate(kern[i]);  

  DoInvFFT(fCompTemp, input);    
  
  return;
}

//Scheme for adding two signals which have an arbitrary
//relative translation.  Shape1 is translated over shape2
//and is replaced with the sum, or the translated result
//if add = false
//--------------------------------------------------
void util::LArFFT::AlignedSum(std::vector<double> & shape1,
			       std::vector<double> & shape2,     
			       bool add)
{  
  double shift = PeakCorrelation(shape1,shape2);    
  
  ShiftData(shape1,shift);  
   
  if(add)for(int i = 0; i < fSize; i++) shape1[i]+=shape2[i];   

  return;
}

//According to the Fourier transform identity
//f(x-a) = Inverse Transform(exp(-2*Pi*i*a*w)F(w))
//--------------------------------------------------
void util::LArFFT::ShiftData(std::vector<TComplex> & input,	
			      double shift)
{ 
  double factor = -2.0*TMath::Pi()*shift/(double)fSize;  

  for(int i = 0; i < fFreqSize; i++)    
    input[i] *= TComplex::Exp(TComplex(0,factor*(double)i));  

  return;
}

//Shifts real vectors using above function
//--------------------------------------------------
void util::LArFFT::ShiftData(std::vector<double> & input,	
			      double shift)
{ 
  DoFFT(input,fCompTemp);
  ShiftData(fCompTemp,shift);
  DoInvFFT(fCompTemp,input); 

  return;
}

//Returns the length of the translation at which the correlation
//of 2 signals is maximal.
//--------------------------------------------------
double util::LArFFT::PeakCorrelation(std::vector<double> & shape1,	
				      std::vector<double> & shape2)
{ 
  fConvHist->Reset("ICE");
  std::vector<double> holder = shape1;  
  Correlate(holder,shape2);  
  
  int maxT = max_element(holder.begin(), holder.end())-holder.begin();  
  double startT = maxT-fFitBins/2;  
  int offset=0;
  
  for(int i = 0; i < fFitBins; i++) { 
    if(startT+i < 0) offset=fSize;
    else if(startT+i > fSize) offset=-fSize;
    else offset = 0;     
    fConvHist->Fill(i,holder[i+startT+offset]);    
  }  
  
  fPeakFit->SetParameters(fConvHist->GetMaximum(),fFitBins/2,fFitBins/2);  
  fConvHist->Fit(fPeakFit,"QWNR","",0,fFitBins);  
  return fPeakFit->GetParameter(1)+startT;
}
