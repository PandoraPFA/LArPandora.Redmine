/// \file    ColorDrawingOptions.cxx
/// \brief   The color scales used by the event display
/// \author  messier@indiana.edu
/// \version $Id: ColorDrawingOptions.cxx,v 1.1.1.1 2010/11/10 19:44:54 p-novaart Exp $
#include "EventDisplay/ColorDrawingOptions.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"

namespace evd{
  
  //......................................................................
  ColorDrawingOptions::ColorDrawingOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
    : fColorOrGray(pset.get< int    >("ColorOrGrayScale"))
    , fRawDiv     (pset.get< int    >("RawDiv"))
    , fRecoDiv    (pset.get< int    >("RecoDiv"))
    , fRawQLow    (pset.get< double >("RawQLow"))	
    , fRawQHigh   (pset.get< double >("RawQHigh"))	
    , fRecoQLow   (pset.get< double >("RecoQLow"))	
    , fRecoQHigh  (pset.get< double >("RecoQHigh"))     
    , fColorScaleRaw(evdb::ColorScale(fRawQLow, fRawQHigh, 
				      evdb::kBlueToRedII, evdb::kLinear,
				      fRawDiv,
				      285.0, 135.0, //angle in the color wheel
				      0.65, 0.25))  //intensity from light to dark, starting with low color wheel value
    , fColorScaleReco(evdb::ColorScale(fRecoQLow, fRecoQHigh, 
				       evdb::kBlueToRedII, evdb::kLinear,
				       fRecoDiv,
				       285.0, 135.0,
				       0.65, 0.25)) 
    , fGrayScaleRaw(evdb::ColorScale(fRawQLow, fRawQHigh, 
				     evdb::kLinGray, evdb::kLinear,
				     fRawDiv,
				     270.0, 0.0, //angle in the color wheel
				     0.5, 0.5))  //intensity from light to dark, starting with low color wheel value
    , fGrayScaleReco(evdb::ColorScale(fRecoQLow, fRecoQHigh, 
				      evdb::kLinGray, evdb::kLinear,
				      fRecoDiv,
				      270.0, 0.0, //angle in the color wheel
				      0.5, 0.5))  //intensity from light to dark, starting with low color wheel value
    
  {

  }

  //......................................................................
  ColorDrawingOptions::~ColorDrawingOptions()
  {
  }

  //......................................................................
  void ColorDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
    fColorOrGray = pset.get< int    >("ColorOrGrayScale");
    fRawDiv      = pset.get< int    >("RawDiv");
    fRecoDiv     = pset.get< int    >("RecoDiv");
    fRawQLow     = pset.get< double >("RawQLow");	
    fRawQHigh    = pset.get< double >("RawQHigh");	
    fRecoQLow    = pset.get< double >("RecoQLow");	
    fRecoQHigh   = pset.get< double >("RecoQHigh");     

    fColorScaleRaw .SetBounds(fRawQLow,  fRawQHigh);
    fColorScaleReco.SetBounds(fRecoQLow, fRecoQHigh);
    fGrayScaleRaw  .SetBounds(fRawQLow,  fRawQHigh);
    fGrayScaleReco .SetBounds(fRecoQLow, fRecoQHigh);

  }

  //......................................................................
  const evdb::ColorScale& ColorDrawingOptions::RawQ() const
  {
    if(fColorOrGray > 0) return fGrayScaleRaw;
    
    return fColorScaleRaw;
  }

  //......................................................................
  const evdb::ColorScale& ColorDrawingOptions::CalQ() const
  {
    if(fColorOrGray > 0) return fGrayScaleReco;
    
    return fColorScaleReco;
  }

  //......................................................................
  const evdb::ColorScale& ColorDrawingOptions::RawT() const
  {
    if(fColorOrGray > 0) return fGrayScaleRaw;
    
    return fColorScaleRaw;
  }

  //......................................................................
  const evdb::ColorScale& ColorDrawingOptions::CalT() const
  {
    if(fColorOrGray > 0) return fGrayScaleReco;
    
    return fColorScaleReco;
  }
}// namespace
////////////////////////////////////////////////////////////////////////
