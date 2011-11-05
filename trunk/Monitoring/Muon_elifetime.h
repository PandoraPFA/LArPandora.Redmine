////////////////////////////////////////////////////////////////////////
/// \file Muon_elifetime.h
///
/// Module to produce a set list of particles for a MC event
///
/// \version $Id: SingleGen.h,v 1.2 2010/02/15 19:10:40 brebel Exp $
/// \author  Josh Spitz
////////////////////////////////////////////////////////////////////////
#ifndef MUON_ELIFETIME_H
#define MUON_ELIFETIME_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"

// ROOT includes.
#include "TH1.h"
#include "TF1.h"

#include <string>

class TTree;

namespace mon{
  class Muon_elifetime : public art::EDAnalyzer {

  public:
    
    explicit Muon_elifetime(fhicl::ParameterSet const& pset);
    virtual ~Muon_elifetime();

    void reconfigure(fhicl::ParameterSet const& p);
    void analyze (const art::Event& evt);
    void beginJob();
    void endRun  (art::Run const&);
    TF1 *langaufit(TH1F *his, 
		   double *fitrange, 
		   double *startvalues, 
		   double *parlimitslo, 
		   double *parlimitshi, 
		   double *fitparams, 
		   double *fiterrors, 
		   double *ChiSqr, 
		   int    *NDF);
    
  private:
    
    TTree*      fTree;
    int         m_run;          
    int         m_subrun;                          
    double      m_run_timestamp;          
    int         m_event;             ///< Event number     		      
    float*      m_mipZ;              ///< MIP deposition of hit (peak height)
    float*      m_drifttimeZ;        ///< drift time of hit		      
    float*      m_upadcZ;            ///< MIP deposition of hit (integrated) 
    int         m_sizeHitZ;          ///< Number of Hits                     
    int         fSlices;
    int         fWirespan;
    int         fPlane;
    float       fChisqrndf;
    int         fADCdef;
    std::string fRawDigitLabel;      ///< module that made the raw digits
    std::string fHitModuleLabel;     ///< module that made the hits you want
    std::string fClusterModuleLabel; ///< module that made the clusters you want
      
};

}

#endif //MUON_ELIFETIME_H
