////////////////////////////////////////////////////////////////////////
/// 
/// \brief Module to tabulate Beam stats
///
/// \author  msoderbe@syr.edu
////////////////////////////////////////////////////////////////////////

#ifndef BEAMANA_H
#define BEAMANA_H

#include <iostream>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

#include "Utilities/LArProperties.h"
#include "T962/T962_Objects/MINOS.h"
#include "RawData/BeamInfo.h"

#include <TH1F.h>



///T962 beam analysis code
namespace beam {
  
  class BeamAna :  public art::EDAnalyzer {
    
  public:
    
    explicit BeamAna(fhicl::ParameterSet const& pset); 
    virtual ~BeamAna();        

    void analyze (const art::Event& evt);
    void beginJob();
    void endJob();
        
  private:
       
    //plots for BeamAna
    TH1F*  fTOR101_Run;
    TH1F*  fTORTGT_Run;
    TH1F*  fTRTGTD_Run;
    TH1F*  fTOR101_Time;
    TH1F*  fTORTGT_Time;
    TH1F*  fTRTGTD_Time;

    TH1F*  fTOR101_WithMINOS_Run;
    TH1F*  fTORTGT_WithMINOS_Run;
    TH1F*  fTRTGTD_WithMINOS_Run;
    TH1F*  fTOR101_WithMINOS_Time;
    TH1F*  fTORTGT_WithMINOS_Time;
    TH1F*  fTRTGTD_WithMINOS_Time;
     
    TH1F*  fLifetime_Run;
    TH1F*  fLifetime_Time;
   
        
  };
    

//-----------------------------------------------------------------------------
  BeamAna::BeamAna(fhicl::ParameterSet const& pset)
  {
  }
    
//-----------------------------------------------------------------------------
  BeamAna::~BeamAna()
  {
  }

//-----------------------------------------------------------------------------
  void BeamAna::beginJob()
  {
    // get access to the TFile service  
    art::ServiceHandle<art::TFileService> tfs;
    
    fTOR101_Run = tfs->make<TH1F>("fTOR101_Run","TOR101 vs Run", 300,600.0,900.0);
    fTORTGT_Run = tfs->make<TH1F>("fTORTGT_Run","TORTGT vs Run", 300,600.0,900.0);
    fTRTGTD_Run = tfs->make<TH1F>("fTRTGTD_Run","TRTGTD vs Run", 300,600.0,900.0);

    fTOR101_Time = tfs->make<TH1F>("fTOR101_Time","TOR101 vs Time", 10000,1253090000.0,1266900000.0);
    fTORTGT_Time = tfs->make<TH1F>("fTORTGT_Time","TORTGT vs Time", 10000,1253090000.0,1266900000.0);
    fTRTGTD_Time = tfs->make<TH1F>("fTRTGTD_Time","TRTGTD vs Time", 10000,1253090000.0,1266900000.0);

    fTOR101_WithMINOS_Run = tfs->make<TH1F>("fTOR101_WithMINOS_Run","TOR101 vs Run, with MINOS", 300,600.0,900.0);
    fTORTGT_WithMINOS_Run = tfs->make<TH1F>("fTORTGT_WithMINOS_Run","TORTGT vs Run, with MINOS", 300,600.0,900.0);
    fTRTGTD_WithMINOS_Run = tfs->make<TH1F>("fTRTGTD_WithMINOS_Run","TRTGTD vs Run, with MINOS", 300,600.0,900.0);

    fTOR101_WithMINOS_Time = tfs->make<TH1F>("fTOR101_WithMINOS_Time","TOR101 vs Time, with MINOS", 10000,1253090000.0,1266900000.0);
    fTORTGT_WithMINOS_Time = tfs->make<TH1F>("fTORTGT_WithMINOS_Time","TORTGT vs Time, with MINOS", 10000,1253090000.0,1266900000.0);
    fTRTGTD_WithMINOS_Time = tfs->make<TH1F>("fTRTGTD_WithMINOS_Time","TRTGTD vs Time, with MINOS", 10000,1253090000.0,1266900000.0);

    fLifetime_Run = tfs->make<TH1F>("fLifetime_Run","Lifetime vs Run", 300,600.0,900.0);
    fLifetime_Time = tfs->make<TH1F>("fLifetime_Time","Lifetime vs Time", 10000,1253090000.0,1266900000.0);
        
  }

//-----------------------------------------------------------------------------
  void BeamAna::endJob()
  {
    mf::LogInfo ("Summary") << "All Events: "
                            << "\nfTOR101_Run->Integral() = " << fTOR101_Run->Integral()
                            << "\nfTORTGT_Run->Integral() = " << fTORTGT_Run->Integral() 
                            << "\nfTRTGTD_Run->Integral() = " << fTRTGTD_Run->Integral()
                            << "\n\nWith MINOS requirement: " 
                            << "\nfTOR101_WithMINOS_Run->Integral() = " << fTOR101_WithMINOS_Run->Integral() 
                            << "\nfTORTGT_WithMINOS_Run->Integral() = " << fTORTGT_WithMINOS_Run->Integral() 
                            << "\nfTRTGTD_WithMINOS_Run->Integral() = " << fTRTGTD_WithMINOS_Run->Integral();
  }
  
//-----------------------------------------------------------------------------
  void BeamAna::analyze(const art::Event& evt) 
  {
    
    art::ServiceHandle<util::LArProperties> LArProp;
    double tau = LArProp->ElectronLifetime();
     
    art::Handle< raw::BeamInfo> beam;
    double beamtime;
    if(evt.getByLabel("beam",beam)){
      beamtime = (double)beam->get_t_ms();
      beamtime/=1000.0;
      double tor101 = beam->get_tor101();
      double tortgt = beam->get_tortgt();
      double trtgtd = beam->get_trtgtd();

      fTOR101_Run->Fill(evt.run(),tor101);
      fTORTGT_Run->Fill(evt.run(),tortgt);
      fTRTGTD_Run->Fill(evt.run(),trtgtd);
      fTOR101_Time->Fill(beamtime,tor101);
      fTORTGT_Time->Fill(beamtime,tortgt);
      fTRTGTD_Time->Fill(beamtime,trtgtd);

      art::Handle< std::vector<t962::MINOS>> minos;
      if(evt.getByLabel("minos",minos)){
        fTOR101_WithMINOS_Run->Fill(evt.run(),tor101);
        fTORTGT_WithMINOS_Run->Fill(evt.run(),tortgt);
        fTRTGTD_WithMINOS_Run->Fill(evt.run(),trtgtd);
        fTOR101_WithMINOS_Time->Fill(beamtime,tor101);
        fTORTGT_WithMINOS_Time->Fill(beamtime,tortgt);
        fTRTGTD_WithMINOS_Time->Fill(beamtime,trtgtd);
      }
      if(evt.event()==1){
        fLifetime_Run->Fill(evt.run(),tau);
        fLifetime_Time->Fill(beamtime,tau);
      }
    }

    return;
   
  }//analyze

//Required
  DEFINE_ART_MODULE(BeamAna);

} // end of namespace

#endif // BEAMANA_H
