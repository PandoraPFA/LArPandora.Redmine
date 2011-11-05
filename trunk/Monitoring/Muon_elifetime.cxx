////////////////////////////////////////////////////////////////////////
//
/// \file Muon_elifetime.cxx
/// \version $Id: NeutrinoAna.h,v 1.3 2010/06/17 12:06:38 antonm Exp $
///
/// \author joshua.spitz@yale.edu
//
//  This class is designed to find electron lifetime from long tracks.
//  The electron lifetime is found using the ADC values and drift times 
//  of hits that have been associated with HoughClusters. 
////////////////////////////////////////////////////////////////////////

#include "Monitoring/Muon_elifetime.h"

// Framework includes
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft inclues
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "RawData/raw.h"

// ROOT includes
#include <TH1D.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TTree.h>
#include "TDatabasePDG.h"
#include "TSystem.h"// 
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"

// C++ language includes
#include <iostream>
#include <sstream>
#include <vector>

namespace mon{

  //......................................................................
  Muon_elifetime::Muon_elifetime(fhicl::ParameterSet const& pset) 
  {
    this->reconfigure(pset);
  }

  //......................................................................
  Muon_elifetime::~Muon_elifetime() 
  {
    delete [] m_mipZ;
    delete [] m_drifttimeZ;
    delete [] m_upadcZ;
  }

  void Muon_elifetime::reconfigure(fhicl::ParameterSet const& p)
  {
    m_run              = 0;
    m_event            = 0;
    m_sizeHitZ         = 10000;
    fSlices            = p.get< int         >("Slices");
    fWirespan          = p.get< int         >("Wirespan");
    fPlane             = p.get< int         >("Plane");
    fChisqrndf         = p.get< float       >("Chisqrndf");
    fADCdef            = p.get< int         >("ADCdef");
    fRawDigitLabel     = p.get< std::string >("RawDigitLabel", "daq");
    fHitModuleLabel    = p.get< std::string >("HitModuleLabel");
    fClusterModuleLabel= p.get< std::string >("ClusterModuleLabel");
      
    return;
  }

  //......................................................................
  void Muon_elifetime::beginJob()
  {

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    fTree        = tfs->make<TTree>("HoughTree","HoughTree");
    m_mipZ       = new Float_t[m_sizeHitZ];
    m_drifttimeZ = new Float_t[m_sizeHitZ];
    m_upadcZ     = new Float_t[m_sizeHitZ];
    fTree->Branch("run", &m_run, "run/I");
    fTree->Branch("run_timestamp", &m_run_timestamp, "run_timestamp/D");
    fTree->Branch("event", &m_event, "event/I");
    fTree->Branch("numberHits",&m_sizeHitZ,"numberHits/I");
    fTree->Branch("mipZ",m_mipZ,"mipZ[numberHits]/F");
    fTree->Branch("drifttimeZ",m_drifttimeZ,"drifttitmeZ[numberHits]/F");
    fTree->Branch("upadcZ",m_upadcZ,"upadcZ[numberHits]/F");

  }

  //......................................................................
  //this job fills the tree (on an event-by-event basis) with the information we need
  void Muon_elifetime::analyze(const art::Event& evt)
  {

    art::ServiceHandle<geo::Geometry> geo;
  
    std::cout << "Event  : " << evt.id() << std::endl;
    m_run          =evt.id().run();
    m_subrun       =0;//evt.luminosityBlock();
    m_event        =evt.id().event();
    m_run_timestamp=evt.time().value();

    // Pull the raw digits out of the event
    art::Handle< std::vector<raw::RawDigit> > digitHandle;
    evt.getByLabel(fRawDigitLabel, digitHandle);
    std::vector< art::Ptr<raw::RawDigit> > rawdigit;
    art::fill_ptr_vector(rawdigit, digitHandle);

    // Pull the hits out of the event
    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fHitModuleLabel,hitHandle);
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitHandle);

    // Pull the clusters out of the event
    art::Handle< std::vector<recob::Cluster> > clusterHandle;
    evt.getByLabel(fClusterModuleLabel,clusterHandle);
    std::vector< art::Ptr<recob::Cluster> > clusters;
    art::fill_ptr_vector(clusters, clusterHandle);
     
  
    int firstwire=0;
    int lastwire=0;
    m_sizeHitZ=0;
    for(unsigned int j=0; j<clusters.size();++j) {      
      art::PtrVector<recob::Hit> _hits=clusters[j]->Hits(fPlane,-1);
      if(_hits.size()!=0){

	//***********have to fix the association from wire->rawdigit
// 	geo->ChannelToWire(_hits[0]->Wire()->RawDigit()->Channel(), p, firstwire);
// 	geo->ChannelToWire(_hits[_hits.size()-1]->Wire()->RawDigit()->Channel(), p, lastwire);
	if((lastwire-firstwire)<fWirespan)
	  continue;       
	  
	m_sizeHitZ=_hits.size();
	for(unsigned int i = 0; i < _hits.size(); ++i) {  
	  m_mipZ[i]=_hits[i]->Charge();
	  m_drifttimeZ[i]= _hits[i]->PeakTime();
	  m_upadcZ[i]=_hits[i]->Charge();
	} 
	fTree->Fill(); 
      } 

    }

    return;
  }


  //......................................................................
  void Muon_elifetime::endRun(art::Run const&){
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1,0);   
    gStyle->SetOptFit(kTRUE);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(111);
    int rangeX,nbins;
    if(fADCdef==0){
      nbins=1000;
      rangeX=100;
    }
    else{
      nbins=10000;
      rangeX=1000;
    }
  
    TH2F* hist = new TH2F("", "",2048,0,2048,nbins,0,rangeX);
    TH1F* hist2 = new TH1F("b", "b",200,0,100);
    TH1F* hist3 = new TH1F("c", "c",fSlices,0,2048);
    TH1F* hist4 = new TH1F("d", "d",fSlices,0,2048); 
    TFile *outfile=new TFile("e_lifetime.root","RECREATE");
    TCanvas *canvas = new TCanvas("plane","slice");
    int           event;
    int           numberHits=10000;
    TBranch        *b_event;  
    TBranch        *b_numberHits;
    TBranch        *b_mipZ;
    TBranch        *b_upadcZ;
    TBranch        *b_drifttimeZ;
    float* mipZ;
    float* upadcZ;
    mipZ = new float[numberHits];
    upadcZ = new float[numberHits];
    float* drifttimeZ;
    drifttimeZ = new float[numberHits]; 
    fTree->SetBranchAddress("event", &event, &b_event);
    fTree->SetBranchAddress("numberHits", &numberHits, &b_numberHits);
    fTree->SetBranchAddress("mipZ", mipZ, &b_mipZ);
    fTree->SetBranchAddress("upadcZ", upadcZ, &b_upadcZ);
    fTree->SetBranchAddress("drifttimeZ", drifttimeZ, &b_drifttimeZ); 
    unsigned int nevt=(int)fTree->GetEntries();        
    for(unsigned int event=0;event<nevt;event++){
      fTree->GetEntry(event);    
      for(int j=0;j<numberHits;j++){
	if(fADCdef==0)
	  hist->Fill(drifttimeZ[j],mipZ[j]); 
   	else  	
	  hist->Fill(drifttimeZ[j],upadcZ[j]);    	
      }    		   
    }
    
    for(int n=0;n<fSlices;n++){
      hist2->Reset();
      for(int i=n*(2048/fSlices);i<(n+1)*(2048/fSlices);i++){
	for(int j=0;j<hist->GetNbinsY();j++){
	  if(hist->GetBinContent(i,j)!=0)
	    hist2->Fill(hist->GetBinCenter(j)/10.,hist->GetBinContent(i,j)); 
	}
      }

      if(hist->GetNbinsY()==0)
	continue;
    
      // Variables for langaufit call:
      //   his             histogram to fit
      //   fitrange[2]     lo and hi boundaries of fit range
      //   startvalues[4]  reasonable start values for the fit
      //   parlimitslo[4]  lower parameter limits
      //   parlimitshi[4]  upper parameter limits
      //   fitparams[4]    returns the final fit parameters
      //   fiterrors[4]    returns the final fit errors
      //   ChiSqr          returns the chi square
      //   NDF             returns ndf
    
      double fr[2];
      double sv[4], pllo[4], plhi[4], fp[4], fpe[4];
      fr[0]=0.3*hist2->GetMean();
      fr[1]=3.0*hist2->GetMean();
      pllo[0]=0.2; pllo[1]=0.0; pllo[2]=1.0; pllo[3]=0.4;
      plhi[0]=5.0; plhi[1]=50.0; plhi[2]=5000.0; plhi[3]=5.0;
      sv[0]=1.0; sv[1]=hist2->GetMean(); sv[2]=500.0; sv[3]=1.5;
    
      double chisqr;
      int    ndf; 
      TF1 *g0 = langaufit(hist2,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
    
      //double SNRPeak, SNRFWHM;
      //langaupro(fp,SNRPeak,SNRFWHM);
    
      hist2->Draw();
      g0->DrawCopy("lsame");
    
      if(g0->GetNDF()>0 && g0->GetParameter(1)>1&&g0->GetChisquare()/g0->GetNDF()<fChisqrndf){
	hist3->Fill(n*(2048/fSlices),g0->GetParameter(1));
	hist3->SetBinError(n+1,g0->GetParError(1));
	hist4->Fill(n*(2048/fSlices),g0->GetChisquare()/g0->GetNDF());
      }
      canvas->Update();
      canvas->Write();
      g0->Clear();
    }

    TCanvas *canvas2 = new TCanvas("lifetime","lifetime");
  
    TF1* g00 = new TF1("g00","expo",80,1600); //probably starts at around 60. choose 80 to be safe
    hist3->Fit(g00,"EMR");
    float frac=-.198/g00->GetParameter(1);
    float unc=(g00->GetParError(1)/g00->GetParameter(1))*.198/g00->GetParameter(1);
    hist3->SetTitle("fbR650.root (100 wire span required)");
    hist3->SetXTitle("Time samples (198ns/sample)");
    hist3->SetYTitle("~ADC counts");
    hist3->Draw();
    g00->DrawCopy("lsame");
    char expression[64];
    sprintf(expression,"e-lifetime=%.0f +/- %.0f microsec", frac,unc);
    TLatex* latex = new TLatex(.2, .8, expression);
    latex->SetNDC();
    latex->SetTextSize(.042);
    latex->SetTextColor(1);
    latex->Draw();
    canvas2->Update();
    canvas2->Write();
  
    TCanvas *canvas3 = new TCanvas("chi2","chi2");
    hist4->SetStats(kFALSE);
    TLine *line = new TLine(0,1,2048,1);
    hist4->SetMarkerStyle(7);
    hist4->SetTitle("#chi^{2}/ndf");
    hist4->SetXTitle("Time samples (198ns/sample)");
    hist4->SetYTitle("#chi^{2}/ndf");
    hist4->DrawCopy("P");
    line->Draw("SAME");
    canvas3->Write();
  
    outfile->Close();

    ///delete all the new'ed pointers
    delete hist;
    delete hist2;
    delete hist3;
    delete hist4;
    delete outfile;
    delete canvas;
    delete [] mipZ;
    delete [] upadcZ;
    delete drifttimeZ;
    delete canvas2;
    delete g00;
    delete latex;
    delete canvas3;
    delete line;

    return;
  }

  ////-----------------------------------------------------------------------
  //
  //	Convoluted Landau and Gaussian Fitting Function
  //         (using ROOT's Landau and Gauss functions)
  //
  //  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
  //  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
  //   Markus Friedl (Markus.Friedl@cern.ch)

  double langaufun(double *x, double *par) {

    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation), 
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.
  
    // Numeric constants
    double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    double mpshift  = -0.22278298;       // Landau maximum location
  
    // Control constants
    double np = 100.0;      // number of convolution steps
    double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
    // Variables
    double xx;
    double mpc;
    double fland;
    double sum = 0.0;
    double xlow,xupp;
    double step;
    double i;
  
  
    // MP shift correction
    mpc = par[1] - mpshift * par[0]; 
  
    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
  
    step = (xupp-xlow) / np;
  
    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
      xx = xupp - (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }
  
    return (par[2] * step * sum * invsq2pi / par[3]);
  }

  TF1 *Muon_elifetime::langaufit(TH1F *his, 
				 double *fitrange, 
				 double *startvalues, 
				 double *parlimitslo, 
				 double *parlimitshi, 
				 double *fitparams, 
				 double *fiterrors, 
				 double *ChiSqr, 
				 int *NDF)
  {
    // Once again, here are the Landau * Gaussian parameters:
    //   par[0]=Width (scale) parameter of Landau density
    //   par[1]=Most Probable (MP, location) parameter of Landau density
    //   par[2]=Total area (integral -inf to inf, normalization constant)
    //   par[3]=Width (sigma) of convoluted Gaussian function
    //
    // Variables for langaufit call:
    //   his             histogram to fit
    //   fitrange[2]     lo and hi boundaries of fit range
    //   startvalues[4]  reasonable start values for the fit
    //   parlimitslo[4]  lower parameter limits
    //   parlimitshi[4]  upper parameter limits
    //   fitparams[4]    returns the final fit parameters
    //   fiterrors[4]    returns the final fit errors
    //   ChiSqr          returns the chi square
    //   NDF             returns ndf

    int i;
    Char_t FunName[100];

    sprintf(FunName,"Fitfcn_%s",his->GetName());

    TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
    if (ffitold) delete ffitold;

    TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
    ffit->SetParameters(startvalues);
    ffit->SetParNames("Width","MP","Area","GSigma");
   
    for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
    }

    his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot

    ffit->GetParameters(fitparams);    // obtain fit parameters
    for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
    }
    ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
    NDF[0] = ffit->GetNDF();           // obtain ndf

    return (ffit);              // return fit function

  }


  int langaupro(double *params, double &maxx, double &FWHM) {

    // Seaches for the location (x value) at the maximum of the 
    // Landau-Gaussian convolute and its full width at half-maximum.
    //
    // The search is probably not very efficient, but it's a first try.
  
    double p,x,fy,fxr,fxl;
    double step;
    double l,lold;
    int i = 0;
    int MAXCALLS = 10000;
  
  
    // Search for maximum
  
    p = params[1] - 0.1 * params[0];
    step = 0.05 * params[0];
    lold = -2.0;
    l    = -1.0;
  
  
    while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
    
      lold = l;
      x = p + step;
      l = langaufun(&x,params);
    
      if (l < lold)
	step = -step/10;
    
      p += step;
    }
  
    if (i == MAXCALLS)
      return (-1);
  
    maxx = x;
  
    fy = l/2;
  
  
    // Search for right x location of fy
  
    p = maxx + params[0];
    step = params[0];
    lold = -2.0;
    l    = -1e300;
    i    = 0;
  
  
    while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
    
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
    
      if (l > lold)
	step = -step/10;
    
      p += step;
    }
  
    if (i == MAXCALLS)
      return (-2);
  
    fxr = x;
  
  
    // Search for left x location of fy
  
    p = maxx - 0.5 * params[0];
    step = -params[0];
    lold = -2.0;
    l    = -1e300;
    i    = 0;
  
    while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
    
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
    
      if (l > lold)
	step = -step/10;
    
      p += step;
    }
  
    if (i == MAXCALLS)
      return (-3);
  
  
    fxl = x;
  
    FWHM = fxr - fxl;
    return (0);
  }
   
   
} //end namespace
