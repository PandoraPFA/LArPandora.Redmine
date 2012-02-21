#define datakin_cxx
#include "datakin.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "cutvariables.h"
void datakin::Loop()
{

std::cout<<denom<<std::endl;

 gROOT->SetStyle("Plain");
 gStyle->SetEndErrorSize(3);

   gStyle->SetPalette(1,0); 
   gStyle->SetLineWidth(2); 
   gStyle->SetHistLineWidth(3);
//    gStyle->SetOptFit(kTRUE);
   gStyle->SetOptStat(kFALSE);
//    gStyle->SetOptStat(kTRUE);
    gStyle->SetOptStat(1001001100);
    gStyle->SetOptFit(0011);
    gStyle->SetLabelSize(.032,"Z");
//    gStyle->SetOptFit(111);
  set_plot_style();
   TH1::AddDirectory(false);
   
    TH2D* histunfold = new TH2D("", "",100,0,1.57,100,0,1.57);
    TH1D* histeff1 = new TH1D("", "",18,0.,36);
    histeff1->Sumw2();
    TH1D* histeff111 = new TH1D("", "",18,0.,36);
    histeff111->Sumw2();
    TH1D* histeff2 = new TH1D("", "",18,0.,36);
    histeff2->Sumw2();
    TH1D* histeff3 = new TH1D("", "",18,0.,36);
    histeff3->Sumw2();
    TH1D* histeff11 = new TH1D("", "",20,0.,25);
    histeff11->Sumw2();
    TH1D* hhisteff11 = new TH1D("", "",20,0.,25);
    hhisteff11->Sumw2();
    TH1D* histeff1111 = new TH1D("", "",20,0.,25);
    histeff1111->Sumw2();
    TH1D* histeff22 = new TH1D("", "",20,0.,25);
    histeff22->Sumw2();
    TH1D* histeff33 = new TH1D("", "",20,0.,25);
    histeff33->Sumw2();
    
    TH1D* mom1 = new TH1D("", "",20,0.,25);
    TH1D* mom2 = new TH1D("", "",20,0.,25);
    TH1D* mom3 = new TH1D("", "",20,0.,25);
    TH1D* cos1 = new TH1D("", "",18,0.,36);
    TH1D* cos2 = new TH1D("", "",18,0.,36);
    TH1D* cos3 = new TH1D("", "",18,0.,36);
    
    TH1D* hhisteff111 = new TH1D("", "",18,0.,36);
    hhisteff111->Sumw2();
    TH1D* hhisteff1111 = new TH1D("", "",20,0.,25);
    hhisteff1111->Sumw2();
    
    TH1D* histerror0 = new TH1D("", "",20,0.,25);
    TH1D* histerror1 = new TH1D("", "",20,0.,25);
    TH1D* histerror2 = new TH1D("", "",20,0.,25);
    TH1D* histerror3 = new TH1D("", "",20,0.,25);
    TH1D* histerror4 = new TH1D("", "",20,0.,25);
    TH1D* histerror5 = new TH1D("", "",20,0.,25);
    TH1D* histerror6 = new TH1D("", "",20,0.,25);
    TH1D* histerror7 = new TH1D("", "",20,0.,25);
    TH1D* histerror0a = new TH1D("", "",20,0.,25);
    TH1D* histerror1a = new TH1D("", "",20,0.,25);
    TH1D* histerror2a = new TH1D("", "",20,0.,25);
    TH1D* histerror3a = new TH1D("", "",20,0.,25);
    TH1D* histerror4a = new TH1D("", "",20,0.,25);
    TH1D* histerror5a = new TH1D("", "",20,0.,25);
    TH1D* histerror6a = new TH1D("", "",20,0.,25);
    TH1D* histerror7a = new TH1D("", "",20,0.,25);
    TH1D* histerror00 = new TH1D("", "",18,0.,36);
    TH1D* histerror11 = new TH1D("", "",18,0.,36);
    TH1D* histerror22 = new TH1D("", "",18,0.,36);
    TH1D* histerror33 = new TH1D("", "",18,0.,36);
    TH1D* histerror44 = new TH1D("", "",18,0.,36);
    TH1D* histerror55 = new TH1D("", "",18,0.,36);
    TH1D* histerror77 = new TH1D("", "",18,0.,36);
    TH1D* histerror00a = new TH1D("", "",18,0.,36);
    TH1D* histerror11a = new TH1D("", "",18,0.,36);
    TH1D* histerror22a = new TH1D("", "",18,0.,36);
    TH1D* histerror33a = new TH1D("", "",18,0.,36);
    TH1D* histerror44a = new TH1D("", "",18,0.,36);
    TH1D* histerror55a = new TH1D("", "",18,0.,36);
    TH1D* histerror77a = new TH1D("", "",18,0.,36);
    
    TH1D* histbkgd = new TH1D("", "",18,0.,36);
     histbkgd->Sumw2();
    TH1D* histbkgd2 = new TH1D("", "",20,0.,25);
     histbkgd2->Sumw2();
     
     TH1D* histbkgd_par = new TH1D("", "",18,0.,36);
     histbkgd_par->Sumw2();
    TH1D* histbkgd2_par = new TH1D("", "",20,0.,25);
     histbkgd2_par->Sumw2();
     
     
     TH1D* histtgmuon = new TH1D("", "",18,0.,36);
     histtgmuon->Sumw2();
    TH1D* histtgmuon2 = new TH1D("", "",20,0.,25);
     histtgmuon2->Sumw2();
     
         
     TH1D* histnutgmuon = new TH1D("", "",18,0.,36);
     histnutgmuon->Sumw2();
    TH1D* histnutgmuon2 = new TH1D("", "",20,0.,25);
     histnutgmuon2->Sumw2(); 
     
       TH1D* histnutgmuon3 = new TH1D("", "",18,0.,36);
     histnutgmuon3->Sumw2();
    TH1D* histnutgmuon4 = new TH1D("", "",20,0.,25);
     histnutgmuon4->Sumw2(); 
     
     
    TH1D* hist1effcorrected = new TH1D("", "",18,0.,36);
    hist1effcorrected->Sumw2();
    TH1D* hist1mueeffcorrected = new TH1D("", "",20,0.,25);
    hist1mueeffcorrected->Sumw2();
    
    TH2D* histmomsyst = new TH2D("", "",20,0.,25,100,0,2);
    TH2D* histthetasyst = new TH2D("", "",18,0.,36,100,0,2);
    TH2D* histmomsyst2 = new TH2D("", "",20,0.,25,100,0,2);
    TH2D* histthetasyst2 = new TH2D("", "",18,0.,36,100,0,2);
    
    
    TH1D* hist2mueeffcorrected = new TH1D("", "",20,0.,25);
    TH1D* hist2effcorrected = new TH1D("", "",18,0.,36);
    
    
    TH2D* hist1_scan= new TH2D("", "",20,0,10,10,0,10);
    TH1D* hist2_scan= new TH1D("", "",3,0,3);
    TH1D* hist3_scan= new TH1D("", "",43,608,651);
    TH2D* hist4_scan= new TH2D("", "",50,0,50,50,0,50);
    TH2D* hist5_scan= new TH2D("", "",40,-20,20,40,-20,20);
    TH2D* hist6_scan= new TH2D("", "",90,0,90,90,0,90);
//    

    TH1F* vtxhist = new TH1F("", "",20,0,50);
    TH1F* vtxhist2 = new TH1F("", "",20,-25,25);
    TH1F* vtxhist3 = new TH1F("", "",45,0,90);
    TH2D* vvtxhist = new TH2D("", "",100,0,50,100,0,50);
    TH2D* vvtxhist2 = new TH2D("", "",100,-25,25,100,-25,25);
    TH2D* vvtxhist3 = new TH2D("", "",180,0,90,180,0,90);
    
    TH2D* vvtxhistsim = new TH2D("", "",100,0,50,100,0,50);
    TH2D* vvtxhist2sim = new TH2D("", "",100,-25,25,100,-25,25);
    TH2D* vvtxhist3sim = new TH2D("", "",180,0,90,180,0,90);
    TH1F* vtxhistsim = new TH1F("", "",20,0,50);
    TH1F* vtxhist2sim = new TH1F("", "",20,-25,25);
    TH1F* vtxhist3sim = new TH1F("", "",45,0,90);
    TH2D* vtxhist4= new TH2D("", "",100,0,50,100,-25,25);
    TH2D* vtxhist5= new TH2D("", "",100,0,50,100,0,100);
    TH2D* vtxhist7= new TH2D("", "",100,0,50,100,-25,25);
    TH2D* vtxhist8= new TH2D("", "",100,0,50,100,0,100);
    TH2D* vtxhist4_sim= new TH2D("", "",100,0,50,100,-25,25);
    TH2D* vtxhist5_sim= new TH2D("", "",100,0,50,100,0,100);
    TH2D* vtxhist7_sim= new TH2D("", "",100,0,50,100,-25,25);
    TH2D* vtxhist8_sim= new TH2D("", "",100,0,50,100,0,100);
      TH1D* histpur1 = new TH1D("", "",18,0.,36);
    histpur1->Sumw2();
    TH1D* histpur2 = new TH1D("", "",18,0.,36);
    histpur2->Sumw2();
    
TH1D* hist1 = new TH1D("", "",18,0.,36);
TH1D* hist2 = new TH1D("", "",18,0.,36);
TH1D* hist3 = new TH1D("h3", "h3",20,0.,40);
TH1D* hist3up = new TH1D("h3up", "h3up",20,0.,40);
TH1D* hist3down = new TH1D("h3down", "h3down",20,0.,40);

TH1D* hist33 = new TH1D("", "",18,0.,36);
TH1D* hist33up = new TH1D("", "",18,0.,36);
TH1D* hist33down = new TH1D("", "",18,0.,36);

TH1D* hist333 = new TH1D("", "",18,0.,36);
TH1D* hist333up = new TH1D("", "",18,0.,36);
TH1D* hist333down = new TH1D("", "",18,0.,36);

TH1D* hist4 = new TH1D("", "",18,0.,36);
TH1D* hist5 = new TH1D("", "",18,0.,36);
TH1D* hist6 = new TH1D("", "",18,0.,36);

TH1D* hist1mue = new TH1D("", "",20,0.,25);
TH1D* hist2mue = new TH1D("", "",20,0.,25);
TH1D* hist3mue = new TH1D("h3m", "h3m",25,0.,31.25);
TH1D* hist3mueup = new TH1D("h3mup", "h3mup",25,0.,31.25);
TH1D* hist3muedown = new TH1D("h3mdown", "h3mdown",25,0.,31.25);

TH1D* hist33mue = new TH1D("", "",20,0.,25);
TH1D* hist33mueup = new TH1D("", "",20,0.,25);
TH1D* hist33muedown = new TH1D("", "",20,0.,25);

TH1D* hist333mue = new TH1D("", "",20,0.,25);
TH1D* hist333mueup = new TH1D("", "",20,0.,25);
TH1D* hist333muedown = new TH1D("", "",20,0.,25);

TH1D* hist8mue = new TH1D("", "",20,0.,25);
TH1D* hist8mueup = new TH1D("", "",20,0.,25);
TH1D* hist8muedown = new TH1D("", "",20,0.,25);

TH1D* hist4mue = new TH1D("", "",20,0.,25);
TH1D* hist5mue = new TH1D("", "",20,0.,25);
TH1D* hist6mue = new TH1D("", "",20,0.,25);

TH1D* histkalsim = new TH1D("", "",100,0.5,1);
TH1D* histkaldata = new TH1D("", "",100,0.5,1);

TH1D* hist11 = new TH1D("", "",18,0.,36);
TH1D* hist11effcorrected = new TH1D("cos", "cos",18,0.,36);
hist11effcorrected->Sumw2();
TH1D* hist11effcorrectedd = new TH1D("", "",18,0.,36);
hist11effcorrectedd->Sumw2();

TH1D* hist11mue = new TH1D("", "",20,0.,25); 




TH1D* hist11mueeffcorrected = new TH1D("mom", "mom",20,0.,25);
hist11mueeffcorrected->Sumw2();
TH1D* hist11mueeffcorrectedd = new TH1D("", "",20,0.,25);
hist11mueeffcorrectedd->Sumw2();


TH1D* hist22mueeffcorrected = new TH1D("", "",20,0.,25);
TH1D* hist22effcorrected = new TH1D("", "",18,0.,36);

TH1D* hist1_minos_sim = new TH1D("", "",60,-2,2);
TH1D* hist2_minos_sim = new TH1D("", "",60,-2,2);
TH1D* hist3_minos_sim = new TH1D("", "",60,-2,2);
TH1D* hist4_minos_sim = new TH1D("", "",25,0,50);
TH1D* hist5_minos_sim = new TH1D("", "",35,-35,35);
TH1D* hist6_minos_sim = new TH1D("", "",35,-35,35);
TH1D* hist7_minos_sim = new TH1D("", "",25,0,50);

TH1D* hist1_minos_data = new TH1D("", "",60,-2,2);
TH1D* hist2_minos_data = new TH1D("", "",60,-2,2);
TH1D* hist3_minos_data = new TH1D("", "",60,-2,2);
TH1D* hist4_minos_data = new TH1D("", "",25,0,50);
TH1D* hist5_minos_data = new TH1D("", "",35,-35,35);
TH1D* hist6_minos_data = new TH1D("", "",35,-35,35);
TH1D* hist7_minos_data = new TH1D("", "",25,0,50);
TH2D* hist8_minos_data = new TH2D("", "",500,-250,250,500,-250,250);
TH2D* hist8_minos_sim= new TH2D("", "",500,-250,250,500,-250,250);
TChain ch("KinTree");


// ch.Add("simkinreco_hist_golden1_0.root/analysistree/anatree");
// ch.Add("simkinreco_hist_golden2_0.root/analysistree/anatree");
// ch.Add("simkinreco_hist_golden3_0.root/analysistree/anatree");
// ch.Add("simkinreco_hist_golden4_0.root/analysistree/anatree");

//default
// ch.Add("simkinreco_hist_golden1_1.root/analysistree/anatree");
//  ch.Add("simkinreco_hist_golden2_1.root/analysistree/anatree");
//  ch.Add("simkinreco_hist_golden3_1.root/analysistree/anatree");
//  ch.Add("simkinreco_hist_golden4_1.root/analysistree/anatree");
 
 //w/ mitch spacepoints 
ch.Add("simkinreco_hist_golden1_2.root/analysistree/anatree");  
ch.Add("simkinreco_hist_golden2_2.root/analysistree/anatree"); 
ch.Add("simkinreco_hist_golden3_2.root/analysistree/anatree"); 
ch.Add("simkinreco_hist_golden4_2.root/analysistree/anatree"); 
//  
 
 datakin aEvent(KinTree);

  Long64_t nentries = KinTree.GetEntries();
  cout << "Number of Entries " << nentries << endl;
  float theta=0.;
  int failed=0;
  int numevents=0;
  int numccevents_beforecuts=0;
  int numccevents_aftercuts=0;
  //obvious cuts
  int numevents_passmatchandcharge=0;
  //not obvious cuts
  int numevents_passntracks=0;
  int numevents_passnclusters=0;
  double totpot=0.;
  double totpot2=0.;
  double totpot3=0.;
  
   TFile *fi=new TFile("profsyst.root","READ");
   TFile *momcosiso=new TFile("momcosiso.root","READ");
   TH1F *profthetasyst2 = (TH1F*)fi->Get("histthetasystfinal");
   TH1F *profmomsyst2 = (TH1F*)fi->Get("histmomsystfinal");
   
   TH1F *cosiso = (TH1F*)momcosiso->Get("isocos");
   TH1F *momiso = (TH1F*)momcosiso->Get("isomom");
   
//    double thetasyst[20]= {profthetasyst2->GetBinContent(1), profthetasyst2->GetBinContent(2), profthetasyst2->GetBinContent(3), profthetasyst2->GetBinContent(4), profthetasyst2->GetBinContent(5), profthetasyst2->GetBinContent(6), profthetasyst2->GetBinContent(7), profthetasyst2->GetBinContent(8), profthetasyst2->GetBinContent(9), profthetasyst2->GetBinContent(10), profthetasyst2->GetBinContent(11), profthetasyst2->GetBinContent(12), profthetasyst2->GetBinContent(13), profthetasyst2->GetBinContent(14),profthetasyst2->GetBinContent(15),profthetasyst2->GetBinContent(16),profthetasyst2->GetBinContent(17),profthetasyst2->GetBinContent(18),profthetasyst2->GetBinContent(19),profthetasyst2->GetBinContent(20)};
  
// std::cout<<profthetasyst2->GetBinContent(20)<<std::endl;

   double thetasyst[30]= {profthetasyst2->GetBinContent(1), profthetasyst2->GetBinContent(2), profthetasyst2->GetBinContent(3), profthetasyst2->GetBinContent(4), profthetasyst2->GetBinContent(5), profthetasyst2->GetBinContent(6), profthetasyst2->GetBinContent(7), profthetasyst2->GetBinContent(8), profthetasyst2->GetBinContent(9), profthetasyst2->GetBinContent(10), profthetasyst2->GetBinContent(11), profthetasyst2->GetBinContent(12), profthetasyst2->GetBinContent(13), profthetasyst2->GetBinContent(14),profthetasyst2->GetBinContent(15),profthetasyst2->GetBinContent(16),profthetasyst2->GetBinContent(17),profthetasyst2->GetBinContent(18),profthetasyst2->GetBinContent(19),profthetasyst2->GetBinContent(20),profthetasyst2->GetBinContent(21),profthetasyst2->GetBinContent(22),profthetasyst2->GetBinContent(23),profthetasyst2->GetBinContent(24),profthetasyst2->GetBinContent(25),profthetasyst2->GetBinContent(26),profthetasyst2->GetBinContent(27),profthetasyst2->GetBinContent(28),profthetasyst2->GetBinContent(29),profthetasyst2->GetBinContent(30)};
  
   double momsyst[30]= {profmomsyst2->GetBinContent(1), profmomsyst2->GetBinContent(2), profmomsyst2->GetBinContent(3), profmomsyst2->GetBinContent(4), profmomsyst2->GetBinContent(5), profmomsyst2->GetBinContent(6), profmomsyst2->GetBinContent(7), profmomsyst2->GetBinContent(8), profmomsyst2->GetBinContent(9), profmomsyst2->GetBinContent(10), profmomsyst2->GetBinContent(11), profmomsyst2->GetBinContent(12), profmomsyst2->GetBinContent(13), profmomsyst2->GetBinContent(14),profmomsyst2->GetBinContent(15),profmomsyst2->GetBinContent(16),profmomsyst2->GetBinContent(17),profmomsyst2->GetBinContent(18),profmomsyst2->GetBinContent(19),profmomsyst2->GetBinContent(20),profmomsyst2->GetBinContent(21),profmomsyst2->GetBinContent(22),profmomsyst2->GetBinContent(23),profmomsyst2->GetBinContent(24),profmomsyst2->GetBinContent(25),profmomsyst2->GetBinContent(26),profmomsyst2->GetBinContent(27),profmomsyst2->GetBinContent(28),profmomsyst2->GetBinContent(29),profmomsyst2->GetBinContent(30)};
  
  
  TFile *f=new TFile("numu_numode_final.root","READ");
   TH1F *reweight = (TH1F*)f->Get("histdiv");

   double weight[14]= {reweight->GetBinContent(1), reweight->GetBinContent(2), reweight->GetBinContent(3), reweight->GetBinContent(4), reweight->GetBinContent(5), reweight->GetBinContent(6), reweight->GetBinContent(7), reweight->GetBinContent(8), reweight->GetBinContent(9), reweight->GetBinContent(10), reweight->GetBinContent(11), reweight->GetBinContent(12), reweight->GetBinContent(13), reweight->GetBinContent(14)};
   
double rweight=0.;
   
for (int event = 0; event<nentries; event++) 
{

 
if(event%10000 == 0)
cout<<event<<" / "<<nentries<<endl;
aEvent.GetEntry(event);  



double dcosy=(aEvent.trackstart_dcosy_reco*TMath::Cos(.0583497))+(aEvent.trackstart_dcosz_reco*TMath::Sin(.0583497));
double dcosz=(aEvent.trackstart_dcosz_reco*TMath::Cos(.0583497))-(aEvent.trackstart_dcosy_reco*TMath::Sin(.0583497));

double dcosytruth=(aEvent.lep_dcosy_truth*TMath::Cos(.0583497))+(aEvent.lep_dcosz_truth*TMath::Sin(.0583497));
double dcosztruth=(aEvent.lep_dcosz_truth*TMath::Cos(.0583497))-(aEvent.lep_dcosy_truth*TMath::Sin(.0583497));



totpot+=aEvent.pot/100.;

if(aEvent.enu_truth>=0&&aEvent.enu_truth<3)
rweight=weight[0];
if(aEvent.enu_truth>=3&&aEvent.enu_truth<4)
rweight=weight[1];
if(aEvent.enu_truth>=4&&aEvent.enu_truth<5)
rweight=weight[2];
if(aEvent.enu_truth>=5&&aEvent.enu_truth<7)
rweight=weight[3];
if(aEvent.enu_truth>=7&&aEvent.enu_truth<9)
rweight=weight[4];
if(aEvent.enu_truth>=9&&aEvent.enu_truth<12)
rweight=weight[5];
if(aEvent.enu_truth>=12&&aEvent.enu_truth<15)
rweight=weight[6];
if(aEvent.enu_truth>=15&&aEvent.enu_truth<18)
rweight=weight[7];
if(aEvent.enu_truth>=18&&aEvent.enu_truth<22)
rweight=weight[8];
if(aEvent.enu_truth>=22&&aEvent.enu_truth<26)
rweight=weight[9];
if(aEvent.enu_truth>=26&&aEvent.enu_truth<30)
rweight=weight[10];
if(aEvent.enu_truth>=30&&aEvent.enu_truth<36)
rweight=weight[11];
if(aEvent.enu_truth>=36&&aEvent.enu_truth<42)
rweight=weight[12];
if(aEvent.enu_truth>=42&&aEvent.enu_truth<50)
rweight=weight[13];


if(
aEvent.nuvtxx_truth>0.+FIDX_CUT 
&& aEvent.nuvtxx_truth<47.-FIDX_CUT 
&& aEvent.nuvtxy_truth>-20.+FIDY_CUT 
&& aEvent.nuvtxy_truth<20.-FIDY_CUT 
&& aEvent.nuvtxz_truth>0.+FIDZup_CUT 
&& aEvent.nuvtxz_truth<90.-FIDZdown_CUT 
)
{
if(aEvent.ccnc_truth==0&&aEvent.nuPDG_truth==14)
{
histeff1->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
histeff11->Fill(aEvent.lep_mom_truth,rweight);
hhisteff11->Fill(aEvent.lep_mom_truth,rweight);
histeff111->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
histeff1111->Fill(aEvent.lep_mom_truth,rweight);
numccevents_beforecuts++;

 if(aEvent.test_charge_minos==-1)
 {
 hhisteff111->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
 hhisteff1111->Fill(aEvent.lep_mom_truth,rweight);
 }
}

}


if(
aEvent.vtxx_reco>0.+FIDX_CUT && aEvent.trackstart_x_reco>0.+FIDX_CUT 
&& aEvent.vtxx_reco<47.-FIDX_CUT && aEvent.trackstart_x_reco<47.-FIDX_CUT  
&& aEvent.vtxy_reco>-20.+FIDY_CUT && aEvent.trackstart_y_reco>-20.+FIDY_CUT 
&& aEvent.vtxy_reco<20.-FIDY_CUT  && aEvent.trackstart_y_reco<20.-FIDY_CUT
&& aEvent.vtxz_reco>0.+FIDZup_CUT && aEvent.trackstart_z_reco>0.+FIDZup_CUT 
&& aEvent.vtxz_reco<90.-FIDZdown_CUT && aEvent.trackstart_z_reco<90.-FIDZdown_CUT 
)
{



if(aEvent.nmatched_reco==1&&aEvent.trk_charge_minos<0&&aEvent.trk_mom_minos>0)
{



if(!aEvent.isdata)
{

double lardirectionEnd[3]={aEvent.trackexit_dcosx_reco, aEvent.trackexit_dcosy_reco, aEvent.trackexit_dcosz_reco};
double larEnd[3]= { aEvent.trackexit_x_reco, aEvent.trackexit_y_reco, aEvent.trackexit_z_reco};
double minosdirectionStart[3]={ aEvent.trk_dcosx_minos, aEvent.trk_dcosy_minos, aEvent.trk_dcosz_minos };

double minosStart[3]={aEvent.trk_vtxx_minos, aEvent.trk_vtxy_minos, aEvent.trk_vtxz_minos };
double xdiff,ydiff,rdiff,totaldiff,thetadiff;


project(lardirectionEnd, larEnd, minosdirectionStart,minosStart,xdiff,ydiff,rdiff,totaldiff,thetadiff);


double minosvtx[3] = {(aEvent.trk_vtxx_minos*100.)-117.4, (aEvent.trk_vtxy_minos*100.)+19.3, (aEvent.trk_vtxz_minos*100.)+147.1};

double distance=sqrt(pow(aEvent.vtxx_reco-minosvtx[0],2)+pow(aEvent.vtxy_reco-minosvtx[1],2)+pow(aEvent.vtxz_reco-minosvtx[2],2));
 
MOM_OFFSET=(0.002198*distance)-.04404; 

if(rdiff>RDIFF_CUT||thetadiff>THETADIFF_CUT)
continue;



histpur1->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);

if(aEvent.ccnc_truth==1)
{
hist2->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);
hist2mue->Fill(aEvent.trk_mom_minos+MOM_OFFSET,rweight);
}
if(aEvent.ccnc_truth==0&&aEvent.nuPDG_truth==14)
{
 histeff2->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
//histeff2->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);
histeff3->Fill((180./3.14159)*TMath::ACos(dcosztruth),rweight);
//histeff3->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);

hist3->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);
int bin=((180./3.14159)*TMath::ACos(dcosz)-0.)*(36./72.)-.0000001;
if(bin<30&&bin>=0)
{


if(((180./3.14159)*TMath::ACos(dcosz)+((180./3.14159)*TMath::ACos(dcosz)*thetasyst[bin]))>0.)
hist3up->Fill((180./3.14159)*TMath::ACos(dcosz)+((180./3.14159)*TMath::ACos(dcosz)*thetasyst[bin]),rweight);
else
{
hist3up->Fill(.000001,rweight);
}
 //std::cout<<aEvent.trackstart_dcosz_reco<<" "<<bin<<" "<<thetasyst[bin]<<" "<<aEvent.trackstart_dcosz_reco+(aEvent.trackstart_dcosz_reco*thetasyst[bin])<<std::endl;


hist3down->Fill((180./3.14159)*TMath::ACos(dcosz)-((180./3.14159)*TMath::ACos(dcosz)*thetasyst[bin]),rweight);



}

int bin2=((aEvent.trk_mom_minos+MOM_OFFSET)-0.)*(40./50.);
if(bin2<30&&bin2>=0)
{
hist3mueup->Fill((aEvent.trk_mom_minos+MOM_OFFSET)+((aEvent.trk_mom_minos+MOM_OFFSET)*momsyst[bin2]),rweight);
// std::cout<<aEvent.trk_mom_minos+MOM_OFFSET<<" "<<bin2<<" "<<momsyst[bin2]<<" "<<aEvent.trk_mom_minos+MOM_OFFSET+((aEvent.trk_mom_minos+MOM_OFFSET)*momsyst[bin2])<<std::endl;

if(((aEvent.trk_mom_minos+MOM_OFFSET)-((aEvent.trk_mom_minos+MOM_OFFSET)*momsyst[bin2]))>0.)
hist3muedown->Fill((aEvent.trk_mom_minos+MOM_OFFSET)-((aEvent.trk_mom_minos+MOM_OFFSET)*momsyst[bin2]),rweight);
else
hist3muedown->Fill(.000001,rweight);
}




histeff22->Fill(aEvent.lep_mom_truth,rweight);
histeff33->Fill(aEvent.lep_mom_truth,rweight);
hist3mue->Fill(aEvent.trk_mom_minos+MOM_OFFSET,rweight);

hist8mue->Fill(aEvent.trk_mom_minos+MOM_OFFSET,rweight);
hist8mueup->Fill(aEvent.trk_mom_minos+MOM_OFFSET+(MOM_OFFSET*MOM_OFFSET_UNC),rweight);

if(aEvent.trk_mom_minos+MOM_OFFSET-(MOM_OFFSET*MOM_OFFSET_UNC)>0)
hist8muedown->Fill(aEvent.trk_mom_minos+MOM_OFFSET-(MOM_OFFSET*MOM_OFFSET_UNC),rweight);
else
hist8muedown->Fill(.000001,rweight);

hist1->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);
hist1effcorrected->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);
hist2effcorrected->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);
hist1mue->Fill(aEvent.trk_mom_minos+MOM_OFFSET,rweight);
hist1mueeffcorrected->Fill(aEvent.trk_mom_minos+MOM_OFFSET,rweight);
hist2mueeffcorrected->Fill(aEvent.trk_mom_minos+MOM_OFFSET,rweight);
histmomsyst->Fill(aEvent.lep_mom_truth, fabs(aEvent.trk_mom_minos+MOM_OFFSET-aEvent.lep_mom_truth)/aEvent.lep_mom_truth,rweight);
histthetasyst->Fill((180./3.14159)*TMath::ACos(dcosztruth),fabs((180./3.14159)*TMath::ACos(dcosz)-(180./3.14159)*TMath::ACos(dcosztruth))/(180./3.14159)*TMath::ACos(dcosztruth),rweight);
histpur2->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);
numccevents_aftercuts=numccevents_aftercuts+(1.*rweight);
}
else
{
histbkgd->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);
histbkgd2->Fill(aEvent.trk_mom_minos+MOM_OFFSET,rweight);

}


if(aEvent.ccnc_truth==0 && aEvent.mode_truth==0 && aEvent.nuPDG_truth==14)
{
hist4->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);
hist4mue->Fill(aEvent.trk_mom_minos+MOM_OFFSET,rweight);
}

if(aEvent.ccnc_truth==0 && aEvent.mode_truth==1 && aEvent.nuPDG_truth==14)
{
hist5->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);
hist5mue->Fill(aEvent.trk_mom_minos+MOM_OFFSET,rweight);
}
if(aEvent.ccnc_truth==0 && aEvent.mode_truth==2 && aEvent.nuPDG_truth==14)
{
hist6->Fill((180./3.14159)*TMath::ACos(dcosz),rweight);
hist6mue->Fill(aEvent.trk_mom_minos+MOM_OFFSET,rweight);
}

}


}

}





}


//TFile *f=new TFile("datakinreco_hist.root");//default
TFile *f=new TFile("datakinreco_hist_space.root");//w/ mitch spacepoint

TTree *dataTree = (TTree*)f->Get("analysistree/anatree");
datakin adataEvent(dataTree);
Long64_t nentries2 = dataTree.GetEntries();

for (int event = 0; event<nentries2; event++) 
{

 
if(event%10000 == 0)
cout<<event<<" / "<<nentries2<<endl;

adataEvent.GetEntry(event);  


double dcosydata=(adataEvent.trackstart_dcosy_reco*TMath::Cos(.0583497))+(adataEvent.trackstart_dcosz_reco*TMath::Sin(.0583497));
double dcoszdata=(adataEvent.trackstart_dcosz_reco*TMath::Cos(.0583497))-(adataEvent.trackstart_dcosy_reco*TMath::Sin(.0583497));

if(
adataEvent.vtxx_reco>0.+FIDX_CUT && adataEvent.trackstart_x_reco>0.+FIDX_CUT 
&& adataEvent.vtxx_reco<47.-FIDX_CUT && adataEvent.trackstart_x_reco<47.-FIDX_CUT  
&& adataEvent.vtxy_reco>-20.+FIDY_CUT && adataEvent.trackstart_y_reco>-20.+FIDY_CUT 
&& adataEvent.vtxy_reco<20.-FIDY_CUT  && adataEvent.trackstart_y_reco<20.-FIDY_CUT
&& adataEvent.vtxz_reco>0.+FIDZup_CUT && adataEvent.trackstart_z_reco>0.+FIDZup_CUT 
&& adataEvent.vtxz_reco<90.-FIDZdown_CUT && adataEvent.trackstart_z_reco<90.-FIDZdown_CUT 
)
{



if(adataEvent.nmatched_reco==1&&adataEvent.trk_charge_minos<0&&adataEvent.trk_mom_minos>0)
{


double lardirectionEnd[3]={adataEvent.trackexit_dcosx_reco, adataEvent.trackexit_dcosy_reco, adataEvent.trackexit_dcosz_reco};
double larEnd[3]= { adataEvent.trackexit_x_reco, adataEvent.trackexit_y_reco, adataEvent.trackexit_z_reco};
double minosdirectionStart[3]={ adataEvent.trk_dcosx_minos, adataEvent.trk_dcosy_minos, adataEvent.trk_dcosz_minos };

double minosStart[3]={adataEvent.trk_vtxx_minos, adataEvent.trk_vtxy_minos, adataEvent.trk_vtxz_minos };
double xdiff,ydiff,rdiff,totaldiff,thetadiff;

double minosvtx[3] = {(adataEvent.trk_vtxx_minos*100.)-117.4, (adataEvent.trk_vtxy_minos*100.)+19.3, (adataEvent.trk_vtxz_minos*100.)+147.1};

double distance=sqrt(pow(adataEvent.vtxx_reco-minosvtx[0],2)+pow(adataEvent.vtxy_reco-minosvtx[1],2)+pow(adataEvent.vtxz_reco-minosvtx[2],2));

MOM_OFFSET=(0.002198*distance)-.044; 


project(lardirectionEnd, larEnd, minosdirectionStart,minosStart,xdiff,ydiff,rdiff,totaldiff,thetadiff);
if(rdiff>RDIFF_CUT||thetadiff>THETADIFF_CUT)
continue;


hist11->Fill((180./3.14159)*TMath::ACos(dcoszdata));
hist11effcorrected->Fill((180./3.14159)*TMath::ACos(dcoszdata));
hist11effcorrectedd->Fill((180./3.14159)*TMath::ACos(dcoszdata));
hist22effcorrected->Fill((180./3.14159)*TMath::ACos(dcoszdata));
hist11mue->Fill(adataEvent.trk_mom_minos+MOM_OFFSET);
hist11mueeffcorrected->Fill(adataEvent.trk_mom_minos+MOM_OFFSET);
hist11mueeffcorrectedd->Fill(adataEvent.trk_mom_minos+MOM_OFFSET);
hist22mueeffcorrected->Fill(adataEvent.trk_mom_minos+MOM_OFFSET);


}
// 
 }
// 
 }
//   


TChain ch2("KinTree2");



//ch2.Add("simkinreco_tgmuon_hist4.root/analysistree/anatree");//default

 ch2.Add("simkinreco_tgmuon_hist_space.root/analysistree/anatree");//mitch's spacepoints


 datakin adataEvent2(KinTree2);


// TFile *ff=new TFile("simkinreco_tgmuon_hist2.root");
// 
// 
// TTree *dataTree2 = (TTree*)ff->Get("analysistree/anatree");
// datakin adataEvent2(dataTree2);
 Long64_t nentries2 = KinTree2.GetEntries();

for (int event = 0; event<nentries2; event++) 
{

totpot2+=adataEvent2.pot/1000.;
 
if(event%10000 == 0)
cout<<event<<" / "<<nentries2<<endl;

adataEvent2.GetEntry(event);  


double dcosydata2=(adataEvent2.trackstart_dcosy_reco*TMath::Cos(.0583497))+(adataEvent2.trackstart_dcosz_reco*TMath::Sin(.0583497));
double dcoszdata2=(adataEvent2.trackstart_dcosz_reco*TMath::Cos(.0583497))-(adataEvent2.trackstart_dcosy_reco*TMath::Sin(.0583497));

if(
adataEvent2.vtxx_reco>0.+FIDX_CUT && adataEvent2.trackstart_x_reco>0.+FIDX_CUT 
&& adataEvent2.vtxx_reco<47.-FIDX_CUT && adataEvent2.trackstart_x_reco<47.-FIDX_CUT  
&& adataEvent2.vtxy_reco>-20.+FIDY_CUT && adataEvent2.trackstart_y_reco>-20.+FIDY_CUT 
&& adataEvent2.vtxy_reco<20.-FIDY_CUT  && adataEvent2.trackstart_y_reco<20.-FIDY_CUT
&& adataEvent2.vtxz_reco>0.+FIDZup_CUT && adataEvent2.trackstart_z_reco>0.+FIDZup_CUT 
&& adataEvent2.vtxz_reco<90.-FIDZdown_CUT && adataEvent2.trackstart_z_reco<90.-FIDZdown_CUT 
)
{


if(adataEvent2.nmatched_reco==1&&adataEvent2.trk_charge_minos<0&&adataEvent2.trk_mom_minos>0)
{
double lardirectionEnd[3]={adataEvent2.trackexit_dcosx_reco, adataEvent2.trackexit_dcosy_reco, adataEvent2.trackexit_dcosz_reco};
double larEnd[3]= { adataEvent2.trackexit_x_reco, adataEvent2.trackexit_y_reco, adataEvent2.trackexit_z_reco};
double minosdirectionStart[3]={ adataEvent2.trk_dcosx_minos, adataEvent2.trk_dcosy_minos, adataEvent2.trk_dcosz_minos };
double minosStart[3]={adataEvent2.trk_vtxx_minos, adataEvent2.trk_vtxy_minos, adataEvent2.trk_vtxz_minos };
double xdiff,ydiff,rdiff,totaldiff,thetadiff;

double minosvtx[3] = {(adataEvent2.trk_vtxx_minos*100.)-117.4, (adataEvent2.trk_vtxy_minos*100.)+19.3, (adataEvent2.trk_vtxz_minos*100.)+147.1};

double distance=sqrt(pow(adataEvent2.vtxx_reco-minosvtx[0],2)+pow(adataEvent2.vtxy_reco-minosvtx[1],2)+pow(adataEvent2.vtxz_reco-minosvtx[2],2));

MOM_OFFSET=(0.0022054*distance)-.04507; 

project(lardirectionEnd, larEnd, minosdirectionStart,minosStart,xdiff,ydiff,rdiff,totaldiff,thetadiff);
if(rdiff>RDIFF_CUT||thetadiff>THETADIFF_CUT)
continue;

histtgmuon->Fill((180./3.14159)*TMath::ACos(dcoszdata2));
histtgmuon2->Fill(adataEvent2.trk_mom_minos+MOM_OFFSET);


///////////////////
hist33->Fill((180./3.14159)*TMath::ACos(dcoszdata2));
hist33mue->Fill(adataEvent2.trk_mom_minos+MOM_OFFSET);
int bin=((180./3.14159)*TMath::ACos(dcoszdata2)-0.)*(36./72.)-.0000001;
if(bin<30&&bin>=0)
{
if(((180./3.14159)*TMath::ACos(dcoszdata2)+((180./3.14159)*TMath::ACos(dcoszdata2)*thetasyst[bin]))>0.)
hist33up->Fill((180./3.14159)*TMath::ACos(dcoszdata2)+((180./3.14159)*TMath::ACos(dcoszdata2)*thetasyst[bin]));
else
{
hist33up->Fill(.0000001);
}
 //std::cout<<adataEvent2.trackstart_dcosz_reco<<" "<<bin<<" "<<thetasyst[bin]<<" "<<adataEvent2.trackstart_dcosz_reco+(adataEvent2.trackstart_dcosz_reco*thetasyst[bin])<<std::endl;
hist33down->Fill((180./3.14159)*TMath::ACos(dcoszdata2)-((180./3.14159)*TMath::ACos(dcoszdata2)*thetasyst[bin]));
}

int bin2=((adataEvent2.trk_mom_minos+MOM_OFFSET)-0.)*(40./50.);
if(bin2<30&&bin2>=0)
{
hist33mueup->Fill((adataEvent2.trk_mom_minos+MOM_OFFSET)+((adataEvent2.trk_mom_minos+MOM_OFFSET)*momsyst[bin2]));
// std::cout<<aEvent.trk_mom_minos+MOM_OFFSET<<" "<<bin2<<" "<<momsyst[bin2]<<" "<<aEvent.trk_mom_minos+MOM_OFFSET+((aEvent.trk_mom_minos+MOM_OFFSET)*momsyst[bin2])<<std::endl;

if(((adataEvent2.trk_mom_minos+MOM_OFFSET)-((adataEvent2.trk_mom_minos+MOM_OFFSET)*momsyst[bin2]))>0)
hist33muedown->Fill((adataEvent2.trk_mom_minos+MOM_OFFSET)-((adataEvent2.trk_mom_minos+MOM_OFFSET)*momsyst[bin2]));
else
hist33muedown->Fill(.000001);
}




}
 
 }

 }



TChain ch3("KinTree3");




//ch3.Add("simkinreco_nutgmuon_hist3.root/analysistree/anatree");//default


ch3.Add("simkinreco_nutgmuon_space.root/analysistree/anatree");//mitch's spacepoints

 datakin adataEvent3(KinTree3);


 Long64_t nentries2 = KinTree3.GetEntries();



// TFile *fff=new TFile("simkinreco_nutgmuon_hist.root");
// TTree *dataTree3 = (TTree*)fff->Get("analysistree/anatree");
// datakin adataEvent3(dataTree3);
// Long64_t nentries2 = dataTree3.GetEntries();


if(adataEvent3.enu_truth>=0&&adataEvent3.enu_truth<3)
rweight=weight[0];
if(adataEvent3.enu_truth>=3&&adataEvent3.enu_truth<4)
rweight=weight[1];
if(adataEvent3.enu_truth>=4&&adataEvent3.enu_truth<5)
rweight=weight[2];
if(adataEvent3.enu_truth>=5&&adataEvent3.enu_truth<7)
rweight=weight[3];
if(adataEvent3.enu_truth>=7&&adataEvent3.enu_truth<9)
rweight=weight[4];
if(adataEvent3.enu_truth>=9&&adataEvent3.enu_truth<12)
rweight=weight[5];
if(adataEvent3.enu_truth>=12&&adataEvent3.enu_truth<15)
rweight=weight[6];
if(adataEvent3.enu_truth>=15&&adataEvent3.enu_truth<18)
rweight=weight[7];
if(adataEvent3.enu_truth>=18&&adataEvent3.enu_truth<22)
rweight=weight[8];
if(adataEvent3.enu_truth>=22&&adataEvent3.enu_truth<26)
rweight=weight[9];
if(adataEvent3.enu_truth>=26&&adataEvent3.enu_truth<30)
rweight=weight[10];
if(adataEvent3.enu_truth>=30&&adataEvent3.enu_truth<36)
rweight=weight[11];
if(adataEvent3.enu_truth>=36&&adataEvent3.enu_truth<42)
rweight=weight[12];
if(adataEvent3.enu_truth>=42&&adataEvent3.enu_truth<50)
rweight=weight[13];

for (int event = 0; event<nentries2; event++) 
{

totpot3+=adataEvent3.pot/100.;
 
if(event%10000 == 0)
cout<<event<<" / "<<nentries2<<endl;

adataEvent3.GetEntry(event);  


double dcosydata3=(adataEvent3.trackstart_dcosy_reco*TMath::Cos(.0583497))+(adataEvent3.trackstart_dcosz_reco*TMath::Sin(.0583497));
double dcoszdata3=(adataEvent3.trackstart_dcosz_reco*TMath::Cos(.0583497))-(adataEvent3.trackstart_dcosy_reco*TMath::Sin(.0583497));

double dcosytruth3=(adataEvent3.lep_dcosy_truth*TMath::Cos(.0583497))+(adataEvent3.lep_dcosz_truth*TMath::Sin(.0583497));
double dcosztruth3=(adataEvent3.lep_dcosz_truth*TMath::Cos(.0583497))-(adataEvent3.lep_dcosy_truth*TMath::Sin(.0583497));

if(
adataEvent3.vtxx_reco>0.+FIDX_CUT && adataEvent3.trackstart_x_reco>0.+FIDX_CUT 
&& adataEvent3.vtxx_reco<47.-FIDX_CUT && adataEvent3.trackstart_x_reco<47.-FIDX_CUT  
&& adataEvent3.vtxy_reco>-20.+FIDY_CUT && adataEvent3.trackstart_y_reco>-20.+FIDY_CUT 
&& adataEvent3.vtxy_reco<20.-FIDY_CUT  && adataEvent3.trackstart_y_reco<20.-FIDY_CUT
&& adataEvent3.vtxz_reco>0.+FIDZup_CUT && adataEvent3.trackstart_z_reco>0.+FIDZup_CUT 
&& adataEvent3.vtxz_reco<90.-FIDZdown_CUT && adataEvent3.trackstart_z_reco<90.-FIDZdown_CUT 
)
{



if(adataEvent3.nmatched_reco==1&&adataEvent3.trk_charge_minos<0&&adataEvent3.trk_mom_minos>0)
{
double lardirectionEnd[3]={adataEvent3.trackexit_dcosx_reco, adataEvent3.trackexit_dcosy_reco, adataEvent3.trackexit_dcosz_reco};
double larEnd[3]= { adataEvent3.trackexit_x_reco, adataEvent3.trackexit_y_reco, adataEvent3.trackexit_z_reco};
double minosdirectionStart[3]={ adataEvent3.trk_dcosx_minos, adataEvent3.trk_dcosy_minos, adataEvent3.trk_dcosz_minos };
double minosStart[3]={adataEvent3.trk_vtxx_minos, adataEvent3.trk_vtxy_minos, adataEvent3.trk_vtxz_minos };
double xdiff,ydiff,rdiff,totaldiff,thetadiff;

double minosvtx[3] = {(adataEvent3.trk_vtxx_minos*100.)-117.4, (adataEvent3.trk_vtxy_minos*100.)+19.3, (adataEvent3.trk_vtxz_minos*100.)+147.1};

double distance=sqrt(pow(adataEvent3.vtxx_reco-minosvtx[0],2)+pow(adataEvent3.vtxy_reco-minosvtx[1],2)+pow(adataEvent3.vtxz_reco-minosvtx[2],2));

MOM_OFFSET=(0.0022054*distance)-.04507; 

project(lardirectionEnd, larEnd, minosdirectionStart,minosStart,xdiff,ydiff,rdiff,totaldiff,thetadiff);
if(rdiff>RDIFF_CUT||thetadiff>THETADIFF_CUT)
continue;

if(adataEvent3.ccnc_truth==0&&adataEvent3.nuPDG_truth==14)
{

    if(totaldiff<15.)//mean is about 8.5. 
	{
    histnutgmuon3->Fill((180./3.14159)*TMath::ACos(dcoszdata3),rweight);
    histnutgmuon4->Fill(adataEvent3.trk_mom_minos+MOM_OFFSET,rweight);
	histmomsyst2->Fill(adataEvent3.lep_mom_truth, fabs(adataEvent3.trk_mom_minos+MOM_OFFSET-adataEvent3.lep_mom_truth)/adataEvent3.lep_mom_truth,rweight);	histthetasyst2->Fill((180./3.14159)*TMath::ACos(dcosztruth3),fabs((180./3.14159)*TMath::ACos(dcoszdata3)-(180./3.14159)*TMath::ACos(dcosztruth3))/(180./3.14159)*TMath::ACos(dcosztruth3),rweight);

	
	
	}
}
else
{
histnutgmuon->Fill((180./3.14159)*TMath::ACos(dcoszdata3),rweight);
histnutgmuon2->Fill(adataEvent3.trk_mom_minos+MOM_OFFSET,rweight);


hist333->Fill((180./3.14159)*TMath::ACos(dcoszdata3),rweight);
hist333mue->Fill(adataEvent3.trk_mom_minos+MOM_OFFSET,rweight);
int bin=((180./3.14159)*TMath::ACos(dcoszdata3)-0.)*(36./72.)-.0000001;
if(bin<30&&bin>=0)
{
if(((180./3.14159)*TMath::ACos(dcoszdata3)+((180./3.14159)*TMath::ACos(dcoszdata3)*thetasyst[bin]))>0.)
hist333up->Fill((180./3.14159)*TMath::ACos(dcoszdata3)+((180./3.14159)*TMath::ACos(dcoszdata3)*thetasyst[bin]),rweight);
else
{
hist333up->Fill(.0000001,rweight);
}
 //std::cout<<adataEvent3.trackstart_dcosz_reco<<" "<<bin<<" "<<thetasyst[bin]<<" "<<adataEvent3.trackstart_dcosz_reco+(adataEvent3.trackstart_dcosz_reco*thetasyst[bin])<<std::endl;
hist333down->Fill((180./3.14159)*TMath::ACos(dcoszdata3)-((180./3.14159)*TMath::ACos(dcoszdata3)*thetasyst[bin]),rweight);
}

int bin2=((aEvent.trk_mom_minos+MOM_OFFSET)-0.)*(40./50.);
if(bin2<30&&bin2>=0)
{
hist333mueup->Fill((adataEvent3.trk_mom_minos+MOM_OFFSET)+((adataEvent3.trk_mom_minos+MOM_OFFSET)*momsyst[bin2]),rweight);
// std::cout<<aEvent.trk_mom_minos+MOM_OFFSET<<" "<<bin2<<" "<<momsyst[bin2]<<" "<<aEvent.trk_mom_minos+MOM_OFFSET+((aEvent.trk_mom_minos+MOM_OFFSET)*momsyst[bin2])<<std::endl;

if(((adataEvent3.trk_mom_minos+MOM_OFFSET)-((adataEvent3.trk_mom_minos+MOM_OFFSET)*momsyst[bin2]))>0)
hist333muedown->Fill((adataEvent3.trk_mom_minos+MOM_OFFSET)-((adataEvent3.trk_mom_minos+MOM_OFFSET)*momsyst[bin2]),rweight);
else
hist333muedown->Fill(.000001,rweight);
}



}


}
 
 }

 }

histmomsyst->Add(histmomsyst2,totpot/totpot3);
histthetasyst->Add(histthetasyst2,totpot/totpot3);


std::cout<<"POT in file "<<totpot<<std::endl;
std::cout<<"POT in tgmuon file "<<totpot2<<std::endl;
std::cout<<"POT in nutgmuon file "<<totpot3<<std::endl;

double potfactor=8.518840e18/totpot;

double potfactor2=8.518840e18/totpot2;

double potfactor3=8.518840e18/totpot3;

std::cout<<" pot factor "<<potfactor<<" "<<potfactor2<<" "<<potfactor3<<std::endl;

std::cout<<numccevents_beforecuts<<" "<<numccevents_aftercuts<<std::endl;
double efficiency=(double)numccevents_aftercuts/(double)numccevents_beforecuts;
double purity=(double)numccevents_aftercuts/(double)numevents;

std::cout<<"Efficiency "<<efficiency<<" , Purity "<<purity<<std::endl;

char expression[64], expressionn[64], expressionnn[64], expressionnnn[64], expressionnnnn[64], expressionnnnnn[64], expressionnnnnnn[64], expressionnnnnnnn[64];
sprintf(expression,"All flavors/interactions");
TLatex* latex = new TLatex(.5, .8, expression);
latex->SetNDC();
latex->SetTextSize(.038);
latex->SetTextColor(1);

sprintf(expressionn,"NC");
TLatex* latexx = new TLatex(.5, .75, expressionn);
latexx->SetNDC();
latexx->SetTextSize(.038);
latexx->SetTextColor(2);


sprintf(expressionnn,"#nu_{#mu} CC");
TLatex* latexxx = new TLatex(.5, .7, expressionnn);
latexxx->SetNDC();
latexxx->SetTextSize(.038);
latexxx->SetTextColor(3);

sprintf(expressionnnn,"#nu_{#mu} CCQE");
TLatex* latexxxx = new TLatex(.5, .65, expressionnnn);
latexxxx->SetNDC();
latexxxx->SetTextSize(.038);
latexxxx->SetTextColor(4);

sprintf(expressionnnnn,"#nu_{#mu} RES");
TLatex* latexxxxx = new TLatex(.5, .6, expressionnnnn);
latexxxxx->SetNDC();
latexxxxx->SetTextSize(.038);
latexxxxx->SetTextColor(7);

sprintf(expressionnnnnn,"#nu_{#mu} DIS");
TLatex* latexxxxxx = new TLatex(.5, .55, expressionnnnnn);
latexxxxxx->SetNDC();
latexxxxxx->SetTextSize(.038);
latexxxxxx->SetTextColor(6);

sprintf(expressionnnnnnn,"Recosim");
TLatex* latexxxxxxx = new TLatex(.6, .8, expressionnnnnnn);
latexxxxxxx->SetNDC();
latexxxxxxx->SetTextSize(.038);
latexxxxxxx->SetTextColor(2);

sprintf(expressionnnnnnnn,"GENIE expectation");
TLatex* latexxxxxxxx = new TLatex(.5, .8, expressionnnnnnnn);
latexxxxxxxx->SetNDC();
latexxxxxxxx->SetTextSize(.038);
latexxxxxxxx->SetTextColor(2);




char expn[64], expnn[64], expnnn[64], expnnnn[64], expnnnnn[64], expnnnnnn[64], expnnnnnnn[64], expnnnnnnnn[64];
sprintf(expn,"Total error (in quad.)");
TLatex* lat = new TLatex(.5, .8, expn);
lat->SetNDC();
lat->SetTextSize(.038);
lat->SetTextColor(1);

sprintf(expnn,"Stat. (data and MC)");
TLatex* latx = new TLatex(.5, .75, expnn);
latx->SetNDC();
latx->SetTextSize(.038);
latx->SetTextColor(2);


sprintf(expnnn,"# targets");
TLatex* latxx = new TLatex(.5, .65, expnnn);
latxx->SetNDC();
latxx->SetTextSize(.038);
latxx->SetTextColor(3);

sprintf(expnnnn,"Flux");
TLatex* latxxx = new TLatex(.5, .7, expnnnn);
latxxx->SetNDC();
latxxx->SetTextSize(.038);
latxxx->SetTextColor(4);

sprintf(expnnnnn,"Measurement scale");
TLatex* latxxxx = new TLatex(.5, .6, expnnnnn);
latxxxx->SetNDC();
latxxxx->SetTextSize(.038);
latxxxx->SetTextColor(6);

sprintf(expnnnnnn,"POT counting");
TLatex* latxxxxx = new TLatex(.5, .55, expnnnnnn);
latxxxxx->SetNDC();
latxxxxx->SetTextSize(.038);
latxxxxx->SetTextColor(7);

// sprintf(expnnnnnnn,"Lost energy estimate");
// TLatex* latxxxxxx = new TLatex(.5, .45, expnnnnnnn);
// latxxxxxx->SetNDC();
// latxxxxxx->SetTextSize(.038);
// latxxxxxx->SetTextColor(2);

sprintf(expnnnnnnnn,"Fiducial volume");
TLatex* latxxxxxxx = new TLatex(.5, .5, expnnnnnnnn);
latxxxxxxx->SetNDC();
latxxxxxxx->SetTextSize(.038);
latxxxxxxx->SetTextColor(9);

char expn2[64], expnn2[64], expnnn2[64], expnnnn2[64], expnnnnn2[64], expnnnnnn2[64], expnnnnnnn2[64], expnnnnnnnn2[64];
sprintf(expn2,"Total error (in quad.)");
TLatex* lat2 = new TLatex(.2, .8, expn2);
lat2->SetNDC();
lat2->SetTextSize(.038);
lat2->SetTextColor(1);

sprintf(expnn2,"Stat. (data and MC)");
TLatex* latx2 = new TLatex(.2, .75, expnn2);
latx2->SetNDC();
latx2->SetTextSize(.038);
latx2->SetTextColor(2);


sprintf(expnnn2,"# targets");
TLatex* latxx2 = new TLatex(.2, .65, expnnn2);
latxx2->SetNDC();
latxx2->SetTextSize(.038);
latxx2->SetTextColor(3);

sprintf(expnnnn2,"Flux");
TLatex* latxxx2 = new TLatex(.2, .7, expnnnn2);
latxxx2->SetNDC();
latxxx2->SetTextSize(.038);
latxxx2->SetTextColor(4);

sprintf(expnnnnn2,"Measurement scale");
TLatex* latxxxx2 = new TLatex(.2, .6, expnnnnn2);
latxxxx2->SetNDC();
latxxxx2->SetTextSize(.038);
latxxxx2->SetTextColor(6);

sprintf(expnnnnnn2,"POT counting");
TLatex* latxxxxx2 = new TLatex(.2, .55, expnnnnnn2);
latxxxxx2->SetNDC();
latxxxxx2->SetTextSize(.038);
latxxxxx2->SetTextColor(7);

// sprintf(expnnnnnnn2,"Lost energy estimate");
// TLatex* latxxxxxx2 = new TLatex(.2, .45, expnnnnnnn2);
// latxxxxxx2->SetNDC();
// latxxxxxx2->SetTextSize(.038);
// latxxxxxx2->SetTextColor(2);

sprintf(expnnnnnnnn2,"Fiducial volume");
TLatex* latxxxxxxx2 = new TLatex(.2, .5, expnnnnnnnn2);
latxxxxxxx2->SetNDC();
latxxxxxxx2->SetTextSize(.038);
latxxxxxxx2->SetTextColor(9);




leg = new TLegend(0.5,0.83,0.6,0.89);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->AddEntry(hist3_minos_data,"Data (w/ stat. and total error)","lep");

leg2 = new TLegend(0.5,0.83,0.6,0.89);
leg2->SetFillColor(0);
leg2->SetBorderSize(0);
leg2->AddEntry(hist3_minos_data,"Data (w/ stat. and total error)","lep");
leg2->SetTextSize(.04);

leg3 = new TLegend(0.5,0.83,0.6,0.89);
leg3->SetFillColor(0);
leg3->SetBorderSize(0);
leg3->AddEntry(hist3_minos_data,"Data (w/ total error)","lep");
leg3->SetTextSize(.04);

TLine l;

char ep[64],epp[64],eppp[64];

sprintf(ep,"Central value");
TLatex* lt = new TLatex(.6, .6, ep);
lt->SetNDC();
lt->SetTextSize(.038);
lt->SetTextColor(1);

sprintf(epp,"+1 #sigma");
TLatex* ltt = new TLatex(.6, .65, epp);
ltt->SetNDC();
ltt->SetTextSize(.038);
ltt->SetTextColor(2);

sprintf(eppp,"-1 #sigma");
TLatex* lttt = new TLatex(.6, .55, eppp);
lttt->SetNDC();
lttt->SetTextSize(.038);
lttt->SetTextColor(3);

// new TCanvas;
//  
//  hist11->SetTitle("Muon cos(#theta_{#mu}) data and recosim (after cuts)");
//  hist11->SetXTitle("cos(#theta_{#mu})");
//  hist11->SetYTitle("events");
//  hist11->Draw("E1"); 
//  hist2->SetLineColor(2);
//  hist2->DrawNormalized("SAME, E1", hist2->Integral()*potfactor);
//  hist3->SetLineColor(3);
//  hist3->DrawNormalized("SAME, E1", hist3->Integral()*potfactor);
//  hist4->SetLineColor(4);
//  hist4->DrawNormalized("SAME, E1", hist4->Integral()*potfactor);
//  hist5->SetLineColor(7);
//  hist5->DrawNormalized("SAME, E1", hist5->Integral()*potfactor);
//  hist6->SetLineColor(6);
//  hist6->DrawNormalized("SAME, E1", hist6->Integral()*potfactor);
//  hist1->DrawNormalized("SAME, E1", hist1->Integral()*potfactor);
// 
// 
// 
// latex->Draw();
// latexx->Draw();
// latexxx->Draw();
// latexxxx->Draw();
// latexxxxx->Draw();
// latexxxxxx->Draw();
// leg->Draw();
// std::cout<<hist3->Integral()*potfactor<<std::endl;
// 
// new TCanvas;
//  hist11mue->SetTitle("Muon momentum data and recosim (after cuts)");
//  hist11mue->SetXTitle("momentum (GeV/c)");
//  hist11mue->SetYTitle("events");
//  hist11mue->Draw("E1");
// 
//  hist2mue->SetLineColor(2);
//  hist2mue->DrawNormalized("SAME, E1", hist2mue->Integral()*potfactor);
//  hist3mue->SetLineColor(3);
//  hist3mue->DrawNormalized("SAME, E1", hist3mue->Integral()*potfactor);
//  hist4mue->SetLineColor(4);
//  hist4mue->DrawNormalized("SAME, E1", hist4mue->Integral()*potfactor);
//  hist5mue->SetLineColor(7);
//  hist5mue->DrawNormalized("SAME, E1", hist5mue->Integral()*potfactor);
//  hist6mue->SetLineColor(6);
//  hist6mue->DrawNormalized("SAME, E1", hist6mue->Integral()*potfactor);
//  hist1mue->DrawNormalized("SAME, E1", hist1mue->Integral()*potfactor);
// 
// latex->Draw();
// latexx->Draw();
// latexxx->Draw();
// latexxxx->Draw();
// latexxxxx->Draw();
// latexxxxxx->Draw();
// leg->Draw();

// new TCanvas;
// 
// std::cout<<"histeff1 15th bin error "<<histeff1->GetBinError(15)<<" "<<histeff1->GetBinContent(15)<<std::endl;
// histeff1->Draw();
// 
// 
// new TCanvas;
// 
// std::cout<<"histeff2 15th bin error "<<histeff2->GetBinError(15)<<" "<<histeff2->GetBinContent(15)<<std::endl;
// histeff2->Draw();
gStyle->SetOptStat(0);
new TCanvas;
histeff111->Scale(potfactor);
histeff111->DrawCopy();

new TCanvas;
histeff1111->Scale(potfactor);
histeff1111->DrawCopy();
gStyle->SetOptStat(1001001100);

 new TCanvas;
  gStyle->SetStatX(.9);
gStyle->SetStatY(.9);
gStyle->SetStatH(.05);
 histbkgd->SetTitle("NC/WS background (after cuts)");
 histbkgd->SetXTitle("#theta_{#mu, recosim} (degrees)");
 histbkgd->SetYTitle("background events");
 histbkgd->Scale(potfactor);
 histbkgd->DrawCopy();
 gPad->GetCanvas()->Print("histbkgd.pdf");
 
 new TCanvas;
  gStyle->SetStatX(.9);
gStyle->SetStatY(.9);
 histbkgd2->SetTitle("NC/WS background (after cuts)");
 histbkgd2->SetXTitle("P_{#mu, recosim} (GeV/c)");
 histbkgd2->SetYTitle("background events");
 histbkgd2->Scale(potfactor);
 histbkgd2->DrawCopy();
 gPad->GetCanvas()->Print("histbkgd2.pdf");
  
 hhisteff111->Divide(histeff1);
 hhisteff1111->Divide(hhisteff11);
 new TCanvas;
 hhisteff111->DrawCopy();
 new TCanvas;
 hhisteff1111->DrawCopy();
 
 new TCanvas;
//   gStyle->SetStatX(.45);
// gStyle->SetStatY(.8);
 histtgmuon->SetTitle("TG muon background (after cuts)");
 histtgmuon->SetXTitle("#theta_{#mu, recosim} (degrees) ");
 histtgmuon->SetYTitle("background events");
 histtgmuon->Scale(potfactor2);
 histtgmuon->Divide(hhisteff111);
 histtgmuon->DrawCopy();
 gPad->GetCanvas()->Print("histtgmuon.pdf");
 
 
 new TCanvas;
  gStyle->SetStatX(.9);
gStyle->SetStatY(.9);
 histtgmuon2->SetTitle("TG muon background (after cuts)");
 histtgmuon2->SetXTitle("P_{#mu, recosim} (GeV/c)");
 histtgmuon2->SetYTitle("background events");
 histtgmuon2->Scale(potfactor2);
 histtgmuon2->Divide(hhisteff1111);
 histtgmuon2->DrawCopy();
 gPad->GetCanvas()->Print("histtgmuon2.pdf");
 new TCanvas;
//   gStyle->SetStatX(.45);
// gStyle->SetStatY(.8);
 histnutgmuon->SetTitle("NC/WS mismatch to TG muon background (after cuts)");
 histnutgmuon->SetXTitle("#theta_{#mu, recosim} (degrees) ");
 histnutgmuon->SetYTitle("background events");
 histnutgmuon->Scale(potfactor3);
 histnutgmuon->Divide(hhisteff111);
 histnutgmuon->DrawCopy();
 gPad->GetCanvas()->Print("histnutgmuon.pdf");
 new TCanvas;
  gStyle->SetStatX(.9);
gStyle->SetStatY(.9);
 
 histnutgmuon2->SetTitle("NC/WS mismatch to TG muon background (after cuts)");
 histnutgmuon2->SetXTitle("P_{#mu, recosim} (GeV/c)");
 histnutgmuon2->SetYTitle("background events");
 histnutgmuon2->Scale(potfactor3);
 histnutgmuon2->Divide(hhisteff1111);
 histnutgmuon2->DrawCopy();
gPad->GetCanvas()->Print("histnutgmuon2.pdf");




 new TCanvas;
//   gStyle->SetStatX(.45);
// gStyle->SetStatY(.8);
 histnutgmuon3->SetTitle("CC #nu_{#mu} mismatch with TG muon (after cuts)");
 histnutgmuon3->SetXTitle("#theta_{#mu, recosim} (degrees) ");
 histnutgmuon3->SetYTitle("background events");
 histnutgmuon3->Scale(potfactor3);
 histnutgmuon3->DrawCopy();
 gPad->GetCanvas()->Print("histnutgmuon3.pdf");
 new TCanvas;
  gStyle->SetStatX(.9);
gStyle->SetStatY(.9);
 histnutgmuon4->SetTitle("CC #nu_{#mu} mismatch with TG muon (after cuts)");
 histnutgmuon4->SetXTitle("P_{#mu, recosim} (GeV/c)");
 histnutgmuon4->SetYTitle("background events");
 histnutgmuon4->Scale(potfactor3);
 histnutgmuon4->DrawCopy();
gPad->GetCanvas()->Print("histnutgmuon4.pdf");



gStyle->SetOptStat(0);
 new TCanvas;
 histpur2->Divide(histpur1);
 histpur2->SetTitle("CC #nu_{#mu} muon reconstruction purity (after cuts)");
 histpur2->SetXTitle("#theta_{#mu, recosim} (degrees) ");
 histpur2->SetYTitle("purity after cuts & matching");
 histpur2->DrawCopy("E1");


new TCanvas;
gPad->SetGridx();
gPad->SetGridy();
histeff22->Divide(histeff11);
 histeff2->Divide(histeff1);
histeff2->SetTitle("CC #nu_{#mu} muon reconstruction probability (after cuts)");
histeff2->SetXTitle("#theta_{#mu, true} (degrees)");
histeff2->SetYTitle("probability after cuts & matching");
std::cout<<"histeff2 15th bin error "<<histeff2->GetBinError(15)<<" "<<histeff2->GetBinContent(15)<<std::endl;
histeff2->DrawCopy("E1");


double thetaintegral_data,thetaintegral_sim,;
 for(int i=1;i<=hist11effcorrected->GetNbinsX();i++)
{
thetaintegral_data+=hist11effcorrected->GetBinContent(i);


// thetaintegral_sim+=hist11effcorrected->GetBinContent(i);

}
 
std::cout<<"Theta integral data "<<thetaintegral_data<<std::endl;

// std::cout<<"Theta integral sim "<<thetaintegral_sim<<std::endl;
 
double momintegral_data,momintegral_sim;
 for(int i=1;i<=hist11mueeffcorrected->GetNbinsX();i++)
{
momintegral_data+=hist11mueeffcorrected->GetBinContent(i);
// momintegral_sim+=hist11mueeffcorrected->GetBinContent(i);

}
 
std::cout<<"Mom integral data "<<momintegral_data<<std::endl;
// std::cout<<"Mom integral sim "<<momintegral_sim<<std::endl;
 

new TCanvas;
gStyle->SetOptStat(1001001100);
   gStyle->SetStatH(.05);
   
hist11effcorrectedd->SetTitle("Muon #theta_{#mu} data (after cuts, not corrected)");
 hist11effcorrectedd->SetXTitle("#theta_{#mu} (degrees)");
 hist11effcorrectedd->SetYTitle("events");
hist11effcorrectedd->Draw();
gPad->GetCanvas()->Print("uncorr1.pdf");



new TCanvas;
gStyle->SetOptStat(1001001100);
hist11mueeffcorrectedd->SetTitle("Muon momentum data (after cuts, not corrected)");
 hist11mueeffcorrectedd->SetXTitle("P_{#mu} (GeV/c)");
hist11mueeffcorrectedd->SetYTitle("events");
hist11mueeffcorrectedd->Draw();
gPad->GetCanvas()->Print("uncorr2.pdf");


new TCanvas;
    gStyle->SetStatX(.9);
   gStyle->SetStatY(.9);
// gStyle->SetStatX(.45);
// gStyle->SetStatY(.8);

histbkgd->SetTitle("CC #nu_{#mu} muon #theta_{#mu} total background (after cuts)");
histbkgd->SetXTitle("#theta_{#mu} (degrees)");
histbkgd->SetYTitle("background events");
histbkgd->Add(histtgmuon,1);
histbkgd->Add(histnutgmuon,1);
// histbkgd->Divide(histeff3);
histbkgd->DrawCopy();
gPad->GetCanvas()->Print("totbackground.pdf");


 gStyle->SetStatX(.9);
   gStyle->SetStatY(.9);

new TCanvas;

histbkgd2->SetTitle("CC #nu_{#mu} muon momentum total background (after cuts)");
histbkgd2->SetXTitle("P_{#mu} (GeV/c)");
histbkgd2->SetYTitle("background events");
histbkgd2->Add(histtgmuon2,1);
histbkgd2->Add(histnutgmuon2,1);
// histbkgd2->Divide(histeff33);
histbkgd2->DrawCopy();
gPad->GetCanvas()->Print("totbackground2.pdf"); 
// gStyle->SetStatX(.45);
// gStyle->SetStatY(.8);
//  
 
////////////////
TH1D *histbkd_par=new TH1D("","",18,0.,36);
//TF1 *gfit = new TF1("pol10","pol10",0.8,1.0);

TF1 *efit = new TF1("f1","landau",0.,36.);//three parametr exponential association fit (3 param)
       gStyle->SetStatX(.9);
   gStyle->SetStatY(.9);
gStyle->SetStatH(0.3);  
efit->SetParNames("constant","MP","sigma");
histbkgd->Fit(efit, "WEMR");
TFitResultPtr r = histbkgd->Fit(efit, "WSEMR");
int nn=18;
double x[18];
for(int i=0;i<nn;i++)
{
x[i]=(double)i*(36/18.)+(1.);
}

double ci[18];
double cl = 0.683;  // for 1 sigma error
r->GetConfidenceIntervals(nn,1,1,x,ci,cl);


for(int i=0;i<histbkgd_par->GetNbinsX();i++)
{
std::cout<<i<<" "<<efit->Eval(x[i])<<" "<<ci[i]<<std::endl;
histbkgd_par->SetBinContent(i+1,efit->Eval(x[i]));
histbkgd_par->SetBinError(i+1,ci[i]);
}

new TCanvas;
   gStyle->SetOptFit(111);
   gStyle->SetStatH(0.2); 
  histbkgd->SetMaximum(3.6); 
histbkgd->DrawCopy();
gPad->GetCanvas()->Print("totbackground_wfit.pdf");



/////////////////////  

//   gStyle->SetStatX(.45);
//    gStyle->SetStatY(.8);
new TCanvas;
   gStyle->SetStatX(.9);
   gStyle->SetStatY(.9);
gStyle->SetStatH(0.15);
histbkgd_par->SetTitle("CC #nu_{#mu} muon #theta_{#mu} total background (after cuts, parameterized)");
histbkgd_par->SetXTitle("#theta_{#mu} (degrees)");
histbkgd_par->SetYTitle("background events (parameterized)");
histbkgd_par->Draw();
gPad->GetCanvas()->Print("histbkgd_par.pdf"); 
   gStyle->SetStatX(.9);
   gStyle->SetStatY(.9);
 new TCanvas;
// gPad->SetGridx();
// gPad->SetGridy();
histeff3->Divide(histeff1);
hist11effcorrected->Add(histbkgd_par,-1);

 TFile *unf2=new TFile("unfold2.root", "RECREATE"); 
 hist11effcorrected->Write();
 unf2->Close();
 

//insert here

//    TFile *cvfi=new TFile("unfoldedCV.root","READ");
//    hist11effcorrected = (TH1D*)cvfi->Get("unfoldedmeascos");
//    hist11effcorrected->SetLineColor(1);



hist11effcorrected->Divide(histeff3);
hist1effcorrected->SetTitle("Muon #theta_{#mu} data and recosim (after cuts, reco-probability corrected)");
 hist1effcorrected->SetXTitle("#theta_{#mu} (degrees)");
 hist1effcorrected->SetYTitle("events (reco-probability corrected)");
hist1effcorrected->Scale(potfactor);
hist1effcorrected->Divide(histeff2);
hist1effcorrected->SetLineColor(2);
 std::cout<<"Theta sim integral "<<hist1effcorrected->Integral()<<std::endl;
// hist1effcorrected->DrawCopy("E1");

// hist11effcorrected->DrawCopy("SAME, E1");
// latexxxxxxx->Draw();
// leg2->Draw();


 new TCanvas;
   TProfile *profthetasyst = histthetasyst->ProfileX("", 0, 100, ""); 
//    profthetasyst->SetTitle("CC #nu_{#mu} muon cos(#theta_{#mu}) fractional error (comparing recosim & truth)");
//    profthetasyst->SetXTitle("cos(#theta_{true})");
//    profthetasyst->SetYTitle("|#Delt(180./3.14159)*TMath::ACos(#theta)|/cos(#theta_{true})");
//    profthetasyst->Draw("E1");

   TFile *fi2=new TFile("cosfverror.root","READ");
   TH1F *cosfverror = (TH1F*)fi2->Get("cosfverror");

double thetaintegral,thetaintegralerror,thetaintegralerror_thetasyst,thetaintegralerror_flux,thetaintegralerror_stat,thetaintegralerror_ntarg;


 


hist3->Scale(potfactor);
hist3up->Scale(potfactor);
hist3down->Scale(potfactor);




hist33->Scale(potfactor2);
hist33up->Scale(potfactor2);
hist33down->Scale(potfactor2);

hist333->Scale(potfactor3);
hist333up->Scale(potfactor3);
hist333down->Scale(potfactor3);

hist3mue->Scale(potfactor);
hist3mueup->Scale(potfactor);
hist3muedown->Scale(potfactor);


TFile *unf3=new TFile("unfold3.root", "RECREATE");
hist3->Write();
hist3up->Write();
hist3down->Write();
hist3mue->Write();
hist3mueup->Write();
hist3muedown->Write();
unf3->Close();

hist33mue->Scale(potfactor2);
hist33mueup->Scale(potfactor2);
hist33muedown->Scale(potfactor2);

hist333mue->Scale(potfactor3);
hist333mueup->Scale(potfactor3);
hist333muedown->Scale(potfactor3);

hist33->Divide(hhisteff111);
hist33up->Divide(hhisteff111);
hist33down->Divide(hhisteff111);
hist333->Divide(hhisteff111);
hist333up->Divide(hhisteff111);
hist333down->Divide(hhisteff111);

hist33mue->Divide(hhisteff1111);
hist33mueup->Divide(hhisteff1111);
hist33muedown->Divide(hhisteff1111);
hist333mue->Divide(hhisteff1111);
hist333mueup->Divide(hhisteff1111);
hist333muedown->Divide(hhisteff1111);

//I don't need to background subtract the neutrino MC. there is no background in the neutrino MC
// hist3->Add(hist33,-1);
// hist3up->Add(hist33up,-1);
// hist3down->Add(hist33down,-1);
// hist3->Add(hist333,-1);
// hist3up->Add(hist333up,-1);
// hist3down->Add(hist333down,-1);
// 
// hist3mue->Add(hist33mue,-1);
// hist3mueup->Add(hist33mueup,-1);
// hist3muedown->Add(hist33muedown,-1);
// hist3mue->Add(hist333mue,-1);
// hist3mueup->Add(hist333mueup,-1);
// hist3muedown->Add(hist333muedown,-1);


TFile *g=new TFile("cos_normalfv.root", "RECREATE");
hist3->Write();
g->Close();

 TFile *gg=new TFile("mom_normalfv.root", "RECREATE");
hist3mue->Write();
gg->Close();


hist3->Divide(histeff2);
hist3up->Divide(histeff2);
hist3down->Divide(histeff2);

hist3mue->Divide(histeff22);
hist3mueup->Divide(histeff22);
hist3muedown->Divide(histeff22);


TFile *cverrortheta=new TFile("unfoldingerror.root","READ");
TH1D *thetaerror = (TH1D*)cverrortheta->Get("thetaerror");

 for(int i=1;i<=hist11effcorrected->GetNbinsX();i++)
{

// double fracthetasyst=((fabs(hist3up->GetBinContent(i)-hist3->GetBinContent(i))+fabs(hist3down->GetBinContent(i)-hist3->GetBinContent(i)))/2.)/hist3->GetBinContent(i);


double fracthetasyst=0.;
if(fabs(hist3up->GetBinContent(i)-hist3->GetBinContent(i))>fabs(hist3down->GetBinContent(i)-hist3->GetBinContent(i)))
fracthetasyst=fabs(hist3up->GetBinContent(i)-hist3->GetBinContent(i))/hist3->GetBinContent(i);
else
fracthetasyst=fabs(hist3down->GetBinContent(i)-hist3->GetBinContent(i))/hist3->GetBinContent(i);


 double totalerror=sqrt(pow(hist11effcorrected->GetBinError(i),2)+pow(delta_INTEGRATEDFLUX_0_50*hist11effcorrected->GetBinContent(i),2)+pow(delta_pot*hist11effcorrected->GetBinContent(i),2)+pow(delta_NTARGETS*hist11effcorrected->GetBinContent(i),2)+pow(thetaerror->GetBinContent(i)*hist11effcorrected->GetBinContent(i),2)+pow(cosfverror->GetBinContent(i)*hist11effcorrected->GetBinContent(i),2));
 
 cos2->SetBinError(i,hist11effcorrected->GetBinError(i));
 histerror00->Fill(hist11effcorrected->GetBinCenter(i),totalerror/hist11effcorrected->GetBinContent(i));
 histerror11->Fill(hist11effcorrected->GetBinCenter(i),hist11effcorrected->GetBinError(i)/hist11effcorrected->GetBinContent(i));
 histerror22->Fill(hist11effcorrected->GetBinCenter(i),(delta_NTARGETS*hist11effcorrected->GetBinContent(i))/hist11effcorrected->GetBinContent(i));
 histerror33->Fill(hist11effcorrected->GetBinCenter(i),(delta_INTEGRATEDFLUX_0_50*hist11effcorrected->GetBinContent(i))/hist11effcorrected->GetBinContent(i));
 histerror44->Fill(hist11effcorrected->GetBinCenter(i),(thetaerror->GetBinContent(i)*hist11effcorrected->GetBinContent(i))/hist11effcorrected->GetBinContent(i));
 histerror55->Fill(hist11effcorrected->GetBinCenter(i),(delta_pot*hist11effcorrected->GetBinContent(i))/hist11effcorrected->GetBinContent(i));
 
 histerror77->Fill(hist11effcorrected->GetBinCenter(i),(cosfverror->GetBinContent(i)*hist11effcorrected->GetBinContent(i))/hist11effcorrected->GetBinContent(i));
 
 histerror00a->Fill(hist11effcorrected->GetBinCenter(i),totalerror);
 
 histerror11a->Fill(hist11effcorrected->GetBinCenter(i),hist11effcorrected->GetBinError(i));
 histerror22a->Fill(hist11effcorrected->GetBinCenter(i),(delta_NTARGETS*hist11effcorrected->GetBinContent(i)));
 histerror33a->Fill(hist11effcorrected->GetBinCenter(i),(delta_INTEGRATEDFLUX_0_50*hist11effcorrected->GetBinContent(i)));
 histerror44a->Fill(hist11effcorrected->GetBinCenter(i),(thetaerror->GetBinContent(i)*hist11effcorrected->GetBinContent(i)));
 histerror55a->Fill(hist11effcorrected->GetBinCenter(i),(delta_pot*hist11effcorrected->GetBinContent(i)));
 
 histerror77a->Fill(hist11effcorrected->GetBinCenter(i),(cosfverror->GetBinContent(i)*hist11effcorrected->GetBinContent(i)));
 
 thetaintegral+=hist11effcorrected->GetBinContent(i);
 
 thetaintegralerror+=pow(totalerror,2);
 
 thetaintegralerror_ntarg+=delta_NTARGETS*hist11effcorrected->GetBinContent(i);
 
 thetaintegralerror_thetasyst+=pow((thetaerror->GetBinContent(i)*hist11effcorrected->GetBinContent(i)),2);
 
 thetaintegralerror_flux+=delta_INTEGRATEDFLUX_0_50*hist11effcorrected->GetBinContent(i);
 
 thetaintegralerror_stat+=pow(hist11effcorrected->GetBinError(i),2);
 
 hist11effcorrected->SetBinError(i,totalerror);
}


std::cout<<histerror00a->Integral()<<" theta total error "<<histerror00a->Integral()/hist11effcorrected->Integral()<<std::endl;

std::cout<<histerror11a->Integral()<<" theta stat error "<<histerror11a->Integral()/hist11effcorrected->Integral()<<std::endl;

std::cout<<histerror22a->Integral()<<" theta ntargets error "<<histerror22a->Integral()/hist11effcorrected->Integral()<<std::endl;

std::cout<<histerror33a->Integral()<<" theta flux error "<<histerror33a->Integral()/hist11effcorrected->Integral()<<std::endl;

std::cout<<histerror44a->Integral()<<" theta measurement error "<<histerror44a->Integral()/hist11effcorrected->Integral()<<std::endl;

std::cout<<histerror55a->Integral()<<" theta pot error "<<histerror55a->Integral()/hist11effcorrected->Integral()<<std::endl;

std::cout<<histerror77a->Integral()<<" theta FV error "<<histerror77a->Integral()/hist11effcorrected->Integral()<<std::endl;


// TH1D* profmomsystup = new TH1D("", "",20,0.,25);
// 
//  for(int i=1;i<=profmomsyst->GetNbinsX();i++)
// {
// std::cout<<profmomsyst->GetBinContent(i)<<std::endl;
// 
// 
// }



// new TCanvas;
// hist1effcorrected->DrawCopy("E1");
// hist11effcorrected->DrawCopy("SAME, E1");

  


  
////////////////
//TF1 *gfit = new TF1("pol7","pol7",0.,25.);

TF1 *gfit = new TF1("f1","[0]*(x-[2])^([1])",0.,25.);//shifted power fit (3 params)
gfit->SetParNames("a","b","c");
       gStyle->SetStatX(.9);
   gStyle->SetStatY(.9);
gStyle->SetStatH(0.4); 
gStyle->SetOptFit(0011);
histbkgd2->Fit(gfit,"WEMR");
TFitResultPtr rr = histbkgd2->Fit( gfit, "WSEMR");
int nnn=20;
double xx[20];

for(int i=0;i<nnn;i++)
{
xx[i]=(double)i*(25./20.)+(25./40.);
}

double cii[20];
double cll = 0.683;  // for 1 sigma error
rr->GetConfidenceIntervals(nnn,1,1,xx,cii,cll);


for(int i=0;i<histbkgd2->GetNbinsX();i++)
{
std::cout<<i<<" "<<gfit->Eval(xx[i])<<std::endl;
histbkgd2_par->SetBinContent(i+1,gfit->Eval(xx[i]));
histbkgd2_par->SetBinError(i+1,cii[i]);
}
/////////////////////  
  
new TCanvas;
gStyle->SetStatH(0.2); 
gStyle->SetOptFit(1);
histbkgd2->DrawCopy();
 gfit->Draw("SAME");
gPad->GetCanvas()->Print("totbackground_wfit2.pdf");
  
new TCanvas;
gStyle->SetStatH(0.2); 
gStyle->SetOptFit(0);
gStyle->SetStatX(.9);
gStyle->SetStatY(.9);
histbkgd2_par->SetTitle("CC #nu_{#mu} muon momentum total background (after cuts, parameterized)");
histbkgd2_par->SetXTitle("P_{#mu} (GeV/c)");
histbkgd2_par->SetYTitle("background events (parameterized)");
histbkgd2_par->Draw();
gPad->GetCanvas()->Print("histbkgd2_par.pdf");

gStyle->SetOptStat(0);


 new TCanvas;
// gPad->SetGridx();
// gPad->SetGridy();


// histeff22->SetTitle("CC #nu_{#mu} muon reconstruction probability (after cuts)");
// histeff22->SetXTitle("momentum (GeV/c)");
// histeff22->SetYTitle("probability after cuts & matching");
// histeff22->Draw("E1");
 
// new TCanvas;
// gPad->SetGridx();
// gPad->SetGridy();
hist11mueeffcorrected->Add(histbkgd2_par,-1);


 TFile *unfmom2=new TFile("unfoldmom2.root", "RECREATE"); 
 hist11mueeffcorrected->Write();
 unfmom2->Close();


// TFile *cvfi2=new TFile("unfoldedCV.root","READ");
// hist11mueeffcorrected = (TH1D*)cvfi2->Get("unfoldedmeasmom");
// hist11mueeffcorrected->SetLineColor(1);

histeff33->Divide(histeff11);
hist11mueeffcorrected->Divide(histeff33);
 hist1mueeffcorrected->SetTitle("Muon momentum data and recosim (after cuts, reco-probability corrected)");
 hist1mueeffcorrected->SetXTitle("P_{#mu} (GeV/c)");
 hist1mueeffcorrected->SetYTitle("events (probability corrected)");
hist1mueeffcorrected->Scale(potfactor);
hist1mueeffcorrected->Divide(histeff22);
hist1mueeffcorrected->SetLineColor(2);
 std::cout<<"Mom integral sim "<<hist1mueeffcorrected->Integral()<<std::endl;
// hist1mueeffcorrected->Draw("E1"); 
// hist11mueeffcorrected->Draw("SAME, E1");
// latexxxxxxxx->Draw();
// leg3->Draw();


TProfile::Approximate(kTRUE);  
//    new TCanvas;
   TProfile *profmomsyst = histmomsyst->ProfileX("", 0, 100, "");
   
//    profmomsyst->SetTitle("CC #nu_{#mu} muon momentum fractional error");
//    profmomsyst->SetXTitle("P_{true}");
//    profmomsyst->SetYTitle("|#DeltaP|/P_{true}");
//    profmomsyst->GetYaxis()->SetTitleOffset(1.25);
//    profmomsyst->Draw("E1");

   TFile *fi=new TFile("momfverror.root","READ");
   TH1F *momfverror = (TH1F*)fi->Get("momfverror");
 
 
   TFile *cverrormom=new TFile("unfoldingerror.root","READ");
   TH1D *momerror = (TH1D*)cverrormom->Get("momerror");
 
 
double momintegral,momintegralerror,momintegralerror_momsyst,momintegralerror_flux,momintegralerror_stat,momintegralerror_ntarg;
 for(int i=1;i<=hist11mueeffcorrected->GetNbinsX();i++)
{

// double fracmomsyst=((fabs(hist3mueup->GetBinContent(i)-hist3mue->GetBinContent(i))+fabs(hist3muedown->GetBinContent(i)-hist3mue->GetBinContent(i)))/2.)/hist3mue->GetBinContent(i);
double fracmomsyst=0.;
if(fabs(hist3mueup->GetBinContent(i)-hist3mue->GetBinContent(i))>fabs(hist3muedown->GetBinContent(i)-hist3mue->GetBinContent(i)))
fracmomsyst=fabs(hist3mueup->GetBinContent(i)-hist3mue->GetBinContent(i))/hist3mue->GetBinContent(i);
else
fracmomsyst=fabs(hist3muedown->GetBinContent(i)-hist3mue->GetBinContent(i))/hist3mue->GetBinContent(i);

double fracmomsyst_lost=0.;
if(fabs(hist8mueup->GetBinContent(i)-hist8mue->GetBinContent(i))>fabs(hist8muedown->GetBinContent(i)-hist8mue->GetBinContent(i)))
fracmomsyst_lost=fabs(hist8mueup->GetBinContent(i)-hist8mue->GetBinContent(i))/hist8mue->GetBinContent(i);
else
fracmomsyst_lost=fabs(hist8muedown->GetBinContent(i)-hist8mue->GetBinContent(i))/hist8mue->GetBinContent(i);


 double totalerror=sqrt(pow(hist11mueeffcorrected->GetBinError(i),2)+pow(delta_pot*hist11mueeffcorrected->GetBinContent(i),2)+pow(delta_INTEGRATEDFLUX_0_50*hist11mueeffcorrected->GetBinContent(i),2)+pow(delta_NTARGETS*hist11mueeffcorrected->GetBinContent(i),2)+pow(momerror->GetBinContent(i)*hist11mueeffcorrected->GetBinContent(i),2)+pow(momfverror->GetBinContent(i)*hist11mueeffcorrected->GetBinContent(i),2));
 
 mom2->SetBinError(i,hist11mueeffcorrected->GetBinError(i));
 histerror0->Fill(hist11mueeffcorrected->GetBinCenter(i),totalerror/hist11mueeffcorrected->GetBinContent(i));
 histerror1->Fill(hist11mueeffcorrected->GetBinCenter(i),hist11mueeffcorrected->GetBinError(i)/hist11mueeffcorrected->GetBinContent(i));
 histerror2->Fill(hist11mueeffcorrected->GetBinCenter(i),(delta_NTARGETS*hist11mueeffcorrected->GetBinContent(i))/hist11mueeffcorrected->GetBinContent(i));
 histerror3->Fill(hist11mueeffcorrected->GetBinCenter(i),(delta_INTEGRATEDFLUX_0_50*hist11mueeffcorrected->GetBinContent(i))/hist11mueeffcorrected->GetBinContent(i));
 histerror4->Fill(hist11mueeffcorrected->GetBinCenter(i),(momerror->GetBinContent(i)*hist11mueeffcorrected->GetBinContent(i))/hist11mueeffcorrected->GetBinContent(i));
 
 histerror5->Fill(hist11mueeffcorrected->GetBinCenter(i),(delta_pot*hist11mueeffcorrected->GetBinContent(i))/hist11mueeffcorrected->GetBinContent(i));
 
 histerror6->Fill(hist11mueeffcorrected->GetBinCenter(i),(fracmomsyst_lost*hist11mueeffcorrected->GetBinContent(i))/hist11mueeffcorrected->GetBinContent(i));
 histerror7->Fill(hist11mueeffcorrected->GetBinCenter(i),(momfverror->GetBinContent(i)*hist11mueeffcorrected->GetBinContent(i))/hist11mueeffcorrected->GetBinContent(i));
 
 histerror0a->Fill(hist11mueeffcorrected->GetBinCenter(i),totalerror);
 
 histerror1a->Fill(hist11mueeffcorrected->GetBinCenter(i),hist11mueeffcorrected->GetBinError(i));
 histerror2a->Fill(hist11mueeffcorrected->GetBinCenter(i),(delta_NTARGETS*hist11mueeffcorrected->GetBinContent(i)));
 histerror3a->Fill(hist11mueeffcorrected->GetBinCenter(i),(delta_INTEGRATEDFLUX_0_50*hist11mueeffcorrected->GetBinContent(i)));
 histerror4a->Fill(hist11mueeffcorrected->GetBinCenter(i),(momerror->GetBinContent(i)*hist11mueeffcorrected->GetBinContent(i)));
 
 histerror5a->Fill(hist11mueeffcorrected->GetBinCenter(i),(delta_pot*hist11mueeffcorrected->GetBinContent(i)));
 
 histerror6a->Fill(hist11mueeffcorrected->GetBinCenter(i),(fracmomsyst_lost*hist11mueeffcorrected->GetBinContent(i)));
 
 histerror7a->Fill(hist11mueeffcorrected->GetBinCenter(i),(momfverror->GetBinContent(i)*hist11mueeffcorrected->GetBinContent(i)));
 
  momintegral+=hist11mueeffcorrected->GetBinContent(i);
  
 momintegralerror+=pow(totalerror,2);
 
 momintegralerror_ntarg+=delta_NTARGETS*hist11mueeffcorrected->GetBinContent(i);
 
 momintegralerror_momsyst+=pow(momerror->GetBinContent(i)*hist11mueeffcorrected->GetBinContent(i),2);
 
 momintegralerror_flux+=delta_INTEGRATEDFLUX_0_50*hist11mueeffcorrected->GetBinContent(i);
 
 momintegralerror_stat+=pow(hist11mueeffcorrected->GetBinError(i),2);
 
 hist11mueeffcorrected->SetBinError(i,totalerror);
}

std::cout<<histerror0a->Integral()<<" mom total error "<<histerror0a->Integral()/hist11mueeffcorrected->Integral()<<std::endl;

std::cout<<histerror1a->Integral()<<" mom stat error "<<histerror1a->Integral()/hist11mueeffcorrected->Integral()<<std::endl;

std::cout<<histerror2a->Integral()<<" mom ntargets error "<<histerror2a->Integral()/hist11mueeffcorrected->Integral()<<std::endl;

std::cout<<histerror3a->Integral()<<" mom flux error "<<histerror3a->Integral()/hist11mueeffcorrected->Integral()<<std::endl;

std::cout<<histerror4a->Integral()<<" mom measurement error "<<histerror4a->Integral()/hist11mueeffcorrected->Integral()<<std::endl;

std::cout<<histerror5a->Integral()<<" mom pot error "<<histerror5a->Integral()/hist11mueeffcorrected->Integral()<<std::endl;

// std::cout<<histerror5a->Integral()<<" mom Elost error "<<histerror6a->Integral()/hist11mueeffcorrected->Integral()<<std::endl;

std::cout<<histerror7a->Integral()<<" mom FV error "<<histerror7a->Integral()/hist11mueeffcorrected->Integral()<<std::endl;


// new TCanvas;
// 
// 
// hist1mueeffcorrected->DrawCopy("E1");
// hist11mueeffcorrected->DrawCopy("SAME, E1");


//TFile *gg=new TFile("mom.root", "RECREATE");
new TCanvas;
histeff1111->SetTitle("CC #nu_{#mu} P_{#mu} (after cuts, reco-probability/background corrected)");
histeff1111->SetXTitle("P_{#mu} (GeV/c)");
histeff1111->SetYTitle("corrected events");
histeff1111->SetLineColor(2);
histeff1111->SetMaximum(215);
histeff1111->SetMinimum(0);
histeff1111->DrawCopy("hist");
hist11mueeffcorrected->DrawCopy("SAME, E1");

// histeff1111->Write();

latexxxxxxxx->Draw();
leg3->Draw();
gPad->GetCanvas()->Print("histmomcorrrate.pdf");




// TFile *g=new TFile("cos.root", "RECREATE");

new TCanvas;
histeff111->SetTitle("CC #nu_{#mu} #theta_{#mu} (after cuts, reco-probability/background corrected)");
histeff111->SetXTitle("#theta_{#mu} (degrees)");
histeff111->SetYTitle("corrected events");
histeff111->SetLineColor(2);
histeff111->SetMaximum(100);
histeff111->SetMinimum(0);
histeff111->DrawCopy("hist");
hist11effcorrected->DrawCopy("SAME, E1");
//histeff111->Write();

latexxxxxxxx->Draw();
leg3->Draw();
gPad->GetCanvas()->Print("histdircoscorrrate.pdf");







thetaintegral_data=0.;

 for(int i=1;i<=hist11effcorrected->GetNbinsX();i++)
{
thetaintegral_data+=hist11effcorrected->GetBinContent(i);

}
 
std::cout<<"Theta integral data2 "<<thetaintegral_data<<std::endl;

 
 
double momintegral_data=0.;
 for(int i=1;i<=hist11mueeffcorrected->GetNbinsX();i++)
{
momintegral_data+=hist11mueeffcorrected->GetBinContent(i);

}
 
std::cout<<"Mom integral data2 "<<momintegral_data<<std::endl;

new TCanvas;

histerror00->SetTitle("Contributions to uncertainty");
histerror00->SetXTitle("#theta_{#mu} (degrees)");
histerror00->SetYTitle("fractional uncertainty");
histerror00->SetMinimum(0);
histerror00->SetMaximum(1.1);
histerror00->Draw();
histerror11->SetLineColor(2);
histerror11->Draw("SAME");
histerror22->SetLineColor(3);
histerror22->Draw("SAME");
histerror33->SetLineColor(4);
histerror33->Draw("SAME");
histerror44->SetLineColor(6);
histerror44->Draw("SAME");
histerror55->SetLineColor(7);
histerror55->Draw("SAME");
histerror77->SetLineColor(9);
histerror77->Draw("SAME");

lat->Draw();
latx->Draw();
latxx->Draw();
latxxx->Draw();
latxxxx->Draw();
latxxxxx->Draw();
latxxxxxxx->Draw();
gPad->GetCanvas()->Print("errorcontr2.pdf");





new TCanvas;

histerror0a->SetTitle("Contributions to uncertainty");
histerror0a->SetXTitle("P_{#mu} (GeV/c)");
histerror0a->SetYTitle("uncertainty (corrected events)");
histerror0a->SetMinimum(0);
histerror0a->Draw();
histerror1a->SetLineColor(2);
histerror1a->Draw("SAME");
histerror2a->SetLineColor(3);
histerror2a->Draw("SAME");
histerror3a->SetLineColor(4);
histerror3a->Draw("SAME");
histerror4a->SetLineColor(6);
histerror4a->Draw("SAME");
histerror5a->SetLineColor(7);
histerror5a->Draw("SAME");
// histerror6a->SetLineColor(2);
// histerror6a->Draw("SAME");
histerror7a->SetLineColor(9);
histerror7a->Draw("SAME");
lat->Draw();
latx->Draw();
latxx->Draw();
latxxx->Draw();
latxxxx->Draw();
latxxxxx->Draw();
// latxxxxxx->Draw();
latxxxxxxx->Draw();
gPad->GetCanvas()->Print("errorcontr3.pdf");





new TCanvas;
histerror00a->SetTitle("Contributions to uncertainty");
histerror00a->SetXTitle("#theta_{#mu} (degrees)");
histerror00a->SetYTitle("uncertainty (corrected events)");
histerror00a->SetMinimum(0);
histerror00a->SetMaximum(26);
histerror00a->Draw();
histerror11a->SetLineColor(2);
histerror11a->Draw("SAME");
histerror22a->SetLineColor(3);
histerror22a->Draw("SAME");
histerror33a->SetLineColor(4);
histerror33a->Draw("SAME");
histerror44a->SetLineColor(6);
histerror44a->Draw("SAME");
histerror55a->SetLineColor(7);
histerror55a->Draw("SAME");
histerror77a->SetLineColor(9);
histerror77a->Draw("SAME");
lat->Draw();
latx->Draw();
latxx->Draw();
latxxx->Draw();
latxxxx->Draw();
latxxxxx->Draw();
latxxxxxxx->Draw();
gPad->GetCanvas()->Print("errorcontr4.pdf");


new TCanvas;
histerror0->SetMaximum(1.05);
histerror0->SetTitle("Contributions to uncertainty");
histerror0->SetXTitle("P_{#mu} (GeV/c)");
histerror0->SetYTitle("fractional uncertainty");
histerror0->SetMinimum(0);
histerror0->Draw();
histerror1->SetLineColor(2);
histerror1->Draw("SAME");
histerror2->SetLineColor(3);
histerror2->Draw("SAME");
histerror3->SetLineColor(4);
histerror3->Draw("SAME");
histerror4->SetLineColor(6);
histerror4->Draw("SAME");
histerror5->SetLineColor(7);
histerror5->Draw("SAME");
// histerror6->SetLineColor(2);
// histerror6->Draw("SAME");
histerror7->SetLineColor(9);
histerror7->Draw("SAME");
lat2->Draw();
latx2->Draw();
latxx2->Draw();
latxxx2->Draw();
latxxxx2->Draw();
latxxxxx2->Draw();
// latxxxxxx2->Draw();
latxxxxxxx2->Draw();
gPad->GetCanvas()->Print("errorcontr1.pdf");

new TCanvas;
TH1D *cos1=(TH1D*)hist11effcorrected->Clone();
for(int i=1;i<=hist11effcorrected->GetNbinsX();i++)
{
cos2->SetBinContent(i,hist11effcorrected->GetBinContent(i));
cos3->SetBinContent(i,hist11effcorrected->GetBinContent(i));
cos3->SetBinError(i,cos2->GetBinError(i));
}



cos1->SetLineColor(1);
cos1->SetFillColor(1);
 cos1->SetFillStyle(0);
cos1->Draw("E1");
cos3->Draw("E1,SAME");

new TCanvas;
TH1D *mom1=(TH1D*)hist11mueeffcorrected->Clone();
for(int i=1;i<=hist11mueeffcorrected->GetNbinsX();i++)
{
mom2->SetBinContent(i,hist11mueeffcorrected->GetBinContent(i));
mom3->SetBinContent(i,hist11mueeffcorrected->GetBinContent(i));
mom3->SetBinError(i,mom2->GetBinError(i));
}


mom1->SetLineColor(1);
mom1->SetFillColor(1);
mom1->SetFillStyle(0);
mom1->Draw("E1");
mom3->Draw("E1,SAME");

new TCanvas;

 hist3->SetTitle("CC #nu_{#mu} muon #theta_{#mu} (measurement scale error)");
 hist3->SetXTitle("#theta_{#mu} (degrees)");
 hist3->SetYTitle("corrected events");
 hist3->SetMaximum(80);
  hist3->SetMinimum(0);
 hist3->SetLineColor(1);
 hist3->DrawNormalized("hist", hist3->Integral());
 hist3up->SetLineColor(2);
 hist3up->DrawNormalized("SAME,hist", hist3up->Integral());
 hist3down->SetLineColor(3);
 hist3down->DrawNormalized("SAME,hist", hist3down->Integral());
 lt->Draw();
 ltt->Draw();
 lttt->Draw();
 gPad->GetCanvas()->Print("cossmear_meas.pdf");
 
new TCanvas;
 
 hist33->SetTitle("Muon #theta_{#mu} (w/ backgorund???) adataevent2");
 hist33->SetXTitle("#theta_{#mu} (degrees)");
 hist33->SetYTitle("events");
 
 hist33->SetLineColor(1);
 hist33->DrawNormalized("", hist33->Integral());
 hist33up->SetLineColor(2);
 hist33up->DrawNormalized("SAME", hist33up->Integral());
 hist33down->SetLineColor(3);
 hist33down->DrawNormalized("SAME", hist33down->Integral());

new TCanvas;
 
 hist333->SetTitle("Muon #theta_{#mu} (w/ background???) adataevent3");
 hist333->SetXTitle("#theta_{#mu} (degrees)");
 hist333->SetYTitle("events");

 hist333->SetLineColor(1);
 hist333->DrawNormalized("", hist333->Integral());
 hist333up->SetLineColor(2);
 hist333up->DrawNormalized("SAME", hist333up->Integral());
 hist333down->SetLineColor(3);
 hist333down->DrawNormalized("SAME", hist333down->Integral());

new TCanvas;
 
 hist33mue->SetTitle("Muon mom (w/ backgorund???) adataevent2");
 hist33mue->SetXTitle("mom");
 hist33mue->SetYTitle("events");
 
 hist33mue->SetLineColor(1);
 hist33mue->DrawNormalized("", hist33mue->Integral());
 hist33mueup->SetLineColor(2);
 hist33mueup->DrawNormalized("SAME", hist33mueup->Integral());
 hist33muedown->SetLineColor(3);
 hist33muedown->DrawNormalized("SAME", hist33muedown->Integral());

new TCanvas;
 
 hist333mue->SetTitle("Muon mom (w/ background???) adataevent3");
 hist333mue->SetXTitle("mom");
 hist333mue->SetYTitle("events");

 hist333mue->SetLineColor(1);
 hist333mue->DrawNormalized("", hist333mue->Integral());
 hist333mueup->SetLineColor(2);
 hist333mueup->DrawNormalized("SAME", hist333mueup->Integral());
 hist333muedown->SetLineColor(3);
 hist333muedown->DrawNormalized("SAME", hist333muedown->Integral());


new TCanvas;

 hist3mue->SetMaximum(215);
 hist3mue->SetTitle("CC #nu_{#mu} muon momentum (measurement scale error)");
 hist3mue->SetXTitle("P_{#mu} (GeV/c)");
 hist3mue->SetYTitle("corrected events");

 hist3mue->SetLineColor(1);
 hist3mue->DrawNormalized("hist", hist3mue->Integral());
 hist3mueup->SetLineColor(2);
 hist3mueup->DrawNormalized("SAME,hist", hist3mueup->Integral());
 hist3muedown->SetLineColor(3);
 hist3muedown->DrawNormalized("SAME,hist", hist3muedown->Integral());
  lt->Draw();
 ltt->Draw();
 lttt->Draw();
gPad->GetCanvas()->Print("momsmear_meas.pdf");
new TCanvas;
 
 hist8mue->SetTitle("CC #nu_{#mu} muon momentum (energy lost error)");
 hist8mue->SetXTitle("P_{#mu} (GeV/c)");
 hist8mue->SetYTitle("corrected events");

 hist8mue->SetLineColor(1);
 hist8mue->DrawNormalized("", hist8mue->Integral()*potfactor);
 hist8mueup->SetLineColor(2);
 hist8mueup->DrawNormalized("SAME", hist8mueup->Integral()*potfactor);
 hist8muedown->SetLineColor(3);
 hist8muedown->DrawNormalized("SAME", hist8muedown->Integral()*potfactor);
   lt->Draw();
 ltt->Draw();
 lttt->Draw();
gPad->GetCanvas()->Print("elostsmear_meas.pdf");




new TCanvas;
gStyle->SetTitleX(.3);
TH1D *histeff1111j=(TH1D*)histeff1111->Clone();
TH1D *hist11mueeffcorrectedj=(TH1D*)hist11mueeffcorrected->Clone();

histeff1111j->Scale(1./denom);
hist11mueeffcorrectedj->Scale(1/denom);

histeff1111j->Scale(1./1.25);//bin width 
hist11mueeffcorrectedj->Scale(1./1.25);//bin width 

histeff1111j->SetTitle("#nu_{#mu} CC d#sigma/dP_{#mu} on Ar");
histeff1111j->SetYTitle("d#sigma/dP_{#mu} [cm^{2}/(GeV/c)]");
histeff1111j->SetMinimum(0);
// histeff1111j->SetMaximum(80*pow(10,-38));
histeff1111j->SetMaximum(32*pow(10,-38));
histeff1111j->GetYaxis()->SetTitleOffset(.9);
histeff1111j->GetXaxis()->SetTitleOffset(.9);
histeff1111j->GetYaxis()->SetTitleSize(.05);
histeff1111j->GetXaxis()->SetTitleSize(.05);
histeff1111j->DrawCopy("hist");
hist11mueeffcorrectedj->DrawCopy("SAME, E1"); 
latexxxxxxxx->Draw();
leg2->Draw();
// gPad->GetCanvas()->Print("momArdiff.pdf");

std::cout<<"KS test "<<histeff1111j->KolmogorovTest(hist11mueeffcorrectedj,"N")<<std::endl;

for(int i=1;i<=hist11mueeffcorrectedj->GetNbinsX();i++)
std::cout<<hist11mueeffcorrectedj->GetBinCenter(i)<<" "<<hist11mueeffcorrectedj->GetBinContent(i)<<" "<<hist11mueeffcorrectedj->GetBinError(i)<<" "<<histeff1111j->GetBinContent(i)<<std::endl;


new TCanvas;

histeff1111j->DrawCopy("hist");
mom1->SetLineColor(1);
mom1->SetFillColor(1);
mom1->SetFillStyle(0);
mom1->DrawNormalized("E1,SAME",hist11mueeffcorrectedj->Integral());
mom3->DrawNormalized("E1,SAME",hist11mueeffcorrectedj->Integral());
latexxxxxxxx->Draw();
leg2->Draw();

gPad->GetCanvas()->Print("momArdiff.pdf");

for(int i=1;i<=hist11mueeffcorrectedj->GetNbinsX();i++)
std::cout<<"mom Ar "<<hist11mueeffcorrectedj->GetBinCenter(i)<<" "<<hist11mueeffcorrectedj->GetBinContent(i)<<" "<<hist11mueeffcorrectedj->GetBinError(i)<<std::endl;

new TCanvas;

TH1D *histeff111j=(TH1D*)histeff111->Clone();
TH1D *hist11effcorrectedj=(TH1D*)hist11effcorrected->Clone();

histeff111j->Scale(1./denom);
hist11effcorrectedj->Scale(1./denom);

histeff111j->Scale(1./2.0);//bin width 
hist11effcorrectedj->Scale(1./2.0);//bin width 

histeff111j->SetTitle("on Ar");
histeff111j->SetMinimum(0);
histeff111j->SetTitle("#nu_{#mu} CC d#sigma/d#theta_{#mu} on Ar");
histeff111j->SetYTitle("d#sigma/d#theta_{#mu} (cm^{2}/degree)");
histeff111j->GetYaxis()->SetTitleOffset(.9);
histeff111j->GetXaxis()->SetTitleOffset(.9);
histeff111j->GetYaxis()->SetTitleSize(.05);
histeff111j->GetXaxis()->SetTitleSize(.05);
// histeff111j->SetMaximum(.6*pow(10,-36));
histeff111j->SetMaximum(.09*pow(10,-36));
histeff111j->DrawCopy("hist");
latexxxxxxxx->Draw();
leg2->Draw();
hist11effcorrectedj->DrawCopy("SAME, E1");

for(int i=1;i<=hist11effcorrectedj->GetNbinsX();i++)
std::cout<<"theta Ar "<<hist11effcorrectedj->GetBinCenter(i)<<" "<<hist11effcorrectedj->GetBinContent(i)<<" "<<hist11effcorrectedj->GetBinError(i)<<std::endl;


std::cout<<"KS test "<<histeff111j->KolmogorovTest(hist11effcorrectedj,"N")<<std::endl;


for(int i=1;i<=hist11effcorrectedj->GetNbinsX();i++)
std::cout<<hist11effcorrectedj->GetBinCenter(i)<<" "<<hist11effcorrectedj->GetBinContent(i)<<" "<<hist11effcorrectedj->GetBinError(i)<<" "<<histeff111j->GetBinContent(i)<<std::endl;

// gPad->GetCanvas()->Print("thetaArdiff.pdf");


new TCanvas;

histeff111j->DrawCopy("hist");
cos1->SetLineColor(1);
cos1->SetFillColor(1);
 cos1->SetFillStyle(0);
cos1->DrawNormalized("E1,SAME",hist11effcorrectedj->Integral());
cos3->DrawNormalized("E1,SAME",hist11effcorrectedj->Integral());
latexxxxxxxx->Draw();
leg2->Draw();
gPad->GetCanvas()->Print("thetaArdiff.pdf");


new TCanvas;

TH1D *histeff1111jj=(TH1D*)histeff1111->Clone();
TH1D *hist11mueeffcorrectedjj=(TH1D*)hist11mueeffcorrected->Clone();

histeff1111jj->Scale(1./denom);
hist11mueeffcorrectedjj->Scale(1/denom);

histeff1111jj->Scale(1./1.25);//bin width 
hist11mueeffcorrectedjj->Scale(1./1.25);//bin width 

histeff1111jj->Scale(1./39.948);//per nucleon 
hist11mueeffcorrectedjj->Scale(1./39.948);//per nucleon 

// histeff1111jj->Multiply(momiso);
// hist11mueeffcorrectedjj->Multiply(momiso);


for(int i=1;i<=histeff1111jj->GetNbinsX();i++)
{
histeff1111jj->SetBinContent(i,histeff1111jj->GetBinContent(i)*((histeff1111jj->GetBinCenter(i)*(-.000014))+.958));

double frac=hist11mueeffcorrectedjj->GetBinError(i)/hist11mueeffcorrectedjj->GetBinContent(i);

hist11mueeffcorrectedjj->SetBinContent(i,hist11mueeffcorrectedjj->GetBinContent(i)*((hist11mueeffcorrectedjj->GetBinCenter(i)*(-.000014))+.958));

hist11mueeffcorrectedjj->SetBinError(i,frac*hist11mueeffcorrectedjj->GetBinContent(i));


}


histeff1111jj->SetTitle("CC #nu_{#mu} d#sigma/dP_{#mu} on isoscalar");
histeff1111jj->SetYTitle("d#sigma/dP_{#mu} [cm^{2}/(GeV/c)]");
histeff1111jj->GetYaxis()->SetTitleOffset(1.2);
histeff1111jj->SetMinimum(0);
// histeff1111jj->SetMaximum(2.*pow(10,-38));
histeff1111jj->SetMaximum(0.77*pow(10,-38));
histeff1111jj->DrawCopy("hist");
hist11mueeffcorrectedjj->DrawCopy("SAME, E1");

for(int i=1;i<=hist11mueeffcorrectedjj->GetNbinsX();i++)
std::cout<<"mom iso "<<hist11mueeffcorrectedjj->GetBinCenter(i)<<" "<<hist11mueeffcorrectedjj->GetBinContent(i)<<" "<<hist11mueeffcorrectedjj->GetBinError(i)<<std::endl;
latexxxxxxxx->Draw();
leg3->Draw();
gPad->GetCanvas()->Print("momisodiff.pdf");

new TCanvas;

TH1D *histeff111jj=(TH1D*)histeff111->Clone();
TH1D *hist11effcorrectedjj=(TH1D*)hist11effcorrected->Clone();

histeff111jj->Scale(1./denom);
hist11effcorrectedjj->Scale(1./denom);

histeff111jj->Scale(1./2.0);//bin width 
hist11effcorrectedjj->Scale(1./2.0);//bin width 

histeff111jj->Scale(1./39.948);//per nucleon 
hist11effcorrectedjj->Scale(1./39.948);//per nucleon 

// histeff111jj->Multiply(cosiso);
// hist11effcorrectedjj->Multiply(cosiso);

for(int i=1;i<=histeff1111jj->GetNbinsX();i++)
{
histeff111jj->SetBinContent(i,histeff111jj->GetBinContent(i)*((histeff111jj->GetBinCenter(i)*(-.00023))+.962));


double frac=hist11effcorrectedjj->GetBinError(i)/hist11effcorrectedjj->GetBinContent(i);

hist11effcorrectedjj->SetBinContent(i,hist11effcorrectedjj->GetBinContent(i)*((hist11effcorrectedjj->GetBinCenter(i)*(-.00023))+.962));

hist11effcorrectedjj->SetBinError(i,frac*hist11effcorrectedjj->GetBinContent(i));
}

histeff111jj->SetTitle("on iso");

histeff111jj->SetMinimum(0);
histeff111jj->SetTitle("CC #nu_{#mu} d#sigma/d#theta_{#mu} on isoscalar");

histeff111jj->SetYTitle("d#sigma/d#theta_{#mu} (cm^{2}/degree)");
// gStyle->SetTitleW(0.5);
gStyle->SetTitleX(.3);
histeff111jj->GetYaxis()->SetTitleOffset(1.2);
// histeff111jj->SetMaximum(14.*pow(10,-39));
histeff111jj->SetMaximum(2.25*pow(10,-39));
histeff111jj->DrawCopy("hist");
hist11effcorrectedjj->DrawCopy("SAME, E1");


for(int i=1;i<=hist11effcorrectedjj->GetNbinsX();i++)
std::cout<<"theta iso "<<hist11effcorrectedjj->GetBinCenter(i)<<" "<<hist11effcorrectedjj->GetBinContent(i)<<" "<<hist11effcorrectedjj->GetBinError(i)<<std::endl;

latexxxxxxxx->Draw();
leg3->Draw();
gPad->GetCanvas()->Print("thetaisodiff.pdf");
// gStyle->SetEndErrorSize(10);
// gStyle->SetErrorX(0.5);





// new TCanvas;
// THStack hs("hs","test stacked histograms");
// histerror00->SetFillColor(1);
// hs.Add(histerror00);
// histerror11->SetFillColor(8);
// hs.Add(histerror11);
// histerror22->SetFillColor(3);
// hs.Add(histerror22);
// histerror33->SetFillColor(4);
// hs.Add(histerror33);
// histerror44->SetFillColor(6);
// hs.Add(histerror44);
// histerror55->SetFillColor(7);
// hs.Add(histerror55);
// histerror77->SetFillColor(9);
// hs.Add(histerror77);
// hs->Draw();

new TCanvas;
histeff111->SetLineColor(2);
histeff111->Draw();
hist3->SetLineColor(1);
hist3->Draw("SAME");




}

void project(double lardirectionEnd[3], double larEnd[3], double minosdirectionStart[3], double minosStart[3], double &xpred, double &ypred, double &rdiff, double &totaldiff, double &thetadiff)
{
    double D=(90*0.5)+(42.4*2.54)-5.588; //distance from the front (upstream) of the TPC to the 1st Minos plane 
      double x_offset=117.4; // previously 116.9;
      double y_offset=19.3; //previously  20.28;
      double dz = D - larEnd[2]+(100.0 * minosStart[2]);//z-difference between end of T962 track and
                                                                  //begin of MINOS track...in centimeters
      double l = dz/(lardirectionEnd[2]);//3-d distance between end of T962 track and begin of MINOS track
      double x_pred = l*lardirectionEnd[0]+larEnd[0];//predicted x-pos. of T962 track at z-position equal to                                                       //start of MINOS track
      double y_pred = l*lardirectionEnd[1]+larEnd[1];//predicted y-pos. of T962 track at z-position equal to                                                       //start of MINOS track
      
      double dx = 100.0*minosStart[0]- x_offset - x_pred;
      double dy = 100.0*minosStart[1] + y_offset - y_pred;
      

      xpred = dx;
      ypred = dy;
      rdiff = sqrt(dx*dx + dy*dy);
      //totaldiff is a measure of the agreement between the ArgoNeuT projected track and the MINOS track based on radial distance and angle. totaldiff= rdiff/cos(theta)  
      totaldiff=fabs(rdiff/((lardirectionEnd[0]*minosdirectionStart[0])+(lardirectionEnd[1]*minosdirectionStart[1])+(lardirectionEnd[2]*minosdirectionStart[2])));
      thetadiff =TMath::ACos((lardirectionEnd[0]*minosdirectionStart[0])+(lardirectionEnd[1]*minosdirectionStart[1])+(lardirectionEnd[2]*minosdirectionStart[2]));
      
return;
}


void
set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}