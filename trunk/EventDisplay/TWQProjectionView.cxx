///
/// \file    TWQProjectionView.cxx
/// \brief   The "main" event display view that most people will want to use
/// \author  brebel@fnal.gov
/// \version $Id: XZYZProjectionsView.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
#include "TCanvas.h"
#include "TGFrame.h" // For TGMainFrame, TGHorizontalFrame
#include "TGLayout.h" // For TGLayoutHints
#include "TGButton.h" // For TGCheckButton
#include "TGDimension.h"
#include "TGNumberEntry.h"
#include "TGLabel.h"
#include "TMath.h"
#include "TString.h"
#include "TRootEmbeddedCanvas.h"
#include "TLine.h"
#include "Buttons.h"
#include "TROOT.h"

#include "EventDisplay/TWQProjectionView.h"
#include "EventDisplay/HeaderPad.h"
#include "EventDisplay/ColorDrawingOptions.h"
#include "EventDisplay/EvdLayoutOptions.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "EventDisplay/SimulationDrawingOptions.h"
#include "EventDisplay/TWireProjPad.h"
#include "EventDisplay/TQPad.h"
#include "EventDisplay/MCBriefPad.h"
#include "EventDisplay/RawDataDrawer.h"
#include "EventDisplay/RecoBaseDrawer.h"
#include "Geometry/geo.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace evd{
  
  static unsigned int kPlane;
  static unsigned int kWire;
 
  
  TWQProjectionView::TWQProjectionView(TGMainFrame* mf) : 
    evdb::Canvas(mf)
  {  

    art::ServiceHandle<geo::Geometry> geo;

    // first make pads for things that don't depend on the number of 
    // planes in the detector
    // bottom left corner is (0.,0.), top right is  (1., 1.)

    evdb::Canvas::fCanvas->cd();  
    fHeaderPad = new HeaderPad("fHeaderPad","Header",0.0,0.0,0.15,0.13,"");  
    fHeaderPad->Draw();  
  
    evdb::Canvas::fCanvas->cd();  
    fMC = new MCBriefPad("fMCPad","MC Info.",0.15,0.13,1.0,0.17,"");  
    fMC->Draw();  

    evdb::Canvas::fCanvas->cd();  
    fWireQ = new TQPad("fWireQPad", "ADCvsTime",0.15,0.0,1.0,0.13,"TQ", 0, 0);  
    fWireQ->Draw();  

    
     // add new "meta frame" to hold the GUI Canvas and a side frame (vframe)
    fMetaFrame  = new TGCompositeFrame(mf, 60, 60, kHorizontalFrame);
   
   //new frame organizing the buttons on the left of the canvas.
    fVFrame  = new TGCompositeFrame(fMetaFrame, 60, 60, kVerticalFrame); 
     // Define a layout for placing the canvas within the frame.
    fLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX |
                              kLHintsExpandY, 5, 5, 5, 5);
    
    mf->RemoveFrame((TGFrame *)fEmbCanvas);
    mf->RemoveFrame(fFrame);
    
    fEmbCanvas->ReparentWindow( fMetaFrame, fXsize, fYsize);
    	  
    
     fMetaFrame->AddFrame(fVFrame,new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandY));
     fMetaFrame->AddFrame(fEmbCanvas, fLayout);
  
    
     mf->AddFrame(fMetaFrame,fLayout);
     mf->AddFrame(fFrame);
    
    // plane number entry  
    fPlaneEntry = new TGNumberEntry(fFrame,
				    0,2,-1,
				    TGNumberFormat::kNESInteger, 
				    TGNumberFormat::kNEAAnyNumber, 
				    TGNumberFormat::kNELLimitMinMax, 
				    0 , geo->Nplanes()-1 );

    kPlane = 0;
    kWire = TMath::Nint(0.5*geo->Nwires(0));
    fWireQ->SetPlaneWire(kPlane, kWire);

    // Initial value
    fPlaneEntry->SetNumber( kPlane );

    // There are two "signals" to which a TGNumberEntry may respond:
    // when the user clicks on the arrows, or when the user types in a
    // new number in the text field.
    fPlaneEntry->Connect("ValueSet(Long_t)", "evd::TWQProjectionView", this, "SetPlane()");
    fPlaneEntry->GetNumberEntry()->Connect("ReturnPressed()", "evd::TWQProjectionView", this, "SetPlane()");
    // Text label for this numeric field.
    fPlaneLabel= new TGLabel(fFrame,"Plane");
 
    // wire number entry 
    fWireEntry = new TGNumberEntry(fFrame,0,6,-1,
				   TGNumberFormat::kNESInteger, 
				   TGNumberFormat::kNEAAnyNumber, 
				   TGNumberFormat::kNELLimitMinMax, 
				   0 , geo->Nwires(0)-1 );
    // Initial value
    fWireEntry->SetNumber( kWire );

    // There are two "signals" to which a TGNumberEntry may respond:
    // when the user clicks on the arrows, or when the user types in a
    // new number in the text field.
    fWireEntry->Connect("ValueSet(Long_t)", "evd::TWQProjectionView", this, "SetWire()");
    fWireEntry->GetNumberEntry()->Connect("ReturnPressed()", "evd::TWQProjectionView", this, "SetWire()");

    // Text label for this numeric field.
    fWireLabel= new TGLabel(fFrame,"Wire");

    // adc threshold number entry 
    fThresEntry = new TGNumberEntry(fFrame,0,6,-1,
				   TGNumberFormat::kNESInteger, 
				   TGNumberFormat::kNEAAnyNumber, 
				   TGNumberFormat::kNELLimitMinMax, 
				   0 , geo->Nwires(0)-1 );
    // Initial value
    art::ServiceHandle<evd::ColorDrawingOptions>      cst;
    art::ServiceHandle<evd::SimulationDrawingOptions> sdo;
    art::ServiceHandle<evd::RawDrawingOptions>        rawopt;
    art::ServiceHandle<evd::EvdLayoutOptions>        evdlayoutopt;
   
    
    fThresEntry->SetNumber( rawopt->fMinSignal );

    // There are two "signals" to which a TGNumberEntry may respond:
    // when the user clicks on the arrows, or when the user types in a
    // new number in the text field.
    fThresEntry->Connect("ValueSet(Long_t)", "evd::TWQProjectionView", this, "SetThreshold()");
    fThresEntry->GetNumberEntry()->Connect("ReturnPressed()", "evd::TWQProjectionView", this, "SetThreshold()");

    // Text label for this numeric field.
    fThresLabel= new TGLabel(fFrame,"ADC Threshold");

    // check button to toggle color vs grey
    fGreyScale = new TGCheckButton(fFrame,"Grayscale",1);
    fGreyScale->Connect("Clicked()", "evd::TWQProjectionView", this, "SetGreyscale()");
    if(cst->fColorOrGray == 1) fGreyScale->SetState(kButtonDown);

    // check button to toggle MC information
    fMCOn = new TGCheckButton(fFrame,"MC Truth",5);
    fMCOn->Connect("Clicked()", "evd::TWQProjectionView", this, "SetMCInfo()");
    if(sdo->fShowMCTruthText == 1) fMCOn->SetState(kButtonDown);

    // radio buttons to toggle drawing raw vs calibrated information
    fRawCalibDraw = new TGRadioButton(fFrame,"Both",          2);
    fCalibDraw    = new TGRadioButton(fFrame,"Reconstructed", 3);
    fRawDraw      = new TGRadioButton(fFrame,"Raw",           4);
    fRawDraw     ->Connect("Clicked()", "evd::TWQProjectionView", this, "SetRawCalib()");
    fCalibDraw   ->Connect("Clicked()", "evd::TWQProjectionView", this, "SetRawCalib()");
    fRawCalibDraw->Connect("Clicked()", "evd::TWQProjectionView", this, "SetRawCalib()");
    if(rawopt->fDrawRawDataOrCalibWires == 0)      fRawDraw->SetState(kButtonDown);
    else if(rawopt->fDrawRawDataOrCalibWires == 1) fCalibDraw->SetState(kButtonDown);
    else if(rawopt->fDrawRawDataOrCalibWires == 2) fRawCalibDraw->SetState(kButtonDown);

    // Put all these widgets into the frame.  The last
    // four numbers in each TGLayoutHint are padleft, padright,
    // padtop, padbottom.

    fFrame->AddFrame(fMCOn,          new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
    fFrame->AddFrame(fGreyScale,     new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
    fFrame->AddFrame(fRawCalibDraw,  new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
    fFrame->AddFrame(fCalibDraw,     new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
    fFrame->AddFrame(fRawDraw,       new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
    fFrame->AddFrame(fPlaneEntry,    new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 2, 1 ) );
    fFrame->AddFrame(fPlaneLabel,    new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
    fFrame->AddFrame(fWireEntry,     new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 2, 1 ) );
    fFrame->AddFrame(fWireLabel,     new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );
    fFrame->AddFrame(fThresEntry,    new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 2, 1 ) );
    fFrame->AddFrame(fThresLabel,    new TGLayoutHints(kLHintsBottom | kLHintsRight, 0,  0, 5, 1 ) );

    // geometry to figure out the number of planes
    unsigned int nplanes = geo->Nplanes();

    if(evdlayoutopt-> fShowSideBar)
	SetUpSideBar();    
    
    // now determine the positions of all the time vs wire number 
    // and charge histograms for the planes
    for(unsigned int i = 0; i < nplanes; ++i){
      double twx1 = 0.;
      double twx2 = 0.97;
      double twx3 = 1.0;
      double twy1 = 0.17 +   (i)*(1.0-0.171)/(1.*nplanes);
      double twy2 = 0.17 + (i+1)*(1.0-0.171)/(1.*nplanes);
      
     
      
      TString padname = "fWireProjP";
      padname += i;

      TString padtitle = "Plane";
      padtitle += i;

      evdb::Canvas::fCanvas->cd();
      fPlanes.push_back(new TWireProjPad(padname, padtitle, twx1, twy1, twx2, twy2, i));
      fPlanes[i]->Draw();

      fPlanes[i]->Pad()->AddExec("selectwire",Form("evd::TWQProjectionView::ChangeWire(%d, (void*)%d)", i, this));
      
      fPlanes[i]->Pad()->AddExec("getmousezoom",Form("evd::TWQProjectionView::SetMouseZoomRegion(%d, (void*)%d)", i,this));
     
      
      padname = "fQPadPlane";
      padname += i;

      padtitle = "QPlane";
      padtitle += i;

      evdb::Canvas::fCanvas->cd();
      fPlaneQ.push_back(new TQPad(padname, padtitle, twx2, twy1, twx3, twy2, "Q", i, 0));
      fPlaneQ[i]->Draw();
    
    } 


    
    evdb::Canvas::fCanvas->Update();

  }

  //......................................................................
  TWQProjectionView::~TWQProjectionView() 
  {  
    if (fHeaderPad) { delete fHeaderPad;  fHeaderPad  = 0; }  
    if (fMC)        { delete fMC;         fMC         = 0; }  
    if (fWireQ)     { delete fWireQ;      fWireQ      = 0; }  
    if (fPlaneEntry){ delete fPlaneEntry; fPlaneEntry = 0; }
    if (fWireEntry) { delete fWireEntry;  fWireEntry  = 0; }
    if (fPlaneLabel){ delete fPlaneLabel; fPlaneLabel = 0; }
    if (fWireLabel) { delete fWireLabel;  fWireLabel  = 0; }
    for(unsigned int i = 0; i < fPlanes.size(); ++i){
      if(fPlanes[i]){ delete fPlanes[i];  fPlanes[i]  = 0; }
      if(fPlaneQ[i]){ delete fPlaneQ[i];  fPlaneQ[i]  = 0; }
    }
    fPlanes.clear();
    fPlaneQ.clear();
  }

//......................................................................
  void TWQProjectionView::Draw(const char* opt) 
  {  
    evdb::Canvas::fCanvas->cd();    

    fHeaderPad->Draw();    
    fMC       ->Draw();  
    fWireQ    ->Draw();
    
    for(unsigned int i = 0; i < fPlanes.size(); ++i){
      fPlanes[i]->Draw(opt);
      fPlaneQ[i]->Draw();
      }

     art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutopt;
     if(evdlayoutopt->fPrintTotalCharge)
        PrintCharge();
	 
    evdb::Canvas::fCanvas->Update();
  }

// comment out this method as for now we don't want to change every
  // plane to have the same range in wire number because wire numbers
  // don't necessarily overlap from plane to plane, ie the same range
  // isn't appropriate for every plane
  //......................................................................
//   void TWQProjectionView::RangeChanged() 
//   {  
//     static int ilolast = -1;  
//     static int ihilast = -1;    
// 
//     int ilo; 
//     int ihi;
//     std::vector<int> lo;
//     std::vector<int> hi;
//     std::vector<bool> axischanged;
//     for(unsigned int i = 0; i < fPlanes.size(); ++i){
//       fPlanes[i]->GetWireRange(&ilo, &ihi);  
//       lo.push_back(ilo);
//       hi.push_back(ihi);
//       axischanged.push_back((ilo != ilolast) || (ihi != ihilast));
//     }
//       
//     TVirtualPad* ori = gPad;  
// 
//     // loop over the bools to see which axes need to change
//     for(unsigned int i = 0; i < axischanged.size(); ++i){
//       if (axischanged[i]) {    
// 	fPlanes[i]->SetWireRange(ilo, ihi);    
// 	fPlanes[i]->Pad()->cd();    
// 	fPlanes[i]->Pad()->Modified();    
// 	fPlanes[i]->Pad()->Update();    
// 
// 	ilolast = ilo;    
// 	ihilast = ihi;  
//       }  
//     }
// 
//     evdb::Canvas::fCanvas->cd();  
//     evdb::Canvas::fCanvas->Modified();  
//     evdb::Canvas::fCanvas->Update();  
//     ori->cd();
//   }
//......................................................................
  void 	TWQProjectionView::PrintCharge()
  {
    
     art::ServiceHandle<geo::Geometry> geo;
     art::ServiceHandle<evd::RawDrawingOptions> rawopt;
     
    for(unsigned int iplane=0;iplane<fPlanes.size(); iplane++)
    { 
      if(geo->Plane(iplane).SignalType()==geo::kCollection)
	  {
	 double ch=0,convch=0;
         if(rawopt->fDrawRawDataOrCalibWires == 0)
	 {fPlanes[iplane]->RawDataDraw()->GetChargeSum(iplane,ch,convch);
	  std::cout << "Warning! Calculating for RawData! " << std::endl;
	 }
         else 
	 {  
	 fPlanes[iplane]->RecoBaseDraw()->GetChargeSum(iplane,ch,convch);  
	 }    
      
	 std::cout << std::endl << "charge collected at collection plane: " << iplane << " " << ch << " " << convch << std::endl;
	  }
     }


  }
  
  
  //......................................................................
  void TWQProjectionView::ChangeWire(int plane, void * wqpv)
  {
    //initial check for a mouse click on a TBox object
    int event = gPad->GetEvent();
    int px = gPad->GetEventX();
    if(event!=11) return;
    TObject *select = gPad->GetSelected();
    if(!select) return;
    if(!select->InheritsFrom("TBox")) return;
  
    //now find wire that was clicked on
    float xx = gPad->AbsPixeltoX(px);
    float x = gPad->PadtoX(xx);
          

    kPlane = plane;
    kWire  = (unsigned int)TMath::Nint(x);

    
    evd::TWQProjectionView *wqpp = (evd::TWQProjectionView*)wqpv;
    wqpp->SetPlaneWire();

    return;
   
  }
//.......................................................................
 void TWQProjectionView::SetMouseZoomRegion(int plane,void *wqpv)
{
//*-*-*-*-*-*-*-*-*-*-*Create a new arrow in this pad*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*                  ==============================
//
  TObject *select = gPad->GetSelected();
  if(!select) return;
  if(!select->InheritsFrom("TBox")) return;
  
  static Float_t x0=-1, y0=-1, x1=-1, y1=-1;
  
  static Int_t pxold, pyold;
  static Int_t px0, py0;
  static Int_t linedrawn;
  //TLine *line;
  
  static int wstart,wend;
  static float tstart,tend;
  
  int event = gPad->GetEvent();
  int px = gPad->GetEventX();
  int py = gPad->GetEventY();
  
  switch (event){
    
  case kButton1Down:{
     gVirtualX->SetLineColor(-1);
     x0 = gPad->AbsPixeltoX(px);
     y0 = gPad->AbsPixeltoY(py);
     px0   = px; py0   = py;
     pxold = px; pyold = py;
     linedrawn = 0;
     float x = gPad->PadtoX(x0);
     tstart = gPad->PadtoY(y0);
     
     wstart  = (unsigned int)TMath::Nint(x);
     
     break;
   }
   case kButton1Motion:{ 
     int lx,hx,ly,hy;
     if (px0 < pxold){
       lx=px0; 
       hx=pxold; 
     }
     else{
       lx=pxold;
       hx=px0;
     }
     
     if (py0 < pyold){
       ly=py0; 
       hy=pyold;
     }
      else{
	ly=pyold;
	hy=py0;
      }
      
     if (linedrawn) gVirtualX->DrawBox(lx, ly, hx, hy,TVirtualX::kHollow);
     pxold = px;
     pyold = py;
     linedrawn = 1;
     
     if (px0 < pxold){
       lx=px0; 
       hx=pxold; 
     }
     else{
       lx=pxold;
       hx=px0;
     }
     
     if (py0 < pyold){
       ly=py0; 
       hy=pyold;
     }
     else{
       ly=pyold;
       hy=py0;
     }
      
     gVirtualX->DrawBox(lx, ly, hx, hy,TVirtualX::kHollow);
     break;
   }
   case kButton1Up:{
     if (px == px0 && py == py0) break;
     x1 = gPad->AbsPixeltoX(px);
     y1 = gPad->AbsPixeltoY(py);
     gPad->Modified(kTRUE);
      
     //   line = new TLine(x0,y0,x1,y1);
     //   line->Draw();
     
     float x = gPad->PadtoX(x1);
     tend = gPad->PadtoY(y1);
     wend  = (unsigned int)TMath::Nint(x);
     
     gROOT->SetEditorMode();
     if(wstart != -1 && tstart != -1){
       
       evd::TWQProjectionView *wqpp = (evd::TWQProjectionView*)wqpv;
       wqpp->SetZoom(plane,wstart,wend,tstart,tend);
       wstart=-1;
       tstart=-1;
     }
     break;
   }
  }// end switch
}
  
  
//......................................................................
// if flag is true then zoom. If flag is false then unzoom.

 void 	TWQProjectionView::ZoomInterest(bool flag)
 {
 art::ServiceHandle<geo::Geometry> geo;
 art::ServiceHandle<evd::RawDrawingOptions> rawopt; 
 
 for(unsigned int iplane=0;iplane<fPlanes.size();iplane++)
 {
   int minw,maxw,mint,maxt;
   if(flag)
   {int test=0;
    if(rawopt->fDrawRawDataOrCalibWires == 0)
      fPlanes[iplane]->RawDataDraw()->GetRegionOfInterest(iplane,minw,maxw,mint,maxt);
    else
      fPlanes[iplane]->RecoBaseDraw()->GetRegionOfInterest(iplane,minw,maxw,mint,maxt);
      
   if(test==-1)
      continue;}
   else
   {
    minw = -0.005*(geo->Nwires(iplane)-1);
    maxw =  1.005*(geo->Nwires(iplane)-1);
    mint = -0.005*fPlanes[iplane]->RawDataDraw()->TotalClockTicks();
    maxt =  1.01*fPlanes[iplane]->RawDataDraw()->TotalClockTicks();
   }
   
  SetZoom(iplane,minw,maxw,mint,maxt);
  

  }
 
 }
 
 
 //......................................................................
  void TWQProjectionView::SetUpSideBar()
  {  
  SetUpZoomButtons();
 
  }
 
 
 
 //......................................................................
  void TWQProjectionView::SetZoomInterest()
  {  
  art::ServiceHandle<evd::EvdLayoutOptions>   evdlayoutopt;
  evdlayoutopt->fAutoZoomInterest = fToggleAutoZoom->GetState();
  }
 
//......................................................................
  void TWQProjectionView::SetUpZoomButtons()
  {
     // enter zoom buttons
     art::ServiceHandle<evd::EvdLayoutOptions>        evdlayoutopt;  
      
    fZoomInterest=new TGTextButton(fVFrame,"&Zoom Interest",150);
    fZoomInterest->Connect("Clicked()", "evd::TWQProjectionView", this, "ZoomInterest()");
    
    fUnZoomInterest=new TGTextButton(fVFrame,"&UnZoom Interest",150);
    fUnZoomInterest->Connect("Clicked()", "evd::TWQProjectionView", this, "ZoomInterest(=false)");
    
    fToggleAutoZoom=new TGCheckButton(fVFrame,"AutoZoom",0);;       ///< Toggle the autozoom setting 
    fToggleAutoZoom->Connect("Clicked()", "evd::TWQProjectionView", this, "SetZoomInterest()");
    if(evdlayoutopt->fAutoZoomInterest == 1) fToggleAutoZoom->SetState(kButtonDown);
    
    fVFrame->AddFrame(fZoomInterest,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
    fVFrame->AddFrame(fUnZoomInterest,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
    fVFrame->AddFrame(fToggleAutoZoom,   new TGLayoutHints(kLHintsTop | kLHintsLeft, 0,  0, 5, 1 ) );
  }

 //------------------------------------
  void    TWQProjectionView::SetZoom(int plane,int wirelow,int wirehi,int timelow,int timehi)
  {

    TVirtualPad *ori = gPad;
      
    // error checking - useful for the mouse zoom.
    if(wirehi<wirelow){
      int temp=wirelow;
      wirelow=wirehi;
      wirehi=temp;
    }
   
    if(timehi<timelow){
      int temp=timelow;
      timelow=timehi;
      timehi=temp;
    }
    
    fPlanes[plane]->SetZoomRange(wirelow, wirehi,timelow,timehi);
    fPlanes[plane]->Draw("1");
    fPlanes[plane]->Pad()->cd();
    fPlanes[plane]->Pad()->Modified();
    fPlanes[plane]->Pad()->Update();
    evdb::Canvas::fCanvas->cd();
    evdb::Canvas::fCanvas->Modified();
    evdb::Canvas::fCanvas->Update();
    
    ori->cd();

    return;  
  }

  //-----------------------------------------------------------------
  void TWQProjectionView::SetPlaneWire()
  {
    TVirtualPad *ori = gPad;

    fWireQ->SetPlaneWire(kPlane, kWire);
    fWireQ->Draw();
    fWireQ->Pad()->cd();
    fWireQ->Pad()->Modified();
    fWireQ->Pad()->Update();

    fPlaneEntry->SetNumber(kPlane);
    fWireEntry->SetNumber(kWire);

    evdb::Canvas::fCanvas->cd();
    evdb::Canvas::fCanvas->Modified();
    evdb::Canvas::fCanvas->Update();

    ori->cd();
  }



  //-----------------------------------------------------------------
  void TWQProjectionView::SetPlane()
  {
    kPlane = (unsigned int)fPlaneEntry->GetNumberEntry()->GetNumber();

    this->SetPlaneWire();
  }

  //-----------------------------------------------------------------
  void TWQProjectionView::SetWire()
  {
    kWire = (unsigned int)fWireEntry->GetNumberEntry()->GetNumber();

    this->SetPlaneWire();
  }

  //-----------------------------------------------------------------
  void TWQProjectionView::SetThreshold()
  {
    double threshold = fThresEntry->GetNumberEntry()->GetNumber();
    
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    rawopt->fMinSignal = threshold;

    TVirtualPad *ori = gPad;
    for(size_t i = 0; i < fPlanes.size(); ++i){
      fPlanes[i]->Draw();
      fPlanes[i]->Pad()->cd();
      fPlanes[i]->Pad()->Modified();
      fPlanes[i]->Pad()->Update();
    }
    evdb::Canvas::fCanvas->cd();
    evdb::Canvas::fCanvas->Modified();
    evdb::Canvas::fCanvas->Update();

    ori->cd();

    return;
  }

  //-----------------------------------------------------------------
  void TWQProjectionView::SetGreyscale()
  {
    art::ServiceHandle<evd::ColorDrawingOptions> cst;

    TGButton *b = (TGButton *)gTQSender;
    if(b->GetState() == kButtonDown){
      cst->fColorOrGray = 1;
    }
    else{
      cst->fColorOrGray = 0;
    }

    TVirtualPad *ori = gPad;
    for(size_t i = 0; i < fPlanes.size(); ++i){
      fPlanes[i]->Draw();
      fPlanes[i]->Pad()->cd();
      fPlanes[i]->Pad()->Modified();
      fPlanes[i]->Pad()->Update();

      fPlaneQ[i]->Draw();
      fPlaneQ[i]->Pad()->cd();
      fPlaneQ[i]->Pad()->Modified();
      fPlaneQ[i]->Pad()->Update();
    }
    evdb::Canvas::fCanvas->cd();
    evdb::Canvas::fCanvas->Modified();
    evdb::Canvas::fCanvas->Update();

    ori->cd();

    return;
  }

  //-----------------------------------------------------------------
  void TWQProjectionView::SetRawCalib()
  {
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;

    TGButton *b = (TGButton *)gTQSender;
    int id = b->WidgetId();

    // id values are set in lines 125 - 127
    if(id == 4){
      rawopt->fDrawRawDataOrCalibWires = 0;
      fRawDraw->SetState(kButtonDown);
      fCalibDraw->SetState(kButtonUp);
      fRawCalibDraw->SetState(kButtonUp);
    }
    else if(id == 3){
      rawopt->fDrawRawDataOrCalibWires = 1;
      fRawDraw->SetState(kButtonUp);
      fCalibDraw->SetState(kButtonDown);
      fRawCalibDraw->SetState(kButtonUp);
    }
    else if(id == 2){
      rawopt->fDrawRawDataOrCalibWires = 2;
      fRawDraw->SetState(kButtonUp);
      fCalibDraw->SetState(kButtonUp);
      fRawCalibDraw->SetState(kButtonDown);
    }

    TVirtualPad *ori = gPad;
    fWireQ->Draw();
    fWireQ->Pad()->cd();
    fWireQ->Pad()->Modified();
    fWireQ->Pad()->Update();

    for(size_t i = 0; i < fPlanes.size(); ++i){
      fPlanes[i]->Draw();
      fPlanes[i]->Pad()->cd();
      fPlanes[i]->Pad()->Modified();
      fPlanes[i]->Pad()->Update();

      fPlaneQ[i]->Draw();
      fPlaneQ[i]->Pad()->cd();
      fPlaneQ[i]->Pad()->Modified();
      fPlaneQ[i]->Pad()->Update();
    }
    evdb::Canvas::fCanvas->cd();
    evdb::Canvas::fCanvas->Modified();
    evdb::Canvas::fCanvas->Update();

    ori->cd();

    return;
  }

  //-----------------------------------------------------------------
  void TWQProjectionView::SetMCInfo()
  {
    art::ServiceHandle<evd::SimulationDrawingOptions> sdo;

    TGButton *b = (TGButton *)gTQSender;
    if(b->GetState() == kButtonDown){
      sdo->fShowMCTruthText    = 1;
      sdo->fShowMCTruthVectors = 1;
    }
    else{
      sdo->fShowMCTruthText    = 0;
      sdo->fShowMCTruthVectors = 0;
    }

    TVirtualPad *ori = gPad;

    fMC->Draw();

    for(size_t i = 0; i < fPlanes.size(); ++i){
      fPlanes[i]->Draw();
      fPlanes[i]->Pad()->cd();
      fPlanes[i]->Pad()->Modified();
      fPlanes[i]->Pad()->Update();

      fPlaneQ[i]->Draw();
      fPlaneQ[i]->Pad()->cd();
      fPlaneQ[i]->Pad()->Modified();
      fPlaneQ[i]->Pad()->Update();
    }
    evdb::Canvas::fCanvas->cd();
    evdb::Canvas::fCanvas->Modified();
    evdb::Canvas::fCanvas->Update();

    ori->cd();
  }
  
}// namespace
