///
/// \file    TQPad.cxx
/// \brief   Drawing pad for time or charge histograms
/// \author  messier@indiana.edu
/// \version $Id: TQPad.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
#include "EventDisplay/TQPad.h"
#include <cassert>
#include "TBox.h"
#include "TH1F.h"
#include "TF1.h"
#include "TPad.h"
#include "EventDisplayBase/evdb.h"
#include "EventDisplayBase/EventHolder.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "EventDisplay/ColorDrawingOptions.h"
#include "EventDisplay/RawDataDrawer.h"
#include "EventDisplay/RecoBaseDrawer.h"

namespace evd{

   static const int kRAW      = 0;
   static const int kCALIB    = 1;
   static const int kRAWCALIB = 2;
   static const int kQ        = 0;
   static const int kTQ       = 1;

   //......................................................................

   TQPad::TQPad(const char* nm, const char* ti,
                double x1, double y1,
                double x2, double y2,
                const char *opt,
                unsigned int plane,
                unsigned int wire) :
      DrawingPad(nm, ti, x1, y1, x2, y2),
      fWire(wire),
      fPlane(plane),
      fRawHisto(0),
      fRecoHisto(0)
   {

      art::ServiceHandle<geo::Geometry> geo;
      unsigned int planes = geo->Nplanes();

      this->Pad()->cd();

      this->Pad()->SetLeftMargin  (0.050);
      this->Pad()->SetRightMargin (0.050);

      this->Pad()->SetTopMargin   (0.005);
      this->Pad()->SetBottomMargin(0.110);

      // there has to be a better way of doing this that does
      // not have a case for each number of planes in a detector
      if(planes == 2 && fPlane > 0){
         this->Pad()->SetTopMargin   (0.110);
         this->Pad()->SetBottomMargin(0.005);
      }
      else if(planes > 2){
         if(fPlane == 1){
            this->Pad()->SetTopMargin   (0.005);
            this->Pad()->SetBottomMargin(0.005);
         }
         else if(fPlane == 2){
            this->Pad()->SetTopMargin   (0.110);
            this->Pad()->SetBottomMargin(0.005);
         }
      }


      std::string opts(opt);
      if(opts == "TQ") fTQ = kTQ;
      if(opts == "Q" ){
         fTQ = kQ;
      }

      this->BookHistogram();
      fView = new evdb::View2D();
   }

   //......................................................................

   TQPad::~TQPad() 
   {
      if (fView)      { delete fView;      fView      = 0; }
      if (fRawHisto)  { delete fRawHisto;  fRawHisto  = 0; }
      if (fRecoHisto) { delete fRecoHisto; fRecoHisto = 0; }
      //if (fHitGaussian) { delete fHitGaussian;}
   }

   //......................................................................

   void TQPad::Draw() 
   {
      art::ServiceHandle<evd::RawDrawingOptions> drawopt;

      if (!fRawHisto || !fRecoHisto) return;

      //grab the singleton with the event
      const art::Event* evt = evdb::EventHolder::Instance()->GetEvent();
      if(!evt) return;

      fPad->Clear();
      fPad->cd();

      std::vector<double> hstart;
      std::vector<double> hend;
      std::vector<double> hamplitudes;
      std::vector<double> hpeaktimes;

      if(fTQ == kTQ){
         fRawHisto->Reset("ICEM");
         fRecoHisto->Reset("ICEM");

         this->RawDataDraw()->FillTQHisto(*evt,
                                          fPlane,
                                          fWire, 
                                          fRawHisto);

         this->RecoBaseDraw()->FillTQHisto(*evt,
                                           fPlane,
                                           fWire, 
                                           fRecoHisto,
                                           hstart,
                                           hend,
                                           hamplitudes,
                                           hpeaktimes);

         if(drawopt->fDrawRawDataOrCalibWires==kRAW)      fRawHisto->Draw();
         if(drawopt->fDrawRawDataOrCalibWires==kCALIB)    fRecoHisto->Draw();
         if(drawopt->fDrawRawDataOrCalibWires==kRAWCALIB){
            fRawHisto->SetMaximum(1.1*TMath::Max(fRawHisto->GetMaximum(), fRecoHisto->GetMaximum()));
            fRawHisto->SetMinimum(1.1*TMath::Min(fRawHisto->GetMinimum(), fRecoHisto->GetMinimum()));
            fRawHisto->Draw();
            fRecoHisto->Draw("same");
         }

         // this loop draws Gaussian shapes for identified hits in the reco histo
         for (size_t i=0; i<hstart.size() && drawopt->fDrawRawDataOrCalibWires!=kRAW; ++i) {
            //hend and hstart are 1-sigma away from peak
            double width = (hend[i]-hstart[i])/2;

            //create a function corresponding to the Gaussian shape
            TF1 *f1 = new TF1("hitshape","gaus(0)",hstart[i]-1.5*width,hend[i]+1.5*width);//5-sigma wide window
            f1->SetParameters(hamplitudes[i],hpeaktimes[i],width);

            //create TPolyLine that actually gets drawn
            TPolyLine& p1 = fView->AddPolyLine(1001, 
                                               kOrange+7,
                                               3,
                                               1);

            //set coordinates of TPolyLine based on Gaussian function
            for(int j = 0; j<1001; ++j){ 
               double x = hstart[i]-1.5*width+j*5*width/1000;
               double y = f1->Eval(x); 
               p1.SetPoint(j, x, y);
            }
            p1.Draw("same");
            if(f1) delete f1;
         }

         if(drawopt->fDrawRawDataOrCalibWires==kCALIB) fRecoHisto->Draw("same");
         else if(drawopt->fDrawRawDataOrCalibWires==kRAWCALIB){
            fRawHisto->Draw("same");
            fRecoHisto->Draw("same");
         }

         fRawHisto->SetTitleOffset(0.2, "Y");
         fRecoHisto->SetLabelSize(0.2, "Y");

      } // end if fTQ == kTQ
    
      if(fTQ == kQ){
         art::ServiceHandle<evd::ColorDrawingOptions> cst;

         TH1F *hist;

         int ndiv = 0;
         if(drawopt->fDrawRawDataOrCalibWires!=kCALIB){
            hist = fRawHisto;
            ndiv = cst->fRawDiv;
         }
         if(drawopt->fDrawRawDataOrCalibWires==kCALIB){
            hist = fRecoHisto;
            ndiv = cst->fRecoDiv;
         }

         hist->SetLabelSize(0, "X");
         hist->SetLabelSize(0, "Y");
         hist->SetTickLength(0, "X");
         hist->SetTickLength(0, "Y");
         hist->Draw("pY+");

         //
         // Use this to fill the histogram with colors from the color scale
         //
         double x1, x2, y1, y2;
         x1 = 0.;
         x2 = 1.;

         for(int i = 0; i < ndiv; ++i){
            y1 = hist->GetMinimum() + i*(hist->GetMaximum()-hist->GetMinimum())/(1.*ndiv);
            y2 = hist->GetMinimum() + (i + 1)*(hist->GetMaximum()-hist->GetMinimum())/(1.*ndiv);

            int c=1;
            if (drawopt->fDrawRawDataOrCalibWires==kRAW) {
               c = cst->RawQ().GetColor(0.5*(y1+y2));
            }
            if (drawopt->fDrawRawDataOrCalibWires!=kRAW) {
               c= cst->CalQ().GetColor(0.5*(y1+y2));
            }
	
            TBox& b = fView->AddBox(x1,y1,x2,y2);
            b.SetFillStyle(1001);
            b.SetFillColor(c);      
            b.Draw();
            //}
         } // end loop over Q histogram bins

         hist->Draw("same");

      } // end if fTQ == kQ

   }


   //......................................................................
  
   void TQPad::BookHistogram() 
   {
      if (fRawHisto) {
         delete fRawHisto; 
         fRawHisto = 0;
      }
      if (fRecoHisto) {
         delete fRecoHisto; 
         fRecoHisto = 0;
      }
  
      art::ServiceHandle<evd::ColorDrawingOptions> cst;
      art::ServiceHandle<evd::RawDrawingOptions> drawopt;
      art::ServiceHandle<geo::Geometry> geo;

      /// \todo decide if ndivraw and ndivreco are useful
//     int    ndivraw   = cst->fRawDiv;
//     int    ndivreco  = cst->fRecoDiv;
      double qxloraw   = cst->fRawQLow;
      double qxhiraw   = cst->fRawQHigh;
      double qxloreco  = cst->fRecoQLow;
      double qxhireco  = cst->fRecoQHigh;
      double tqxlo     = 0.;
      double tqxhi     = 1.*this->RawDataDraw()->TotalClockTicks();

      if(fTQ == kQ){
         fRawHisto = new TH1F("fRAWQHisto",
                              ";;",
                              2,0.,1.);

         fRawHisto->SetMaximum(qxhiraw);
         fRawHisto->SetMinimum(qxloraw);

         fRecoHisto = new TH1F("fCALQHisto",
                               ";;",
                               1,0.,1.);

         fRecoHisto->SetMaximum(qxhireco);
         fRecoHisto->SetMinimum(qxloreco);

      } // end if fTQ == kQ

      if(fTQ == kTQ){
         fRawHisto = new TH1F("fRAWTQHisto",
                              ";t [ticks];q [ADC]",
                              (int)tqxhi,tqxlo,tqxhi);

         fRecoHisto = new TH1F("fCALTQHisto",
                               ";t [ticks];q [ADC]",
                               (int)tqxhi,tqxlo,tqxhi);
         fRecoHisto->SetLineColor(kBlue);
      }//end if fTQ == kTQ

      // By this time I must have a histogram booked
      assert(fRecoHisto);
      assert(fRawHisto);

      fRawHisto->SetLabelSize  (0.07,"X");
      fRawHisto->SetLabelOffset(0.00,"X");
      fRawHisto->SetTitleSize  (0.10,"X");
      fRawHisto->SetTitleOffset(0.50,"X");

      fRawHisto->SetLabelSize  (0.07,"Y");
      fRawHisto->SetLabelOffset(0.01,"Y");
      fRawHisto->SetTitleSize  (0.10,"Y");
      fRawHisto->SetTitleOffset(0.30,"Y");

      fRecoHisto->SetLabelSize  (0.07,"X");
      fRecoHisto->SetLabelOffset(0.00,"X");
      fRecoHisto->SetTitleSize  (0.10,"X");
      fRecoHisto->SetTitleOffset(0.50,"X");

      fRecoHisto->SetLabelSize  (0.07,"Y");
      fRecoHisto->SetLabelOffset(0.01,"Y");
      fRecoHisto->SetTitleSize  (0.10,"Y");
      fRecoHisto->SetTitleOffset(0.30,"Y");
   }

}
//////////////////////////////////////////////////////////////////////////
