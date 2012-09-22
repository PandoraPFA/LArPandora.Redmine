#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TLine.h"
#include "TBox.h"
#include "TText.h"
#include "TPolyMarker3D.h"

#include "Geometry/Geometry.h"
#include "EventDisplayBase/View2D.h"
#include "EventDisplayBase/View3D.h"
#include "T962/ArgoneutEventDisplay/ArgoneutDrawingOptions.h"
#include "T962/ArgoneutEventDisplay/MinosDrawer.h"
#include "EventDisplay/RecoDrawingOptions.h"

#include "T962/T962_Objects/MINOS.h"
#include "RecoBase/recobase.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Utilities/AssociationUtil.h"
#include "art/Framework/Core/FindMany.h"



namespace argoevd{

  //......................................................................
  MinosDrawer::MinosDrawer()
  {
  }

  //......................................................................
  MinosDrawer::~MinosDrawer()
  {
  }

   //......................................................................
   int MinosDrawer::GetMinos(const art::Event&            evt,
                             const std::string&           which,
                             art::PtrVector<t962::MINOS>& minos) 
  {
    minos.clear();

    art::Handle< std::vector<t962::MINOS> > minoshandle;
    art::PtrVector<t962::MINOS> temp;
    
    evt.getByLabel(which, minoshandle);
    if(minoshandle.failedToGet()){
      mf::LogWarning("GetMinos") << "No MINOS information found with label = "
                                 << which << ".  Skipping.\n";
      return 0;
    }
       
    for(unsigned int i = 0; i < minoshandle->size(); ++i){
      art::Ptr<t962::MINOS> w(minoshandle, i);
      temp.push_back(w);
    }
    temp.swap(minos);
    
    return minos.size();
  }

 

  //......................................................................
   void MinosDrawer::MinosTrack(const art::Ptr<t962::MINOS> minos,
                                int                 id, 
                                evdb::View3D*       view,
                                bool                matched)
   {
      //if(!matched) return;
      int s = 1;//style
      int w = 2;//width
      int c = 5;//kColor[id%kNCOLS];
      if(minos->fcharge==-1) c = kBlue;//negative
      if(minos->fcharge==1)  c = kRed; //positive
      if(!matched) w = 1;
      
      if(matched){
         TPolyLine3D& pm = view->AddPolyLine3D(minos->ftrkstpX.size(), c, w, s);
         for(size_t s = 0; s < minos->ftrkstpX.size(); ++s){
            pm.SetPoint(s, 100.0*minos->ftrkstpX[s]-117.4, 100.0*minos->ftrkstpY[s]+19.3, 100.0*minos->ftrkstpZ[s]+147.108);
         }
      }else{
         TPolyMarker3D& pm = view->AddPolyMarker3D(minos->ftrkstpX.size(), c, 8, 3);
         for(size_t s = 0; s < minos->ftrkstpX.size(); ++s){
            pm.SetPoint(s, 100.0*minos->ftrkstpX[s]-117.4, 100.0*minos->ftrkstpY[s]+19.3, 100.0*minos->ftrkstpZ[s]+147.108);
         }
      }
      
      
      return;
   }

   //......................................................................
   void MinosDrawer::Minos3D(const art::Event& evt,
                            evdb::View3D*     view)
  {
    art::ServiceHandle<argoevd::ArgoneutDrawingOptions> argoopt;
    
     art::PtrVector<t962::MINOS> minos;
     int minossize = this->GetMinos(evt, argoopt->fMINOSLabel, minos);
     if(minossize==0) return;

     art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
     art::Handle< std::vector<recob::Track> > LarTrackHandle;
     
     for(size_t p = 0; p < minos.size(); ++p){
        bool matched = false;
       
        if(recoOpt->fDrawTracks!=0){
           for (unsigned int imod=0; imod<recoOpt->fTrackLabels.size(); ++imod){

              bool get = evt.getByLabel(recoOpt->fTrackLabels[imod],LarTrackHandle);
              if(!get) continue;

              //find matched MINOS information 
              art::FindOne<t962::MINOS> fomatch(LarTrackHandle, evt, argoopt->fMatchLabel);
              for(size_t q= 0; q < fomatch.size(); ++q)
              {
                 if(fomatch.at(q).isValid()){//Found matching MINOS track
                    if(minos[p]->ftrkIndex == fomatch.at(q).ref().ftrkIndex) matched = true;
                 }
              }  
           }
        }
        this->MinosTrack(minos[p], (int)minos[p]->ftrkIndex, view, matched);
     }
  }


   
   //......................................................................
  void MinosDrawer::DetOutline3D(evdb::View3D*        view)
  {
    art::ServiceHandle<geo::Geometry> geo;

    double xlo =  0.;
    double xhi =  2.*geo->DetHalfWidth();
    double ylo = -geo->DetHalfHeight();
    double yhi =  geo->DetHalfHeight();
    double zlo =  0.0;
    double zhi =  geo->DetLength();
  
    int c = kGray;
    int s = 1;
    int w = 1;
    TPolyLine3D& top = view->AddPolyLine3D(5, c, w, s);
    top.SetPoint(0, xlo, yhi, zlo);
    top.SetPoint(1, xhi, yhi, zlo);
    top.SetPoint(2, xhi, yhi, zhi);
    top.SetPoint(3, xlo, yhi, zhi);
    top.SetPoint(4, xlo, yhi, zlo);

    TPolyLine3D& side = view->AddPolyLine3D(5, c, w, s);
    side.SetPoint(0, xhi, yhi, zlo);
    side.SetPoint(1, xhi, ylo, zlo);
    side.SetPoint(2, xhi, ylo, zhi);
    side.SetPoint(3, xhi, yhi, zhi);
    side.SetPoint(4, xhi, yhi, zlo);

    TPolyLine3D& side2 = view->AddPolyLine3D(5, c, w, s);
    side2.SetPoint(0, xlo, yhi, zlo);
    side2.SetPoint(1, xlo, ylo, zlo);
    side2.SetPoint(2, xlo, ylo, zhi);
    side2.SetPoint(3, xlo, yhi, zhi);
    side2.SetPoint(4, xlo, yhi, zlo);

    TPolyLine3D& bottom = view->AddPolyLine3D(5, c, w, s);
    bottom.SetPoint(0, xlo, ylo, zlo);
    bottom.SetPoint(1, xhi, ylo, zlo);
    bottom.SetPoint(2, xhi, ylo, zhi);
    bottom.SetPoint(3, xlo, ylo, zhi);
    bottom.SetPoint(4, xlo, ylo, zlo);

    
    c = kGray+2;
    s = 1;
    w = 1;
  
  }

   //......................................................................
   void MinosDrawer::MinosOutline3D(evdb::View3D*        view)
   {
      Float_t x[] = 
         {   -177.292,   -308.432,   -308.432,   -305.435,   -292.456,    -280.01
             ,    -241.91,    -241.91,   -177.292,   -177.292,    177.292,    177.292
             ,     241.91,     241.91,     280.06,    297.942,    305.435,    308.432
             ,    308.432,    177.292,    177.292,   -177.292, -177.292 };
      Float_t y[] = 
         {    154.711,    23.5712,     1.1938,     1.1938,     8.6868,     8.6868
              ,    -3.7592,   -90.0938,   -154.711,   -190.602,   -190.602,   -154.711
              ,   -90.0938,    -3.7592,     8.6868,     8.6868,     1.1938,     1.1938
              ,    23.5712,    154.711,    190.602,    190.602, 154.711 };

      int c = kGray;//color
      int s = 1;//style
      int w = 1;//width
      TPolyLine3D& front = view->AddPolyLine3D(23, c, w, s);
      TPolyLine3D& back  = view->AddPolyLine3D(23, c, w, s);
      Float_t zfront = 147.108;
      Float_t zback  = zfront+1660.0;
  
      for(int i = 0; i<23; ++i){
         front.SetPoint(i, x[i], y[i], zfront);
         back. SetPoint(i, x[i], y[i], zback);
      }

      for(int i = 0; i<22; ++i){
         TPolyLine3D& connect = view->AddPolyLine3D(2, c, w, s);
         connect.SetPoint(0,x[i],y[i],zfront);
         connect.SetPoint(1,x[i],y[i],zback);
      }

      

   }


}// namespace
////////////////////////////////////////////////////////////////////////
