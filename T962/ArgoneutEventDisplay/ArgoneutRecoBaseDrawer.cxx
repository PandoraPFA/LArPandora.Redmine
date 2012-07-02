/// \file    ArgoneutRecoBaseDrawer.cxx
/// \brief   Class to aid in the rendering of ArgoNeuT reconstruction objects
/// \author  msoderbe@syr.edu
#include <cmath>
#include <map>

#include "TMarker.h"
#include "TBox.h"
#include "TH1.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"

#include "EventDisplayBase/evdb.h"
#include "Geometry/geo.h"
#include "T962/ArgoneutEventDisplay/ArgoneutRecoBaseDrawer.h"
#include "EventDisplay/RecoDrawingOptions.h"
#include "EventDisplay/ColorDrawingOptions.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "Simulation/sim.h"
#include "RecoBase/recobase.h"
#include "Filters/ChannelFilter.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "T962/T962_Objects/MINOS.h"
#include "T962/T962_Objects/MINOSTrackMatch.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace {
   // Utility function to make uniform error messages.
   void writeErrMsg(const char* fcn,
                    cet::exception const& e)
   {
      mf::LogWarning("ArgoneutRecoBaseDrawer") << "ArgoneutRecoBaseDrawer::" << fcn
                                               << " failed with message:\n"
                                               << e;
   }
}

static const int kNCOLS = 14;
static const int kColor[kNCOLS] = { 2, 3, 4, 5, 6, 7, 8, 29, 30, 38, 40, 41, 42, 46 };

namespace argoevd{

   //......................................................................
   ArgoneutRecoBaseDrawer::ArgoneutRecoBaseDrawer() 
   {

   }

   //......................................................................
   ArgoneutRecoBaseDrawer::~ArgoneutRecoBaseDrawer() 
   {
 
   }

   //......................................................................
  void ArgoneutRecoBaseDrawer::ArgoSpacePoint(std::vector<const recob::SpacePoint*> sps,
                                               int                 id, 
                                               evdb::View3D*       view,
                                               bool                matched,
                                               float               charge)
   {
      // loop over all space points in the prong and draw them
      int c = 3;
      if(matched && charge==-1) c = kBlue;//negative
      if(matched && charge==1)  c = kRed; //positive
      TPolyMarker3D& pm = view->AddPolyMarker3D(sps.size(), c, 1, 3);

      for(size_t s = 0; s < sps.size(); ++s){
         const double *xyz = sps[s]->XYZ();
         pm.SetPoint(s, xyz[0], xyz[1], xyz[2]);
      }
    
      return;
   }

   //......................................................................
   void ArgoneutRecoBaseDrawer::ArgoProng3D(const art::Event& evt,
                                            evdb::View3D*     view)
   {
      art::ServiceHandle<evd::RawDrawingOptions>  rawopt;
      if (rawopt->fDrawRawDataOrCalibWires < 1) return;
      art::ServiceHandle<evd::RecoDrawingOptions> drawopt;
//       if (drawopt->fDrawProngs!=0) {

//          for (unsigned int imod=0; imod<drawopt->fProngLabels.size(); ++imod) {
//             std::string const which = drawopt->fProngLabels[imod];
      
//             std::vector<const recob::Prong*> prong;
//             this->GetArgoProngs(evt, which, prong);

//             for(size_t p = 0; p < prong.size(); ++p)
// 	      this->DrawProng3D(*(prong[p]), prong[p]->ID(), view);
//          }
//       }
  
      if(drawopt->fDrawTracks!=0){
         for (unsigned int imod=0; imod<drawopt->fTrackLabels.size(); ++imod){
            std::string const which = drawopt->fTrackLabels[imod];
      
            art::PtrVector<recob::Track> track;
            this->GetArgoTracks(evt, which, track);

	    art::FindMany<recob::SpacePoint> fmsp(track, evt, which);

            art::PtrVector<t962::MINOS> minos;
            this->GetMinos(evt, "minos", minos);

            art::PtrVector<t962::MINOSTrackMatch> minosmatch;
            this->GetMinosTrackMatch(evt, "matchtracks", minosmatch); 
          
            for(size_t p = 0; p < track.size(); ++p){
               bool matched = false;
               float charge = 0;
               for(size_t q= 0; q < minosmatch.size(); ++q)
               {
                  if(minosmatch[q]->fArgoNeuTtrackid == track[p]->ID()){
                     matched = true;
                     for(size_t r = 0; r<minos.size(); ++r){
                        if(minosmatch[q]->fMINOStrackid == minos[r]->ftrkIndex) charge = minos[r]->fcharge;
                     }
                  }
               }  
	       std::vector<const recob::SpacePoint*> sps = fmsp.at(p);
               this->ArgoSpacePoint(sps, track[p]->ID(), view, matched,charge);
            }
         }
      }
    
      if (drawopt->fDrawShowers!=0) {
         for (unsigned int imod=0; imod<drawopt->fShowerLabels.size(); ++imod){
            std::string const which = drawopt->fShowerLabels[imod];
      
            art::PtrVector<recob::Shower> shower;
            this->GetArgoShowers(evt, which, shower);

            for(size_t p = 0; p < shower.size(); ++p)
	      this->DrawShower3D(*(shower[p]), shower[p]->ID(), view);
         }
      }
    
      return;
   }




   //.....................................................................
   int ArgoneutRecoBaseDrawer::GetArgoTracks(const art::Event&                 evt,
                                             const std::string&                which,
                                             art::PtrVector<recob::Track>&     track)
   {
      art::PtrVector<recob::Track> temp;
      art::Handle< std::vector<recob::Track> > handle;
      try{
         evt.getByLabel(which,handle);
         for(unsigned int i = 0; i < handle->size(); ++i){
	   art::Ptr<recob::Track> w(handle, i);
            temp.push_back(w);
         }
	 
         temp.swap(track);
      }
      catch(cet::exception& e){
         writeErrMsg("GetArgoTracks", e);
      }
      
      return track.size();

   }

   //.....................................................................
   int ArgoneutRecoBaseDrawer::GetArgoShowers(const art::Event&                 evt,
					      const std::string&                which,
					      art::PtrVector<recob::Shower>&    shower)
   {
      art::PtrVector<recob::Shower> temp;
      art::Handle< std::vector<recob::Shower> > handle;
      try{
         evt.getByLabel(which,handle);
         for(unsigned int i = 0; i < handle->size(); ++i){
	   art::Ptr<recob::Shower> w(handle, i);
            temp.push_back(w);
         }
	 
         temp.swap(shower);
      }
      catch(cet::exception& e){
         writeErrMsg("GetArgoShowers", e);
      }
      
      return shower.size();

   }

   //......................................................................
   int ArgoneutRecoBaseDrawer::GetMinos(const art::Event&            evt,
                                const std::string&           which,
                                art::PtrVector<t962::MINOS>& minos) 
   {
      minos.clear();

      art::Handle< std::vector<t962::MINOS> > minoshandle;
      art::PtrVector<t962::MINOS> temp;

      try{
         evt.getByLabel(which, minoshandle);
       
         for(unsigned int i = 0; i < minoshandle->size(); ++i){
            art::Ptr<t962::MINOS> w(minoshandle, i);
            temp.push_back(w);
         }
         temp.swap(minos);
      }
      catch(cet::exception& e){
         writeErrMsg("GetMinos", e);
      }
    
      return minos.size();
   }
   //......................................................................
   int ArgoneutRecoBaseDrawer::GetMinosTrackMatch(const art::Event&            evt,
                                          const std::string&           which,
                                          art::PtrVector<t962::MINOSTrackMatch>& minosmatch) 
   {
      minosmatch.clear();

      art::Handle< std::vector<t962::MINOSTrackMatch> > minosmatchhandle;
      art::PtrVector<t962::MINOSTrackMatch> temp;

      try{
         evt.getByLabel(which, minosmatchhandle);
       
         for(unsigned int i = 0; i < minosmatchhandle->size(); ++i){
            art::Ptr<t962::MINOSTrackMatch> w(minosmatchhandle, i);
            temp.push_back(w);
         }
         temp.swap(minosmatch);
      }
      catch(cet::exception& e){
         writeErrMsg("GetMinosTrackMatch", e);
      }
    
      return minosmatch.size();
   }


 

}// namespace
////////////////////////////////////////////////////////////////////////
