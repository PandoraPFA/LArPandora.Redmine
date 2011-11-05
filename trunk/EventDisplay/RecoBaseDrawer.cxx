/// \file    RecoBaseDrawer.cxx
/// \brief   Class to aid in the rendering of RecoBase objects
/// \author  brebel@fnal.gov
/// \version $Id: RecoBaseDrawer.cxx,v 1.3 2010/11/11 22:47:19 p-novaart Exp $
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
#include "EventDisplay/RecoBaseDrawer.h"
#include "EventDisplay/RecoDrawingOptions.h"
#include "EventDisplay/ColorDrawingOptions.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "Simulation/sim.h"
#include "RecoBase/recobase.h"
#include "Filters/ChannelFilter.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"

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
    mf::LogWarning("RecoBaseDrawer") << "RecoBaseDrawer::" << fcn
				     << " failed with message:\n"
				     << e;
  }
}

static const int kNCOLS = 14;
static const int kColor[kNCOLS] = { 2, 3, 4, 5, 6, 7, 8, 29, 30, 38, 40, 41, 42, 46 };

namespace evd{

  //......................................................................
  RecoBaseDrawer::RecoBaseDrawer() 
  {
    art::ServiceHandle<geo::Geometry> geo;

    // set the list of bad channels in this detector
    filter::ChannelFilter cf;
    for(size_t t = 0; t < geo->NTPC(); ++t){
      for(size_t p = 0; p < geo->Nplanes(t); ++p){
	unsigned int nplanes=geo->Nplanes(t);
	fWireMin.resize(nplanes,-1);   
	fWireMax.resize(nplanes,-1);    
	fTimeMin.resize(nplanes,-1);    
	fTimeMax.resize(nplanes,-1);    
	fRawCharge.resize(nplanes,0);   
	fConvertedCharge.resize(nplanes,0);	
	for(size_t w = 0; w < geo->TPC(t).Plane(p).Nwires(); ++w){
	  unsigned int channel = geo->PlaneWireToChannel(p, w, t);
	  if( cf.BadChannel(channel) ) fBadChannels.push_back(channel);
	}
      }
    }

  

  }

  //......................................................................
  RecoBaseDrawer::~RecoBaseDrawer() 
  {
 
  }

  //......................................................................
  void RecoBaseDrawer::Wire2D(const art::Event& evt,
			      evdb::View2D*     view,
			      unsigned int      plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions>  rawopt;
    art::ServiceHandle<evd::RecoDrawingOptions> recoopt;
    if (rawopt->fDrawRawDataOrCalibWires < 1) return;
  
    art::ServiceHandle<geo::Geometry> geo;
    art::ServiceHandle<evd::ColorDrawingOptions> cst;

    unsigned int p    = 0;
    unsigned int w    = 0;
    unsigned int t    = 0;
    int ticksPerPoint = rawopt->fTicksPerPoint;
    int ticks         = 0; // fill this below

    for (unsigned int imod=0; imod<recoopt->fWireLabels.size(); ++imod) {
      std::string const which = recoopt->fWireLabels[imod];
    
      art::PtrVector<recob::Wire> wires;
      this->GetWires(evt, which, wires);

      if(wires.size() < 1) return;

      ticks = wires[0]->fSignal.size();
      
      for (unsigned int i=0; i<wires.size(); ++i) {
      
	unsigned int channel = wires[i]->RawDigit()->Channel();
	geo->ChannelToWire(channel, t, p, w);

	// check plane and tpc are correct
	if(p != plane || t != rawopt->fTPC) continue;
      
	double wire = 1.*w;
	double tick = 0;
	// get an iterator over the adc values
	std::vector<double>::const_iterator itr = wires[i]->fSignal.begin();
	while( itr != wires[i]->fSignal.end() ){
	  int ticksUsed = 0;
	  double tdcsum = 0.;
	  double adcsum = 0.;
	  while(ticksUsed < ticksPerPoint && itr != wires[i]->fSignal.end()){
	    tdcsum  += tick;
	    adcsum  += (1.*(*itr));
	    ++ticksUsed;
	    tick += 1.;
	    itr++; // this advance of the iterator is sufficient for the external loop too
	  }
	  double adc = adcsum/ticksPerPoint;
	  double tdc = tdcsum/ticksPerPoint;
	
	  if(TMath::Abs(adc) < rawopt->fMinSignal) continue;
	
	  int    co = 0;
	  double sf = 1.;
	  double q0 = 1000.0;
	
	  co = cst->CalQ().GetColor(adc);
	  if (rawopt->fScaleDigitsByCharge) {
	    sf = sqrt(adc/q0);
	    if (sf>1.0) sf = 1.0;
	  }

	  if(rawopt->fAxisOrientation < 1){	
	    TBox& b1 = view->AddBox(wire-sf*0.5,tdc-sf*0.5*ticksPerPoint,wire+sf*0.5,tdc+sf*0.5*ticksPerPoint);
	    b1.SetFillStyle(1001);
	    b1.SetFillColor(co);    
	    b1.SetBit(kCannotPick);
	  }
	  else{
	    TBox& b1 = view->AddBox(tdc-sf*0.5*ticksPerPoint,wire-sf*0.5,tdc+sf*0.5*ticksPerPoint,wire+sf*0.5);
	    b1.SetFillStyle(1001);
	    b1.SetFillColor(co);    
	    b1.SetBit(kCannotPick);
	  }
	  
	// 	  TBox& b2 = view->AddBox(wire-0.1,tdc-0.1,wire+0.1,tdc+0.1);
	// 	  b2.SetFillStyle(0);
	// 	  b2.SetLineColor(15);
	// 	  b2.SetBit(kCannotPick);
	 
	}// end loop over samples 
      }//end loop over wires
    }// end loop over wire module labels

    // now loop over all the bad channels and set them to 0 adc
    for(size_t bc = 0; bc < fBadChannels.size(); ++bc){
      unsigned int t = 0;
      unsigned int p = 0;
      unsigned int w = 0;
      geo->ChannelToWire(fBadChannels[bc], t, p, w);

      // check that we have correct plane and tpc
      if(p != plane || t != rawopt->fTPC) continue;

      if(rawopt->fMinSignal > 0) continue;
	
      int      co = cst->CalQ().GetColor(0);
      double wire = 1.*w;

      for(int i = 0; i < ticks; i += ticksPerPoint){
	double tdc = i + 0.5*ticksPerPoint;
	
	if(rawopt->fAxisOrientation < 1){	
	  TBox& b1 = view->AddBox(wire-0.5,tdc-0.5*ticksPerPoint,wire+0.5,tdc+0.5*ticksPerPoint);
	  b1.SetFillStyle(1001);
	  b1.SetFillColor(co);    
	  b1.SetBit(kCannotPick);
	}
	else{
	  TBox& b1 = view->AddBox(tdc-0.5*ticksPerPoint,wire-0.5,tdc+0.5*ticksPerPoint,wire+0.5);
	  b1.SetFillStyle(1001);
	  b1.SetFillColor(co);    
	  b1.SetBit(kCannotPick);
	}
      }
    }// end loop over bad channels    


  }

  //......................................................................

  ///
  /// Render Hit objects on a 2D viewing canvas
  ///
  /// @param evt    : Event handle to get data objects from
  /// @param view   : Pointer to view to draw on
  /// @param plane  : plane number of view
  ///
  void RecoBaseDrawer::Hit2D(const art::Event& evt,
			     evdb::View2D*     view,
			     unsigned int      plane) 
  {
    art::ServiceHandle<evd::RawDrawingOptions>  rawopt;
    art::ServiceHandle<evd::RecoDrawingOptions> recoopt;
    if (recoopt->fDrawHits == 0)               return;
    if (rawopt->fDrawRawDataOrCalibWires < 1) return;
    art::ServiceHandle<util::LArProperties> larp;
    art::ServiceHandle<util::DetectorProperties> detp;
    
    
    fRawCharge[plane]=0;
    fConvertedCharge[plane]=0;
    art::ServiceHandle<geo::Geometry> geo;

    unsigned int p = 0;
    unsigned int w  = 0;
    unsigned int t  = 0;
    // to make det independent later:
    double mint=5000,maxt=0,minw=5000,maxw=0;
    
    for (unsigned int imod=0; imod<recoopt->fHitLabels.size(); ++imod) {
      std::string const which = recoopt->fHitLabels[imod];
    
      art::PtrVector<recob::Hit> hits;
      this->GetHits(evt, which, hits);
  
      // Display all hits on the two 2D views provided
      for (unsigned int i=0; i<hits.size(); ++i) {
	art::Ptr<recob::Hit> h = hits[i];
	unsigned short channel = h->Wire()->RawDigit()->Channel();

	geo->ChannelToWire(channel, t, p, w);

	// check that we are on the correct plane and in the correct TPC
	if(p != plane || t != rawopt->fTPC) continue;

	evdb::View2D*  v = view;
      
	// Try to get the "best" charge measurement, ie. the one last in
	// the calibration chain
      
	float time = h->PeakTime();
      
	double sf =   1.0;
	
	if(w<minw)
	    minw=w;
	if(w>maxw)
	    maxw=w;
	if(time<mint)
	    mint=time;
	if(time>maxt)
	    maxt=time;
	
	fRawCharge[p]+=h->Charge(true);
	fConvertedCharge[p]+=larp->BirksCorrectionAmplitude(h->Charge(true)/geo->WirePitch(),1./detp->ElectronsToADC());
	
	if(rawopt->fAxisOrientation < 1){
	  TBox& b1 = v->AddBox(w-sf*0.5, time-sf*0.5, w+sf*0.5, time+sf*0.5);
	  b1.SetFillStyle(0);
	  b1.SetLineColor(kBlack);
	  b1.SetBit(kCannotPick);
	}
	else{
	  TBox& b1 = v->AddBox(time-sf*0.5, w-sf*0.5, time+sf*0.5, w+sf*0.5);
	  b1.SetFillStyle(0);
	  b1.SetLineColor(kBlack);
	  b1.SetBit(kCannotPick);
	}

// 	TBox& b2 = v->AddBox(w-0.5, t-0.5, w+0.5, t+0.5);
// 	b2.SetFillStyle(0);
// 	b2.SetFillColor(17);
// 	b2.SetBit(kCannotPick);

      } // loop on i cell hits
    } // loop on imod folders
      fWireMin[plane]=minw;   
      fWireMax[plane]=maxw;    
      fTimeMin[plane]=mint;    
      fTimeMax[plane]=maxt; 
    
  }
  //......................................................................

  ///
  /// Render Hit objects on a 2D viewing canvas
  ///
  /// @param hits   : vector of hits for the veiw
  /// @param color  : color of associated cluster/prong
  /// @param view   : Pointer to view to draw on
  ///
  /// assumes the hits are all from the correct plane for the given view
  void RecoBaseDrawer::Hit2D(art::PtrVector<recob::Hit> hits,
			     int                        color,
			     evdb::View2D*              view) 
  {
    art::ServiceHandle<geo::Geometry> geo;
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    
    unsigned int p = 0;
    unsigned int w  = 0;
    unsigned int t  = 0;

    
    
    for(unsigned int c = 0; c < hits.size(); ++c){
      art::Ptr<recob::Hit> h = hits[c];

      unsigned short channel = hits[c]->Wire()->RawDigit()->Channel();

      geo->ChannelToWire(channel, t, p, w);

      // check that we are in the correct TPC
      // the view should tell use we are in the correct plane
      if(t != rawopt->fTPC) continue;
      
      // Try to get the "best" charge measurement, ie. the one last in
      // the calibration chain
      float t = h->PeakTime();

      
      
      if(rawopt->fAxisOrientation < 1){
	TBox& b1 = view->AddBox(w-0.5, t-0.5, w+0.5, t+0.5);
	b1.SetFillStyle(0);
	b1.SetLineColor(color);
	b1.SetBit(kCannotPick);
      }
      else{
	TBox& b1 = view->AddBox(t-0.5, w-0.5, t+0.5, w+0.5);
	b1.SetFillStyle(0);
	b1.SetLineColor(color);
	b1.SetBit(kCannotPick);
      }

    } // loop on c cell hits
  }


//........................................................................
int RecoBaseDrawer::GetRegionOfInterest(int plane,int& minw,int& maxw,int& mint,int& maxt)
{
  art::ServiceHandle<geo::Geometry> geo;
 
  if((unsigned int)plane>fWireMin.size())
    {mf::LogWarning  ("RecoBaseDrawer") << " Requested plane " << plane <<" is larger than those available " << std::endl;
    return -1;
    }
  
  minw=fWireMin[plane];
  maxw=fWireMax[plane];
  mint=fTimeMin[plane];
  maxt=fTimeMax[plane];
  
  //make values a bit larger, but make sure they don't go out of bounds 
  minw= (minw-30<0) ? 0 : minw-30;
  mint= (mint-10<0) ? 0 : mint-10;

  art::ServiceHandle<evd::RawDrawingOptions> rawopt;
  int fTicks = rawopt->fTicks;
  
  maxw= (maxw+10>(int)geo->Nwires(plane)) ? geo->Nwires(plane) : maxw+10;
  maxt= (maxt+10>fTicks) ? fTicks : maxt+10;
    
  return 0;
}

//......................................................................
  void RecoBaseDrawer::GetChargeSum(int plane,double& charge,double& convcharge)
  {
   charge=fRawCharge[plane]; 
   convcharge=fConvertedCharge[plane];    
    
  }



  //......................................................................

  void RecoBaseDrawer::Cluster2D(const art::Event& evt,
				 evdb::View2D*     view,
				 unsigned int      plane) 
  {
    art::ServiceHandle<evd::RecoDrawingOptions> drawopt;
    art::ServiceHandle<evd::RawDrawingOptions>   rawopt;
    if (rawopt->fDrawRawDataOrCalibWires < 1) return;
    if (drawopt->fDrawClusters == 0)          return;

    art::ServiceHandle<geo::Geometry> geo;

    // if user sets "DrawClusters" to 2, draw the clusters differently:
    bool drawAsMarkers = (drawopt->fDrawClusters==2);
    int count = 0; // how many clusters have we drawn?

    for (unsigned int imod=0; imod<drawopt->fClusterLabels.size(); ++imod) {
      std::string const which = drawopt->fClusterLabels[imod];
    
      art::PtrVector<recob::Cluster> clust;
      this->GetClusters(evt, which, clust);
      for (unsigned int ic=0; ic<clust.size(); ++ic) {

	art::Ptr<recob::Cluster> c = clust[ic];

	art::PtrVector<recob::Hit> hits = c->Hits();
	unsigned int w = 0;
	unsigned int t = 0;
	unsigned int p = 0;
	geo->ChannelToWire(hits[0]->Wire()->RawDigit()->Channel(), t, p, w);

	// check for correct plane and tpc
	if(p != plane || t != rawopt->fTPC) continue;

	if (drawAsMarkers) {
	  // draw cluster with unique marker  
	  // Place this cluster's unique marker at the hit's location
	  int color  = kColor[count%kNCOLS];
	  this->Hit2D(hits, color, view);
	}
	else {

	  // default "outline" method:
	  std::vector<double> tpts, wpts;
      
	  this->GetClusterOutlines(c, tpts, wpts, plane);
      
	  int lcolor = 9; // line color
	  int fcolor = 9; // fill color
	  int width = 2;  // line width
	  int style = 1;  // 1=solid line style
	  if (view!=0) {
	    TPolyLine& p1 = view->AddPolyLine(wpts.size(), 
					      lcolor,
					      width,
					      style);
	    TPolyLine& p2 = view->AddPolyLine(wpts.size(),
					      lcolor,
					      width,
					      style);
	    p1.SetOption("f");
	    p1.SetFillStyle(3003);
	    p1.SetFillColor(fcolor);
	    for (unsigned int i=0; i<wpts.size(); ++i) {
	      if(rawopt->fAxisOrientation < 1){
		p1.SetPoint(i, wpts[i], tpts[i]);
		p2.SetPoint(i, wpts[i], tpts[i]);
	      }
	      else{
		p1.SetPoint(i, tpts[i], wpts[i]);
		p2.SetPoint(i, tpts[i], wpts[i]);
	      }
	    } // loop on i points in ZX view
	  } // if we have a cluster in the ZX view
	}// end if outline mode

	// draw the direction cosine of the cluster as well as it's starting point
	this->Draw2DEndPointAndSlope(c->StartPos()[0], c->StartPos()[1], c->dTdW(), 
				     kColor[count%kNCOLS], view);

	count++;
      } // loop on ic clusters
    } // loop on imod folders
  }

  //......................................................................
  void RecoBaseDrawer::Draw2DEndPointAndSlope(double        x,
					      double        y,
					      double        slope,
					      int           color,
					      evdb::View2D* view)
  {
    art::ServiceHandle<evd::RecoDrawingOptions>  recoopt;
    if(recoopt->fDraw2DEndPoints < 1) return;

    double x1 = x;
    double y1 = y;
    double slope1 = slope;

    art::ServiceHandle<evd::RawDrawingOptions>   rawopt;
    if(rawopt->fDrawRawDataOrCalibWires < 1) return;
    if(rawopt->fAxisOrientation > 0){
      x1 = y;
      y1 = x;
      if(fabs(slope) > 0.) slope1 = 1./slope;
      else slope1 = 1.e6;
    }		

    TMarker& strt = view->AddMarker(x1, y1, color, 29, 2.0);
    strt.SetMarkerColor(color); // stupid line to shut up compiler warning
    TLine& l = view->AddLine(x1, y1, x1+50, y1 + slope1*50);
    l.SetLineColor(color);
    l.SetLineWidth(2);

    return;
  }

  //......................................................................
  ///
  /// Make a set of points which outline a cluster
  ///
  /// @param c      : Reco base cluster to outline
  /// @param wpts   : wire values of the outlines
  /// @param tpts   : tdc values of the outlines
  /// @param plane  : plane number 
  ///
  void RecoBaseDrawer::GetClusterOutlines(art::Ptr<recob::Cluster> &c,
					  std::vector<double>&  wpts,
					  std::vector<double>&  tpts,
					  unsigned int          plane)
  {
    art::ServiceHandle<geo::Geometry>          geo;
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    
    art::PtrVector<recob::Hit> hits = c->Hits(plane);

    unsigned int p = 0;
    unsigned int tpc = 0;
    unsigned int w  = 0;

    // Map wire numbers to highest and lowest in the plane
    std::map<unsigned int, double> wlo, whi;
    // On first pass, initialize
    for (size_t j=0; j<hits.size(); ++j) {
      unsigned short channel = hits[j]->Wire()->RawDigit()->Channel();

      geo->ChannelToWire(channel, tpc, p, w);

      // check that we are on the correct plane and TPC
      if(p != plane || tpc != rawopt->fTPC) continue;

      wlo[w] = hits[j]->PeakTime();
      whi[w] = hits[j]->PeakTime();
    }

    double t = 0.;

    // Finalize on second pass
    for (size_t j=0; j<hits.size(); ++j) {
      unsigned short channel = hits[j]->Wire()->RawDigit()->Channel();

      geo->ChannelToWire(channel, tpc, p, w);

      t = hits[j]->PeakTime();

      if (t < wlo[w]) wlo[w] = t;
      if (t > whi[w]) whi[w] = t;
    }
    
    // Loop over wires and low times to make lines along bottom
    // edge. Work from upstream edge to downstream edge
    std::map<unsigned int, double>::iterator itr   (wlo.begin());
    std::map<unsigned int, double>::iterator itrEnd(wlo.end());
    for (; itr!=itrEnd; ++itr) {
	w = itr->first;
	t = itr->second;
	
	wpts.push_back(1.*w-0.1); tpts.push_back(t-0.1);
	wpts.push_back(1.*w+0.1); tpts.push_back(t-0.1);
      }
    
      // Loop over planes and high cells to make lines along top
      // edge. Work from downstream edge toward upstream edge
      std::map<unsigned int, double>::reverse_iterator ritr   (whi.rbegin());
      std::map<unsigned int, double>::reverse_iterator ritrEnd(whi.rend());
      for (; ritr!=ritrEnd; ++ritr) {
	w = ritr->first;
	t = ritr->second;

	wpts.push_back(1.*w+0.1); tpts.push_back(t+0.1);
	wpts.push_back(1.*w-0.1); tpts.push_back(t+0.1);
      }
    
      // Add link to starting point to close the box
      wpts.push_back(wpts[0]); tpts.push_back(tpts[0]);
    
  }

  //......................................................................

  void RecoBaseDrawer::DrawProng2D(std::vector<const recob::Prong*>& prong,
				   evdb::View2D*                     view,
				   unsigned int                      plane,
				   int                               id) 
  {
    art::ServiceHandle<evd::RawDrawingOptions>  rawopt;
    art::ServiceHandle<evd::RecoDrawingOptions> drawopt;
    art::ServiceHandle<util::LArProperties>     larp;

    art::ServiceHandle<geo::Geometry> geo;
    geo::View_t gview = geo->TPC(rawopt->fTPC).Plane(plane).View();

    art::ServiceHandle<util::DetectorProperties> detprop;
    
    for (unsigned int i=0; i<prong.size(); ++i) {

      int colID = id;
      if(id < 0) colID = prong[i]->ID();

      art::PtrVector<recob::Cluster> clusts = prong[i]->Clusters(gview);
      
      for(size_t c = 0; c < clusts.size(); ++c){
	art::PtrVector<recob::Hit> hits = clusts[c]->Hits(gview);
	this->Hit2D(hits, kColor[colID%kNCOLS], view);
	this->Draw2DEndPointAndSlope(clusts[c]->StartPos()[0], clusts[c]->StartPos()[1], 
				     clusts[c]->dTdW(), kColor[colID%kNCOLS], view);	
      }// end loop over clusters

//       // make a line that has the correct direction cosines for the prong at 
//       // its start and 
//       double dcstart[3] = {0.};
//       double dcend[3]   = {0.};
//     double sampleRate = detprop->SamplingRate();

//       prong[i]->Direction(dcstart, dcend);

//       std::vector<double> xyz(3, 0.);
//       std::vector<double> xyz2(3, 0.);

//       prong[i]->Extent(xyz, xyz2);
//       double length = sqrt((xyz2[0]-xyz[0])*(xyz2[0]-xyz[0]) +
// 			   (xyz2[1]-xyz[1])*(xyz2[1]-xyz[1]) + 
// 			   (xyz2[2]-xyz[2])*(xyz2[2]-xyz[2]) );

//       xyz2[0] = xyz[0] + length*dcstart[0];
//       xyz2[1] = xyz[1] + length*dcstart[1];
//       xyz2[2] = xyz[2] + length*dcstart[2];
	
//       if(xyz2[2] < 0.)                          xyz2[2] = 0.;
//       if(xyz2[2] > geo->DetLength() )           xyz2[2] = geo->DetLength();
//       if(fabs(xyz2[1]) > geo->DetHalfHeight() ) xyz2[1] = geo->DetHalfHeight();

//       dcstart[0] = xyz[0];  dcstart[1] = xyz[1];  dcstart[2] = xyz[2];
//       dcend[0]   = xyz2[0]; dcend[1]   = xyz2[1]; dcend[2]   = xyz2[2];

//       unsigned int chan1 = geo->NearestChannel(dcstart); 
//       unsigned int chan2 = geo->NearestChannel(dcend); 
//       unsigned int tpc   = 0;
//       unsigned int pln   = 0;
//       unsigned int w1    = 0;
//       unsigned int w2    = 0;
//       geo->ChannelToWire(chan1, tpc, pln, w1);
//       geo->ChannelToWire(chan2, tpc, pln, w2);

//       double time = xyz[0]*1e3/larp->DriftVelocity(larp->Efield(), larp->Temperature())/sampleRate;
//       double time2 = xyz2[0]*1e3/larp->DriftVelocity(larp->Efield(), larp->Temperature())/sampleRate;

//       time += detprop->TriggerOffset();
//       time2 += detprop->TriggerOffset();

//       std::cout << "draw line length: " << length
// 		<< " with start: " << w1 << "," << time 
// 		<< " and end: " << w2 << "," << time2 << std::endl;

//       if(rawopt->fAxisOrientation < 1){
// 	TLine& p = view->AddLine(w1, time, w2, time2);
// 	p.SetLineColor(kColor[colID%kNCOLS]);
// 	p.SetLineWidth(2);
//       }
//       else{
// 	TLine& p = view->AddLine(time, w1, time2, w2);
// 	p.SetLineColor(kColor[colID%kNCOLS]);
// 	p.SetLineWidth(2);
//       }
    } // end loop over prongs in event

    return;
  }

  //......................................................................
  void RecoBaseDrawer::Prong2D(const art::Event& evt,
			       evdb::View2D*     view,
			       unsigned int      plane) 
  {    

    art::ServiceHandle<evd::RawDrawingOptions>  rawopt;
    if (rawopt->fDrawRawDataOrCalibWires < 1) return;

    art::ServiceHandle<evd::RecoDrawingOptions> drawopt;
    if (drawopt->fDrawProngs!=0) {

      for (unsigned int imod=0; imod<drawopt->fProngLabels.size(); ++imod) {
	std::string const which = drawopt->fProngLabels[imod];
      
	std::vector<const recob::Prong*> prong;
	this->GetProngs(evt, which, prong);
	this->DrawProng2D(prong, view, plane);
      }
    }
  
    if(drawopt->fDrawTracks!=0){
      for (unsigned int imod=0; imod<drawopt->fTrackLabels.size(); ++imod){
	std::string const which = drawopt->fTrackLabels[imod];
      
	std::vector<const recob::Prong*> track;
	this->GetProngs(evt, which, track);

	this->DrawProng2D(track, view, plane);
      }
    }
  
    if (drawopt->fDrawShowers!=0) {
      for (unsigned int imod=0; imod<drawopt->fShowerLabels.size(); ++imod){
	std::string const which = drawopt->fShowerLabels[imod];
      
	std::vector<const recob::Prong*> shower;
	this->GetProngs(evt, which, shower);

	this->DrawProng2D(shower, view, plane);
      }
    }
  
  }

  //......................................................................
  void RecoBaseDrawer::Vertex2D(const art::Event& evt,
				evdb::View2D*     view,
				unsigned int      plane) 
  {    
    art::ServiceHandle<evd::RawDrawingOptions>  rawopt;
    if (rawopt->fDrawRawDataOrCalibWires < 1) return;
    art::ServiceHandle<evd::RecoDrawingOptions> drawopt;
    if(drawopt->fDrawVertices != 0){

      for (unsigned int imod = 0; imod < drawopt->fVertexLabels.size(); ++imod) {
	std::string const which = drawopt->fVertexLabels[imod];

	art::PtrVector<recob::Vertex> vertex;
	this->GetVertices(evt, which, vertex);

	for(size_t v = 0; v < vertex.size(); ++v){

	  // grab the Prongs from the vertex and draw those
	  std::vector<const recob::Prong*> prong;
	  for(size_t t = 0; t < vertex[v]->Tracks().size(); ++t)
	    prong.push_back(vertex[v]->Tracks()[t].get());

	  for(size_t s = 0; s < vertex[v]->Showers().size(); ++s)
	    prong.push_back(vertex[v]->Showers()[s].get());

	  this->DrawProng2D(prong, view, plane, vertex[v]->ID());
	} // end loop over vertices to draw from this label
      } // end loop over vertex module lables
    } // end if we are drawing vertices

    return;
  }

  //......................................................................
  void RecoBaseDrawer::Event2D(const art::Event& evt,
			       evdb::View2D*     view,
			       unsigned int      plane) 
  {    
    art::ServiceHandle<evd::RawDrawingOptions>  rawopt;
    if (rawopt->fDrawRawDataOrCalibWires < 1) return;
    art::ServiceHandle<evd::RecoDrawingOptions> drawopt;
    if(drawopt->fDrawEvents != 0){

      for (unsigned int imod = 0; imod < drawopt->fEventLabels.size(); ++imod) {
	std::string const which = drawopt->fEventLabels[imod];

	art::PtrVector<recob::Event> event;
	this->GetEvents(evt, which, event);

	for(size_t e = 0; e < event.size(); ++e){

	  // grab the vertices for this event
	  art::PtrVector<recob::Vertex> vertex = event[e]->Vertices();

	  for(size_t v = 0; v < vertex.size(); ++v){

	    // grab the Prongs from the vertex and draw those
	    std::vector<const recob::Prong*> prong;
	    for(size_t t = 0; t < vertex[v]->Tracks().size(); ++t)
	      prong.push_back(vertex[v]->Tracks()[t].get());

	    for(size_t s = 0; s < vertex[v]->Showers().size(); ++s)
	      prong.push_back(vertex[v]->Showers()[s].get());

	    this->DrawProng2D(prong, view, plane, event[e]->ID());
	  } // end loop over vertices from this event
	} // end loop over events 
      } // end loop over event module lables
    } // end if we are drawing events

    return;
  }

  //......................................................................
  void RecoBaseDrawer::SpacePoint(const recob::Prong* prong,
				  int                 id, 
				  evdb::View3D*       view)
  {
    // loop over all space points in the prong and draw them
    const std::vector<recob::SpacePoint> sps = prong->SpacePoints();

    TPolyMarker3D& pm = view->AddPolyMarker3D(sps.size(), kColor[id%kNCOLS], 1, 3);

    for(size_t s = 0; s < sps.size(); ++s){
      const double *xyz = sps[s].XYZ();
      pm.SetPoint(s, xyz[0], xyz[1], xyz[2]);
    }
    
    return;
  }

  //......................................................................
  void RecoBaseDrawer::Prong3D(const art::Event& evt,
			       evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RawDrawingOptions>  rawopt;
    if (rawopt->fDrawRawDataOrCalibWires < 1) return;
    art::ServiceHandle<evd::RecoDrawingOptions> drawopt;
    if (drawopt->fDrawProngs!=0) {

      for (unsigned int imod=0; imod<drawopt->fProngLabels.size(); ++imod) {
	std::string const which = drawopt->fProngLabels[imod];
      
	std::vector<const recob::Prong*> prong;
	this->GetProngs(evt, which, prong);

	for(size_t p = 0; p < prong.size(); ++p)
	  this->SpacePoint(prong[p], prong[p]->ID(), view);
      }
    }
  
    if(drawopt->fDrawTracks!=0){
      for (unsigned int imod=0; imod<drawopt->fTrackLabels.size(); ++imod){
	std::string const which = drawopt->fTrackLabels[imod];
      
	std::vector<const recob::Prong*> track;
	this->GetProngs(evt, which, track);

	for(size_t p = 0; p < track.size(); ++p)
	  this->SpacePoint(track[p], track[p]->ID(), view);
      }
    }
  
    if (drawopt->fDrawShowers!=0) {
      for (unsigned int imod=0; imod<drawopt->fShowerLabels.size(); ++imod){
	std::string const which = drawopt->fShowerLabels[imod];
      
	std::vector<const recob::Prong*> shower;
	this->GetProngs(evt, which, shower);

	for(size_t p = 0; p < shower.size(); ++p)
	  this->SpacePoint(shower[p], shower[p]->ID(), view);
      }
    }
    
    return;
  }

  //......................................................................
  void RecoBaseDrawer::Vertex3D(const art::Event& evt,
				evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RawDrawingOptions>  rawopt;
    if (rawopt->fDrawRawDataOrCalibWires < 1) return;
    art::ServiceHandle<evd::RecoDrawingOptions> drawopt;
    if(drawopt->fDrawVertices != 0){

      for (unsigned int imod = 0; imod < drawopt->fVertexLabels.size(); ++imod) {
	std::string const which = drawopt->fVertexLabels[imod];

	art::PtrVector<recob::Vertex> vertex;
	this->GetVertices(evt, which, vertex);

	for(size_t v = 0; v < vertex.size(); ++v){

	  // grab the Prongs from the vertex and draw those
	  std::vector<const recob::Prong*> prong;
	  for(size_t t = 0; t < vertex[v]->Tracks().size(); ++t)
	    this->SpacePoint(vertex[v]->Tracks()[t].get(), vertex[v]->ID(), view);

	  for(size_t s = 0; s < vertex[v]->Showers().size(); ++s)
	    this->SpacePoint(vertex[v]->Showers()[s].get(), vertex[v]->ID(), view);

	  double xyz[3] = {0.};
	  vertex[v]->XYZ(xyz);
	  TPolyMarker3D& pm = view->AddPolyMarker3D(1, kColor[vertex[v]->ID()%kNCOLS], 29, 6);
	  pm.SetPoint(0, xyz[0], xyz[1], xyz[2]);

	} // end loop over vertices to draw from this label
      } // end loop over vertex module lables
    } // end if we are drawing vertices

    return;
  }

  //......................................................................
  void RecoBaseDrawer::Event3D(const art::Event& evt,
			       evdb::View3D*     view) 
  {    
    art::ServiceHandle<evd::RawDrawingOptions>  rawopt;
    if (rawopt->fDrawRawDataOrCalibWires < 1) return;
    art::ServiceHandle<evd::RecoDrawingOptions> drawopt;
    if(drawopt->fDrawEvents != 0){

      for (unsigned int imod = 0; imod < drawopt->fEventLabels.size(); ++imod) {
	std::string const which = drawopt->fEventLabels[imod];

	art::PtrVector<recob::Event> event;
	this->GetEvents(evt, which, event);

	for(size_t e = 0; e < event.size(); ++e){

	  // grab the vertices for this event
	  art::PtrVector<recob::Vertex> vertex = event[e]->Vertices();

	  for(size_t v = 0; v < vertex.size(); ++v){

	    // grab the Prongs from the vertex and draw those
	    std::vector<const recob::Prong*> prong;
	    for(size_t t = 0; t < vertex[v]->Tracks().size(); ++t)
	      this->SpacePoint(vertex[v]->Tracks()[t].get(), event[e]->ID(), view);

	    for(size_t s = 0; s < vertex[v]->Showers().size(); ++s)
	      this->SpacePoint(vertex[v]->Showers()[s].get(), event[e]->ID(), view);

	  } // end loop over vertices from this event

	  double xyz[3] = {0.};
	  event[e]->PrimaryVertex()->XYZ(xyz);
	  TPolyMarker3D& pm = view->AddPolyMarker3D(1, kColor[event[e]->ID()%kNCOLS], 29, 6);
	  pm.SetPoint(0, xyz[0], xyz[1], xyz[2]);

	} // end loop over events 
      } // end loop over event module lables
    } // end if we are drawing events

    return;
  }

  //......................................................................
  int RecoBaseDrawer::GetWires(const art::Event&            evt,
			       const std::string&           which,
			       art::PtrVector<recob::Wire>& wires) 
  {
    wires.clear();

    art::Handle< std::vector<recob::Wire> > wcol;
    art::PtrVector<recob::Wire> temp;

    try{
      evt.getByLabel(which, wcol);
      
      for(unsigned int i = 0; i < wcol->size(); ++i){
	art::Ptr<recob::Wire> w(wcol, i);
	temp.push_back(w);
      }
      temp.swap(wires);
    }
    catch(cet::exception& e){
      writeErrMsg("GetWires", e);
    }

    return wires.size();
  }

  //......................................................................
  int RecoBaseDrawer::GetHits(const art::Event&           evt,
			      const std::string&          which,
			      art::PtrVector<recob::Hit>& hits) 
  {
    hits.clear();

    art::Handle< std::vector<recob::Hit> > hcol;
    art::PtrVector<recob::Hit> temp;

    try{
      evt.getByLabel(which, hcol);
      
      for(unsigned int i = 0; i < hcol->size(); ++i){
	art::Ptr<recob::Hit> h(hcol, i);
	temp.push_back(h);
      }
      temp.swap(hits);
    }
    catch(cet::exception& e){
      writeErrMsg("GetHits", e);
    }

    return hits.size();
  }

  //......................................................................
  int RecoBaseDrawer::GetClusters(const art::Event&               evt, 
				  const std::string&              which,
				  art::PtrVector<recob::Cluster>& clust)
  {
    clust.clear();
    art::PtrVector<recob::Cluster> temp;

    art::Handle< std::vector<recob::Cluster> > clcol;

    try{
      evt.getByLabel(which, clcol);
      for(unsigned int i = 0; i < clcol->size(); ++i){
	art::Ptr<recob::Cluster> cl(clcol, i);
	temp.push_back(cl);
      }
      temp.swap(clust);
    }
    catch(cet::exception& e){
      writeErrMsg("GetClusters", e);
    }

    return clust.size();
  }

  //......................................................................
  int RecoBaseDrawer::GetProngs(const art::Event&                 evt, 
				const std::string&                which,
				std::vector<const recob::Prong*>& prong)
  {
    std::vector<const recob::Prong*> temp(prong);
    try{
      evt.getView(which,temp);
      temp.swap(prong);
    }
    catch(cet::exception& e){
      writeErrMsg("GetProngs", e);
    }

    return prong.size();
  }

  //......................................................................
  int RecoBaseDrawer::GetVertices(const art::Event&              evt, 
				  const std::string&             which,
				  art::PtrVector<recob::Vertex>& vertex)
  {
    vertex.clear();
    art::PtrVector<recob::Vertex> temp;

    art::Handle< std::vector<recob::Vertex> > vcol;

    try{
      evt.getByLabel(which, vcol);
      for(unsigned int i = 0; i < vcol->size(); ++i){
	art::Ptr<recob::Vertex> v(vcol, i);
	temp.push_back(v);
      }
      temp.swap(vertex);
    }
    catch(cet::exception& e){
      writeErrMsg("GetVertices", e);
    }

    return vertex.size();
  }

  //......................................................................
  int RecoBaseDrawer::GetEvents(const art::Event&             evt, 
				const std::string&            which,
				art::PtrVector<recob::Event>& event)
  {
    event.clear();
    art::PtrVector<recob::Event> temp;

    art::Handle< std::vector<recob::Event> > ecol;

    try{
      evt.getByLabel(which, ecol);
      for(unsigned int i = 0; i < ecol->size(); ++i){
	art::Ptr<recob::Event> e(ecol, i);
	temp.push_back(e);
      }
      temp.swap(event);
    }
    catch(cet::exception& e){
      writeErrMsg("GetEvents", e);
    }

    return event.size();
  }

  //......................................................................
  void RecoBaseDrawer::FillTQHisto(const art::Event&    evt,
				   unsigned int         plane,
				   unsigned int         wire,
				   TH1F*                histo,
				   std::vector<double>& hstart,
                   std::vector<double>& hend,
                   std::vector<double>& hitamplitudes,
                   std::vector<double>& hpeaktimes)
  {

    // Check if we're supposed to draw raw hits at all
    art::ServiceHandle<evd::RecoDrawingOptions> drawopt;
    art::ServiceHandle<evd::RawDrawingOptions>  rawopt;
    if (rawopt->fDrawRawDataOrCalibWires==0) return;

    art::ServiceHandle<geo::Geometry> geo;
  
    unsigned int p = 0;
    unsigned int w  = 0;
    unsigned int t  = 0;

    for (unsigned int imod=0; imod<drawopt->fWireLabels.size(); ++imod) {
      std::string const which = drawopt->fWireLabels[imod];

      art::PtrVector<recob::Wire> wires;
      this->GetWires(evt, which, wires);

      for (unsigned int i=0; i<wires.size(); ++i) {

	unsigned short channel = wires[i]->RawDigit()->Channel();

	geo->ChannelToWire(channel, t, p, w);

	// check for correct plane, wire and tpc
	if(p != plane || w != wire || t != rawopt->fTPC) continue;
	
	for(size_t s = 0; s < wires[i]->fSignal.size(); ++s)
	  histo->Fill(1.*s, wires[i]->fSignal[s]);

      }//end loop over wires
    }//end loop over wire modules

    for (unsigned int imod=0; imod<drawopt->fHitLabels.size(); ++imod) {
      std::string const which = drawopt->fHitLabels[imod];

      art::PtrVector<recob::Hit> hits;
      this->GetHits(evt, which, hits);


      for (unsigned int i=0; i<hits.size(); ++i) {
	unsigned short channel = hits[i]->Wire()->RawDigit()->Channel();
	
	geo->ChannelToWire(channel, t, p, w);

	// check for correct plane, wire and tpc
	if(p != plane || w != wire || t != rawopt->fTPC) continue;

	hstart.push_back(hits[i]->StartTime());
	hend.push_back(hits[i]->EndTime());
    hitamplitudes.push_back(hits[i]->Charge(true));
    hpeaktimes.push_back(hits[i]->PeakTime());

      }//end loop over reco hits
    }//end loop over HitFinding modules

  }

  //......................................................................
  void RecoBaseDrawer::FillQHisto(const art::Event& evt,
				  unsigned int      plane,
				  TH1F*             histo)
  {

    // Check if we're supposed to draw raw hits at all
    art::ServiceHandle<evd::RawDrawingOptions>   rawopt;
    art::ServiceHandle<evd::RecoDrawingOptions> drawopt;
    if (rawopt->fDrawRawDataOrCalibWires==0) return;

    art::ServiceHandle<geo::Geometry> geo;
  
    unsigned int p = 0;
    unsigned int w  = 0;
    unsigned int t  = 0;

    for (unsigned int imod=0; imod<drawopt->fWireLabels.size(); ++imod) {
      std::string const which = drawopt->fWireLabels[imod];

      art::PtrVector<recob::Wire> wires;
      this->GetWires(evt, which, wires);

      for (unsigned int i=0; i<wires.size(); ++i) {

	unsigned short channel = wires[i]->RawDigit()->Channel();

	geo->ChannelToWire(channel, t, p, w);

	// check for correct plane and tpc
	if(p != plane || t != rawopt->fTPC) continue;
	
	for(size_t s = 0; s < wires[i]->fSignal.size(); ++s)
	  histo->Fill(wires[i]->fSignal[s]);

      }//end loop over raw hits
    }//end loop over Wire modules

  }

}// namespace
////////////////////////////////////////////////////////////////////////
