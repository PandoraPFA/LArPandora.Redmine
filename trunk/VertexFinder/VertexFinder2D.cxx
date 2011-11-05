////////////////////////////////////////////////////////////////////////
//
// VertexFinder2D class
//
// tjyang@fnal.gov
//
// This algorithm is designed to reconstruct the vertices using the
// 2D cluster information
// 
// This is Preliminary Work and needs modifications
// ////////////////////////////////////////////////////////////////////////


#include <iostream>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "VertexFinder/VertexFinder2D.h"
#include <iomanip>
#include <ios>
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <TVector3.h>
#include <vector>

#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"

#include "TH1D.h"
#include "TVectorD.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"

namespace vertex{

//-----------------------------------------------------------------------------
  VertexFinder2D::VertexFinder2D(fhicl::ParameterSet const& pset)
  {  
    this->reconfigure(pset);    
    produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::EndPoint2D> >();
  }
//-----------------------------------------------------------------------------
  VertexFinder2D::~VertexFinder2D()
  {
  }

  //---------------------------------------------------------------------------
  void VertexFinder2D::reconfigure(fhicl::ParameterSet const& p) 
  {
    fClusterModuleLabel  = p.get< std::string >("ClusterModuleLabel");
    return;
  }
  //-------------------------------------------------------------------------
  void VertexFinder2D::beginJob(){
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    dtIC = tfs->make<TH1D>("dtIC","It0-Ct0",100,-5,5);
    dtIC->Sumw2();
  }

// //-----------------------------------------------------------------------------
  void VertexFinder2D::produce(art::Event& evt)
  {

    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;

    // define TPC parameters
    TString tpcName = geom->GetLArTPCVolumeName();
    
    //  double YC =  (m_TPCHalfZ-5.)*2.; // TPC height in cm
    double YC =  (geom->DetHalfHeight())*2.; // *ArgoNeuT* TPC active-volume height in cm
    double Angle = geom->Plane(1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to the vertical direction
    // Parameters temporary defined here, but possibly to be retrieved somewhere in the code
    double timetick = 0.198;    //time sample in us
    double presamplings = 75.;
    //const double wireShift=50.; // half the number of wires from the Induction(Collection) plane intersecting with a wire from the Collection(Induction) plane.
    //double plane_pitch = geom->PlanePitch(0,1);   //wire plane pitch in cm 
    double wire_pitch = geom->WirePitch(0,1,0);    //wire pitch in cm
    double Efield_drift = 0.5;  // Electric Field in the drift region in kV/cm
    //double Efield_SI = 0.7;     // Electric Field between Shield and Induction planes in kV/cm
    //double Efield_IC = 0.9;     // Electric Field between Induction and Collection planes in kV/cm
    double Temperature = 90.;  // LAr Temperature in K
    
    double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature); //drift velocity in the drift 
    //region (cm/us)
    //double driftvelocity_SI = larprop->DriftVelocity(Efield_SI,Temperature); //drift velocity between shield 
    //and induction (cm/us)
    //double driftvelocity_IC = larprop->DriftVelocity(Efield_IC,Temperature); //drift velocity between induction 
    //and collection (cm/us)
    double timepitch = driftvelocity*timetick;                               //time sample (cm) 
    //double tSI = plane_pitch/driftvelocity_SI/timetick;                      //drift time between Shield and 
    //Collection planes (time samples)
    //double tIC = plane_pitch/driftvelocity_IC/timetick;                      //drift time between Induction and 
    //Collection planes (time samples)
    
    
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    evt.getByLabel(fClusterModuleLabel,clusterListHandle);

    art::PtrVector<recob::Cluster> clusters;
    for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii){
      art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      clusters.push_back(clusterHolder);
    }


    //Point to a collection of vertices to output.
    std::auto_ptr<std::vector<recob::Vertex> > vcol(new std::vector<recob::Vertex>);          //3D vertex
    std::auto_ptr<std::vector<recob::EndPoint2D> >epcol(new std::vector<recob::EndPoint2D>);  //2D vertex

    for(size_t tpc = 0; tpc < geom->NTPC(); ++tpc){
      int nplanes = geom->Nplanes(tpc);
    
      std::vector<int> Cls[nplanes]; //index to clusters in each view
      std::vector<double> dtdwstart;
      
      //loop over clusters
      for(unsigned int iclu=0; iclu<clusters.size();++iclu){
	
	switch(clusters[iclu]->View()){
	  
	case geo::kU :
	  Cls[0].push_back(iclu);
	  break;
	case geo::kV :
	  Cls[1].push_back(iclu);
	  break;
	case geo::kW :
	  Cls[2].push_back(iclu);
	  break;
	default :
	  break;
	}
	
	unsigned int channel;
	unsigned int wire;
	unsigned int t;
	unsigned int plane;
	double time;
	
	std::vector<double> wires;
	std::vector<double> times;
	
	art::PtrVector<recob::Hit> hit;
	hit = clusters[iclu]->Hits();
	int n = 0;
	for (unsigned int i = 0; i<hit.size(); i++){
	  time = hit[i]->PeakTime();
	  channel = hit[i]->Channel();
	  geom->ChannelToWire(channel,t,plane,wire);
	  wires.push_back(wire);
	  times.push_back(time);
	  n++;
	  //if (n==20) break;
	}
	if (n>=2){
	  TGraph *the2Dtrack = new TGraph(n,&wires[0],&times[0]);           
	  try
	    {
	      the2Dtrack->Fit("pol1","Q");
	      TF1 *pol1=(TF1*) the2Dtrack->GetFunction("pol1");
	      double par[2];
	      pol1->GetParameters(par);
	      //std::cout<<iclu<<" "<<par[1]<<" "<<clusters[iclu]->dTdW()<<std::endl;
	      dtdwstart.push_back(par[1]);
	    }
	  catch(...){std::cout<<"Fitter failed"<<std::endl;
	    delete the2Dtrack;
	    dtdwstart.push_back(clusters[iclu]->dTdW());
	    continue;}
	  delete the2Dtrack;
	}
	else dtdwstart.push_back(clusters[iclu]->dTdW());
      }
      
      std::vector<int> cluvtx[nplanes];
      std::vector<double> vtx_w;
      std::vector<double> vtx_t;
      
      for (int i = 0; i<nplanes; i++){
	if (Cls[i].size()>=1){//at least one cluster
	  //find the longest two clusters
	  int c1 = -1;
	  int c2 = -1;
	  double ww0 = -999;
	  double wb1 = -999;
	  double we1 = -999;
	  double wb2 = -999;
	  double we2 = -999;
	  double tt0 = -999;
	  double tt1 = -999;
	  double tt2 = -999;
	  double dtdw1 = -999;
	  double dtdw2 = -999;
	  double lclu1 = -999;
	  double lclu2 = -999;
	  for (unsigned j = 0; j<Cls[i].size(); j++){
	    double lclu = sqrt(pow((clusters[Cls[i][j]]->StartPos()[0]-clusters[Cls[i][j]]->EndPos()[0])*13.5,2)+pow(clusters[Cls[i][j]]->StartPos()[1]-clusters[Cls[i][j]]->EndPos()[1],2));
	    bool rev = false;
	    bool deltaraylike = false;
	    bool enoughhits = false;
	    if (c1!=-1){
	      double wb = clusters[Cls[i][j]]->StartPos()[0];
	      double we = clusters[Cls[i][j]]->EndPos()[0];
	      double tt = clusters[Cls[i][j]]->StartPos()[1];
	      double dtdw = dtdwstart[Cls[i][j]];
	      int nhits = clusters[Cls[i][j]]->Hits().size();
	      tt0 = (dtdw1*dtdw*(wb1-wb)+dtdw1*tt-dtdw*tt1)/(dtdw1-dtdw);
	      ww0 = (tt-tt1+dtdw1*wb1-dtdw*wb)/(dtdw1-dtdw);
	      if (fabs(wb1-ww0)>fabs(we1-ww0)) rev = true;//reverse cluster dir
	      if ((!rev&&ww0>wb1+15)||(rev&&ww0<we1-15)) deltaraylike = true;
	      if (((!rev&&ww0>wb1+10)||(rev&&ww0<we1-10))&&nhits<5) deltaraylike = true;
	      if (wb>wb1+20&&nhits<10) deltaraylike = true;
	      if (wb>wb1+50&&nhits<20) deltaraylike = true;
	      if (wb>wb1+8&&TMath::Abs(dtdw1-dtdw)<0.15) deltaraylike = true;
	      if (TMath::Abs(wb-wb1)>30&&TMath::Abs(we-we1)>30) deltaraylike = true;
	      //std::cout<<rev<<" "<<ww0<<" "<<wb1<<" "<<wb<<" "<<we<<" "<<we1<<" "<<nhits<<std::endl;
	      //make sure there are enough hits in the cluster
	      //at leaset 2 hits if goes horizentally, at leaset 4 hits if goes vertically
	      double alpha = TMath::ATan(dtdw);
	      if (nhits>=int(2+3*(1-TMath::Abs(TMath::Cos(alpha))))) enoughhits = true;
	      if (nhits<5&&(ww0<wb1-20||ww0>we1+20)) enoughhits = false;
	      //std::cout<<nhits<<" "<<int(2+3*(1-TMath::Abs(TMath::Cos(alpha))))<<" "<<alpha<<std::endl;
	      //std::cout<<Cls[i][j]<<" "<<deltaraylike<<" "<<enoughhits<<" "<<dtdw<<" "<<dtdw1<<std::endl;
	    }
	    //do not replace the second cluster if the 3rd cluster is not consistent with the existing 2
	    bool replace = true;
	    if (c1!=-1&&c2!=-1){
	      double wb = clusters[Cls[i][j]]->StartPos()[0];
	      double we = clusters[Cls[i][j]]->EndPos()[0];
	      tt0 = (dtdw1*dtdw2*(wb1-wb2)+dtdw1*tt2-dtdw2*tt1)/(dtdw1-dtdw2);
	      ww0 = (tt2-tt1+dtdw1*wb1-dtdw2*wb2)/(dtdw1-dtdw2);
	      if ((TMath::Abs(ww0-wb1)<10||TMath::Abs(ww0-we1)<10)&&
		  (TMath::Abs(ww0-wb2)<10||TMath::Abs(ww0-we2)<10)){
		if (TMath::Abs(ww0-wb)>15&&TMath::Abs(ww0-we)>15) replace = false;
	      }
	      //std::cout<<c1<<" "<<c2<<" "<<ww0<<" "<<wb1<<" "<<wb2<<" "<<wb<<" "<<we<<std::endl;
	    }
	    if (lclu1<lclu){
	      if (c1!=-1&&!deltaraylike&&enoughhits){
		lclu2 = lclu1;
		c2 = c1;
		wb2 = wb1;
		we2 = we1;
		tt2 = tt1;
		dtdw2 = dtdw1;
	      }
	      lclu1 = lclu;
	      c1 = Cls[i][j];
	      wb1 = clusters[Cls[i][j]]->StartPos()[0];
	      we1 = clusters[Cls[i][j]]->EndPos()[0];
	      tt1 = clusters[Cls[i][j]]->StartPos()[1];
	      dtdw1 = dtdwstart[Cls[i][j]];
	    }
	    else if (lclu2<lclu){
	      if (!deltaraylike&&enoughhits&&replace){
		lclu2 = lclu;
		c2 = Cls[i][j];
		wb2 = clusters[Cls[i][j]]->StartPos()[0];
		we2 = clusters[Cls[i][j]]->EndPos()[0];
		tt2 = clusters[Cls[i][j]]->StartPos()[1];
		dtdw2 = dtdwstart[Cls[i][j]];
	      }
	    }
	  }
	  if (c1!=-1&&c2!=-1){
	    cluvtx[i].push_back(c1);
	    cluvtx[i].push_back(c2);
	    //std::cout<<i<<" "<<c1<<" "<<c2<<std::endl;
	    double w1 = clusters[c1]->StartPos()[0];
	    double t1 = clusters[c1]->StartPos()[1];
	    double k1 = dtdwstart[c1];
	    double w2 = clusters[c2]->StartPos()[0];
	    double t2 = clusters[c2]->StartPos()[1];
	    double k2 = dtdwstart[c2];
	    //std::cout<<k1<<" "<<k2<<" "<<k1-k2<<std::endl;
	    //calculate the vertex
	    if (fabs(k1-k2)<1e-10) continue;
	    double t0 = (k1*k2*(w1-w2)+k1*t2-k2*t1)/(k1-k2);
	    double w0 = (t2-t1+k1*w1-k2*w2)/(k1-k2);
	    vtx_w.push_back(w0);
	    vtx_t.push_back(t0);
	  }
	  else if (Cls[i].size()>=1){
	    if (c1!=-1){
	      cluvtx[i].push_back(c1);
	      vtx_w.push_back(clusters[c1]->StartPos()[0]);
	      vtx_t.push_back(clusters[c1]->StartPos()[1]);
	    }
	    else{
	      cluvtx[i].push_back(Cls[i][0]);
	      vtx_w.push_back(clusters[Cls[i][0]]->StartPos()[0]);
	      vtx_t.push_back(clusters[Cls[i][0]]->StartPos()[1]);
	    }
	  }
	  //save 2D vertex
	  recob::EndPoint2D vertex;
	  vertex.SetWireNum(int(vtx_w.back()));
	  vertex.SetDriftTime(vtx_t.back());
	  vertex.SetView(clusters[Cls[i][0]]->View());
	  epcol->push_back(vertex);
	}
	else {//no cluster found
	  vtx_w.push_back(-1);
	  vtx_t.push_back(-1);
	}
      }
      //std::cout<<vtx_w[0]<<" "<<vtx_t[0]<<" "<<vtx_w[1]<<" "<<vtx_t[1]<<std::endl;
      
      Double_t vtxcoord[3];
      if (Cls[0].size()>0&&Cls[1].size()>0){//ignore w view
	double Iw0 = (vtx_w[0]+3.95)*wire_pitch;
	double Cw0 = (vtx_w[1]+1.84)*wire_pitch;
	
	double It0 = vtx_t[0] - presamplings;
	It0 *= timepitch;
	double Ct0 = vtx_t[1] - presamplings ;
	Ct0 *= timepitch;
	vtxcoord[0] = Ct0;
	vtxcoord[1] = (Cw0-Iw0)/(2.*TMath::Sin(Angle));
	vtxcoord[2] = (Cw0+Iw0)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle);
	
	double yy,zz;       
	if(vtx_w[0]>=0&&vtx_w[0]<=239&&vtx_w[1]>=0&&vtx_w[1]<=239){
	  if(geom->ChannelsIntersect(geom->PlaneWireToChannel(0,(int)((Iw0/wire_pitch)-3.95)),      
				     geom->PlaneWireToChannel(1,(int)((Cw0/wire_pitch)-1.84)),
				     yy,zz)){
	      //channelsintersect provides a slightly more accurate set of y and z coordinates. use channelsintersect in case the wires in question do cross.
	      vtxcoord[1] = yy;
	      vtxcoord[2] = zz;
	  }
	  else{
	    vtxcoord[0] = -99999;
	    vtxcoord[1] = -99999;
	    vtxcoord[2] = -99999;
	  }
	}	
	dtIC->Fill(It0-Ct0);
      }
      else{
	vtxcoord[0] = -99999;
	vtxcoord[1] = -99999;
	vtxcoord[2] = -99999;
      }
      
      //std::cout<<vtxcoord[0]<<" "<<vtxcoord[1]<<" "<<vtxcoord[2]<<std::endl;
      //need to implement this
      art::PtrVector<recob::Track> vTracks_vec;
      art::PtrVector<recob::Shower> vShowers_vec;
      
      recob::Vertex the3Dvertex(vTracks_vec, vShowers_vec, vtxcoord);
      vcol->push_back(the3Dvertex);

    }//end loop over tpc

    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "VertexFinder2D Summary:";
    for(unsigned int i = 0; i<epcol->size(); ++i) mf::LogVerbatim("Summary") << epcol->at(i) ;
    for(unsigned int i = 0; i<vcol->size(); ++i) mf::LogVerbatim("Summary") << vcol->at(i) ;

    evt.put(epcol);
    evt.put(vcol);
    
  } // end of produce
} // end of vertex namespace

// //-----------------------------------------------------------------------------
