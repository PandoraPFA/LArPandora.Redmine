////////////////////////////////////////////////////////////////////////
//
// PrimaryVertexFinder class
//
// saima@ksu.edu
//
// This algorithm is designed to reconstruct the vertices by finding the
// distance of closest approach between the tracks
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

#include "VertexFinder/PrimaryVertexFinder.h"

#include <iomanip>
#include <ios>
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "TMath.h"
#include <TVector3.h>
#include <vector>



#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TGeoManager.h"


// // //------------------------------------------------------------------------------ // not needed for now
// static bool sp_sort_z(const recob::SpacePoint& sp1, const recob::SpacePoint& sp2)
// {
//   const double* xyz1 = sp1.XYZ();
//   const double* xyz2 = sp2.XYZ();
//   return xyz1[2] < xyz2[2];
// }
// // //------------------------------------------------------------------------------ // not needed for now
// static bool sp_sort_zz(const std::vector<recob::SpacePoint>& sp_vec1, const std::vector<recob::SpacePoint>& sp_vec2)
// {
//   const double* xyz1 = sp_vec1[0].XYZ();
//   const double* xyz2 = sp_vec2[0].XYZ();
//   return xyz1[2] < xyz2[2];
// }
// // //------------------------------------------------------------------------------
bool sort_pred2(const std::pair<art::Ptr<recob::Track>,double>& left, const std::pair<art::Ptr<recob::Track>,double>& right)
{
  return left.second < right.second;
}

namespace vertex{

//-----------------------------------------------------------------------------
  PrimaryVertexFinder::PrimaryVertexFinder(fhicl::ParameterSet const& pset)
  {  
    this->reconfigure(pset);    
    produces< std::vector<recob::Vertex> >();
  }
//-----------------------------------------------------------------------------
  PrimaryVertexFinder::~PrimaryVertexFinder()
  {
  }

  //---------------------------------------------------------------------------
  void PrimaryVertexFinder::reconfigure(fhicl::ParameterSet const& p) 
  {
    fTrackModuleLabel  = p.get< std::string >("TrackModuleLabel");
    fVertexWindow      = p.get<double     >  ("VertexWindow");
    return;
  }
  //-------------------------------------------------------------------------
  void PrimaryVertexFinder::beginJob(){
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    
    //    fNoVertices= tfs->make<TH2F>("fNoVertices", ";Event No; No of vertices", 100,0, 100, 30, 0, 30);
     fNoTracks= tfs->make<TH2F>("fNoTracks", ";Event No; No of Tracks", 10,0, 10, 10, 0, 10);
     fLength_1stTrack = tfs->make<TH1F>("fLength_Track1", "Muon Track Length", 100,0,100);
     fLength_2ndTrack = tfs->make<TH1F>("fLength_Track2", "2nd Track Length", 100,0,100);
     fLength_3rdTrack = tfs->make<TH1F>("fLength_Track3", "3rd Track Length", 100,0,100);
     fLength_4thTrack = tfs->make<TH1F>("fLength_Track4", "4th Track Length", 100,0,100);
     fLength_5thTrack = tfs->make<TH1F>("fLength_Track5", "5th Track Length", 100,0,100);
  }

// //-----------------------------------------------------------------------------
  void PrimaryVertexFinder::produce(art::Event& evt)
  {

    std::cout<<std::endl;
    std::cout<<"------------------------------------------------------------------------------"<<std::endl;
    //   std::cout << "run    : " << evt.Header().Run() << std::endl;
    //   std::cout << "subrun : " << evt.Header().Subrun() << std::endl;
    //std::cout << "event  : " << evt.Header().Event() << std::endl;
    
    std::cout << "event  : " << evt.id().event() << std::endl;
    
    
    art::ServiceHandle<geo::Geometry> geom;
    
    //std::cout << "I am in Primary vertex finder " << std::endl;
        
    art::Handle< std::vector<recob::Track> > trackListHandle;
    evt.getByLabel(fTrackModuleLabel,trackListHandle);
    
    //Point to a collection of vertices to output.
    std::auto_ptr<std::vector<recob::Vertex> > vcol(new std::vector<recob::Vertex>);
    

    art::PtrVector<recob::Track> trkIn;
    for(unsigned int ii = 0; ii < trackListHandle->size(); ++ii)
      {
	art::Ptr<recob::Track> track(trackListHandle, ii);
	trkIn.push_back(track);
      }

    std::cout << "number of tracks in this event = " << trkIn.size() << std::endl;
    fNoTracks->Fill(evt.id().event(),trkIn.size());

    std::vector<recob::SpacePoint> spacepoints;  // space points associated to each track 
    std::vector<recob::SpacePoint> startpoints_vec; // first space point of each track

    std::vector<double> start;
    std::vector<double> end;
    double startcos[3]={0,0,0};
    double endcos[3]={0,0,0};


    std::vector <TVector3> startvec;
    TVector3 startXYZ;
    
    std::vector <TVector3> endvec;
    TVector3 endXYZ;

    std::vector <TVector3> dircosvec;
    TVector3 dircosXYZ;

    std::vector< std::pair<art::Ptr<recob::Track>, double> > trackpair;
    
    for(unsigned int i = 0; i<trkIn.size(); ++i){
      trkIn[i]->Extent(start, end);
      startXYZ.SetXYZ(start[0],start[1],start[2]);
      endXYZ.SetXYZ(end[0],end[1],end[2]);


      double length = (endXYZ-startXYZ).Mag();// (endvec[i]-startvec[i]).Mag();
      //std::cout << "Track length calculated = " << length << std::endl;
      trackpair.push_back(std::pair<art::Ptr<recob::Track>,double>(trkIn[i],length));
    }
    
    for(unsigned int i = 0; i<trackpair.size(); ++i){
      std::cout << "track id is  = " << (trackpair[i].first)->ID() << " track length = " << (trackpair[i].second) << std::endl;
    }
    
    std::sort(trackpair.rbegin(), trackpair.rend(),sort_pred2);
    
    std::cout << "AFTER SORTING " << std::endl;
    for(unsigned int i = 0; i<trackpair.size(); ++i){
      std::cout << "track id is  = " << (trackpair[i].first)->ID() << " track length = " << (trackpair[i].second) << std::endl;
    }
    
    if(trackpair.size()>0)
    fLength_1stTrack->Fill(trackpair[0].second);

    if(trackpair.size()>1)
    fLength_2ndTrack->Fill(trackpair[1].second);

    if(trackpair.size()>2)
    fLength_3rdTrack->Fill(trackpair[2].second);

    if(trackpair.size()>3)
    fLength_4thTrack->Fill(trackpair[3].second);

    if(trackpair.size()>4)
    fLength_5thTrack->Fill(trackpair[4].second);

    for(unsigned int j=0; j<trackpair.size();++j) { //loop over tracks
      spacepoints=trackpair[j].first->SpacePoints();

     
      trackpair[j].first->Extent(start, end);
      trackpair[j].first->Direction(startcos, endcos);

      startXYZ.SetXYZ(start[0],start[1],start[2]);
      endXYZ.SetXYZ(end[0],end[1],end[2]);
      dircosXYZ.SetXYZ(startcos[0],startcos[1],startcos[2]);
 
      startvec.push_back(startXYZ);
      endvec.push_back(endXYZ);
      dircosvec.push_back(dircosXYZ);
      
      start.clear();
      end.clear();
      startcos[3] = 0.;
      

      std::cout<<"PrimaryVertexFinder got "<< spacepoints.size() <<" 3D spacepoint(s) from Track3Dreco.cxx"<<std::endl;
            
      for(unsigned int i=0; i<1; ++i) { // save the first SpacePoint of each Track... from now the SpacePoint ID represents the Track ID!!
	spacepoints[i].SetID(startpoints_vec.size());
	startpoints_vec.push_back(spacepoints[i]);
      }
    }// loop over tracks

    for(unsigned int i=0; i<startvec.size(); i++){ //trackpair.size()
      std::cout << "Tvector3 start point SORTED = " << std::endl; 
      startvec[i].Print();
    }
    for(unsigned int i=0; i<dircosvec.size(); i++){ //trackpair.size()
      std::cout << "Tvector3 dir cos SORTED = " << std::endl; 
      dircosvec[i].Print();
    }

    std::vector<std::vector<int> > vertex_collection_int;
    std::vector <std::vector <TVector3> > vertexcand_vec;

    for (unsigned int i=0; i<trackpair.size(); ++i){
      for (unsigned int j=i+1; j<trackpair.size(); ++j){
	std::cout << "distance between " << i << " and " << j << " = " << StartPointSeperation(startpoints_vec[i], startpoints_vec[j]) << std::endl;
	double GAMMA = gammavalue(startvec[i], startvec[j], dircosvec[i], dircosvec[j]);
	double ALPHA = alphavalue(GAMMA, startvec[i], startvec[j], dircosvec[i], dircosvec[j]);
	double MINDIST = MinDist(ALPHA, GAMMA, startvec[i], startvec[j], dircosvec[i], dircosvec[j]);
	std::cout << "alpha = " << ALPHA << " gamma = " << GAMMA << " MINIMUM DISTANCE = " << MINDIST << std::endl;
	TVector3 TRACK1POINT = PointOnExtendedTrack(ALPHA, startvec[i], dircosvec[i]);
	TVector3 TRACK2POINT = PointOnExtendedTrack(GAMMA, startvec[j], dircosvec[j]);
	std::cout << "POINTS ON THE TRACKS ARE:: " << std::endl;
	TRACK1POINT.Print();
	TRACK2POINT.Print();
	std::cout << std::endl;

	//if(StartPointSeperation(startpoints_vec[i], startpoints_vec[j])<fVertexWindow){ ///// correct this
	//if(MINDIST<2 && trackpair[i].second >30 && trackpair[j].second >30){
	  if(MINDIST<fVertexWindow && ((TRACK1POINT-startvec[i]).Mag())<fVertexWindow){

	  if((!IsInVertexCollection(i, vertex_collection_int)) && (!IsInVertexCollection(j, vertex_collection_int))){
	    std::vector<int> newvertex_int;
	    std::vector <TVector3> vertexcand;
	    newvertex_int.push_back(i);
	    newvertex_int.push_back(j);
	    vertex_collection_int.push_back(newvertex_int);
	    //newvertex.clear();
	    vertexcand.push_back(TRACK1POINT);
	    vertexcand.push_back(TRACK2POINT);
	    vertexcand_vec.push_back(vertexcand);
	  }
	  else
	    {
	      int index = IndexInVertexCollection(i, j, vertex_collection_int);
	      //std::cout << "index where a new vertex will be added = " << index << std::endl;
	      if(!IsInNewVertex(i, vertex_collection_int[index])){
	      vertex_collection_int[index].push_back(i);
	      vertexcand_vec[index].push_back(TRACK1POINT); //need to fix for delta rays
	      }
	     if(!IsInNewVertex(j, vertex_collection_int[index])){
	      vertex_collection_int[index].push_back(j);
	      vertexcand_vec[index].push_back(TRACK2POINT); //need to fix for delta rays
	     }
	    }
	}
      }
    }

    
    //now add the unmatched track IDs to the collection
    for(unsigned int i=0; i<trackpair.size(); i++){
      if(!IsInVertexCollection(i, vertex_collection_int)){
	//if(trackpair[i].second>30){
	std::vector<int> temp;
	std::vector <TVector3> temp1;
	temp.push_back(i);
	temp1.push_back(startvec[i]);	
	vertex_collection_int.push_back(temp);
	vertexcand_vec.push_back(temp1);
	//}
      }
    }

    art::PtrVector<recob::Track> vTracks_vec;
    art::PtrVector<recob::Shower> vShowers_vec;

    for(unsigned int i=0; i<vertex_collection_int.size(); i++){
      double x = 0.;
      double y = 0.;
      double z = 0.;
      int elemsize = 0.;
      for(std::vector<int>::iterator itr = vertex_collection_int[i].begin(); itr < vertex_collection_int[i].end(); ++itr){
	std::cout << "vector elements at index " << i << " are " << *itr <<std::endl;
	//	std::cout << "the number of elements at index " << i << " are " << vertex_collection_int[i].size() <<std::endl;
	std::cout << "track original ID = " << (trackpair[*itr].first)->ID() << std::endl;
	vTracks_vec.push_back(trackpair[*itr].first);
      }
      std::cout << "------------" <<std::endl;

      
      for(std::vector<TVector3>::iterator itr = vertexcand_vec[i].begin(); itr < vertexcand_vec[i].end(); ++itr){
	//calculate sum of x, y and z of a vertex
	x += (*itr).X();
	y += (*itr).Y();
	z += (*itr).Z();
	elemsize = vertexcand_vec[i].size();
      }

      double avgx = x/elemsize;
      double avgy = y/elemsize;
      double avgz = z/elemsize;
      
      Double_t vtxcoord[3];
	  vtxcoord[0] = avgx;
	  vtxcoord[1] = avgy;
	  vtxcoord[2] = avgz;    
            
	  recob::Vertex the3Dvertex(vTracks_vec, vShowers_vec, vtxcoord);
	  vcol->push_back(the3Dvertex);
	  vTracks_vec.clear();

    }

    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "PrimaryVertexFinder Summary:";
    for(unsigned int i = 0; i<vcol->size(); ++i) mf::LogVerbatim("Summary") << vcol->at(i) ;

    evt.put(vcol);
    
  } // end of produce
} // end of vertex namespace

// //-----------------------------------------------------------------------------
double vertex::PrimaryVertexFinder::StartPointSeperation(recob::SpacePoint sp1, recob::SpacePoint sp2)
{
  double x= (sp2.XYZ()[0])-(sp1.XYZ()[0]);
  double y= (sp2.XYZ()[1])-(sp1.XYZ()[1]);
  double z= (sp2.XYZ()[2])-(sp1.XYZ()[2]);
  double distance = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
  return distance;
}
// //---------------------------------------------------------------------------------
bool vertex::PrimaryVertexFinder::IsInVertexCollection(int a, std::vector<std::vector<int> > vertex_collection)
{
  int flag = 0;
  
  for(unsigned int i = 0; i < vertex_collection.size() ; i++){
    for(std::vector<int>::iterator itr = vertex_collection[i].begin(); itr < vertex_collection[i].end(); ++itr){
      if (a == *itr){
	flag = 1;
	break;
      }
    }
  }
  if(flag==1)
    return true;
  return false; 
}
// //------------------------------------------------------------------------------
int vertex::PrimaryVertexFinder::IndexInVertexCollection(int a, int b, std::vector<std::vector<int> > vertex_collection)
{
  int index;
  for(unsigned int i = 0; i < vertex_collection.size() ; i++){
    for(std::vector<int>::iterator itr = vertex_collection[i].begin(); itr < vertex_collection[i].end(); ++itr){
      if (a == *itr || b == *itr)
	index = i; 
    }
  }
  return index;
}
// //------------------------------------------------------------------------------
bool vertex::PrimaryVertexFinder::IsInNewVertex(int a, std::vector<int> newvertex)
{
  int flag = 0;
  for(unsigned int i = 0; i < newvertex.size() ; i++){
    if (a == newvertex[i]){
      flag = 1;
      break;
    }
  }
  
  if(flag==1)
    return true;
  return false; 
}
// //------------------------------------------------------------------------------
double vertex::PrimaryVertexFinder::gammavalue(TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2)
{
  double gamma = ((startpoint1*dircos2)-(startpoint2*dircos2)+((dircos1*dircos2)*(startpoint2*dircos1))-((dircos1*dircos2)*(startpoint1*dircos1)))/(1-((dircos1*dircos2)*(dircos1*dircos2)));

  return gamma;
}
// //------------------------------------------------------------------------------
double vertex::PrimaryVertexFinder::alphavalue(double gamma, TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2)
{
  double alpha = (gamma*(dircos1*dircos2)) + (startpoint2*dircos1) - (startpoint1*dircos1);

  return alpha;
}
// //------------------------------------------------------------------------------
double vertex::PrimaryVertexFinder::MinDist(double alpha, double gamma, TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2)
{
  TVector3 mindis_vector = startpoint1 - startpoint2 + alpha*dircos1 - gamma*dircos2;
  double mindis = mindis_vector.Mag();
  return mindis;
}
// //------------------------------------------------------------------------------
TVector3 vertex::PrimaryVertexFinder::PointOnExtendedTrack(double alphagamma, TVector3 startpoint,  TVector3 dircos)
{
  TVector3 PointOnExtendedTrack = startpoint + (alphagamma * dircos);
  return PointOnExtendedTrack;
}
// //------------------------------------------------------------------------------

