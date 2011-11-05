////////////////////////////////////////////////////////////////////////
//
// ShowerFinder class
//
// roxanne.guenette@yale.edu
//
//  This algorithm is designed to find all the showers in an event.
////////////////////////////////////////////////////////////////////////

#include "ShowerFinder/ShowerFinder.h"

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <sstream>
#include <math.h>
#include <algorithm>

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

// LArSoft Includes
#include "RawData/RawDigit.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"


// ROOT 
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"


namespace shwf{

  //-------------------------------------------------
  ShowerFinder::ShowerFinder(fhicl::ParameterSet const& pset) :
    fVertexModuleLabel        (pset.get<std::string > ("VertexModuleLabel")),
    fClusterModuleLabel       (pset.get<std::string > ("ClusterModuleLabel")),
    fHoughLineModuleLabel     (pset.get<std::string > ("HoughLineModuleLabel")),
    fVertexStrengthModuleLabel(pset.get<std::string > ("VertexStrengthModuleLabel")),
    fRcone                    (pset.get<double      > ("Rcone")), //Radius of the cone
    fLcone                    (pset.get<double      > ("Lcone"))  //Length (perpendicular to the base) of the cone    
  {
    produces< std::vector<recob::Shower> >();
  }
  
  
  //-------------------------------------------------

  ShowerFinder::~ShowerFinder()
  {
  }
  
  //
  //-------------------------------------------------
  void ShowerFinder::produce(art::Event& evt)
  { 
    
    
    //////////////////////////////////////////////////////
    // Make a std::auto_ptr<> for the thing you want to put into the event
    // because that handles the memory management for you
    //////////////////////////////////////////////////////
   
    //std::auto_ptr<std::vector<recob::Shower> > show_ID(new std::vector<recob::Shower>);
    //std::auto_ptr<std::vector<recob::Shower> > show_vcoord(new std::vector<recob::Shower>);
    std::auto_ptr<std::vector<recob::Shower> > showercol(new std::vector<recob::Shower>);


    // Read in the vertex List object(s).
    art::Handle< std::vector<recob::EndPoint2D> > vertexListHandle;
    evt.getByLabel(fVertexModuleLabel,vertexListHandle);
    // Read in the hough line List object(s).
    art::Handle< std::vector<recob::Cluster> > houghListHandle;
    evt.getByLabel(fHoughLineModuleLabel,houghListHandle);
    // Read in the cluster List object(s).
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    evt.getByLabel(fClusterModuleLabel,clusterListHandle);
    // Read in the vertex Strength List object(s).
    art::Handle< std::vector<recob::EndPoint2D> > vertexStrengthListHandle;
    evt.getByLabel(fVertexStrengthModuleLabel,vertexStrengthListHandle);
    
    // art::PtrVector<recob::Cluster> clust;
    //for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii)
    //{
    // art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
    // clust.push_back(clusterHolder);
    //}

    art::PtrVector<recob::Cluster> protoShowers; //vector of clusters associated to a cone

    art::PtrVector<recob::Hit> clusterhits; //hits in the cluster
    art::PtrVector<recob::Hit> hlhits; //hits in the hough Lines


    art::ServiceHandle<geo::Geometry> geom;
    
    unsigned int channel,plane,wire,tpc;
    
    //This vector will contain all strong and strongest vertices
    art::PtrVector<recob::EndPoint2D> vertSel;
    
    //This loop is going over all the vertices in the event 
    //and is interested in ONLY strong and strongest vertices.
    std::cout << "Vertex STRENGTH list size = " << vertexStrengthListHandle->size() << " AND vertices:" << vertexListHandle->size()<< std::endl;
    std::cout << "CLUSTER list size = " << clusterListHandle->size() << " AND Hough: :" << houghListHandle->size()<< std::endl;
    
    for(unsigned int iv = 0; iv < vertexListHandle->size(); ++iv)
      {
	art::Ptr<recob::EndPoint2D> vertex(vertexListHandle, iv);
	//std::cout << "Vertex " << iv << " :  str = " << vertex->ID() << std::endl;
	//if(vertex->Strength() == 4 || vertex->Strength() == 3){
	if(vertex->ID() == 1 || vertex->Strength() == 3){ //only use Strongest and strong
	  vertSel.push_back(vertex);
	}
	else continue;
      }
    
    //Definition of the geometry of the cone (which is basically a triangle)
    double scan_angle = 0; //angle of the scan steps
    double xa_cone = 0; // x coordinate of the cone's apex (wire number)
    double ya_cone = 0; // y coordinate of the cone's apex (drift time)
    double x1_cone = 0; // x coordinate of the cone's top right point (wire number)
    double y1_cone = 0; // y coordinate of the cone's top right point (drift time)
    double x2_cone = 0; // x coordinate of the cone's top left point (wire number)
    double y2_cone = 0; // y coordinate of the cone's top left point  (drift time)
    
    double fScone = sqrt( (fRcone*fRcone) + (fLcone*fLcone)); //The length of the side of the cone
    
    double cone_angle = (TMath::ATan(fRcone / fLcone)) / 2.0; // Opening angle of the cone (defined from input parameters)
    std::cout << "Cone Opening Angle: " << (180.0*cone_angle)/TMath::Pi() << std::endl;
    double compl_angle = 0;  
    
    unsigned int n_scan =1 + (int)(TMath::Pi() / (2.0*cone_angle)); 
    std::cout << "N scan: " << n_scan << std::endl;
    
    double x_hit = 0; //x coordinate of hit
    double y_hit = 0; //y coordinate of hit
    
    int hits_cluster_counter = 0; //count the number of hits in a cluster that is inside a cone
    //int hits_cluster_Total = 0; //The total number of hits in a cluster  
    
    //For EVERY vertex, the algorithm is going to scan the plane to find clusters contained in the scanning cones
    
    for(unsigned int p = 0; p < geom->Nplanes(); p++) {  

      //  std::cout << "AT PLANE: " << p << std::endl;
      //for(unsigned int ivert = 0; ivert < vertSel.size(); ++ivert){
      for(unsigned int ivert = 7; ivert < 8; ++ivert){
	//for(unsigned int ivert = 0; ivert < 1; ++ivert){
	std::cout << "Number of STRONG vertices = " << vertSel.size()<< std::endl;
	//plane = -1;
	//if(p == plane){
	// std::cout << "Checking PLANE: " << p << std::endl;
	//get the coordinates of the vertex for the summit of the cone
	xa_cone = vertSel[ivert]->WireNum();
	ya_cone = vertSel[ivert]->DriftTime();
	std::cout << "Vertex at: (" << xa_cone << ", " << ya_cone << ")" <<std::endl;
	
	//Beginning of the scan!
	for(unsigned int iscan = 0; iscan < n_scan; iscan++){
	  std::cout << ">>>> Start SCAN: " << iscan << std::endl;  
	  
	  //define the scan anlge
	  scan_angle = (TMath::Pi()/2.0) - (iscan*(2.0*cone_angle));
	  std::cout << "Scan Angle: " << (180.*scan_angle)/TMath::Pi() << std::endl;
	  
	  //get the complementary angle for geometry puurposes
	  compl_angle = scan_angle - cone_angle;  
	  //std::cout << "Complementary Angle: " << (180.*compl_angle)/TMath::Pi() << std::endl;
	  
	  //Calculate the coordinates of the top right corner of the cone
	  x1_cone = xa_cone + fScone*(TMath::Cos(compl_angle));
	  y1_cone = ya_cone + fScone*(TMath::Sin(compl_angle));
	  
	  //Calculate the coordinates of the top left corner of the cone
	  x2_cone = xa_cone + fScone*(TMath::Cos(scan_angle + cone_angle));
	  y2_cone = ya_cone + fScone*(TMath::Sin(scan_angle + cone_angle));
	  
	  //std::cout << "cone vertex: ("<< xa_cone << ", " << ya_cone << ")"<< std::endl;
	  //std::cout << "cone bottom: ("<< x1_cone << ", " << y1_cone << ")"<< std::endl;
	  //std::cout << "cone top: ("<< x2_cone << ", " << y2_cone << ")"<< std::endl;
	  
	  //Looking if a cluster is in this cone (loop over all hits of all clusters)
	  for(int unsigned iclust = 0; iclust < clusterListHandle->size(); iclust++){
	  
	    art::Ptr<recob::Cluster> clust(clusterListHandle, iclust);
	    //std::cout << "Number of clusters: " << clusterListHandle->size() << std::endl;
	    
	    //Get the hits vector from the cluster
	    plane = p;
	    clusterhits = clust->Hits(plane);
	    if(clusterhits.size() == 0)continue;
	    //std::cout << "Cluster ID : " << clust->ID()<< std::endl;
	    //std::cout << "Got the clusterhits vector : iclust = " << iclust << "  : size = " << clusterhits.size() << std::endl;
	    //Loop over ALL hits in the cluster. Looking if the cluster's hit is comprised in the cone
	    for(unsigned int ihits = 0; ihits < clusterhits.size(); ihits++){
	      
	      //std::cout << "Nhits in cluster: " << clusterhits.size() << std::endl;
	      channel = clusterhits[ihits]->Wire()->RawDigit()->Channel();
	      geom->ChannelToWire(channel,tpc,plane,wire);
	      x_hit = channel;
	      y_hit = clusterhits[ihits]->PeakTime();
	      
	      //std::cout << " --- " << ihits << " : hits: (" << x_hit << ", " << y_hit << ")"<< std::endl;
	      
	      // Check in hits is INSIDE cone
	      
	      //define the 2 line equations:
	      //y2 = ((y2_cone - ya_cone)/(x2_cone - xa_cone))*x + ya_cone
	      //y1 = ((y1_cone - ya_cone)/(x1_cone - xa_cone))*x + ya_cone
	      if(y_hit <= ((y2_cone - ya_cone)/(x2_cone - xa_cone))*x_hit + ya_cone && y_hit >= ((y1_cone - ya_cone)/(x1_cone - xa_cone))*x_hit + ya_cone){
		hits_cluster_counter++;
		//std::cout << "***shower hit  found: " << hits_cluster_counter<< std::endl;
	      }		    
	      
	    }//end hits loop
	    
	    //std::cout << "hit_counter: "<< hits_cluster_counter << "  size: " << clusterhits.size() <<std::endl;
	    //std::cout << "division : " << (double)hits_cluster_counter / (double)clusterhits.size() << std::endl;
	    
	    //If there is more than 50% if the cluster INSIDE the cone, this is a protoshower
	    if(clusterhits.size() == 0) continue;
	    if(((double)hits_cluster_counter / (double)clusterhits.size()) >= 0.5){
	      std::cout << "GOT A SHOWER!!!  in scan " << iscan << "  cluster: " << iclust << " : " << clust->ID() << std::endl;
	      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      //NEED TO TAKE OUT THE HOUGH LINES FROM THE PROTOSHOWERS!!!!!  
	      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      //if(clust->ID() == showerID) continue;
	      //showerID = clust->ID();
	      protoShowers.push_back(clust);
	      //pshower_clust.push_back(clust); //this vector contains all the cluster in the cone proto-shower
		}
	    clusterhits.clear();
	    hits_cluster_counter = 0;
	    
	  } //end cluster loop

	  if(protoShowers.size() == 0) continue;
	  //std::cout << "@@@@@@ ASSOCIATING A SHOWER @@@@@@"<< std::endl;
	  recob::Shower shower(protoShowers);
	  showercol->push_back(shower);
	  protoShowers.clear();

	} //end scan loop
      } //end vertices loop
    }//end of plane loop 


    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //NEED TO SEPARATE THE SHOWERS FROM THE DIFFERENT VERTEX!!!!  
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
      
    std::cout << "--->  Recorded shower =  "<< showercol->size() << std::endl;
      //check if protoshower from further vertex is also contained in vertex nearer... TODO!!!!
      //if the shower is stand alone ok, else, erase the next one
      //shower.SetID(is);
      //shower.SetVertexCoord(xa_cone, ya_cone);

    vertSel.clear();

    evt.put(showercol);

  } // end of produce
} // end of namespace
