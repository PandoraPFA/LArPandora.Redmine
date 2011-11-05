////////////////////////////////////////////////////////////////////////
//
// \file LArTracker.cxx
//
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

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

// LArSoft includes
#include "LArTracker.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"
#include "Utilities/LArProperties.h"

// ROOT includes
#include "TMath.h"


//-------------------------------------------------
trkf::LArTracker::LArTracker(fhicl::ParameterSet const& pset) 
{
   this->reconfigure(pset);
   produces< std::vector<recob::Track> >();
}

//-------------------------------------------------
trkf::LArTracker::~LArTracker()
{
}

//-------------------------------------------------
void trkf::LArTracker::reconfigure(fhicl::ParameterSet const& pset)
{
   fClusterModuleLabel     = pset.get< std::string >("ClusterModuleLabel");
   ftmatch                 = pset.get< double >("TMatch");
}

//-------------------------------------------------
void trkf::LArTracker::beginJob()
{
   //someday this might do something
}

void trkf::LArTracker::endJob()
{
}

//------------------------------------------------------------------------------------//
void trkf::LArTracker::produce(art::Event& evt)
{ 


  // get services
  art::ServiceHandle<geo::Geometry> geom;

  // std::auto_ptr<> for the thing you want to put into the event
  std::auto_ptr<std::vector<recob::Track> > tcol(new std::vector<recob::Track>);

  // get input Cluster object(s).
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);
  

  for(unsigned int i = 0; i < clusterListHandle->size(); ++i)//grab info. for cluster #1
  {
     art::Ptr<recob::Cluster> cl1(clusterListHandle, i);
     
     // Figure out which View/Plane cluster #1 belongs to 
     int plane1 = cl1->View()-1;//there has to be a better way to deal with this...
 
     // Get all hits for this cluster
     art::PtrVector<recob::Hit> hits1 = cl1->Hits();

     // counter to use in setting ID of newly created spacepoints
     int num_spacepoints = 0;
     
     art::PtrVector<recob::Cluster> clusters;//create a vector for clusters associated with the track we create
     std::vector<recob::SpacePoint> spacepoints;//create a vector for spacepoints we may find within this cluster.

     for(art::PtrVector<recob::Hit>::const_iterator hit1 = hits1.begin(); hit1 != hits1.end();  hit1++) //loop over cluster #1 hits
     {
        if(FindHit(spacepoints,*hit1)) continue;//this hit already associated with a spacepoint

        //now make sure the hit wasn't already included in some other track we created
        bool result1 = false;
        for(unsigned int t = 0; t<tcol->size(); ++t){
           if(FindHit(tcol->at(t),*hit1)) result1 = true;
        }
        if(result1) continue;

        unsigned short channel1 = (*hit1)->Channel();
        unsigned int wire1,p1, t1;
        geom->ChannelToWire(channel1,t1, p1,wire1);
        double time1 = (*hit1)->PeakTime() ;
        
        art::PtrVector<recob::Hit> matched_hits;
        
        for(unsigned int j = 0; j < clusterListHandle->size(); ++j)//grab info. for cluster #2
        {
           art::Ptr<recob::Cluster> cl2(clusterListHandle,j);
           // Figure out which View/Plane cluster #2 belongs to 
           int plane2 = cl2->View()-1;
           
           if(plane1!=plane2 && !ClusterEndPointsMatch(cl1,cl2)) continue;//clusters from different planes that don't coincide in time
           
           art::PtrVector<recob::Hit> hits2 = cl2->Hits();
           
           for(art::PtrVector<recob::Hit>::const_iterator hit2 = hits2.begin(); hit2 != hits2.end();  hit2++) //loop over cluster #2 hits
           {

              if(FindHit(spacepoints,*hit2)) continue;
                
              bool result2 = false;
              for(unsigned int t = 0; t<tcol->size(); ++t){
                 result2 = FindHit(tcol->at(t),*hit2);
              }
              if(result2) continue;

              unsigned short channel2 = (*hit2)->Channel();
              double y,z;
              if(!geom->ChannelsIntersect(channel1,channel2,y,z) && plane1!=plane2) continue;//wires for the 2 hits come from different planes.
              //and they don't intersect.
              
              unsigned int wire2,p2,t2;
              geom->ChannelToWire(channel2,t2,p2,wire2);
              double time2 = (*hit2)->PeakTime() ;
              if(channel1==channel2 && time1==time2) continue;//don't compare a hit to itself
              
              double x1 = DriftCoordinate(plane1,time1);
              double x2 = DriftCoordinate(plane2,time2);
              
              //If hits coincide in the drift coordinate, and they are on wires that intersect (or from the same plane)
              //they belong in the same spacepoint.
              if(fabs(x1-x2)<ftmatch){
                 //Add the cluster/hit info. for the matched hits to the vectors which we'll use to create the track.
                 //Check for uniqueness in all cases, since we don't want to double-count the same hit/cluster.
                 if(!FindCluster(clusters,cl1)) clusters.push_back(cl1);
                 if(!FindCluster(clusters,cl2)) clusters.push_back(cl2);
                 if(!FindHit(matched_hits,*hit1)) matched_hits.push_back(*hit1);
                 if(!FindHit(matched_hits,*hit2)) matched_hits.push_back(*hit2);
              }
              
           }//loop over cluster #2 hits

        }//loop over other clusters

        if(matched_hits.size()>1 && MultiPlane(matched_hits)){
           matched_hits.sort();//sort the hits that will go into the SpacePoint
           recob::SpacePoint mysp(matched_hits);//create a new spacepoint with matched hit info.
           mysp.SetID(num_spacepoints++);
           spacepoints.push_back(mysp);//add it to the working list of spacepoints for this track
        }

     }//loop over cluster #1 hits     

     if(spacepoints.size()>0 && MultiPlane(clusters) && clusters.size()>1){//require at least one spacepoint, containing hits from multiple planes 
        clusters.sort();//sort the clusters that will go into the Track.
        std::sort(spacepoints.begin(),spacepoints.end());//sort the spacepoints that will go into the Track.
        recob::Track  the3DTrack(clusters,spacepoints);
        double dircos_start[3] = {0.};
        double dircos_end[3] = {0.};
        //some method to define track direction?
        the3DTrack.SetDirection(dircos_start,dircos_end);
        the3DTrack.SetID(tcol->size());
        tcol->push_back(the3DTrack);
     }

  }//1st loop over input clusters
    
  //print out some summary information for the new track collection
  mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
  mf::LogVerbatim("Summary") << "LArTracker Summary:";
  for(unsigned int i = 0; i<tcol->size(); ++i){
     std::cout << tcol->at(i) << std::endl;
     for(unsigned int j = 0; j<tcol->at(i).SpacePoints().size();++j){
        std::cout << tcol->at(i).SpacePoints().at(j) << std::endl;
        tcol->at(i).SpacePoints().at(j).PrintHits();
        // for(int k = 0; k<tcol->at(i).SpacePoints().at(j).Hits(-1).size(); ++k){
        //    mf::LogVerbatim("Summary") << tcol->at(i).SpacePoints().at(j).PrintHits();
        // }
     }
  }

  
  //save the new track collection
  evt.put(tcol);
  
}// produce

//------------------------------------------------------------------------------------//
double trkf::LArTracker::DriftCoordinate(int plane, double time)
{ 

   art::ServiceHandle<geo::Geometry> geom;
   // Parameters hard-coded here, but someday to be retrieved from a database...
   double timetick = 0.198;    //time sample in us
   double presamplings = 60.;
 
   double plane_pitch = geom->PlanePitch(0,1);   //wire plane pitch in cm 
   double Efield_drift = 0.5;  // Electric Field in the drift region in kV/cm
   double Efield_SI = 0.7;     // Electric Field between Shield and Induction planes in kV/cm
   double Efield_IC = 0.9;     // Electric Field between Induction and Collection planes in kV/cm
   double Temperature = 87.6;  // LAr Temperature in K

   art::ServiceHandle<util::LArProperties> larprop;
   double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature);    //drift velocity in the drift region (cm/us)
   double driftvelocity_SI = larprop->DriftVelocity(Efield_SI,Temperature);    //drift velocity between shield and induction (cm/us)
   double driftvelocity_IC = larprop->DriftVelocity(Efield_IC,Temperature);    //drift velocity between induction and collection (cm/us)
   double timepitch = driftvelocity*timetick;                                  //time sample (cm) 
   double tSI = plane_pitch/driftvelocity_SI/timetick;                         //drift time between Shield and Collection planes (time samples)
   double tIC = plane_pitch/driftvelocity_IC/timetick;                         //drift time between Induction and Collection planes (time samples)

   time -= presamplings;
	  	  
   //correct for the distance between wire planes
   if(plane==0) time -= tSI;         // Induction
   if(plane==1) time -= (tSI+tIC);   // Collection
	
   //transform hit wire and time into cm
   double time_cm = (double)(time * timepitch);

   return time_cm;

}

//------------------------------------------------------------------------------------//
bool trkf::LArTracker::ClusterEndPointsMatch(art::Ptr<recob::Cluster> c1, art::Ptr<recob::Cluster> c2)
{
   double x1_1 = DriftCoordinate(c1->View(),c1->StartPos()[1]); //start drift coordinate of cluster #1
   double x1_2 = DriftCoordinate(c1->View(),c1->EndPos()[1]); //end drift coordinate of cluster #1

   double x2_1 = DriftCoordinate(c2->View(),c2->StartPos()[1]); //start drift coordinate of cluster #2
   double x2_2 = DriftCoordinate(c2->View(),c2->EndPos()[1]); //end drift coordinate of cluster #2

   if(fabs(x1_1 - x2_1)<ftmatch && fabs(x1_2 - x2_2)) return true;
   else return false;
} 


//------------------------------------------------------------------------------------//
bool trkf::LArTracker::MultiPlane(art::PtrVector<recob::Hit> hitList)
{ 
  for(art::PtrVector<recob::Hit>::const_iterator hit = hitList.begin(); hit != hitList.end(); hit++){
    for(art::PtrVector<recob::Hit>::const_iterator hit1 = hitList.begin()+1; hit1 != hitList.end(); hit1++){
         if((*hit)->View() != (*hit1)->View() ) return true;
      }
      return false;
   }

   return true;
}

//------------------------------------------------------------------------------------//
bool trkf::LArTracker::MultiPlane(art::PtrVector<recob::Cluster> clusterList)
{ 
  for(art::PtrVector<recob::Cluster>::const_iterator cluster = clusterList.begin(); cluster != clusterList.end(); cluster++){
    for(art::PtrVector<recob::Cluster>::const_iterator cluster1 = clusterList.begin()+1; cluster1 != clusterList.end(); cluster1++){
         if((*cluster)->View() != (*cluster1)->View()) return true;
      }
      return false;
   }

   return true;
}

//------------------------------------------------------------------------------------//
bool trkf::LArTracker::FindCluster(art::PtrVector<recob::Cluster> clusterList, art::Ptr<recob::Cluster> c)
{ 
  for(art::PtrVector<recob::Cluster>::const_iterator cluster = clusterList.begin(); cluster != clusterList.end(); cluster++){
      if( ((*cluster)->View()==c->View()) && ((*cluster)->ID()==c->ID()) ) return true;
   }

   return false;
}

//------------------------------------------------------------------------------------//
bool trkf::LArTracker::FindHit(art::PtrVector<recob::Hit> hitList, art::Ptr<recob::Hit> h)
{ 
  for(art::PtrVector<recob::Hit>::const_iterator hit = hitList.begin(); hit != hitList.end(); hit++){
      if( ((*hit)->Channel()==h->Channel()) && ((*hit)->PeakTime() == h->PeakTime())) return true;
   }

   return false;
}

//------------------------------------------------------------------------------------//
bool trkf::LArTracker::FindHit(std::vector<recob::SpacePoint> spList, art::Ptr<recob::Hit> h)
{   
   for(unsigned int i=0;i<spList.size();++i){
      art::PtrVector<recob::Hit> hitList = spList[i].Hits(-1);//get all Hits in the spacepoint
      bool result = FindHit(hitList,h);
      if(result) return true;
   }

   return false;
}


//-----------------------------------------------------------------------------------//
bool trkf::LArTracker::FindHit(recob::Track track, art::Ptr<recob::Hit> h)
{
   std::vector<recob::SpacePoint> spList = track.SpacePoints();//get all spacepoints in the track
   return FindHit(spList,h);
}
