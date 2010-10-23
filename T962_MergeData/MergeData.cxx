////////////////////////////////////////////////////////////////////////
//
//   MergeData Class
//
//   soderber@fnal.gov
//   kinga.partyka@yale.edu
//
////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <math.h>
#include "T962_MergeData/MergeData.h"


#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include<dirent.h>

namespace merge{


int no=0;
int no_events_daq=0,no_events_daq_pot=0,no_e_pot_and_time_matched=0,no_e_pot_and_z_matched=0,no_e_pot_matched=0, e_w_predicted_tracks_outside_minos=0;
int no_distance_cut=0,no_degree_cut=0,no_all_cuts=0;
 int no_matched_tracks=0;
 std::vector<raw::MINOS > minos_tracks;
//-------------------------------------------------
MergeData::MergeData(edm::ParameterSet const& pset) : 
  
  fdaq_modulelabel(pset.getParameter< std::string >("daq")),
  //file(pset.getParameter< std::string >("file")),
  
  foundbeaminfo(false),
  foundpaddleinfo(false),
  foundminosinfo(false)
{
  }
  
  void MergeData::beginJob(edm::EventSetup const&)
  {
  // get access to the TFile service
    edm::Service<edm::TFileService> tfs;
    
  fdiff_x = tfs->make<TH1F>("fdiff_x","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
  fdiff_y = tfs->make<TH1F>("fdiff_y","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_z = tfs->make<TH1F>("fdiff_z"," Minos z coordinate value of tracks (that are candidates for matching)", 8000,-400 ,400);

fdiff_x_MIN= tfs->make<TH1F>("fdiff_x_MIN","Difference between predicted (from Argoneut) and closest Minos x coordinate value of track", 6500,-250 ,400);
  fdiff_y_MIN= tfs->make<TH1F>("fdiff_y_MIN","Difference between predicted (from Argoneut) and closest Minos y coordinate value of track", 6500,-250 ,400);

   fcos_diff_x = tfs->make<TH1F>("fcos_diff_x","Difference between  Argoneut and Minos cosx  value of track", 200,-1 ,1  );
   fcos_diff_y = tfs->make<TH1F>("fcos_diff_y","Difference between  Argoneut and Minos cosy value of track", 200,-1 ,1  );
   fcos_diff_z = tfs->make<TH1F>("fcos_diff_z","Difference between  Argoneut and Minos cosz  value of track", 200,-1 ,1  );

 fcos_diff_x_MIN = tfs->make<TH1F>("fcos_diff_x_MIN","Difference between  Argoneut and Minos cosx  value of the closest track", 200,-1 ,1  );
   fcos_diff_y_MIN = tfs->make<TH1F>("fcos_diff_y_MIN","Difference between  Argoneut and Minos cosy value of the closest track", 200,-1 ,1  );
fdiff_x_minos_coord= tfs->make<TH1F>("fdiff_x_minos_coord","Difference between  Argoneut and Minos candidate x tracks after adjusting to minos coord", 4000,-200 ,200  );
fdiff_y_minos_coord= tfs->make<TH1F>("fdiff_y_minos_coord","Difference between  Argoneut and Minos candidate y tracks after adjusting to minos coord", 4000,-200 ,200  );

 fdiff_x_m26 = tfs->make<TH1F>("fdiff_x_m26","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
 fdiff_y_m26 = tfs->make<TH1F>("fdiff_y_m26","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_m24 = tfs->make<TH1F>("fdiff_x_m24","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
 fdiff_y_m24 = tfs->make<TH1F>("fdiff_y_m24","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_m22 = tfs->make<TH1F>("fdiff_x_m22","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
 fdiff_y_m22 = tfs->make<TH1F>("fdiff_y_m22","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_m20 = tfs->make<TH1F>("fdiff_x_m20","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
 fdiff_y_m20 = tfs->make<TH1F>("fdiff_y_m20","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_m18 = tfs->make<TH1F>("fdiff_x_m18","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
 fdiff_y_m18 = tfs->make<TH1F>("fdiff_y_m18","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_m16 = tfs->make<TH1F>("fdiff_x_m16","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
 fdiff_y_m16 = tfs->make<TH1F>("fdiff_y_m16","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_m14 = tfs->make<TH1F>("fdiff_x_m14","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
 fdiff_y_m14 = tfs->make<TH1F>("fdiff_y_m14","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_m12 = tfs->make<TH1F>("fdiff_x_m12","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
 fdiff_y_m12 = tfs->make<TH1F>("fdiff_y_m12","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_m10 = tfs->make<TH1F>("fdiff_x_m10","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
 fdiff_y_m10 = tfs->make<TH1F>("fdiff_y_m10","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_m8 = tfs->make<TH1F>("fdiff_x_m8","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
 fdiff_y_m8 = tfs->make<TH1F>("fdiff_y_m8","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_m6 = tfs->make<TH1F>("fdiff_x_m6","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_m6 = tfs->make<TH1F>("fdiff_y_m6","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_m4 = tfs->make<TH1F>("fdiff_x_m4","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_m4 = tfs->make<TH1F>("fdiff_y_m4","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
   fdiff_x_m2 = tfs->make<TH1F>("fdiff_x_m2","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_m2 = tfs->make<TH1F>("fdiff_y_m2","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);



fdiff_x_p26 = tfs->make<TH1F>("fdiff_x_p26","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p26 = tfs->make<TH1F>("fdiff_y_p26","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
fdiff_x_p24 = tfs->make<TH1F>("fdiff_x_p24","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p24 = tfs->make<TH1F>("fdiff_y_p24","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
fdiff_x_p22 = tfs->make<TH1F>("fdiff_x_p22","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p22 = tfs->make<TH1F>("fdiff_y_p22","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
fdiff_x_p20 = tfs->make<TH1F>("fdiff_x_p20","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p20 = tfs->make<TH1F>("fdiff_y_p20","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
fdiff_x_p18 = tfs->make<TH1F>("fdiff_x_p18","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p18 = tfs->make<TH1F>("fdiff_y_p18","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
fdiff_x_p16 = tfs->make<TH1F>("fdiff_x_p16","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p16 = tfs->make<TH1F>("fdiff_y_p16","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
fdiff_x_p14 = tfs->make<TH1F>("fdiff_x_p14","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p14 = tfs->make<TH1F>("fdiff_y_p14","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_p12 = tfs->make<TH1F>("fdiff_x_p12","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p12 = tfs->make<TH1F>("fdiff_y_p12","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_p10 = tfs->make<TH1F>("fdiff_x_p10","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p10 = tfs->make<TH1F>("fdiff_y_p10","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_p8 = tfs->make<TH1F>("fdiff_x_p8","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p8 = tfs->make<TH1F>("fdiff_y_p8","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_p6 = tfs->make<TH1F>("fdiff_x_p6","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p6 = tfs->make<TH1F>("fdiff_y_p6","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
 fdiff_x_p4 = tfs->make<TH1F>("fdiff_x_p4","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p4 = tfs->make<TH1F>("fdiff_y_p4","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
   fdiff_x_p2 = tfs->make<TH1F>("fdiff_x_p2","Difference between predicted (from Argoneut) and Minos x coordinate value of track", 6500,-250 ,400);
   fdiff_y_p2 = tfs->make<TH1F>("fdiff_y_p2","Difference between predicted (from Argoneut) and Minos y coordinate value of track", 8000,-400 ,400);
   //.....................................................................................
   // Combinations:

fcosx_vs_cosy = tfs->make<TH2F>("fcosx_vs_cosy", ";Difference between  Argoneut and Minos cosx  value of track; Difference between  Argoneut and Minos cosy value of track", 200,     -1, 1, 200, -1, 1);
fcosx_vs_cosz = tfs->make<TH2F>("fcosx_vs_cosz", ";Difference between  Argoneut and Minos cosx  value of track; Difference between  Argoneut and Minos cosz value of track", 200,     -1, 1, 200, -1, 1);
fcosy_vs_cosz = tfs->make<TH2F>("fcosy_vs_cosz", ";Difference between  Argoneut and Minos cosy  value of track; Difference between  Argoneut and Minos cosz value of track", 200,     -1, 1, 200, -1, 1);

fdiff_x_vs_diff_y= tfs->make<TH2F>("fdiff_x_vs_diff_y", ";Difference between predicted (from Argoneut) and Minos x coordinate value of track; Difference between predicted (from Argoneut) and Minos y coordinate value of track",650 ,-250, 400, 800, -400, 400);



	ftrkqp= tfs->make<TH1F>("ftrkqp","Trkqp from Minos ", 60,-3 ,3);
ftime_diff_tms= tfs->make<TH1F>("ftime_diff_tms","utc (Minos)-tms/1000 ", 40000,-20000 ,20000);
ftime_diff_spilltime= tfs->make<TH1F>("ftime_diff_spilltime","utc (Minos)-spilltime ", 40000,-20000 ,20000);
ftrkErange=tfs->make<TH1F>("ftrkErange","trkErange ", 60,0 ,15);
ftrkE=tfs->make<TH1F>("ftrkE","trkE ", 180,0 ,45);

ftrkqp_single_track= tfs->make<TH1F>("ftrkqp_single_track","Trkqp for single best matched track to argoneut ", 60,-3 ,3);
ftrkErange_single_track=tfs->make<TH1F>("ftrkErange_single_track","trkErange for single best matched track to argoneut", 60,0 ,15);
ftrkE_single_track=tfs->make<TH1F>("ftrkE_single_track","trkE for single best matched track to argoneut", 180,0 ,45);

ftrkE_single_track_neg=tfs->make<TH1F>("ftrkE_single_track_neg","trkE for single best matched track to argoneut for negative charge", 300,0 ,30);

ftrkE_single_track_pos=tfs->make<TH1F>("ftrkE_single_track_pos","trkE for single best matched track to argoneut for positive charge", 300,0 ,30);


	f_distance_MIN=tfs->make<TH1F>("f_distance_MIN","sqrt(delta_x2 + delta_y2) for a minimum distance between argoneut's projected track and the closest minos track ", 3000,0 ,300);
	f_distance_MIN_minos_coord=tfs->make<TH1F>("f_distance_MIN_minos_coord","sqrt(delta_x2 + delta_y2) for a minimum distance between argoneut's projected track and the closest minos track after coordinate transformation ", 3000,0 ,300);

f_distance_MIN_minos_coord_squared=tfs->make<TH1F>("f_distance_MIN_minos_coord_squared","(sqrt(delta_x2 + delta_y2))^2 for a minimum distance between argoneut's projected track and the closest minos track after coordinate transformation", 1000,0 ,1000);



 f_distance_less_than_11=tfs->make<TH1F>(" f_distance_less_than_11","tracks that are less than 6cm from minos' closest track ", 3000,0 ,300);

fdiff_x_combined_minos_coord = tfs->make<TH1F>("fdiff_x_combined_minos_coord","Difference between projected and Minos x coordinate value of the best track", 6500,-250 ,400);
fdiff_y_combined_minos_coord = tfs->make<TH1F>("fdiff_y_combined_minos_coord","Difference between projected and Minos x coordinate value of the best track", 6500,-250 ,400);


 fcos_diff_x_combined_minos_coord = tfs->make<TH1F>("fcos_diff_x_combined_minos_coord","Difference between  Argoneut and Minos cosx value of track chosen by the shortest distance", 200,-1 ,1  );
 fcos_diff_y_combined_minos_coord = tfs->make<TH1F>("fcos_diff_y_combined_minos_coord","Difference between  Argoneut and Minos cosx value of track chosen by the shortest distance", 200,-1 ,1  );

fdegree_diff_y = tfs->make<TH1F>("fdegree_diff_y","Difference between Minos and Argoneut value of acos(cos(theta)) track chosen by the shortest distance", 20000,-100 ,100  );
fdegree_diff_x = tfs->make<TH1F>("fdegree_diff_x","Difference between Minos and Argoneut value of acos(cos(theta)) track chosen by the shortest distance", 20000,-100 ,100  );

ftime_diff_spilltime_all_cuts= tfs->make<TH1F>("ftime_diff_spilltime_all_cuts","utc (Minos)-spilltime ", 40000,-20000 ,20000);
 f_distance_all_cuts=tfs->make<TH1F>(" f_distance_all_cuts","distance in the x-y plane between us and minos' best track  ", 3000,0 ,300);
fdegree_diff_x_all_cuts= tfs->make<TH1F>("fdegree_diff_x_all_cuts","Difference between Minos and Argoneut value of acos(cos(theta)) track chosen by the shortest distance", 20000,-100 ,100  );
fdegree_diff_y_all_cuts= tfs->make<TH1F>("fdegree_diff_y_all_cuts","Difference between Minos and Argoneut value of acos(cos(theta)) track chosen by the shortest distance", 20000,-100 ,100  );

fno_minos_candidates= tfs->make<TH1F>("fno_minos_candidates","No of minos candidate tracks", 15,0 ,15  );
fno_tracks_less_than_11_apart= tfs->make<TH1F>("fno_tracks_less_than_11_apart","No of minos candidate tracks that are less than 11 cm apart for each argoneut track", 10,0 ,10  );
fno_tracks_less_than_15_apart= tfs->make<TH1F>("fno_tracks_less_than_15_apart","No of minos candidate tracks that are less than 15 cm apart for each argoneut track", 10,0 ,10  );
fno_tracks_less_than_20_apart= tfs->make<TH1F>("fno_tracks_less_than_20_apart","No of minos candidate tracks that are less than 20 cm apart for each argoneut track", 10,0 ,10  );
fno_tracks_less_than_30_apart= tfs->make<TH1F>("fno_tracks_less_than_30_apart","No of minos candidate tracks that are less than 30 cm apart for each argoneut track", 10,0 ,10  );
fno_tracks_less_than_40_apart= tfs->make<TH1F>("fno_tracks_less_than_40_apart","No of minos candidate tracks that are less than 40 cm apart for each argoneut track", 10,0 ,10  );


ftheta_x_argoneut = tfs->make<TH1F>("ftheta_x_argoneut","Angle of Argoneut tracks wrt x direction", 46000,-100 ,360  );
ftheta_y_argoneut = tfs->make<TH1F>("ftheta_y_argoneut","Angle of Argoneut tracks wrt x direction", 46000,-100 ,360  );
ftheta_x_minos = tfs->make<TH1F>("ftheta_x_minos","Angle of Minos tracks wrt x direction", 46000,-100 ,360  );
ftheta_y_minos = tfs->make<TH1F>("ftheta_y_minos","Angle of Minos tracks wrt y direction", 46000,-100 ,360  );


fdiff_x_less_than11 = tfs->make<TH1F>("fdiff_x_less_than11","Difference between Argoneut's( projected) and Minos x coordinate value of track for Minos tracks less than 11cm away from argoneut's tracks ", 6500,-250 ,400);
fdiff_y_less_than11 = tfs->make<TH1F>("fdiff_y_less_than11","Difference between Argoneut's( projected) and Minos y coordinate value of track for Minos tracks less than 11cm away from argoneut's tracks ", 8000,-400 ,400);

f_r_vs_abs_q_div_p= tfs->make<TH2F>("f_r_vs_abs_q_div_p", ";r;abs(q/p)", 2200,     0, 220, 80, 0, 8);
f_r_vs_q_div_p= tfs->make<TH2F>("f_r_vs_q_div_p", ";r;q/p", 2200,0 , 220, 100, -5, 5);

f_r_vs_abs_q_div_p_cut_11= tfs->make<TH2F>("f_r_vs_abs_q_div_p_cut_11", ";r;abs(q/p)", 2200,     0, 220, 80, 0, 8);
f_r_vs_q_div_p_cut_11= tfs->make<TH2F>("f_r_vs_q_div_p_cut_11", ";r;q/p", 2200,0 , 220, 100, -5, 5);

f_r_vs_abs_p= tfs->make<TH2F>("f_r_vs_abs_p", ";r;abs(p)", 2200, 0, 220, 1000, 0, 100);
f_r_vs_p= tfs->make<TH2F>("f_r_vs_p", ";r;p", 2200, 0, 220, 1000,-50, 50);

f_p_vs_r_cut11= tfs->make<TH2F>("f_p_vs_r_cut11", ";r;p", 2200, 0, 220, 3000, -150, 150);

f_p_for_positive_cut11= tfs->make<TH1F>("f_p_for_positive_cut11", "Momentum distribution for mu+", 500, 0, 50);
f_p_for_negative_cut11= tfs->make<TH1F>("f_p_for_negative_cut11", "Momentum distribution for mu-", 500, 0, 50);


}

//-------------------------------------------------
MergeData::~MergeData()
{
 
}

//-------------------------------------------------

void MergeData::produce(edm::Event& evt, edm::EventSetup const&)
{

std::auto_ptr<std::vector<raw::BeamInfo> > Beam_coll(new std::vector<raw::BeamInfo> );
std::auto_ptr<std::vector<raw::Paddles> > Paddles_coll(new std::vector<raw::Paddles> );
std::auto_ptr<std::vector<raw::MINOS > > Minos_coll(new std::vector<raw::MINOS > );
  
  edm::Handle< std::vector<raw::DAQHeader> > daqHandle;
    evt.getByLabel(fdaq_modulelabel,daqHandle);
    
  
    MergeBeam(Beam_coll); 
    MergePMT(Paddles_coll);
    MergeMINOS(Minos_coll);
    
  if(foundpaddleinfo){  evt.put(Beam_coll);}
    evt.put(Paddles_coll);
    evt.put(Minos_coll);
 return;
}



void MergeData::MergeBeam(std::auto_ptr<std::vector<raw::BeamInfo> > Beam_coll)
{
  std::cout<<"in mergebeam"<<std::endl;
  time_t spilltime = fDAQHeader[fDAQHeader.size()-1]->GetTimeStamp();//time info. from DAQ480 software
  std::cout<<"1***"<<std::endl;
  // std::cout<<"DAQ480(in MergeBeam()) tells us that the spilltime to match is: "<<spilltime<<std::endl;
  tm *timeinfo = localtime(&spilltime);
  // std::cout << "Run " << fDAQHeader[fDAQHeader.size()-1]->GetRun() << " Event = " << fDAQHeader[fDAQHeader.size()-1]->GetEvent() 
  //<< " time = " << fDAQHeader[fDAQHeader.size()-1]->GetTimeStamp() << " pretty = " << ctime(&spilltime) << std::endl;
  // printf("Date is (in MergeBeam) %d/%02d/%02d\n",timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday);
  //  printf("Time is (in MergeBeam) %02d:%02d:%02d\n",timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec);
  
  char beamfilename[20];
  
  
  sprintf(beamfilename,"matched_%02d_%02d_%d.txt",timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_year+1900);
  
  std::ifstream beamfile(beamfilename);
  
  if(!beamfile.is_open()){
    std::cerr << "MergeBeam:  Could not open file named " << beamfilename << std::endl;
    return;
  }
  
  long long int tms;
  double tor101;
  double tortgt;
  double trtgtd;
  std::string first,k;
  int event,run;
  
  
  
  while(getline(beamfile,k)){
    
    std::istringstream ins3;
    ins3.clear();
    ins3.str(k);
    ins3>>first>>run>>event>>tms>> tor101 >> tortgt >> trtgtd;
    //  std::cout<<std::setprecision(9)<<first<<"  "<<run<<"  "<<event<<"  "<<tms<<" "<<tor101<<"  "<< tortgt <<"  "<< trtgtd <<std::endl;
    
    
    if((fDAQHeader[fDAQHeader.size()-1]->GetRun()==run)&&(fDAQHeader[fDAQHeader.size()-1]->GetEvent()==event ) ){ foundbeaminfo=true;}
    
    
    if(foundbeaminfo){
      
      fBeamInfo.SetTOR101(tor101);
      fBeamInfo.SetTORTGT(tortgt);
      fBeamInfo.SetTRTGTD(trtgtd);
      fBeamInfo.SetT_MS(tms);
      Beam_coll->push_back(fBeamInfo);
      
      //  fBeamInfo.SetDATE(date);
      //       fBeamInfo.SetTIME(time);
      //  std::cout<<"-------------------------------"<<std::endl;
      // std::cout<<fBeamInfo<<std::endl;
      // std::cout<<"-------------------------------"<<std::endl;
      break;
    }
  }
  
  beamfile.close();
  
  return;
  
}

//-------------------------------------------------
void MergeData::MergePMT(std::auto_ptr<std::vector<raw::Paddles> >Paddles_coll)

{
  time_t spilltime = fDAQHeader[fDAQHeader.size()-1]->GetTimeStamp();//time info. from DAQ480 software
  // std::cout<<"DAQ480(in MergePMT()) tells us that the spilltime to match is: "<<spilltime<<std::endl;
  tm *timeinfo = localtime(&spilltime);
  // std::cout << "Run " << fDAQHeader[fDAQHeader.size()-1]->GetRun() << " Event = " << fDAQHeader[fDAQHeader.size()-1]->GetEvent() 
  //    << " time = " << fDAQHeader[fDAQHeader.size()-1]->GetTimeStamp() << " pretty = " << ctime(&spilltime) << std::endl;
  // printf("Date is (inMergePMT) %d/%02d/%02d\n",timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday);
  // printf("Time is (in MergePMT) %02d:%02d:%02d\n",timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec);
  
  char pmtfilename[20];
  sprintf(pmtfilename,"pmt_%02d_%02d_%d.txt",timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_year+1900);
  
  std::ifstream pmtfile(pmtfilename);
  std::string line;
  
  if(!pmtfile.is_open()){
    std::cerr << "Merge:  Could not open file named " << pmtfilename << std::endl;
    return;
  }
  
  std::string g;
  time_t t;
  std::string word;
  std::string pflag;
  std::string footer;
  int hit[4];
  int array1[4];
  int array2[4];
  int array3[4];
  int array4[4];
  
  
  while(getline(pmtfile,g)){
    if(g.find("Couldn't")){
      if(!g.find("time")){//lines that contain time
	
	std::istringstream ins2;
	ins2.clear();
	ins2.str(g);
	ins2>>word>>t;
	
	
	if(fabs(t-spilltime)==0){ foundpaddleinfo=true;}
	
      }      
      else if (!g.find("end")){
	std::istringstream ins3;
	ins3.clear();
	ins3.str(g);
	ins3>>footer;
	g.clear();
	
	if(foundpaddleinfo){
	  
	  fPaddles.SetTime(t);
	  fPaddles.SetPMT(0,array1);
	  fPaddles.SetPMT(1,array2);
	  fPaddles.SetPMT(2,array3);
	  fPaddles.SetPMT(3,array4);
	  Paddles_coll->push_back(fPaddles);
	  // std::cout<<"-------------------"<<std::endl;
	  // std::cout<<fPaddles<<std::endl;
	  //       std::cout<<"-------------------"<<std::endl;
	  break;
	}
      }
      else{
	
	std::istringstream ins4;  
	ins4.clear();
	ins4.str(g);
	
	ins4>>pflag;
	
	for(int i=0; i<4; ++i){
	  ins4>>hit[i];
	}
	
	
	if(!g.find("pmt1")){
	  for(int j=0;j<4;++j){
	    array1[j]=hit[j];
	  }
	  
	}
	
	if(!g.find("pmt2")){
	  for(int j=0;j<4;++j){
	    array2[j]=hit[j];
	  }
	  
	}
	
	if(!g.find("pmt3")){
	  for(int j=0;j<4;++j){
	    array3[j]=hit[j];
	  }
	  
	}
	
	if(!g.find("pmt4")){
	  for(int j=0;j<4;++j){
	    array4[j]=hit[j];
	  }
	  
	}
	
      }
      
      
      
    }
  }
  
  
  
  pmtfile.close();
  return ;
  
  
}


//-------------------------------------------------
void MergeData::MergeMINOS(std::auto_ptr<std::vector<raw::MINOS> >Minos_coll)

{
  
  
  time_t spilltime = fDAQHeader[fDAQHeader.size()-1]->GetTimeStamp();//time info. from DAQ480 software
  ////////////////////////////////////////////////////////////////////////////////////
  
  tm *timeinfo = localtime(&spilltime);
  
  char matchedfilename[20];

  sprintf(matchedfilename,"matched_%02d_%d_%d.txt",timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_year+1900);
 
  std::ifstream matchedfile(matchedfilename);
  std::string line;
  //std::cout<<"****Needs to open :"<<matchedfilename<<std::endl;
 

  if(!matchedfile.is_open()){
    std::cerr << "Merge:  Could not open file named " << matchedfilename << std::endl;
    return;
  }
  //std::cout<<"In "<<matchedfilename<<std::endl;
  
  //..........................
  
  
//take a look here later **************


  //jobc::CmdLine blah= jobc::CmdLine::Instance();
  //std::list<std::string> file= blah.InFileList(); //this is for taking a file that was entered on the cmd line and opening the appropriate file from Maddalena
  
  //std::cout<<"!!!!!!!!!!!file="<<file.front()<<std::endl;
  
  //std::string::size_type pos_start=file.rfind("/");
  // std::cout<<"pos_start= "<<pos_start<<std::endl;
  
  //std::string::size_type pos_end=file.rfind("root");
  //std::cout<<"pos_end= "<<pos_end<<std::endl;
  
 // std::string finalname;
 // finalname=file.substr(pos_start+2,pos_end-pos_start-3);//changing pos_start+1 to pos_start+2 for mR729.root (extra m in front) and "pos_end-pos_start-2" to pos_end-pos_start-3
  
  //std::cout<<"WILL BE READING:"<<finalname<<std::endl;
  
  
  // for(int i=(int)pos;i<file.front().length;i++)
  //      {
  
  //        std::cout<<"and i get the filename of *** "<<file.front()[i]<<std::endl;
  //      }
  
  
  
  
  
  // std::string run_cmdline;
  //   std::stringstream out;
  //   out<<fEvt;
  //   run_cmdline=out.str();
  
  std::ifstream  dircosfile;
  
  //std::string fevt=fEvt;
  
  //i will hardcode dircos file to read for now (since we just have one):
  //  std::string name="dircos3_"; //changed dircos_ to dircos2_
  //   std::string txt=".txt";
  //   std::string fileName=name+finalname+txt;
  //   const char *dircosfilename=fileName.c_str();
  //   std::cout<<"WILL BE READING: "<<dircosfilename<<std::endl;
  
  char dircosfilename[38]="dircos_bis_mR728-9_sel.txt";
  std::cout<<"WILL BE READING: "<<dircosfilename<<std::endl;
  
  dircosfile.open(dircosfilename, std::ios::in);
  
  if(!dircosfile.is_open()){
    std::cerr << "Merge:  Could not open file named " << dircosfilename << std::endl;
    return;
  }
  //...........................
  
  
  


  
  

 


  ///////////////////////////////////////////////////////////////////////////////////
  //std::string path ="/argoneut/data/users/soderber/MINOS/";
    
    std::string path= "/argoneut/app/users/kpartyka/larsoft/MINOS/";//gotta replace New_MINOS with MINOS later. I am testing new file 09.15.2010
  // int no_files=0;
   
  //TFile *f = new TFile("/argoneut/data/users/soderber/MINOS/RunN00017288.root");
   
 
  int no_argoneut_events=0;
  //***************************************
  //MINOS VARIABLES:

  int run;
  int subRun;
  int snarl;
  double utc=18;
  double day;
  float trkIndex;
  float trkE;
  float shwE;
  double crateT0;
  double tmframe;
  double year;
  float vtxX;
  float vtxY;
  float vtxZ;
  float trkErange;
  double sgate53;
  float trkqp;     
  float trkVtxX;
  float trkVtxY;
  float trkVtxZ;
  float trkdcosx;
  float trkdcosy;
  float trkdcosz;
  double month;
  float trkChi2;
  float trkTimeT0;
  float trkVtxT;
  float tor101;
  float tortgt;
  float trtgtd;
  float tr101d;

  float trkmom;
  int charge;
  float trkstpU;
  float trkstpV;
  int ntrkstp;
  float trkstpX;
  float trkstpY;
  float trkstpZ;
  float trkeqp;
  float trkVtxeX;
  float trkVtxeY;
  
  std::vector<float> vtx;
    std::vector<float> trkVtx;
    std::vector<float> trkdcos;
     std::vector<float> trkstp;
     std::vector<float> trkVtxe;
  //***************************************

 
  // long long int tms;
  double tms;
  double tor101_m;
  double tortgt_m;
  double trtgtd_m;
  std::string first,k;
  int event,run_;
  std::string l,w_event,semicolon,sem,w_run,dir,cos,start,pt,cm;
  double cosx,cosy,cosz;
  int run_no,event_no;
   
  int pass1=0,pass2=0,pass3=0;
  int Run=0;
  int Event=0;
  double Cosx=0,Cosy=0,Cosz=0;
   
  double x_s,y_s,z_s;
  double x_start_a=0;
  double y_start_a=0;
  double z_start_a=0;
  int in_argoneut=0;
  double D=(90*0.5)+(42.4*2.54)-5.588; //distance from the front (upstream) of the TPC to the 1st Minos plane (this minus number is the one we measured with Mitch)
  double x_predicted=0, y_predicted=0;
  std::vector<double> minos_trk_x;
  std::vector<double> minos_trk_y, minos_trk_z,time_diff_tms,time_diff_spilltime;
  std::vector<float> minos_cos_x,minos_cos_y,minos_cos_z, minos_trk_index;
  std::vector<double> 	diff_x;
  std::vector<int> events;
  double d_z=0;
  double vector_length= 0;
  double cosx_ar=0,cosy_ar=0,cosz_ar=0;

  std::vector<double> v_x_start_a, v_y_start_a,v_z_start_a, v_Cosx, v_Cosy, v_Cosz;
  v_x_start_a.clear();
 v_y_start_a.clear();
v_z_start_a.clear();
 v_Cosx.clear();
 v_Cosy.clear();
 v_Cosz.clear();
  std::vector<float> minos_trkqp,minos_trkE,minos_trkErange;
  //for determination of the right distance z from minos:
  //....................................
  // double vector_length_m26= 0,vector_length_m24= 0,vector_length_m22= 0,vector_length_m20= 0,vector_length_m18= 0,vector_length_m16= 0,vector_length_m14= 0,vector_length_m12= 0,vector_length_m10= 0 ,vector_length_m8= 0, vector_length_m6= 0, vector_length_m4= 0, vector_length_m2= 0, vector_length_p2= 0, vector_length_p4= 0, vector_length_p6= 0, vector_length_p8= 0, vector_length_p10= 0,vector_length_p12= 0,vector_length_p14= 0,vector_length_p16= 0,vector_length_p18= 0,vector_length_p20= 0,vector_length_p22= 0,vector_length_p24= 0,vector_length_p26= 0;
  // double d_z_m26= 0,d_z_m24= 0,d_z_m22= 0,d_z_m20= 0,d_z_m18= 0,d_z_m16= 0,d_z_m14= 0,d_z_m12= 0, d_z_m10= 0,d_z_m8= 0, d_z_m6= 0, d_z_m4= 0, d_z_m2= 0, d_z_p2= 0, d_z_p4= 0, d_z_p6= 0, d_z_p8= 0, d_z_p10= 0,d_z_p12= 0,d_z_p14= 0,d_z_p16= 0,d_z_p18= 0,d_z_p20= 0,d_z_p22= 0,d_z_p24= 0,d_z_p26= 0;
  //double x_predicted_m26= 0,x_predicted_m24= 0,x_predicted_m22= 0,x_predicted_m20= 0,x_predicted_m18= 0,x_predicted_m16= 0 ,x_predicted_m14= 0,x_predicted_m12= 0,x_predicted_m10= 0,x_predicted_m8= 0,x_predicted_m6= 0, x_predicted_m4= 0, x_predicted_m2= 0, x_predicted_p2= 0, x_predicted_p4= 0, x_predicted_p6= 0, x_predicted_p8= 0, x_predicted_p10= 0, x_predicted_p12= 0,x_predicted_p14= 0,x_predicted_p16= 0,x_predicted_p18= 0,x_predicted_p20= 0,x_predicted_p22= 0,x_predicted_p24= 0,x_predicted_p26= 0;
  // double y_predicted_m26= 0,y_predicted_m24= 0,y_predicted_m22= 0,y_predicted_m20= 0,y_predicted_m18= 0, y_predicted_m16= 0, y_predicted_m14= 0,y_predicted_m12= 0,y_predicted_m10= 0 ,y_predicted_m8= 0,y_predicted_m6= 0, y_predicted_m4= 0, y_predicted_m2= 0, y_predicted_p2= 0, y_predicted_p4= 0, y_predicted_p6= 0,y_predicted_p8= 0,y_predicted_p10= 0,y_predicted_p12= 0,y_predicted_p14= 0,y_predicted_p16= 0,y_predicted_p18= 0,y_predicted_p20= 0,y_predicted_p22= 0,y_predicted_p24= 0,y_predicted_p26= 0;
   
  //........................................................................
  double x_offset=116.9; // previously 118;
  double y_offset=20.28; //previously  19;
  int in_matching=0,pot_matched=0,pot_and_z_matched=0,pot_and_time_matched=0;
  double const pi=4.0*atan(1.0);
   
   
  //std::cout <<" utc = " <<std::setprecision(10)<< utc<<std::endl; 
  while(getline(matchedfile,k))
    {
      //  std::cout<<"here3"<<std::endl;
      std::istringstream ins3;
      ins3.clear();
      ins3.str(k);
      ins3>>first>>run_>>event>>tms>> tor101_m >> tortgt_m >> trtgtd_m;
      //  std::cout<<std::setprecision(9)<<first<<"  "<<run<<"  "<<event<<"  "<<tms<<" "<<tor101<<"  "<< tortgt <<"  "<< trtgtd <<std::endl;
       
      //  std::cout <<"run= "<<run_<< "event = " << event <<" tms= "<<std::setprecision(13)<<tms/1000<< " utc = " <<std::setprecision(10)<< utc << " trkIndex = " << trkIndex << " Vertex(X,Y,Z) = (" 	      << trkVtxX << "," << trkVtxY << "," << trkVtxZ << ") " << " cos(x,y,z) = (" << trkdcosx   << "," << trkdcosy << "," << trkdcosz << ")" << std::endl;
       
      //if we can match DAQ with matched file
      if((run_== fDAQHeader[fDAQHeader.size()-1]->GetRun()) && (event==fDAQHeader[fDAQHeader.size()-1]->GetEvent()))
	{
	  no_events_daq++;
	   
	  //std::cout<<"Run: "<<run_<<" Event : "<<event<<std::endl;
	  //loop through MINOS root files
	   
	  int no_files=0;
	  DIR *pDIR;
	  struct dirent *entry;
	 
	   if( pDIR=opendir("/argoneut/app/users/kpartyka/larsoft/MINOS") )
	  //if( pDIR=opendir("/argoneut/app/users/kpartyka/larsoft/New_MINOS") )//checking the new file from rashid 09.15.2010
	    {
	       
	       
	      while(entry = readdir(pDIR))
		{
		  int utc_same_as_t=0;
		   
		  if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 ){
		    std::cout <<"file is "<< entry->d_name << "\n";
		    no_files ++;}
		  else continue;
		   
		  std::string fileName = entry->d_name;
		  std::string file=path+fileName;
		   
		  TFile *f = new TFile(file.c_str());
		  TTree *minitree = (TTree*)f->Get("minitree");
		  minitree->SetBranchAddress("run",&run);
		  minitree->SetBranchAddress("subRun",&subRun);
		  minitree->SetBranchAddress("snarl",&snarl);
		  minitree->SetBranchAddress("utc",&utc);
		  minitree->SetBranchAddress("day",&day);
		  minitree->SetBranchAddress("trkIndex",&trkIndex);
		  minitree->SetBranchAddress("trkE",&trkE);
		  minitree->SetBranchAddress("shwE",&shwE);
		  minitree->SetBranchAddress("crateT0",&crateT0);
		  minitree->SetBranchAddress("tmframe",&tmframe);
		  minitree->SetBranchAddress("year",&year);
		  minitree->SetBranchAddress("vtxX",&vtxX);
		  minitree->SetBranchAddress("vtxY",&vtxY);
		  minitree->SetBranchAddress("vtxZ",&vtxZ);
		  minitree->SetBranchAddress("trkErange",&trkErange);
		  minitree->SetBranchAddress("sgate53",&sgate53);
		  minitree->SetBranchAddress("trkqp",&trkqp);
		  minitree->SetBranchAddress("trkVtxX",&trkVtxX);
		  minitree->SetBranchAddress("trkVtxY",&trkVtxY);
		  minitree->SetBranchAddress("trkVtxZ",&trkVtxZ);
		  minitree->SetBranchAddress("trkdcosx",&trkdcosx);
		  minitree->SetBranchAddress("trkdcosy",&trkdcosy);
		  minitree->SetBranchAddress("trkdcosz",&trkdcosz);
		  minitree->SetBranchAddress("month",&month);
		  minitree->SetBranchAddress("trkChi2",&trkChi2);
		  minitree->SetBranchAddress("trkTimeT0",&trkTimeT0);
		  minitree->SetBranchAddress("trkVtxT",&trkVtxT);
		  minitree->SetBranchAddress("tr101d",&tr101d);
		  minitree->SetBranchAddress("trtgtd",&trtgtd);
		  minitree->SetBranchAddress("tor101",&tor101);
		  minitree->SetBranchAddress("tortgt",&tortgt);

		  minitree->SetBranchAddress("charge",&charge);
 		  minitree->SetBranchAddress("trkmom",&trkmom);
		  minitree->SetBranchAddress("trkstpU",&trkstpU);
		  minitree->SetBranchAddress("trkstpV",&trkstpV);
		  minitree->SetBranchAddress("ntrkstp",&ntrkstp);
		  minitree->SetBranchAddress("trkstpX",&trkstpX);
		  minitree->SetBranchAddress("trkstpY",&trkstpY);
		  minitree->SetBranchAddress("trkstpZ",&trkstpZ);
		  minitree->SetBranchAddress("trkeqp",&trkeqp);
		  minitree->SetBranchAddress("trkVtxeX",&trkVtxeX);
		  minitree->SetBranchAddress("trkVtxeY",&trkVtxeY);
		  

 
		  //...................................................

	

		  //------------------------------------------
		  Long64_t nentries = minitree->GetEntries();
		  Long64_t nbytes = 0;
		  for(Long64_t i=0;i<nentries;i++){
		    nbytes+=minitree->GetEntry(i);
		    //----------------------------------------------
		    //for print out of minos data:
		     
		    //std::cout<< " utc= "<<std::setprecision(13)<<utc<<" tor101=  "<<std::setprecision(13)<<tor101<<" tortgt=  "<<std::setprecision(13)<<tortgt<<" trtgtd= "<<std::setprecision(13)<<trtgtd<<" tr101d= "<<tr101d<<" trkIndex = " << trkIndex <<" Vertex(X,Y,Z) = (" 	      << trkVtxX << "," << trkVtxY << "," << trkVtxZ << ") " << " cos(x,y,z) = (" << trkdcosx   << "," << trkdcosy << "," << trkdcosz << ")" << std::endl;
		     

		     
		     
		    //........................................................
		     
 

		    // if((utc-spilltime)<0) break;//no further match for ArgoNeuT event in this file.
		     
		     
		     
		    //  if(fabs(utc-spilltime)<1000000 && fabs(trtgtd_m-trtgtd)<0.000001 && fabs(tor101_m-tor101)<0.00001 && fabs(tortgt_m-tortgt)<0.00001){
		    // if(fabs(utc-(tms/1000))<60 && fabs(trtgtd_m-trtgtd)<0.001 && fabs(tor101_m-tor101)<0.001 && fabs(tortgt_m-tortgt)<0.001){
		     
		    if( fabs(trtgtd_m-trtgtd)<0.001 && fabs(tor101_m-tor101)<0.001 && fabs(tortgt_m-tortgt)<0.001){pot_matched=1; }
		    if( fabs(trtgtd_m-trtgtd)<0.001 && fabs(tor101_m-tor101)<0.001 && fabs(tortgt_m-tortgt)<0.001 && fabs(trkVtxZ*100)<10){pot_and_z_matched=1;}
		    if( fabs(utc-(tms/1000))<60 && fabs(trtgtd_m-trtgtd)<0.001 && fabs(tor101_m-tor101)<0.001 && fabs(tortgt_m-tortgt)<0.001 ){pot_and_time_matched=1;}
		    if(fabs(utc-(tms/1000))<2 && fabs(trtgtd_m-trtgtd)<0.001 && fabs(tor101_m-tor101)<0.001 && fabs(tortgt_m-tortgt)<0.001 && fabs(trkVtxZ*100)<10){

		    

		    

		     
		      in_matching=1;
		       
		      events.push_back(event);
		      minos_trk_x.push_back(trkVtxX*100);
		      minos_trk_y.push_back(trkVtxY*100);
		      minos_trk_z.push_back(trkVtxZ*100);
		       
		      minos_cos_x.push_back(trkdcosx);
		      minos_cos_y.push_back(trkdcosy);
		      minos_cos_z.push_back(trkdcosz);
		       
		       
		      minos_trkqp.push_back(trkqp);
		      minos_trkE.push_back(trkE);
		      minos_trkErange.push_back(trkErange);
		      // std::cout<<"pushing with "<<fabs(utc-(tms/1000))<<std::endl;
		      time_diff_tms.push_back(utc-(tms/1000));
		      time_diff_spilltime.push_back(utc-spilltime);
		      minos_trk_index.push_back(trkIndex);
		       
		      std::cout<<"************************************YAY*****************"<<std::endl;

		


		       
		      std::cout<<"spilltime= "<<std::setprecision(13)<<spilltime<< " utc= "<<std::setprecision(13)<<utc<<" trtgtd_m= "<<std::setprecision(13)<<trtgtd_m<<"   "<<" trtgtd= "<<std::setprecision(13)<<trtgtd<<" tor101_m= "<<std::setprecision(13)<<tor101_m<<" tor101=  "<<std::setprecision(13)<<tor101<<" tortgt_m= "<<std::setprecision(13)<<tortgt_m<<" tortgt=  "<<std::setprecision(13)<<tortgt<<" tr101d= "<<tr101d<<" trkIndex = " << trkIndex <<" Run = "<<fDAQHeader[fDAQHeader.size()-1]->GetRun()  <<" Event = "<< fDAQHeader[fDAQHeader.size()-1]->GetEvent() <<"tms="<<tms<<std::endl;
		       
		       
		      std::cout<< " Vertex(X,Y,Z) = (" 	      << trkVtxX << "," << trkVtxY << "," << trkVtxZ << ") " << " cos(x,y,z) = (" << trkdcosx   << "," << trkdcosy << "," << trkdcosz << ")" << std::endl;


		      std::cout<<"utc-spilltime= "<<std::setprecision(10)<<(utc-spilltime)<<" s"<<" , "<<(utc-spilltime)/3600<<" h"<<std::endl;
		      std::cout<<"utc-tms= "<<std::setprecision(10)<<(utc-(tms/1000))<<" s"<<" , "<<(utc-(tms/1000))/3600<<" h"<<std::endl;
		      //std::cout<<"trtgtd (matched)= "<<std::setprecision(13)<<trtgtd<<"   "<<" trkTimeT0= "<<std::setprecision(13)<<trkTimeT0<<" trkVtxT= "<<std::setprecision(13)<<trkVtxT<<" trkChi2=  "<<trkChi2<<std::endl;
		      utc_same_as_t=1;
		      //--------------------------------
		       
		       
		       
		      //ok, so we found a match between daq and minos times,and now I want a match up to some number also between direction cosines.Open up and read dircos files
		      //................................................
		       

		      // DIR *mDIR;
		      // 	  struct dirent *m_entry;
		       
		      // 	   if( mDIR=opendir("/argoneut/app/users/kpartyka/larsoft/Maddalena") )
		      // 	    {
		       
		       
		      // 	      while(entry = readdir(mDIR))
		      // 		{
		       

		      // 		  if( strcmp(m_entry->d_name, ".") != 0 && strcmp(m_entry->d_name, "..") != 0 ){
		      // 		    std::cout <<"file is "<< m_entry->d_name << "\n";
		      // 		    }
		      // 		  else continue;
		       

		      // We now open argoneut 3-d reconstruction file....it is put inside minos tracks because we are only interested in argoneut's events that CAN be matched with minos ...

		  
		     

		   while(getline(dircosfile,l))
		     {
		       Run=0;
		       Event=0;
		       Cosx=0;
		       Cosy=0;
		       Cosz=0;
		       x_start_a=0;
		       y_start_a=0;
		       z_start_a=0;
		      

			   
			   
			   
			  pass3=0;
			  if(!l.find("run")){//lines that contain run#
			    
			    std::istringstream ins2;
			    ins2.clear();
			    ins2.str(l);
			    ins2>>w_run>>semicolon>>run_no;
			    // l.clear();
			    Run=run_no;
			    pass1=1;
			    //std::cout<<"**Run  "<<run_no<<std::endl;
			  }
			  else{ continue;} //
			  getline(dircosfile,l);
			   
			  if(!l.find("event")){//lines that contain event#
			    no_argoneut_events++;
			    std::istringstream ins3;
			    ins3.clear();
			    ins3.str(l);
			    ins3>>w_event>>semicolon>>event_no;
			    Event=event_no;
			    // std::cout<<"**Event  "<<event_no<<" "<<Event<<std::endl;
			    // l.clear();
			    pass2=1;
			     
			  }
			  getline(dircosfile,l);
			  getline(dircosfile,l);
			   
			   
			   
			   
			  if(!l.find("Direction")){//lines that contain event#
			     
			    std::istringstream ins4;
			    ins4.clear();
			    ins4.str(l);
			    ins4>>dir>>cos>>sem>>cosx>>cosy>>cosz;
			    // std::cout<<"**Direction cosine  "<<cosx<<"  "<<cosy<<"  "<<cosz<<std::endl;
			    Cosx=cosx;
			    Cosy=cosy;
			    Cosz=cosz;
			    pass3=1;
			    //l.clear();
			  }
			  //std::cout<<"pass3="<<pass3<<std::endl;
			  //	if(pass3!=1) continue;
			  //if(pass3==1){
			  //...................................................
			  // std::cout<<"$$$$$$$Direction cosine  "<<Cosx<<"  "<<Cosy<<"  "<<Cosz<<std::endl;
			   
			  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
			  getline(dircosfile,l);
			  getline(dircosfile,l);
			  getline(dircosfile,l);
		
			  // std::cout<<"x_start_a= "<<x_start_a<<std::endl;
			   
			  if(!l.find("Start")){//lines that contain event#
			     
			    std::istringstream ins5;
			    ins5.clear();
			    ins5.str(l);
			    ins5>>start>>pt>>cm>>semicolon>>x_s>>y_s>>z_s;
			    // std::cout<<"**  "<<start<<"  "<<pt<<"  "<<x_s<<" "<<y_s<<" "<<z_s<<std::endl;
			    x_start_a=x_s;
			    y_start_a=y_s;
			    z_start_a=z_s;
			    
			    // l.clear();
			  }
			   
			  //std::cout<<"2: x_start_a= "<<x_start_a<<std::endl;
			   
			   
			  // if((Run== fDAQHeader[fDAQHeader.size()-1]->GetRun()) && (Event==fDAQHeader[fDAQHeader.size()-1]->GetEvent()) && (fabs(Cosx-trkdcosx )<1) && (fabs(Cosy-trkdcosy )<1) && (fabs(Cosz-trkdcosz )<1)){
			   
			  if((Run== fDAQHeader[fDAQHeader.size()-1]->GetRun()) && (Event==fDAQHeader[fDAQHeader.size()-1]->GetEvent())){
			     
			    std::cout<<"found a match in our 3d reco file for event "<<Event<<std::endl;
			    in_argoneut=1;
			     
			    no_matched_tracks++;
			     
			    //  d_z=D-z_start_a;
			    // 			     vector_length= d_z/Cosz;
			    // 			     x_predicted=vector_length*Cosx;
			    // 			     y_predicted=vector_length*Cosy;
			     
			     
			    cosx_ar=Cosx;
			    cosy_ar=Cosy;
			    cosz_ar=Cosz;
			     
			     
			  

			    v_x_start_a.push_back(x_start_a);
			    v_y_start_a.push_back(y_start_a);
			    v_z_start_a.push_back(z_start_a);
			    v_Cosx.push_back(Cosx);
			    v_Cosy.push_back(Cosy);
			    v_Cosz.push_back(Cosz);

			   
			    continue;
 






			    
			    
			  }//if run==run and event==event
			  getline(dircosfile,l);
			  getline(dircosfile,l);
			  
			  if (in_argoneut==1 )break; 
			}// while getline
		       
		   //std::cout<<"check0. no of argo tracks= "<<v_x_start_a.size()<<std::endl;
		       
		       
		     
		    //std::cout<<"No of events reconstructed by Maddalena= "<<no_argoneut_events<<std::endl;
		       
		       
		      dircosfile.close();
		       
		       
		       
		       
		      //-------------------end of maddalena file-----------------------------------
		      //.........................................................................
		      // if using maddalena's file is not necessary ---->uncomment this:
		      foundminosinfo=true;
		      
		      vtx.push_back(vtxX);
		      vtx.push_back(vtxY);
		      vtx.push_back(vtxZ);
		      
		      trkVtx.push_back(trkVtxX);
		      trkVtx.push_back(trkVtxY);
		      trkVtx.push_back(trkVtxZ);
		      
		      trkdcos.push_back(trkdcosx);
		      trkdcos.push_back(trkdcosy);
		      trkdcos.push_back(trkdcosz);
		      
		      trkstp.push_back(trkstpU);
		      trkstp.push_back(trkstpV);
		      trkstp.push_back(ntrkstp);
		      
		      trkstp.push_back(trkstpX);
		      trkstp.push_back(trkstpY);
		      trkstp.push_back(trkstpZ);
		      
		      trkVtxe.push_back(trkVtxeX);
		      trkVtxe.push_back(trkVtxeY);
		      
		      
		      raw::MINOS my_minos(run,subRun,snarl,utc,day,trkIndex,trkE,shwE,crateT0,tmframe,year,
		      vtx,trkErange,sgate53,trkqp,trkVtx, trkdcos,month,trkmom, charge, trkstp,trkeqp,trkVtxe,0);
		      Minos_coll->push_back(my_minos);
		       
		       
		      
		       
		       
		       
		      minos_tracks.push_back(my_minos);
		      std::cout<<"length of minos_tracks="<<minos_tracks.size()<<std::endl;
		      //...................................................................
		       
		       
		      //  std::cout<<"exited while"<<std::endl;
		    }//found a match
		    //	std::cout<<"after matched"<<std::endl;
		  }//loop over Minos root file
		  // std::cout<<"finished looping minosfile"<<std::endl;
		   
		   
		  minitree->Delete();
		  f->Close();
		  f->Delete();
		  // if( utc_same_as_t==1)break;
		}//while
	      closedir(pDIR);
	      // std::cout<<no_files<<std::endl;
	       
	       
	       
	      //Take diff between argoneut x coordinate and the one predicted from argoneuut (same for y)
	       	std::vector<double> all_cosx,degree_diff_x;
	       
	      if (in_argoneut==1){
		int  no_tracks_less_than_11_apart=0,no_tracks_less_than_15_apart=0,no_tracks_less_than_20_apart=0,no_tracks_less_than_30_apart=0,no_tracks_less_than_40_apart=0;	
	
		fno_minos_candidates->Fill(minos_tracks.size());
		 
		 
		no++;
		//std::cout<<"Predicted X= "<<x_predicted<<std::endl;
		 
		std::vector<double> all_diff_x,all_diff_y,distance,distance_minos_coord,all_diff_x_minos_coord,all_diff_y_minos_coord;
		 
		 
		//loop thru the number of tracks in argoneut for each event:
		std::cout<<"argoneut has "<<v_x_start_a.size()<<" tracks in event# "<<event<<std::endl;
		for(unsigned int t=0; t<v_x_start_a.size(); t++){

		  std::cout<<"ANALYSIS FOR ARGONEUT TRACK #"<<t<<std::endl;


		for(unsigned int i=0;i<minos_trk_x.size();i++){
		  std::cout<<"CHECKING MINOS TRACK #"<<i<<std::endl;
		  d_z=D-v_z_start_a[t]+minos_trk_z[i];
		  std::cout<<"d_z= "<<d_z<<" for trkvtxZ= "<<minos_trk_z[i]<<" (i.e minos_trk_z["<<i<<"],)"<<"  v_z_start_a["<<i<<" = "<<v_z_start_a[t]<<std::endl;
		  vector_length= d_z/v_Cosz[t];
		  
		  x_predicted=vector_length*v_Cosx[t]+v_x_start_a[t];
		  y_predicted=vector_length*v_Cosy[t]+v_y_start_a[t];
		  std::cout<<"x_predicted= "<<vector_length*v_Cosx[t]+v_x_start_a[t]<<" ,	y_predicted= "<<vector_length*v_Cosy[t]+v_y_start_a[t]<<std::endl;
		  
		  
		  //.......................................................................
		  //let's check if the predicted x&y coordinate actually lie on the minos detector. If they don't then there is no way of matching them to minos data and it would explain why we don't match some of them
		  if((x_predicted+x_offset)>297.7 || (x_predicted+x_offset)<-187.45 || (y_predicted-y_offset)>181.71 || (y_predicted-y_offset)<-199.29)
		    {std::cout<<" PREDICTED COORDINATES DON'T FALL ONTO MINOS !!! CANNOT MATCH THIS EVENT! because";
		      
		      if((x_predicted+x_offset)>297.7 || (x_predicted+x_offset)<-187.45) std::cout<<" x_predicted_in_minos_coord= "<<x_predicted+x_offset<<std::endl;
		      if((y_predicted-y_offset)>181.71 || (y_predicted-y_offset)<-199.29) std::cout<<" y_predicted_in_minos_coord= "<<y_predicted-y_offset<<std::endl;
		      e_w_predicted_tracks_outside_minos++;
		    }
		      
		  



		   
		  //	std::cout<<" minos_trk_x= "<< minos_trk_x[i]<<std::endl;
		  //std::cout<<"so I am filling with the diff= "<<minos_trk_x[i]-x_predicted<<std::endl;
		  fdiff_x->Fill( minos_trk_x[i]-x_predicted);
		  fdiff_x_minos_coord->Fill( minos_trk_x[i]-(x_predicted +x_offset ));
		   
		  all_diff_x.push_back( minos_trk_x[i]-x_predicted );
		  all_diff_x_minos_coord.push_back( minos_trk_x[i]-(x_predicted+x_offset ));

		  //NOW SAME THING FOR Y:
		    fdiff_y->Fill( minos_trk_y[i]-y_predicted);
		  fdiff_y_minos_coord->Fill( minos_trk_y[i]-(y_predicted-y_offset));
		   
		  fdiff_z->Fill( minos_trk_z[i]);
		   
		  fdiff_x_vs_diff_y->Fill( minos_trk_x[i]-x_predicted, minos_trk_y[i]-y_predicted);
		   
		  all_diff_y.push_back( minos_trk_y[i]-y_predicted );
		  all_diff_y_minos_coord.push_back( minos_trk_y[i]-(y_predicted-y_offset));
		   
		  distance.push_back(sqrt(all_diff_x[i]*all_diff_x[i] + all_diff_y[i]*all_diff_y[i]));
		  distance_minos_coord.push_back(sqrt(all_diff_x_minos_coord[i]*all_diff_x_minos_coord[i] + all_diff_y_minos_coord[i]*all_diff_y_minos_coord[i]));
		  std::cout<<"distance_minos_coord["<<i<<"]= "<<sqrt(all_diff_x_minos_coord[i]*all_diff_x_minos_coord[i] + all_diff_y_minos_coord[i]*all_diff_y_minos_coord[i])<<std::endl;
		  //lets test how many of them would fall within 11 cm circle, (i am doing it b/c i want to see if maybe I should not be just taking min # of the distance between argoneut and minos track BUT consider all the ones that fall within the specified distance (11cm)
		   
		  if(sqrt(all_diff_x_minos_coord[i]*all_diff_x_minos_coord[i] + all_diff_y_minos_coord[i]*all_diff_y_minos_coord[i])<11) no_tracks_less_than_11_apart++;
		  if(sqrt(all_diff_x_minos_coord[i]*all_diff_x_minos_coord[i] + all_diff_y_minos_coord[i]*all_diff_y_minos_coord[i])<15) no_tracks_less_than_15_apart++;
		  if(sqrt(all_diff_x_minos_coord[i]*all_diff_x_minos_coord[i] + all_diff_y_minos_coord[i]*all_diff_y_minos_coord[i])<20) no_tracks_less_than_20_apart++;
		  if(sqrt(all_diff_x_minos_coord[i]*all_diff_x_minos_coord[i] + all_diff_y_minos_coord[i]*all_diff_y_minos_coord[i])<30) no_tracks_less_than_30_apart++;
		  if(sqrt(all_diff_x_minos_coord[i]*all_diff_x_minos_coord[i] + all_diff_y_minos_coord[i]*all_diff_y_minos_coord[i])<40) no_tracks_less_than_40_apart++;





	       
		   
		  //std::cout<<" minos_cos_x= "<< minos_cos_x[i]<<" argoneut= "<<cosx_ar<<std::endl;
		  //std::cout<<"so I am filling with the diff= "<<minos_cos_x[i]-cosx_ar<<std::endl;
		  fcos_diff_x->Fill(minos_cos_x[i]-v_Cosx[t] );//??????
		  all_cosx.push_back(minos_cos_x[i]-v_Cosx[t] );
		  degree_diff_x.push_back(acos(minos_cos_x[i])-acos(v_Cosx[t]));
		   
		

	      }
		 
		fdiff_x_MIN->Fill( *( std::min_element(all_diff_x.begin(),all_diff_x.end() ) ));
		//	 delta_x2=(*( std::min_element(all_diff_x.begin(),all_diff_x.end() ) ))*(*( std::min_element(all_diff_x.begin(),all_diff_x.end() ) ));      
		//std::cout<<"no of events matched with maddalena is ***"<<no<<std::endl;
		 
		 
		//std::cout<<"Predicted Y= "<<y_predicted<<std::endl;
		 
		 
	
		 
	
		std::cout<<"No of tracks that are less than 11 cm apart from argoneut projected track = "<< no_tracks_less_than_11_apart<<std::endl;
		std::cout<<"No of tracks that are less than 15 cm apart from argoneut projected track = "<< no_tracks_less_than_15_apart<<std::endl;
		std::cout<<"No of tracks that are less than 20 cm apart from argoneut projected track = "<< no_tracks_less_than_20_apart<<std::endl;
		std::cout<<"No of tracks that are less than 30 cm apart from argoneut projected track = "<< no_tracks_less_than_30_apart<<std::endl;
		std::cout<<"No of tracks that are less than 40 cm apart from argoneut projected track = "<< no_tracks_less_than_40_apart<<std::endl;
		std::cout<<std::endl;
		fno_tracks_less_than_11_apart->Fill( no_tracks_less_than_11_apart);
		fno_tracks_less_than_15_apart->Fill( no_tracks_less_than_15_apart);
		fno_tracks_less_than_20_apart->Fill( no_tracks_less_than_20_apart);
		fno_tracks_less_than_30_apart->Fill( no_tracks_less_than_30_apart);
		fno_tracks_less_than_40_apart->Fill( no_tracks_less_than_40_apart);
		 
		 
		 
		fdiff_y_MIN->Fill( *( std::min_element(all_diff_y.begin(),all_diff_y.end() ) ));
		//	double delta_y2=( *( std::min_element(all_diff_y.begin(),all_diff_y.end() ) ))*( *( std::min_element(all_diff_y.begin(),all_diff_y.end() ) ));
		f_distance_MIN->Fill(*( std::min_element(distance.begin(),distance.end() ) ));
		//std::cout<<"normal d= "<<*( std::min_element(distance.begin(),distance.end() ) )<<"  ";

		f_distance_MIN_minos_coord->Fill(*( std::min_element(distance_minos_coord.begin(),distance_minos_coord.end() ) ));
		f_distance_MIN_minos_coord_squared->Fill((*( std::min_element(distance_minos_coord.begin(),distance_minos_coord.end() ) ))*(*( std::min_element(distance_minos_coord.begin(),distance_minos_coord.end() ) )));
		//std::cout<<"squared= "<<(*( std::min_element(distance_minos_coord.begin(),distance_minos_coord.end() ) ))*(*( std::min_element(distance_minos_coord.begin(),distance_minos_coord.end() ) ))<<"  ";

		double min_dist=*( std::min_element(distance_minos_coord.begin(),distance_minos_coord.end() ) );
		//	std::cout<<"min_dist= "<<min_dist<<std::endl;
		 
		int index=(std::min_element(distance_minos_coord.begin(),distance_minos_coord.end() ))-distance_minos_coord.begin();//get the position of the minimum value in the vector
		//std::cout<<"index= "<<index<<std::endl;
		 
		fdiff_x_combined_minos_coord->Fill(all_diff_x_minos_coord[index]);
		fdiff_y_combined_minos_coord->Fill(all_diff_y_minos_coord[index]);
		 
		//want to plot r vs p where r is the distance in the x-y plane. But I have q/p so i can't do it yet...I will do 2 plots: r vs fabs(q/p) and r vs(q/p):

		f_r_vs_abs_q_div_p->Fill(min_dist,fabs(minos_trkqp[index]));

	      f_r_vs_q_div_p->Fill(min_dist,minos_trkqp[index]);
	      std::cout<<"q/p= "<<minos_trkqp[index]<<std::endl;
	      if(minos_trkqp[index]==0){std::cout<<" q/p=0 !!!!!!!!!!!!!!!!!!!!!"<<std::endl;}

	      if(minos_trkqp[index]!=0){
		f_r_vs_abs_p->Fill(min_dist,1/(fabs(minos_trkqp[index])));
		f_r_vs_p->Fill(min_dist,1/(minos_trkqp[index]));
	      }


		if(min_dist<11)
		  {
		    std::cout<<"for r<11, q/p= "<<minos_trkqp[index]<<std::endl;

		     
		    f_distance_less_than_11->Fill(min_dist);
fdiff_x_less_than11->Fill( minos_trk_x[index]-x_predicted);
fdiff_y_less_than11->Fill( minos_trk_y[index]-y_predicted);

f_r_vs_abs_q_div_p_cut_11->Fill(min_dist,fabs(minos_trkqp[index]));

	      f_r_vs_q_div_p_cut_11->Fill(min_dist,minos_trkqp[index]);

if(minos_trkqp[index]!=0){
		
		f_p_vs_r_cut11->Fill(min_dist,1/(minos_trkqp[index]));

//make momentum distribution for postive and negative muons:

		if(minos_trkqp[index]>0){
		  
		  f_p_for_positive_cut11->Fill(fabs(1/(minos_trkqp[index])));
		}

	if(minos_trkqp[index]<0){
		  
		  f_p_for_negative_cut11->Fill(fabs(1/(minos_trkqp[index])));
		}




 }







		  }
		all_diff_x.clear();
		all_diff_y.clear();
		distance.clear();
		distance_minos_coord.clear();
		all_diff_x_minos_coord.clear();
		all_diff_y_minos_coord.clear();
		//-------------------------------------------------
		 
		ftheta_x_argoneut->Fill((acos(v_Cosx[t]))*180.0/pi);
		ftheta_y_argoneut->Fill((acos(v_Cosy[t]))*180.0/pi);
		ftheta_x_minos->Fill((acos(minos_cos_x[index]))*180.0/pi);
		ftheta_y_minos->Fill((acos(minos_cos_y[index]))*180.0/pi);
	
		//-----------------------------------------------------------------------------------------------------------------------------------------------------

	
		 
		 
		fcos_diff_x_MIN->Fill( *( std::min_element( all_cosx.begin(), all_cosx.end() ) ));
		fcos_diff_x_combined_minos_coord->Fill(all_cosx[index]);
		//fdegree_diff_x->Fill((degree_diff_x[index]*180)/pi);
		 
		all_cosx.clear();
		 
		 
		 
		std::vector<double> all_cosy,degree_diff_y;
		for(unsigned int i=0;i<minos_cos_y.size();i++){
		
		  //std::cout<<" minos_cos_y= "<< minos_cos_y[i]<<" argoneut= "<<cosy_ar<<std::endl;
		  //std::cout<<"so I am filling with the diff= "<<minos_cos_y[i]-cosy_ar<<std::endl;
		  fcos_diff_y->Fill(minos_cos_y[i]-v_Cosy[t] );
		  fcosx_vs_cosy->Fill(minos_cos_x[i]-v_Cosx[t],minos_cos_y[i]-v_Cosy[t]);
		
		  all_cosy.push_back(minos_cos_y[i]-v_Cosx[t] );
		  degree_diff_y.push_back(acos(minos_cos_y[i])-acos(v_Cosy[t]));
		
		}
		fcos_diff_y_MIN->Fill( *( std::min_element( all_cosy.begin(), all_cosy.end() ) ));
		fcos_diff_y_combined_minos_coord->Fill(all_cosy[index]);
		//fdegree_diff_y->Fill((degree_diff_y[index]*180)/pi);
	      
		all_cosy.clear();
		//ok, now let's see how many events we'll get when we put all the cuts together:
		//.................................................................
	      
		if((fabs(((degree_diff_y[index])*180.0)/pi)<7) && (fabs(((degree_diff_x[index])*180.0)/pi)<16)){
		  //std::cout<<"Events that pass degree cuts are: "<<events[index]<<" ";
		  no_degree_cut++;
		
		}
		if(min_dist<11){
		  //std::cout<<"Events that pass distance cut are: "<<events[index]<<" ";
		  no_distance_cut++;
		  fdegree_diff_x->Fill((degree_diff_x[index]*180)/pi);
		  fdegree_diff_y->Fill((degree_diff_y[index]*180)/pi);
		  ftime_diff_tms->Fill(time_diff_tms[index]);
		  ftime_diff_spilltime->Fill(time_diff_spilltime[index]);
		  minos_tracks[index].SetMatched(1);
		  ftrkqp_single_track->Fill(minos_trkqp[index]);
		  ftrkE_single_track->Fill(minos_trkE[index]);
		  ftrkErange_single_track->Fill(minos_trkErange[index]);
		  if(minos_trkqp[index]<0){
		    ftrkE_single_track_neg->Fill(minos_trkE[index]);

		  }


		  if(minos_trkqp[index]>0){
		    ftrkE_single_track_pos->Fill(minos_trkE[index]);

		  }
		
		}
	      
		if(min_dist<11 && (fabs(((degree_diff_y[index])*180.0)/pi)<7) && (fabs(((degree_diff_x[index])*180.0)/pi)<16) )
		  {
		    no_all_cuts++;
		    //std::cout<<"Events that pass all the cuts are: "<<events[index]<<" ";
		    ftime_diff_spilltime_all_cuts->Fill(time_diff_spilltime[index]);
		    f_distance_all_cuts->Fill(min_dist);
		    fdegree_diff_x_all_cuts->Fill((degree_diff_x[index]*180)/pi);
		    fdegree_diff_y_all_cuts->Fill((degree_diff_y[index]*180)/pi);
		  
		  
		  
		  }
	      
	      
		//.................................................................
	      
	      
	      
	      
	      
		for(unsigned int i=0;i<minos_cos_z.size();i++){
		
		  //	std::cout<<" minos_cos_z= "<< minos_cos_z[i]<<" argoneut= "<<cosz_ar<<std::endl;
		  //	std::cout<<"so I am filling with the diff= "<<minos_cos_z[i]-cosz_ar<<std::endl;
		  fcos_diff_z->Fill(minos_cos_z[i]-v_Cosz[t] );
		  fcosx_vs_cosz->Fill(minos_cos_x[i]-v_Cosx[t],minos_cos_z[i]-v_Cosz[t]);
		  fcosy_vs_cosz->Fill(minos_cos_y[i]-v_Cosy[t],minos_cos_z[i]-v_Cosz[t]);
		
		}
		//................................................................................
	      
	// 	for(unsigned int i=0;i<minos_trk_x.size();i++){
// 		  fdiff_x_m26->Fill( minos_trk_x[i]-x_predicted_m26);
// 		  fdiff_x_m24->Fill( minos_trk_x[i]-x_predicted_m24);
// 		  fdiff_x_m22->Fill( minos_trk_x[i]-x_predicted_m22);
// 		  fdiff_x_m20->Fill( minos_trk_x[i]-x_predicted_m20);
// 		  fdiff_x_m18->Fill( minos_trk_x[i]-x_predicted_m18);
// 		  fdiff_x_m16->Fill( minos_trk_x[i]-x_predicted_m16);
// 		  fdiff_x_m14->Fill( minos_trk_x[i]-x_predicted_m14);
// 		  fdiff_x_m12->Fill( minos_trk_x[i]-x_predicted_m12);
// 		  fdiff_x_m10->Fill( minos_trk_x[i]-x_predicted_m10);
// 		  fdiff_x_m8->Fill( minos_trk_x[i]-x_predicted_m8);
// 		  fdiff_x_m6->Fill( minos_trk_x[i]-x_predicted_m6);
// 		  fdiff_x_m4->Fill( minos_trk_x[i]-x_predicted_m4);
// 		  fdiff_x_m2->Fill( minos_trk_x[i]-x_predicted_m2);
// 		  fdiff_x_p2->Fill( minos_trk_x[i]-x_predicted_p2);
// 		  fdiff_x_p4->Fill( minos_trk_x[i]-x_predicted_p4);
// 		  fdiff_x_p6->Fill( minos_trk_x[i]-x_predicted_p6);
// 		  fdiff_x_p8->Fill( minos_trk_x[i]-x_predicted_p8);
// 		  fdiff_x_p10->Fill( minos_trk_x[i]-x_predicted_p10);
// 		  fdiff_x_p12->Fill( minos_trk_x[i]-x_predicted_p12);
// 		  fdiff_x_p14->Fill( minos_trk_x[i]-x_predicted_p14);
// 		  fdiff_x_p16->Fill( minos_trk_x[i]-x_predicted_p16);
// 		  fdiff_x_p18->Fill( minos_trk_x[i]-x_predicted_p18);
// 		  fdiff_x_p20->Fill( minos_trk_x[i]-x_predicted_p20);
// 		  fdiff_x_p22->Fill( minos_trk_x[i]-x_predicted_p22);
// 		  fdiff_x_p24->Fill( minos_trk_x[i]-x_predicted_p24);
// 		  fdiff_x_p26->Fill( minos_trk_x[i]-x_predicted_p26);

// 		}
	      
		for(unsigned int i=0;i<minos_trk_y.size();i++){
		  // fdiff_y_m26->Fill( minos_trk_y[i]-y_predicted_m26);
// 		  fdiff_y_m24->Fill( minos_trk_y[i]-y_predicted_m24);
// 		  fdiff_y_m22->Fill( minos_trk_y[i]-y_predicted_m22);
// 		  fdiff_y_m20->Fill( minos_trk_y[i]-y_predicted_m20);
// 		  fdiff_y_m18->Fill( minos_trk_y[i]-y_predicted_m18);
// 		  fdiff_y_m16->Fill( minos_trk_y[i]-y_predicted_m16);
// 		  fdiff_y_m14->Fill( minos_trk_y[i]-y_predicted_m14);
// 		  fdiff_y_m12->Fill( minos_trk_y[i]-y_predicted_m12);
// 		  fdiff_y_m10->Fill( minos_trk_y[i]-y_predicted_m10);
// 		  fdiff_y_m8->Fill( minos_trk_y[i]-y_predicted_m8);
// 		  fdiff_y_m6->Fill( minos_trk_y[i]-y_predicted_m6);
// 		  fdiff_y_m4->Fill( minos_trk_y[i]-y_predicted_m4);
// 		  fdiff_y_m2->Fill( minos_trk_y[i]-y_predicted_m2);
// 		  fdiff_y_p2->Fill( minos_trk_y[i]-y_predicted_p2);
// 		  fdiff_y_p4->Fill( minos_trk_y[i]-y_predicted_p4);
// 		  fdiff_y_p6->Fill( minos_trk_y[i]-y_predicted_p6);
// 		  fdiff_y_p8->Fill( minos_trk_y[i]-y_predicted_p8);
// 		  fdiff_y_p10->Fill( minos_trk_y[i]-y_predicted_p10);
// 		  fdiff_y_p12->Fill( minos_trk_y[i]-y_predicted_p12);
// 		  fdiff_y_p14->Fill( minos_trk_y[i]-y_predicted_p14);
// 		  fdiff_y_p16->Fill( minos_trk_y[i]-y_predicted_p16);
// 		  fdiff_y_p18->Fill( minos_trk_y[i]-y_predicted_p18);
// 		  fdiff_y_p20->Fill( minos_trk_y[i]-y_predicted_p20);
// 		  fdiff_y_p22->Fill( minos_trk_y[i]-y_predicted_p22);
// 		  fdiff_y_p24->Fill( minos_trk_y[i]-y_predicted_p24);
// 		  fdiff_y_p26->Fill( minos_trk_y[i]-y_predicted_p26);



		  ftrkqp->Fill(minos_trkqp[i]);
		
		  ftrkE->Fill(minos_trkE[i]);
		  ftrkErange->Fill(minos_trkErange[i]);
		
		
		  // 	if(minos_trk_index[i]==0){
		  // 	  //std::cout<<"filling with time diff="<<time_diff_tms[i]<<std::endl;
		  // 	ftime_diff_tms->Fill(time_diff_tms[i]);
		  // 	ftime_diff_spilltime->Fill(time_diff_spilltime[i]);
		
		  // 	}
		}
	      
	      
		
	      
		//minos_tracks[index]->SetMatched(1);
		//.................................................................................
	      
	      
	      
		}
	      
	}//if (in_argoneut==1)
	       
	       
    }// if( pDIR=opendir("/argoneut/app/users/kpartyka/larsoft/MINOS") )
	   
	   
	   
	   
}//match between DAQ and matched file
       
}//while in matched file
   
  matchedfile.close();
   
  minos_trk_x.clear();
  minos_trk_y.clear();
  minos_trkqp.clear();
  minos_trkE.clear();
  minos_trkErange.clear();
  minos_trk_index.clear();
  time_diff_tms.clear();
  time_diff_spilltime.clear();
  minos_tracks.clear();
  v_x_start_a.clear();
  v_y_start_a.clear();
  v_z_start_a.clear();
  v_Cosx.clear();
  v_Cosy.clear();
  v_Cosz.clear();
   
  if(in_matching==1)no_events_daq_pot++;
  if(pot_matched==1)no_e_pot_matched++;
  if(pot_and_z_matched==1)no_e_pot_and_z_matched++;
  if(pot_and_time_matched==1)no_e_pot_and_time_matched++;
   
   
   
   
  //std::cout<<"No of completely matched events(with Maddalena's file)= "<< no_matched_tracks<<std::endl;
  std::cout<<"Total no of events from DAQ and matched with our matched  file is: "<< no_events_daq<<" out of these " <<	no_events_daq_pot<< " events we match by requring same POTs, time cut of 60s and z cut of 10cm. Out of these "<<no<<" is matched also with Argoneut's tracks that we have in dircos file."<<std::endl; 
   
  std::cout<<"No of events matched just by same POTs: "<<no_e_pot_matched  <<std::endl;
  std::cout<<"No of events matched by same POTs and Z cut: "<<no_e_pot_and_z_matched  <<std::endl;
  std::cout<<"No of events matched by same POTs and time cut: "<<no_e_pot_and_time_matched <<std::endl;
  std::cout<<"No of tracks left after requiring them to be within 11 cm radius of each other on x-y plane is "<<no_distance_cut <<" (compare this with "<<no_matched_tracks<<" )"<<std::endl;
  std::cout<<"No of events left after requiring to be within 7 degrees in y and 16 degrees in x is " <<no_degree_cut <<std::endl;
  std::cout<<"No of events left after requiring distance and degree cuts is " <<no_all_cuts <<std::endl;
  std::cout<<"No of events that were calculated as having tracks that are outside of Minos' reach = "<< e_w_predicted_tracks_outside_minos<<std::endl;
   
   
       
   
   
  return;
}
}
