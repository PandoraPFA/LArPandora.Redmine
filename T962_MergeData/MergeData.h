////////////////////////////////////////////////////////////////////////
//
// MergeData Class
// Merge PMT/POT/MINOS information into DAQ480 data,
// if this info. exists for a given event.
// Kinga.partyka@yale.edu
// soderber@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef MERGEDATA_H
#define MERGEDATA_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "RawData/raw.h"


#include <vector>
#include <string>
#include <time.h>


///Merge different data streams 
class TH1F;
class TH2F;



namespace merge {

  class MergeData : public edm::EDProducer {

  public:
          
    explicit MergeData(edm::ParameterSet const& pset); 
    virtual ~MergeData();
    void produce(edm::Event& evt, edm::EventSetup const&);
    void beginJob(edm::EventSetup const&);
 
   

  private:
    TH1F* fdiff_x;
    TH1F* fdiff_y;
    TH1F* fdiff_z;
    TH1F* fcos_diff_x;
    TH1F* fcos_diff_y;
    TH1F* fcos_diff_z;
    TH1F* fcos_diff_x_MIN;
    TH1F* fcos_diff_y_MIN;
    TH1F* 	ftrkqp;
    TH1F *ftrkE;
    TH1F *ftrkErange;
    TH1F* 	ftrkqp_single_track;
    TH1F *ftrkE_single_track;
    TH1F *ftrkErange_single_track;
    TH1F *ftime_diff_tms;
    TH1F *ftime_diff_spilltime;
    TH1F *f_distance_MIN;
    TH1F *f_distance_MIN_minos_coord;
    TH1F *f_distance_MIN_minos_coord_squared;
    TH1F *f_distance_less_than_11;
    TH1F *fdiff_x_minos_coord;
    TH1F *fdiff_y_minos_coord;
    TH1F *fdegree_diff_x;
    TH1F *fdegree_diff_y;
    TH1F *ftime_diff_spilltime_all_cuts;
    TH1F *f_distance_all_cuts;
    TH1F *fdegree_diff_x_all_cuts;
    TH1F *fdegree_diff_y_all_cuts;




    TH1F *fdiff_x_MIN;
    TH1F *fdiff_y_MIN;
    TH1F* fdiff_x_m26;
    TH1F* fdiff_y_m26;
    TH1F* fdiff_x_m24;
    TH1F* fdiff_y_m24;
    TH1F* fdiff_x_m22;
    TH1F* fdiff_y_m22;
    TH1F* fdiff_x_m20;
    TH1F* fdiff_y_m20;
    TH1F* fdiff_x_m18;
    TH1F* fdiff_y_m18;
    TH1F* fdiff_x_m16;
    TH1F* fdiff_y_m16;
    TH1F* fdiff_x_m14;
    TH1F* fdiff_y_m14;
    TH1F* fdiff_x_m12;
    TH1F* fdiff_y_m12;
    TH1F* fdiff_x_m10;
    TH1F* fdiff_y_m10;
    TH1F* fdiff_x_m8;
    TH1F* fdiff_y_m8;
    TH1F* fdiff_x_m6;
    TH1F* fdiff_y_m6;
    TH1F* fdiff_x_m4;
    TH1F* fdiff_y_m4;
    TH1F* fdiff_x_m2;
    TH1F* fdiff_y_m2;

    TH1F* fdiff_x_p26;
    TH1F* fdiff_y_p26;
    TH1F* fdiff_x_p24;
    TH1F* fdiff_y_p24;
    TH1F* fdiff_x_p22;
    TH1F* fdiff_y_p22;
    TH1F* fdiff_x_p20;
    TH1F* fdiff_y_p20;
    TH1F* fdiff_x_p18;
    TH1F* fdiff_y_p18;
    TH1F* fdiff_x_p16;
    TH1F* fdiff_y_p16;
    TH1F* fdiff_x_p14;
    TH1F* fdiff_y_p14;
    TH1F* fdiff_x_p12;
    TH1F* fdiff_y_p12;
    TH1F* fdiff_x_p10;
    TH1F* fdiff_y_p10;
    TH1F* fdiff_x_p8;
    TH1F* fdiff_y_p8;
    TH1F* fdiff_x_p6;
    TH1F* fdiff_y_p6;
    TH1F* fdiff_x_p4;
    TH1F* fdiff_y_p4;
    TH1F* fdiff_x_p2;
    TH1F* fdiff_y_p2;
    TH2F* fcosx_vs_cosy;
    TH2F* fcosx_vs_cosz;
    TH2F* fcosy_vs_cosz;
    TH2F* fdiff_x_vs_diff_y;

    TH1F *fdiff_x_combined_minos_coord;
    TH1F *fdiff_y_combined_minos_coord;
    TH1F *fcos_diff_x_combined_minos_coord;
    TH1F *fcos_diff_y_combined_minos_coord;
    TH1F *fno_minos_candidates;
    TH1F *fno_tracks_less_than_11_apart;
    TH1F *fno_tracks_less_than_15_apart;
    TH1F *fno_tracks_less_than_20_apart;
    TH1F *fno_tracks_less_than_30_apart;
    TH1F *fno_tracks_less_than_40_apart;
    TH1F *ftheta_x_argoneut;
    TH1F *ftheta_y_argoneut;
    TH1F *ftheta_x_minos;
    TH1F *ftheta_y_minos;
    TH1F *fdiff_x_less_than11;
    TH1F *fdiff_y_less_than11;
    TH2F *f_r_vs_abs_q_div_p;
    TH2F *f_r_vs_q_div_p;
    TH2F *f_r_vs_abs_q_div_p_cut_11;
    TH2F *f_r_vs_q_div_p_cut_11;
    TH2F *f_r_vs_abs_p;
    TH2F *f_r_vs_p;
    TH2F *f_p_vs_r_cut11;
    TH1F *ftrkE_single_track_neg;
    TH1F *ftrkE_single_track_pos;
    TH1F *f_p_for_positive_cut11;
    TH1F *f_p_for_negative_cut11;
  protected: 

   
    
   
    void MergeBeam( std::auto_ptr<std::vector <raw::BeamInfo> >Beam_coll);                 ///method to merge beam data
    void MergePMT(std::auto_ptr<std::vector <raw::Paddles> >Paddles_coll);                 ///method to merge pmt data
    void MergeMINOS(std::auto_ptr<std::vector<raw::MINOS> >Minos_coll);                 ///method to merge MINOS data

   
	edm::Ptr<raw::DAQHeader> fdaq;
    raw::BeamInfo  fBeamInfo;
    raw::Paddles   fPaddles;
    
	 
    ///<parameters to set
   std::string  fdaq_modulelabel;               ///< folder for input 
   std::string file;	
   
    bool foundbeaminfo;
    bool foundpaddleinfo;
    bool foundminosinfo;
    
	 
  }; // class MergeData

}

#endif // MERGEDATA_H
