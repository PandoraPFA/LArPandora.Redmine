////////////////////////////////////////////////////////////////////////
//
// RadArgon class
//
// msoderbe@syr.edu
//
// Make plots for through-going muons.
////////////////////////////////////////////////////////////////////////


#include <iostream>


// Framework includes
//   HOW MANY OF THESE DO I NEED? 
//   DO I NEED OTHERS?
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"

#include "T962/RadArgonAna/RadArgonAna.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "TMath.h"

#include "Filters/ChannelFilter.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

namespace radargon {
//-----------------------------------------------------------------------------
  RadArgonAna::RadArgonAna(fhicl::ParameterSet const& pset) :
    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
    fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
    fCalDataModuleLabel       (pset.get< std::string >("CalDataModuleLabel")      ),
    fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
    fHitSensitivity           (pset.get< std::string >("HitSensitivity")          ),
    fHits_label               (pset.get< std::string >("HitsModeulLabel")         ),
    ftmatch                   (pset.get<      int    >("TMatch")                  ),
    fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel")      )
  {

  }
    
//-----------------------------------------------------------------------------
  RadArgonAna::~RadArgonAna()
  {
  }

//-----------------------------------------------------------------------------
  void RadArgonAna::beginJob()
  {
    // get access to the TFile service  
    art::ServiceHandle<art::TFileService> tfs;

    // these are not set up correctly, yet
    // for now, just but place-holders for info we want
//    fPeakTime = tfs->make<TH1F>("fPeakTime","Peak Time of Hit", 200,-25.0,25.0);
//    fCharge = tfs->make<TH1F>("fCharge","Charge of Hit", 200,-25.0,25.0);
        
  }
    
//-----------------------------------------------------------------------------
  void RadArgonAna::analyze(const art::Event& evt) 
  {
    // get services (maybe only need geometry?)
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
////// grabbed from Track3Dreco.cxx 
    double timetick = 0.198;    //time sample in us
    double presamplings = 60.;  // what is this??? number of ticks prior to effective start?
    double plane_pitch = geom->PlanePitch(0,1);    //wire plane pitch in cm 
    double wire_pitch = geom->WirePitch(0,1,0);    //wire pitch in cm
    double Efield_drift = 0.5;    // Electric Field in the drift region in kV/cm
    double Efield_SI = 0.7;       // Electric Field between Shield and Induction planes in kV/cm
    double Efield_IC = 0.9;       // Electric Field between Induction and Collection planes in kV/cm
    double Temperature = 87.6;    // LAr Temperature in K
    double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature); //drift velocity in the drift region (cm/us)
    double driftvelocity_IC = larprop->DriftVelocity(Efield_IC,Temperature); //drift velocity between induction and collection (cm/us)
    double timepitch = driftvelocity*timetick;           //time sample (cm) 
    double tIC = plane_pitch/driftvelocity_IC/timetick;  //drift time between Induction and Collection planes (time samples)
////// end of Track3Dreco.cxx grab

    // declare cluster handle and fill
    art::Handle< std::vector<recob::Cluster> > ClusterHandle;
    evt.getByLabel(fClusterModuleLabel,ClusterHandle);

    // define some vectors to use in plane sorting
    //std::vector<art::PtrVector<recob::Hit>>  ind_Singles;
    art::PtrVector<recob::Hit>  ind_Singles;
    std::vector<unsigned int>                ind_SinglesCount;
    //std::vector<art::PtrVector<recob::Hit>>  col_Singles;
    art::PtrVector<recob::Hit>  col_Singles;
    std::vector<unsigned int>                col_SinglesCount;

    // scan all clusters
    for(unsigned int cluster_itr = 0 ; cluster_itr < ClusterHandle->size() ; cluster_itr++){ 
      art::Ptr<recob::Cluster> cluster(ClusterHandle,cluster_itr);
      art::PtrVector<recob::Hit> hitlist;

      // require clusters contain only one hit 
      if( hitlist.size() == 1 ){ 
        std::cout << " ONE HIT IN CLUSTER #" << cluster_itr << std::endl; 

        // sort single-hit clusters by plane, can adjust for larger number of hits with loop
        for(art::PtrVector<recob::Hit>::const_iterator aHit = hitlist.begin(); aHit != hitlist.end();  aHit++){
          unsigned int channel, plane, tpc, wire;
          channel =  (*aHit)->Channel();
          geom->ChannelToWire(channel,tpc,plane,wire);

          // fill singles vectors
          switch(plane){
            case 0:
              ind_Singles.push_back((*aHit));
              //ind_SinglesCount.push_back(cluster_itr);
              break;
            case 1:
              col_Singles.push_back((*aHit));
              //col_SingleCount.push_back(cluster_itr);
              break;
          } // fill singles vectors
        } // sort by plane
      } // require single hit
    } // scan all clusters

    // now i have ind/col_Singles which contain hit lists
    //   want to see if any of the hits on the collection plane match
    //   in time and have wires that cross in the induction plane
    unsigned int channel_col,channel_ind,wire_col,wire_ind,plane_col,plane_ind,tpc;
    for(art::PtrVector<recob::Hit>::const_iterator ind_hit = ind_Singles.begin(); ind_hit != ind_Singles.end(); ind_hit++){
      channel_ind = (*ind_hit)->Channel();
      double time_ind = (*ind_hit)->PeakTime() - presamplings;
      double y_pos, z_pos;
      for(art::PtrVector<recob::Hit>::const_iterator col_hit = col_Singles.begin(); col_hit != col_Singles.end(); col_hit++){
        channel_col = (*col_hit)->Channel();
        double time_col = (*col_hit)->PeakTime() - presamplings - tIC;

        bool TimesMatch = ( fabs(time_col-time_ind) < ftmatch * timepitch );
        if( geom->ChannelsIntersect(channel_ind,channel_col,y_pos,z_pos) && TimesMatch ) {
          // add to a vector?
          std::cout << "  We have a lone hit! " << std::endl; 
        }
      }

    }
  }//analyze
}//namespace
