////////////////////////////////////////////////////////////////////////
//
// RadArgon class
//
// patch@fnal.gov
//
// Make plots for radioactive argon decay.
////////////////////////////////////////////////////////////////////////

// Framework includes
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

#include "T962/MuonAna/MuonAna.h"
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
#include "T962/T962_Objects/MINOS.h"
#include "T962/T962_Objects/MINOSTrackMatch.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

    
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

    //fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
    //fHitSensitivity           (pset.get< std::string >("HitSensitivity")          ),
    //fHits_label               (pset.get< std::string >("HitsModuleLabel")         ),
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
    fPeakTime = tfs->make<TH1F>("fPeakTime","Peak Time of Hit", 200,-25.0,25.0);
    fCharge = tfs->make<TH1F>("fCharge","Charge of Hit", 200,-25.0,25.0);
    fSingleHitIndCharge = tfs->make<TH1F>("fSingleIndCharge","Charge of Induction Hit", 200,-25.0,25.0);  // these be ???
    fSingleHitColCharge = tfs->make<TH1F>("fSingleColCharge","Charge of Collection Hit", 200,-25.0,25.0); // what size should
    fTimeDifference= tfs->make<TH1F>("fTimeDifference","Time Difference between Hits", 250,0,50.0);

    //fChannel = tfs->make<TH1S>("fChannel","Channel Number of Hit",1000,0,999);
        
    fIndPlaneHits = tfs->make<TH1S>("fIndPlaneHits","Number of Hits in Induction Plane",4,0,3);
    fColPlaneHits = tfs->make<TH1S>("fColPlaneHits","Number of Hits in Collection Plane",4,0,3);
    fPlaneHits = tfs->make<TH2S>("fPlaneHits","Number of Hits in Induction and Collection Planes", 21, 0.0, 20.0, 21, 0.0, 20.0);
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
      hitlist = cluster->Hits();

      // require clusters contain nHits
      unsigned int nHits = 3; 
      if( hitlist.size() <= nHits ){ 
        std::cout << hitlist.size() <<  " HITS IN CLUSTER #" << cluster_itr << std::endl; 

        // sort single-hit clusters by plane, can adjust for larger number of hits with loop
        for(art::PtrVector<recob::Hit>::const_iterator aHit = hitlist.begin(); aHit != hitlist.end();  aHit++){
          unsigned int channel, plane, tpc, wire;
          channel =  (*aHit)->Channel();
          geom->ChannelToWire(channel,tpc,plane,wire);

          double peaktime=(*aHit)->PeakTime();
    //      fPeakTime->Fill(peaktime);
    //      fChannel->Fill(channel);

          // fill singles vectors
          switch(plane){
            case 0:
              ind_Singles.push_back((*aHit));
    //          fIndPlaneHits->Fill(hitlist.size());
              break;
            case 1:
              col_Singles.push_back((*aHit));
    //          fColPlaneHits->Fill(hitlist.size());
              break;
          } // fill singles vectors
        } // sort by plane
      } // require single hit
    } // scan all clusters


    // now i have ind/col_Singles which contain hit lists
    //   want to see if any of the hits on the collection plane match
    //   in time and have wires that cross in the induction plane
    unsigned int channel_col,channel_ind,wire_col,wire_ind,plane_col,plane_ind,tpc;

    // here i think there is a problem due to the fact that i'm looping over single hits... 
    // ... probably not collecting everything correctly
    for(art::PtrVector<recob::Hit>::const_iterator ind_hit = ind_Singles.begin(); ind_hit != ind_Singles.end(); ind_hit++){

      channel_ind = (*ind_hit)->Channel();

      double time_ind = (*ind_hit)->PeakTime() - presamplings;
      double y_pos, z_pos;

    //  double indHitCharge=(*ind_hit)->Charge();
    //  std::cout << " Charge of induction hit = " << indHitCharge << std::endl;
    //  fSingleHitIndCharge->Fill(indHitCharge);

      for(art::PtrVector<recob::Hit>::const_iterator col_hit = col_Singles.begin(); col_hit != col_Singles.end(); col_hit++){

        channel_col = (*col_hit)->Channel();
        double time_col = (*col_hit)->PeakTime() - presamplings - tIC;
        double time_diff = fabs(time_col - time_ind);
        bool TimesMatch = ( time_diff < ( ftmatch * timepitch ) );

        if( geom->ChannelsIntersect(channel_ind,channel_col,y_pos,z_pos) && TimesMatch ) {
          std::cout << "  We have a lone hit! " << std::endl; 
          fSingleHitColCharge->Fill((*col_hit)->Charge());
          fTimeDifference->Fill(time_diff);
        }

      }
    }
  }//analyze
}//namespace
