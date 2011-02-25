////////////////////////////////////////////////////////////////////////
/// \file  T962G4Ana.cxx
/// \brief Use Geant4 to run the LArSoft detector simulation
///
/// \version $Id: T962G4.cxx,v 1.22 2010/07/20 06:08:30 bjpjones Exp $
/// \author  joshua.spitz@yale.edu
/// \author echurch@fnal.gov
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

// LArSoft Includes
#include "Simulation/SimListUtils.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Geometry/geo.h"
#include "T962/T962G4Ana/T962G4Ana.h"

// ROOT includes
#include "TH1.h"

// C++ Includes
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstring>
#include <sys/stat.h>

namespace T962G4 {

  //-----------------------------------------------------------------------
  // Constructor
  T962G4Ana::T962G4Ana(fhicl::ParameterSet const& pset) :
    fG4ModuleLabel(pset.get< std::string >("GeantModuleLabel")),
    fm_run(0), 
    fm_event(0),
    fm_trackID(-1),
    fm_pdgcode(-1),
    fm_mass(-1.),
    fm_energylost(0.),
    fm_deflectionangle(0.),
    fm_tpcexit_px(0.),
    fm_tpcexit_py(0.),
    fm_tpcexit_pz(0.),
    fm_tpcexit_x(0.),
    fm_tpcexit_y(0.),
    fm_tpcexit_z(0.),
    fm_tpcexit_energy(0.),
    fm_minosenter_px(0.),
    fm_minosenter_py(0.),
    fm_minosenter_pz(0.),
    fm_minosenter_x(0.),
    fm_minosenter_y(0.),
    fm_minosenter_z(0.),
    fm_minosenter_energy(0.)
  {
  }

  //-----------------------------------------------------------------------
  // Destructor
  T962G4Ana::~T962G4Ana() 
  {
  }

  //-----------------------------------------------------------------------
  void T962G4Ana::beginJob()
  {
  art::ServiceHandle<art::TFileService> tfs; 
  ftree= tfs->make<TTree>("MINOSinTree","MINOSinTree");
  ftree->Branch("run", &fm_run, "run/I");
  ftree->Branch("event", &fm_event, "event/I");
  ftree->Branch("trackID", &fm_trackID, "trackID/I");
  ftree->Branch("pdgcode", &fm_pdgcode, "pdgcode/I");
  ftree->Branch("mass", &fm_mass, "mass/F");
  ftree->Branch("energylost", &fm_energylost, "energylost/F");
  ftree->Branch("deflectionangle", &fm_deflectionangle, "deflectionangle/F");
  ftree->Branch("tpcexit_px", &fm_tpcexit_px, "tpcexit_px/F");
  ftree->Branch("tpcexit_py", &fm_tpcexit_py, "tpcexit_py/F");
  ftree->Branch("tpcexit_pz", &fm_tpcexit_pz, "tpcexit_pz/F");
  ftree->Branch("tpcexit_x", &fm_tpcexit_x, "tpcexit_x/F");
  ftree->Branch("tpcexit_y", &fm_tpcexit_y, "tpcexit_y/F");
  ftree->Branch("tpcexit_z", &fm_tpcexit_z, "tpcexit_z/F");
  ftree->Branch("tpcexit_energy", &fm_tpcexit_energy, "tpcexit_energy/F");
  ftree->Branch("minosenter_px", &fm_minosenter_px, "minosenter_px/F");
  ftree->Branch("minosenter_py", &fm_minosenter_py, "minosenter_py/F");
  ftree->Branch("minosenter_pz", &fm_minosenter_pz, "minosenter_pz/F");
  ftree->Branch("minosenter_x", &fm_minosenter_x, "minosenter_x/F");
  ftree->Branch("minosenter_y", &fm_minosenter_y, "minosenter_y/F");
  ftree->Branch("minosenter_z", &fm_minosenter_z, "minosenter_z/F");
  ftree->Branch("minosenter_energy", &fm_minosenter_energy, "minosenter_energy/F");  
  }

  //-----------------------------------------------------------------------
  void T962G4Ana::analyze(const art::Event& evt) 
  {

    //get the list of particles from this event
    sim::ParticleList plist = sim::SimListUtils::GetParticleList(evt, fG4ModuleLabel);
    art::ServiceHandle<geo::Geometry> geom;

    // get the particles from the event handle
    art::Handle< std::vector<sim::Particle> > parHandle;
    evt.getByLabel(fG4ModuleLabel, parHandle);

    art::PtrVector<sim::Particle> pvec;
    for(unsigned int i = 0; i < parHandle->size(); ++i){
      art::Ptr<sim::Particle> p(parHandle, i);      
      pvec.push_back(p);
    }

    fm_run=evt.id().run();
    fm_event=evt.id().event(); 
    fm_trackID=-1;
    fm_tpcexit_px=0.;
    fm_tpcexit_py=0.;
    fm_tpcexit_pz=0.;
    fm_tpcexit_x=0.;
    fm_tpcexit_y=0.;
    fm_tpcexit_z=0.;
    fm_tpcexit_energy=0.;
    fm_minosenter_energy=0.;
    fm_minosenter_px=0.;
    fm_minosenter_py=0.;
    fm_minosenter_pz=0.;
    fm_minosenter_x=0.;
    fm_minosenter_y=0.;
    fm_minosenter_z=0.;
    fm_pdgcode=0;
    fm_mass=0.;
    fm_energylost=0.;
    fm_deflectionangle=0.;
    float x_offset=116.9;//cm
    float y_offset=-20.28;//cm
    float z_offset=-147.1;//cm
    
    int numprimaries=0;
    for(unsigned int i = 0; i < pvec.size(); ++i)
    {
     if(pvec[i]->Process()!="primary")
     continue;     
    simb::MCTrajectory trajectory = pvec[i]->Trajectory();
    int numberofpoints= pvec[i]->NumberTrajectoryPoints();    
    	for(int j=1; j<numberofpoints; j++)
    	{
    	TLorentzVector prevposition = trajectory.Position(j-1);
    	TLorentzVector position = trajectory.Position(j);
    
    	if(pvec[i]->Process()=="primary"&&prevposition.Z()<154.22&&position.Z()>154.22)
    	numprimaries++; 
    	}
    
    }
    ofstream myfile;
    myfile.open("MINOSin.txt", std::ios::out | std::ios::app);
    myfile << "________________________________________________________________________________\n";    
    myfile << "***** HEPEVT Common Event#: "<<fm_event<<", "<<numprimaries<<" particles (max 4000) ***** Double Precision\n";
    myfile << "4-byte integers, 8-byte floating point numbers, 4000-allocated entries.\n";
    myfile << "Indx Stat Par- chil-       (  P_x,       P_y,       P_z,    Energy,       M )\n";
    myfile << "      ID  ents dren    Prod (   X,         Y,         Z,        cT)      [mm]\n";
    myfile << "________________________________________________________________________________\n";
    
    
    
    for(unsigned int i = 0; i < pvec.size(); ++i){
    
    if(pvec[i]->Process()!="primary")
    continue;
    
    simb::MCTrajectory trajectory = pvec[i]->Trajectory();
    fm_pdgcode=pvec[i]->PdgCode();       
    fm_mass=pvec[i]->Mass();
    fm_trackID=i;
    
    int numberofpoints= pvec[i]->NumberTrajectoryPoints();    
    for(int j=1; j<numberofpoints; j++)
    {
      TLorentzVector prevposition = trajectory.Position(j-1);
      TLorentzVector position = trajectory.Position(j);
      TLorentzVector prevmomentum = trajectory.Momentum(j-1);
      TLorentzVector momentum = trajectory.Momentum(j);      
      TLorentzVector initmomentum = trajectory.Momentum(0);
    
       if(       (((prevposition.X()<0.&&position.X()>0.)||(prevposition.X()>0.&&position.X()<0.))&&(prevposition.Y()<20.2&&prevposition.Y()>-19.8&&prevposition.Z()>0.&&prevposition.Z()<90.))
       ||       (((prevposition.X()<47.&&position.X()>47.)||(prevposition.X()>47.&&position.X()<47.))&&(prevposition.Y()<20.2&&prevposition.Y()>-19.8&&prevposition.Z()>0.&&prevposition.Z()<90.))
       ||       (((prevposition.Y()<20.2&&position.Y()>20.2)||(prevposition.Y()>20.2&&position.Y()<20.2))&&(prevposition.X()<47.&&prevposition.X()>0.&&prevposition.Z()>0.&&prevposition.Z()<90.))
       ||       (((prevposition.Y()<-19.8&&position.Y()>-19.8)||(prevposition.Y()>-19.8&&position.Y()<-19.8))&&(prevposition.X()<47.&&prevposition.X()>0.&&prevposition.Z()>0.&&prevposition.Z()<90.))
       ||       (((prevposition.Z()<0.&&position.Z()>0.)||(prevposition.Z()>0.&&position.Z()<0.))&&(prevposition.Y()<20.2&&prevposition.Y()>-19.8&&prevposition.X()>0.&&prevposition.X()<47.))
       ||       (((prevposition.Z()<90.&&position.Z()>90.)||(prevposition.Z()>90.&&position.Z()<90.))&&(prevposition.Y()<20.2&&prevposition.Y()>-19.8&&prevposition.X()>0.&&prevposition.X()<47.))
       )
       {    
       fm_tpcexit_x=prevposition.X();
       fm_tpcexit_y=prevposition.Y();
       fm_tpcexit_z=prevposition.Z();
       fm_tpcexit_px=prevmomentum.Px();
       fm_tpcexit_py=prevmomentum.Py();
       fm_tpcexit_pz=prevmomentum.Pz();
       fm_tpcexit_energy=prevmomentum.E(); 
       }
   //  if(prevposition.Z()<154.22&&position.Z()>154.22)
//           {
       fm_minosenter_x=prevposition.X()+x_offset;
       fm_minosenter_y=prevposition.Y()+y_offset;
       fm_minosenter_z=prevposition.Z()+z_offset;
       fm_minosenter_px=prevmomentum.Px();
       fm_minosenter_py=prevmomentum.Py();
       fm_minosenter_pz=prevmomentum.Pz();
       fm_minosenter_energy=prevmomentum.E();
       
       

//     myfile << "   1   +0    0    0    (        0,         0,         0,         0,         0)\n";
//     myfile << "       +0    0    0    (        0,         0,         0,         0)\n";
    myfile <<"   1   +1    0    0    (    "<<fm_minosenter_px<<",  "<<fm_minosenter_py<<",  "<<fm_minosenter_pz<<",  "<<fm_minosenter_energy<<",  "<<fm_mass<<")\n";
    myfile <<"      +"<<fm_pdgcode<<"    0    0    (      "<<fm_minosenter_x*10.<<",       "<<fm_minosenter_y*10.<<",       "<<fm_minosenter_z*10.<<",  2.64e+06)\n";
    myfile << "________________________________________________________________________________\n";         

    
       
//          
        std::cout<<fm_pdgcode<<" at z position "<<position.Z()<<" lost "<<((prevmomentum.E()-momentum.E())*1000.)<<"MeV in "<<(position.Z()-prevposition.Z())<<"cm which is "<<((prevmomentum.E()-momentum.E())*1000.)/(position.Z()-prevposition.Z())<< "MeV/cm VolumeName is " << geom->VolumeName(position.Vect())<<" "<<geom->MaterialName(position.Vect()) << std::endl;      
      
       
       fm_energylost=fm_tpcexit_energy-fm_minosenter_energy;
       fm_deflectionangle= TMath::ACos(((fm_tpcexit_px*fm_minosenter_px)+(fm_tpcexit_py*fm_minosenter_py)+(fm_tpcexit_pz*fm_minosenter_pz))/(sqrt(pow(fm_tpcexit_px,2)+pow(fm_tpcexit_py,2)+pow(fm_tpcexit_pz,2))*sqrt(pow(fm_minosenter_px,2)+pow(fm_minosenter_py,2)+pow(fm_minosenter_pz,2))));

       
      // std::cout<<"exit-enter "<<(fm_tpcexit_energy-fm_minosenter_energy)*1000.<<" "<<fm_deflectionangle<<std::endl;
          // }
    }
    ftree->Fill(); 

  }
myfile.close(); 
return;  
}

} // namespace T962G4

