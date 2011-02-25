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
    fG4ModuleLabel(pset.get< std::string >("GeantModuleLabel"))
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

    int fm_run=evt.id().run();
    int fm_event=evt.id().event();   
    float tpcexit_energy=0.;
    float minosenter_energy=0.;
    float minosenter_px=0.;
    float minosenter_py=0.;
    float minosenter_pz=0.;
    float minosenter_x=0.;
    float minosenter_y=0.;
    float minosenter_z=0.;
    int pdgcode=0;
    float mass=0.;
    float x_offset=116.9;//cm
    float y_offset=-20.28;//cm
    float z_offset=-147.1;//cm
    for(unsigned int i = 0; i < pvec.size(); ++i){    
    std::cout<<pvec[i]->NumberTrajectoryPoints()<<std::endl;
    simb::MCTrajectory trajectory = pvec[i]->Trajectory();
    
    int numberofpoints= pvec[i]->NumberTrajectoryPoints();    
    for(int j=1; j<numberofpoints; j++)
    {
      TLorentzVector prevposition = trajectory.Position(j-1);
      TLorentzVector position = trajectory.Position(j);
      TLorentzVector prevmomentum = trajectory.Momentum(j-1);
      TLorentzVector momentum = trajectory.Momentum(j);      
      TLorentzVector initmomentum = trajectory.Momentum(0);
            
       if(
       (prevposition.X()<0.&&position.X()>0.)
       ||
       (prevposition.X()<47.&&position.X()>47.)
       ||
       (prevposition.Y()<20.2&&position.Y()>20.2)
       ||
       (prevposition.Y()<-19.8&&position.Y()>-19.8)
       ||
       (prevposition.Z()<0.&&position.Z()>0.)
       ||
       (prevposition.Z()<90.&&position.Z()>90.)
       )
       tpcexit_energy=momentum.E();

       if(prevposition.Z()<147.1&&position.Z()>147.1)
       {
       minosenter_energy=momentum.E();
       pdgcode=pvec[i]->PdgCode();       
       mass=pvec[i]->Mass();
       minosenter_x=position.X()+x_offset;
       minosenter_y=position.Y()+y_offset;
       minosenter_z=position.Z()+z_offset;
       minosenter_px=momentum.Px();
       minosenter_py=momentum.Py();
       minosenter_pz=momentum.Pz();
       ofstream myfile;
       
    myfile.open("MINOSin.txt", std::ios::out | std::ios::app);
    myfile << "________________________________________________________________________________\n";    
    myfile << "***** HEPEVT Common Event#: "<<fm_event<<", 1 particles (max 4000) ***** Double Precision\n";
    myfile << "4-byte integers, 8-byte floating point numbers, 4000-allocated entries.\n";
    myfile << "Indx Stat Par- chil-       (  P_x,       P_y,       P_z,    Energy,       M )\n";
    myfile << "      ID  ents dren    Prod (   X,         Y,         Z,        cT)      [mm]\n";
    myfile << "________________________________________________________________________________\n";
//     myfile << "   1   +0    0    0    (        0,         0,         0,         0,         0)\n";
//     myfile << "       +0    0    0    (        0,         0,         0,         0)\n";
    myfile <<"   1   +1    0    0    (    "<<minosenter_px<<",  "<<minosenter_py<<",  "<<minosenter_pz<<",  "<<minosenter_energy<<",  "<<mass<<")\n";
    myfile <<"      +"<<pdgcode<<"    0    0    (      "<<minosenter_x*10.<<",       "<<minosenter_y*10.<<",       "<<minosenter_z*10.<<",  2.64e+06)\n";
    myfile << "________________________________________________________________________________\n";         
    myfile.close();
    
       }
         
      std::cout<<"At z position "<<position.Z()<<" lost "<<((prevmomentum.E()-momentum.E())*1000.)<<"MeV in "<<(position.Z()-prevposition.Z())<<"cm which is "<<((prevmomentum.E()-momentum.E())*1000.)/(position.Z()-prevposition.Z())<< "MeV/cm VolumeName is " << geom->VolumeName(position.Vect())<<" "<<geom->MaterialName(position.Vect()) << std::endl;       
    }

    return;
  }
}

} // namespace T962G4

