////////////////////////////////////////////////////////////////////////
/// \file  T962G4Ana.cxx
/// \brief Use Geant4 to run the LArSoft detector simulation
///
/// \version $Id: T962G4.cxx,v 1.22 2010/07/20 06:08:30 bjpjones Exp $
/// \author  joshua.spitz@yale.edu
/// \author echurch@fnal.gov
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"

// LArSoft Includes
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCNeutrino.h"
#include "Simulation/ParticleList.h"
#include "MCCheater/BackTracker.h"
#include "Geometry/Geometry.h"
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
   static int event=1;
  //-----------------------------------------------------------------------
  // Constructor
  T962G4Ana::T962G4Ana(fhicl::ParameterSet const& pset) :
    fG4ModuleLabel(pset.get< std::string >("GeantModuleLabel")),
    fGenieGenModuleLabel             (pset.get< std::string >("GenieGenModuleLabel")),
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
    fm_minosenter_energy(0.),
    fm_init_px(0.),
    fm_init_py(0.),
    fm_init_pz(0.),
    fm_init_x(0.),
    fm_init_y(0.),
    fm_init_z(0.),
    fm_init_energy(0.),
    fm_offset_x(0.),
    fm_offset_y(0.),
    fm_offset_z(0.)
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
  ftree->Branch("init_px", &fm_init_px, "init_px/F");
  ftree->Branch("init_py", &fm_init_py, "init_py/F");
  ftree->Branch("init_pz", &fm_init_pz, "init_pz/F");
  ftree->Branch("init_x", &fm_init_x, "init_x/F");
  ftree->Branch("init_y", &fm_init_y, "init_y/F");
  ftree->Branch("init_z", &fm_init_z, "init_z/F");
  ftree->Branch("init_energy", &fm_init_energy, "init_energy/F"); 
  ftree->Branch("neutrino_px", &fm_neutrino_px, "neutrino_px/F");
  ftree->Branch("neutrino_py", &fm_neutrino_py, "neutrino_py/F");
  ftree->Branch("neutrino_pz", &fm_neutrino_pz, "neutrino_pz/F");
  ftree->Branch("neutrino_x", &fm_neutrino_x, "neutrino_x/F");
  ftree->Branch("neutrino_y", &fm_neutrino_y, "neutrino_y/F");
  ftree->Branch("neutrino_z", &fm_neutrino_z, "neutrino_z/F");
  ftree->Branch("neutrino_energy", &fm_neutrino_energy, "neutrino_energy/F");
  ftree->Branch("neutrino_pdgcode", &fm_neutrino_pdgcode, "neutrino_pdgcode/I");
  ftree->Branch("offset_x", &fm_offset_x, "offset_x/F");
  ftree->Branch("offset_y", &fm_offset_y, "offset_y/F");
  ftree->Branch("offset_z", &fm_offset_z, "offset_z/F");
  }

  //-----------------------------------------------------------------------
void T962G4Ana::analyze(const art::Event& evt) 
{
  
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);


  art::PtrVector<simb::MCTruth> mclist;
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    } 

  art::Ptr<simb::MCTruth> mc(mclist[0]);
// 

//comment out this section for through going muon simulation

if(mc->NeutrinoSet())
{
	simb::MCNeutrino neut(mc->GetNeutrino());
    fm_neutrino_x =neut.Nu().Vx();
    fm_neutrino_y =neut.Nu().Vy();
    fm_neutrino_z =neut.Nu().Vz();
    fm_neutrino_px =neut.Nu().Px();
    fm_neutrino_py =neut.Nu().Py();
    fm_neutrino_pz =neut.Nu().Pz();
    fm_neutrino_energy=neut.Nu().E();
    fm_neutrino_pdgcode=neut.Nu().PdgCode();
}
//     
        

    //get the list of particles from this event
 art::ServiceHandle<cheat::BackTracker> bt;
  const sim::ParticleList& plist = bt->ParticleList();
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
    fm_offset_x=117.3;//cm
    fm_offset_y=-19.3;//cm
    fm_offset_z=-152.95;//cm
    
    int numprimaries=0;
    int numinminos=0;
    int primindex=1;
    for(unsigned int i = 0; i < pvec.size(); ++i)
    {
     if(pvec[i]->Process()!="primary")
     continue; 
     numinminos++;
    simb::MCTrajectory trajectory = pvec[i]->Trajectory();
    int numberofpoints= pvec[i]->NumberTrajectoryPoints();    
    	for(int j=1; j<numberofpoints; j++)
    	{
    	TLorentzVector prevposition = trajectory.Position(j-1);
    	TLorentzVector position = trajectory.Position(j);
        	if(pvec[i]->Process()=="primary"&&prevposition.Z()<152.95&&position.Z()>=152.95&&!(abs(pvec[i]->PdgCode()==14)||abs(pvec[i]->PdgCode()==12)))
    	numprimaries++;    	
    	}    
    }
    ofstream myfile;
    myfile.open("MINOSin.txt", std::ios::out | std::ios::app);
    myfile << "________________________________________________________________________________\n";    
    myfile << "***** HEPEVT Common Event#: "<<event<<", "<<numprimaries+numinminos+1<<" particles (max 4000) ***** Double Precision\n";
    myfile << "4-byte integers, 8-byte floating point numbers, 4000-allocated entries.\n";
    myfile << "Indx Stat Par- chil-       (  P_x,       P_y,       P_z,    Energy,       M )\n";
    myfile << "      ID  ents dren    Prod (   X,         Y,         Z,        cT)      [mm]\n";
    myfile << "________________________________________________________________________________\n";
    myfile <<"   "<<1<<" +10001    "<<0<<"    0    (    "<<fm_neutrino_px<<",  "<<fm_neutrino_py<<",  "<<fm_neutrino_pz<<",  "<<fm_neutrino_energy<<",  "<<fm_mass<<")\n";
    if(fm_neutrino_pdgcode>0)
    myfile <<"      +"<<fm_neutrino_pdgcode<<"    0    0    (      "<<(fm_neutrino_x+fm_offset_x)*10.<<",       "<<(fm_neutrino_y+fm_offset_y)*10.<<",       "<<(fm_neutrino_z+fm_offset_z)*10.<<",  0.)\n";
    else
    myfile <<"      "<<fm_neutrino_pdgcode<<"    0    0    (      "<<(fm_neutrino_x+fm_offset_x)*10.<<",       "<<(fm_neutrino_y+fm_offset_y)*10.<<",       "<<(fm_neutrino_z+fm_offset_z)*10.<<",  0.)\n";
    for(unsigned int i = 0; i < pvec.size(); ++i){
    
    if(pvec[i]->Process()!="primary")
    continue;
    
    simb::MCTrajectory trajectory = pvec[i]->Trajectory();
    fm_pdgcode=pvec[i]->PdgCode();       
    fm_mass=pvec[i]->Mass();
    fm_trackID=i;
    TLorentzVector initposition = trajectory.Position(0);
    TLorentzVector initmomentum = trajectory.Momentum(0);
    int numberofpoints= pvec[i]->NumberTrajectoryPoints();   
    
       fm_init_x=initposition.X();
       fm_init_y=initposition.Y();
       fm_init_z=initposition.Z();
       fm_init_px=initmomentum.Px();
       fm_init_py=initmomentum.Py();
       fm_init_pz=initmomentum.Pz();
       fm_init_energy=initmomentum.E(); 
      
      primindex++;

      myfile <<"   "<<primindex<<" +10001    1    0    (    "<<fm_init_px<<",  "<<fm_init_py<<",  "<<fm_init_pz<<",  "<<fm_init_energy<<",  "<<fm_mass<<")\n";
      if(fm_pdgcode>0)
      myfile <<"      +"<<fm_pdgcode<<"    1    0    (      "<<(fm_init_x+fm_offset_x)*10.<<",       "<<(fm_init_y+fm_offset_y)*10.<<",       "<<(fm_init_z+fm_offset_z)*10.<<",  2.64e+06)\n";
      else
      myfile <<"      "<<fm_pdgcode<<"    1    0    (      "<<(fm_init_x+fm_offset_x)*10.<<",       "<<(fm_init_y+fm_offset_y)*10.<<",       "<<(fm_init_z+fm_offset_z)*10.<<",  2.64e+06)\n";
    for(int j=1; j<numberofpoints; j++)
    {
      TLorentzVector prevposition = trajectory.Position(j-1);
      TLorentzVector position = trajectory.Position(j);
      TLorentzVector prevmomentum = trajectory.Momentum(j-1);
      TLorentzVector momentum = trajectory.Momentum(j);

       if(       (((prevposition.X()<0.&&position.X()>=0.)||(prevposition.X()>0.&&position.X()<=0.))&&(prevposition.Y()<20.&&prevposition.Y()>=-20.&&prevposition.Z()>0.&&prevposition.Z()<=90.))
       ||       (((prevposition.X()<47.&&position.X()>=47.)||(prevposition.X()>47.&&position.X()<=47.))&&(prevposition.Y()<20.&&prevposition.Y()>=-20.&&prevposition.Z()>0.&&prevposition.Z()<=90.))
       ||       (((prevposition.Y()<20.&&position.Y()>=20.)||(prevposition.Y()>20.&&position.Y()<=20.))&&(prevposition.X()<47.&&prevposition.X()>=0.&&prevposition.Z()>0.&&prevposition.Z()<=90.))
       ||       (((prevposition.Y()<-20.&&position.Y()>=-20.)||(prevposition.Y()>-20.&&position.Y()<=-20.))&&(prevposition.X()<47.&&prevposition.X()>=0.&&prevposition.Z()>0.&&prevposition.Z()<=90.))
       ||       (((prevposition.Z()<0.&&position.Z()>=0.)||(prevposition.Z()>0.&&position.Z()<=0.))&&(prevposition.Y()<20.&&prevposition.Y()>=-20.&&prevposition.X()>0.&&prevposition.X()<=47.))
       ||       (((prevposition.Z()<90.&&position.Z()>=90.)||(prevposition.Z()>90.&&position.Z()<=90.))&&(prevposition.Y()<20.&&prevposition.Y()>=-20.&&prevposition.X()>0.&&prevposition.X()<=47.))
       )
       {    
       fm_tpcexit_x=position.X();
       fm_tpcexit_y=position.Y();
       fm_tpcexit_z=position.Z();
       fm_tpcexit_px=momentum.Px();
       fm_tpcexit_py=momentum.Py();
       fm_tpcexit_pz=momentum.Pz();
       fm_tpcexit_energy=momentum.E(); 
       }
    // if(prevposition.Z()<154.22&&position.Z()>154.22&&abs(pvec[i]->PdgCode()!=14)&&abs(pvec[i]->PdgCode()!=12))
//            {
       if(prevposition.Z()<152.95&&position.Z()>=152.95&&abs(pvec[i]->PdgCode()!=14)&&abs(pvec[i]->PdgCode()!=12))
        {
       fm_minosenter_x=position.X();
       fm_minosenter_y=position.Y();
       fm_minosenter_z=position.Z();
       fm_minosenter_px=momentum.Px();
       fm_minosenter_py=momentum.Py();
       fm_minosenter_pz=momentum.Pz();
       fm_minosenter_energy=momentum.E();

    primindex++;
    myfile <<"   "<<primindex<<"   +1    "<<primindex-1<<"    0    (    "<<fm_minosenter_px<<",  "<<fm_minosenter_py<<",  "<<fm_minosenter_pz<<",  "<<fm_minosenter_energy<<",  "<<fm_mass<<")\n";
    if(fm_pdgcode>0)
    myfile <<"       +"<<fm_pdgcode<<"    "<<primindex-1<<"    0    (      "<<(fm_minosenter_x+fm_offset_x)*10.<<",       "<<(fm_minosenter_y+fm_offset_y)*10.<<",       "<<(fm_minosenter_z+fm_offset_z)*10.<<",  2.64e+06)\n";
    else
     myfile <<"      "<<fm_pdgcode<<"    "<<primindex-1<<"    0    (      "<<(fm_minosenter_x+fm_offset_x)*10.<<",       "<<(fm_minosenter_y+fm_offset_y)*10.<<",       "<<(fm_minosenter_z+fm_offset_z)*10.<<",  2.64e+06)\n";
    
//          
        // std::cout<<fm_pdgcode<<" at z position "<<position.Z()<<" lost "<<((prevmomentum.E()-momentum.E())*1000.)<<"MeV in "<<(position.Z()-prevposition.Z())<<"cm which is "<<((prevmomentum.E()-momentum.E())*1000.)/(position.Z()-prevposition.Z())<< "MeV/cm VolumeName is " << geom->VolumeName(position.Vect())<<" "<<geom->MaterialName(position.Vect()) << std::endl;      

       fm_energylost=fm_init_energy-fm_minosenter_energy;
       fm_deflectionangle= TMath::ACos(((fm_tpcexit_px*fm_minosenter_px)+(fm_tpcexit_py*fm_minosenter_py)+(fm_tpcexit_pz*fm_minosenter_pz))/(sqrt(pow(fm_tpcexit_px,2)+pow(fm_tpcexit_py,2)+pow(fm_tpcexit_pz,2))*sqrt(pow(fm_minosenter_px,2)+pow(fm_minosenter_py,2)+pow(fm_minosenter_pz,2))));
//         std::cout<<"exit-enter "<<(fm_tpcexit_energy-fm_minosenter_energy)*1000.<<" "<<fm_deflectionangle<<" "<<fm_minosenter_energy<<std::endl;
            }
    }
    ftree->Fill();    
  }
  myfile << "________________________________________________________________________________\n";
  myfile.close(); 
  event++;
return;  
}

} // namespace T962G4

