////////////////////////////////////////////////////////////////////////
/// \file  LArG4Ana.cxx
/// \brief Use Geant4 to run the LArSoft detector simulation
///
/// \version $Id: LArG4.cxx,v 1.22 2010/07/20 06:08:30 bjpjones Exp $
/// \author  brebel@fnal.gov
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
#include "Simulation/SimListUtils.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Geometry/geo.h"
#include "LArG4/LArG4Ana.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

// C++ Includes
#include <iostream>
#include <cstring>
#include <sys/stat.h>

namespace larg4 {

  //-----------------------------------------------------------------------
  // Constructor
  LArG4Ana::LArG4Ana(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //-----------------------------------------------------------------------
  // Destructor
  LArG4Ana::~LArG4Ana() 
  {
  }

  //-----------------------------------------------------------------------
  void LArG4Ana::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    
    fDistFromVertex = tfs->make<TH1D>("distFromVertex", ";#Deltar (cm);", 500, 0., 500.);
    fEnergyWithin5cm  = tfs->make<TH2D>("energyWithin5cm",  ";E (GeV); Fraction Deposited", 10, 0., 5., 50, 0., 1.);
    fEnergyWithin10cm = tfs->make<TH2D>("energyWithin10cm", ";E (GeV); Fraction Deposited", 10, 0., 5., 50, 0., 1.);
    fEnergyWithin15cm = tfs->make<TH2D>("energyWithin15cm", ";E (GeV); Fraction Deposited", 10, 0., 5., 50, 0., 1.);
    fEnergyWithin20cm = tfs->make<TH2D>("energyWithin20cm", ";E (GeV); Fraction Deposited", 10, 0., 5., 50, 0., 1.);
    fEnergyWithin25cm = tfs->make<TH2D>("energyWithin25cm", ";E (GeV); Fraction Deposited", 10, 0., 5., 50, 0., 1.);
    fEnergyWithin30cm = tfs->make<TH2D>("energyWithin30cm", ";E (GeV); Fraction Deposited", 10, 0., 5., 50, 0., 1.);
    fEnergyWithin35cm = tfs->make<TH2D>("energyWithin35cm", ";E (GeV); Fraction Deposited", 10, 0., 5., 50, 0., 1.);
    fEnergyWithin40cm = tfs->make<TH2D>("energyWithin40cm", ";E (GeV); Fraction Deposited", 10, 0., 5., 50, 0., 1.);
    fEnergyWithin45cm = tfs->make<TH2D>("energyWithin45cm", ";E (GeV); Fraction Deposited", 10, 0., 5., 50, 0., 1.);
    fEnergyWithin50cm = tfs->make<TH2D>("energyWithin50cm", ";E (GeV); Fraction Deposited", 10, 0., 5., 50, 0., 1.);

    fPDGCodes    = tfs->make<TH1D>("pdgcodes", ";PDG Code;",               5000, -2500, 2500);
    fPi0Momentum = tfs->make<TH1D>("pi0mom",   ";#pi^{0} Momentum (GeV);", 1000, 0.,    1000.);

    fTree = tfs->make<TTree>("MCTTree","MCTTree");
    fnEnergy = tfs->make<TH1D>("nEnergy", ";n,#Lambda^{0},K^{0} Momentum (GeV);", 100, 0., 10.);
    fnDist = tfs->make<TH1D>("nDistance", ";n,#Lambda^{0},K^{0} Distance (m);", 200, -30000.0, +30000.);


    fT4Origin = new Float_t[4];
    fT4DOrigin = new Float_t[fTNds*4];
    fT4Termination = new Float_t[4];
    fT4Momentum = new Float_t[4];
    fT4DMomentum = new Float_t[fTNds*4];
    fTDID = new Int_t[fTNds];
    fTDPdg = new Int_t[fTNds];
    fTDWt  = new Float_t[fTNds];
    fTNds4 = fTNds*4; //  TTree/Branch requirement to store this.

    fTree->Branch("MCEvt", &fTEvt, "MCEvt/I");
    fTree->Branch("MCSub", &fTSub, "MCSub/I");
    fTree->Branch("MCRun", &fTRun, "MCRun/I");
    fTree->Branch("MCWt", &fTWeight, "MCWt/F");
    fTree->Branch("MCPdg", &fTPdg, "MCPdg/I");
    fTree->Branch("MCID", &fTID, "MCID/I");
    fTree->Branch("MCParentID", &fTParentID, "MCParentID/I");
    fTree->Branch("MCNumDs", &fTNds, "MCNumDs/I");
    fTree->Branch("MCNumDs4", &fTNds4, "MCNumDs4/I");
    fTree->Branch("MCDID", fTDID, "MCDID[MCNumDs]/I");
    fTree->Branch("MCDPdg", fTDPdg, "MCDPdg[MCNumDs]/I");
    fTree->Branch("MCDWt", fTDWt, "MCDWt[MCNumDs]/I");
    fTree->Branch("MCProcess", fTProcess, "MCProcess/C");
    fTree->Branch("MCVolume", fTVolume, "MCVolume/C");
    fTree->Branch("MCTVolume", fTTVolume, "MCTVolume/C");
    fTree->Branch("MCMaterial", fTMaterial, "MCMaterial/C");
    fTree->Branch("MCDProcess", fTDProcess, "MCDProcess[MCNumDs]/C");
    fTree->Branch("MCStatus", &fTStatus, "MCStatus/I");
    fTree->Branch("MCOrigin", fT4Origin, "MCOrigin[4]/F");
    fTree->Branch("MCDOrigin", fT4DOrigin, "MCDOrigin[MCNumDs4]/F");
    fTree->Branch("MCTermination", fT4Termination, "MCTermination[4]/F");
    fTree->Branch("MCMomentum", fT4Momentum, "MCMomentum[4]/F");
    fTree->Branch("MCDMomentum", fT4DMomentum, "MCDMomentum[MCNumDs4]/F");
  
  }

  //-----------------------------------------------------------------------
  void LArG4Ana::reconfigure(fhicl::ParameterSet const& p)
  {
    fG4ModuleLabel    = p.get< std::string >("GeantModuleLabel");
    fTNdsOriginal     = p.get< int         >("Ndaughters"      );
    fTruthModuleLabel = p.get< std::string >("TruthModuleLabel"); 
    fTNds = fTNdsOriginal;
    
    return;
  }

  //-----------------------------------------------------------------------
  void LArG4Ana::analyze(const art::Event& evt) 
  {

    //get the list of particles from this event
    sim::ParticleList plist = sim::SimListUtils::GetParticleList(evt,fG4ModuleLabel);
    art::ServiceHandle<geo::Geometry> geom;

    // loop over all sim::SimChannels in the event and make sure there are no
    // sim::IDEs with trackID values that are not in the sim::ParticleList
    std::vector<const sim::SimChannel*> sccol;
    evt.getView(fG4ModuleLabel, sccol);

    std::vector<const simb::MCTruth*> mctcol;
    evt.getView(fTruthModuleLabel, mctcol);
    
    for(size_t sc = 0; sc < sccol.size(); ++sc){
      const std::map<unsigned short, std::vector<sim::IDE> >& tdcidemap = sccol[sc]->TDCIDEMap();
      std::map<unsigned short, std::vector<sim::IDE> >::const_iterator mapitr;
      for(mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
	const std::vector<sim::IDE> idevec = (*mapitr).second;
	for(size_t iv = 0; iv < idevec.size(); ++iv){
	  if(plist.find( idevec[iv].trackID ) == plist.end()
	     && idevec[iv].trackID != sim::NoParticleId) 
	  mf::LogWarning("LArG4Ana") << idevec[iv].trackID << " is not in particle list"; 
	}
      }
    }

    // check the containment of the energy for this event
    //get the list of voxels from this event
    sim::LArVoxelList vlist = sim::SimListUtils::GetLArVoxelList(evt, fG4ModuleLabel);
    this->CheckContainment(vlist, mctcol[0]);

    // loop over the particle list and fill the pdg code histogram
    /*
    for(sim::ParticleList::const_iterator itr = plist.begin(); itr != plist.end(); itr++){
      const sim::Particle *p = (*itr).second;
      fPDGCodes->Fill(p->PdgCode());

      std::cout << "mother = "     << p->Mother() 
		<< " pdgcode = "   << p->PdgCode()
		<< " daughters = " << p->NumberDaughters() << std::endl;
      std::cout << *p;
    }
    */

    // get the particles from the event handle
    art::Handle< std::vector<sim::Particle> > parHandle;
    evt.getByLabel(fG4ModuleLabel, parHandle);

    art::PtrVector<sim::Particle> pvec;
    for(unsigned int i = 0; i < parHandle->size(); ++i){
      art::Ptr<sim::Particle> p(parHandle, i);      
      pvec.push_back(p);
      fPDGCodes->Fill(p->PdgCode());
    }

    // now look for pi0's that decay to 2 gammas
    int pi0loc = -1;
    int numpi0gamma = 0;
    for(unsigned int i = 0; i < pvec.size(); ++i){
      if(pvec[i]->PdgCode() == 111) pi0loc = i;
      if(pvec[i]->Mother() == pi0loc+1 && pi0loc > 0 && pvec[i]->PdgCode() == 22)
	{
	  mf::LogInfo("LArG4Ana") << pvec[i]->E() << " gamma energy ";
	  ++numpi0gamma;
	}
      
      if (pvec[i]->PdgCode() == 2112 || pvec[i]->PdgCode() == 3122 || pvec[i]->PdgCode() == 130|| pvec[i]->PdgCode() == 310|| pvec[i]->PdgCode() == 311 ) // n,Lambda,K0s,K0L,K0
	{
	  fnEnergy->Fill(pvec[i]->E(),pvec[i]->Weight());
	  fnDist->Fill(pvec[i]->Vx(),pvec[i]->Weight());
	}

      fTPdg = pvec[i]->PdgCode();
      fTID = pvec[i]->TrackId();
      // 0 out strings, else there may be cruft in here from prev evt.
      for (unsigned int s=0;s<35;s++) 
	{
	  *(fTProcess+s)=0; *(fTProcess+s)=0;
	  *(fTMaterial+s)=0; *(fTMaterial+s)=0;
	  *(fTVolume+s)=0; *(fTVolume+s)=0;
	  *(fTTVolume+s)=0; *(fTTVolume+s)=0;
	}

      for (unsigned int s=0;s<pvec[i]->Process().length();s++) 
	*(fTProcess+s) = pvec[i]->Process()[s];
      TVector3 dum = pvec[i]->Position().Vect();
      //std::cout << " dum is " << std::endl;
      //dum.Print();
      //std::cout << " VolumeName is " << geom->VolumeName(pvec[i]->Position().Vect()) << std::endl;
      for (unsigned int s=0;s<geom->MaterialName(pvec[i]->Position().Vect()).length();s++) 
	*(fTMaterial+s) = geom->MaterialName(pvec[i]->Position().Vect())[s];
      for (unsigned int s=0;s<geom->VolumeName(pvec[i]->Position().Vect()).length();s++) 
	*(fTVolume+s) = geom->VolumeName(pvec[i]->Position().Vect())[s];
      for (unsigned int s=0;s<geom->VolumeName(pvec[i]->EndPoint()).length();s++) 
	*(fTTVolume+s) = geom->VolumeName(pvec[i]->EndPoint())[s];

      fTEvt = evt.id().event(); 
      fTSub = evt.subRun();
      fTRun = evt.run();
      fTParentID = pvec[i]->Mother();
      fTStatus = pvec[i]->StatusCode();
      int daughter = 9999;
      fTNds = TMath::Min(pvec[i]->NumberDaughters(),fTNdsOriginal);
      for ( int d = 0; d < fTNds; d++ )
        {
	  daughter = pvec[i]->Daughter(d);
	  fTDID[d] = daughter; 
	  // zero it out.
	  for (unsigned int s=0;s<35;s++) *(fTDProcess[d]+s)=0; 

	  for(unsigned int jj = i; jj < pvec.size(); jj++) // Don't look below i.
	    {
	      if (fTDID[d]==pvec[jj]->TrackId())
		{
		  fTDPdg[d] = pvec[jj]->PdgCode(); // get the pointer,  
		  fTDWt[d]  = pvec[jj]->Weight();
		  for (unsigned int s=0;s<pvec[jj]->Process().length();s++) 
		    *(fTDProcess[d]+s) = pvec[jj]->Process()[s];
		  for (unsigned int kk=0;kk<4;kk++)
		    {
		      fT4DOrigin[d*4+kk] = pvec[jj]->Position()[kk];
		      fT4DMomentum[d*4+kk] = pvec[jj]->Momentum()[kk];
		    }
		  break;
		}	      
	    }
        }
      
      for (unsigned int ii=0;ii<4;ii++)
	{
	  fT4Termination[ii] = 1e9;
	  fT4Origin[ii] = pvec[i]->Position()[ii];
	  if (ii!=3) fT4Termination[ii] = pvec[i]->EndPoint()[ii];
	  if (ii==4) fT4Termination[ii] = pvec[i]->Momentum()[ii]; // yes, odd
	  fT4Momentum[ii] = pvec[i]->Momentum()[ii];
	  //	  std::cout << "LArG4Ana: Origin, Momentum, ii " << fT4Origin[ii]  << ", "<< fT4Momentum[ii] << ", " << ii << "." << std::endl;
	}

      fTWeight = pvec[i]->Weight();
      fTree->Fill();
      
      //std::cout << "LArG4Ana:: pid,Vx,Position(1),process is " << pvec[i]->PdgCode() << " " << pvec[i]->Vx() <<  " " << pvec[i]->Position()[0] << " " << pvec[i]->Process() << std::endl;

    } // end loop on particles in list 
    if(numpi0gamma == 2 && pi0loc>0){
      mf::LogInfo("LArG4Ana") << pvec[pi0loc]->E();
      fPi0Momentum->Fill(pvec[pi0loc]->E());
    }

    return;
  }


  //-----------------------------------------------------------------------
  void LArG4Ana::CheckContainment(sim::LArVoxelList const& vlist,
				  const simb::MCTruth*     mct) 
  {
    // loop over the voxels to figure out how much energy is deposited
    // and at what distance from the vertex.
    sim::LArVoxelList::const_iterator vitr;

    // make a map of track id to histograms for the fraction of energy contained
    // within some radius of the central axis
    double idToE10cm;
    double idToE20cm;
    double idToE30cm;
    double idToE40cm;
    double idToE50cm;

    // get the xy vertex of the interaction
    double vtxX = mct->GetNeutrino().Nu().Vx();
    double vtxY = mct->GetNeutrino().Nu().Vy();
    double enu  = mct->GetNeutrino().Nu().E();

    if(abs(mct->GetNeutrino().Nu().PdgCode()) != 12 &&
       mct->GetNeutrino().Mode() != simb::kQE       &&
       mct->GetNeutrino().CCNC() != simb::kCC) return;

    for(vitr = vlist.begin(); vitr != vlist.end(); vitr++){
      // figure out if this voxel is more than n cm away from 
      // the vertex in either the x or y directions
      
      if(vitr->first.X() - vtxX < 10. &&
	 fabs(vitr->first.Y() - vtxY) < 1000.) 
	idToE10cm += vitr->second.Energy();
      if(vitr->first.X() - vtxX < 20. &&
	 fabs(vitr->first.Y() - vtxY) < 2000.) 
	idToE20cm += vitr->second.Energy();
      if(vitr->first.X() - vtxX < 30. &&
	 fabs(vitr->first.Y() - vtxY) < 3000.) 
	idToE30cm += vitr->second.Energy();
      if(vitr->first.X() - vtxX < 40. &&
	 fabs(vitr->first.Y() - vtxY) < 4000.) 
	idToE40cm += vitr->second.Energy();
      if(vitr->first.X() - vtxX < 50. &&
	 fabs(vitr->first.Y() - vtxY) < 5000.) 
	idToE50cm += vitr->second.Energy();

    }// end loop over voxels	  

    fEnergyWithin10cm->Fill(enu, idToE10cm/enu);
    fEnergyWithin20cm->Fill(enu, idToE20cm/enu);
    fEnergyWithin30cm->Fill(enu, idToE30cm/enu);
    fEnergyWithin40cm->Fill(enu, idToE40cm/enu);
    fEnergyWithin50cm->Fill(enu, idToE50cm/enu);

    return;
  }

} // namespace larg4
