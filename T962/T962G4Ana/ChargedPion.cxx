////////////////////////////////////////////////////////////////////////
/// \file  ChargedPion.cxx
/// \brief Extract PDG codes, energy, etc from mc events with pi+
///
/// \author Ellen.Klein@Yale.edu
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <math.h>
#include <iomanip>
#include "TH2.h"
#include "TSystem.h"


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
#include "Geometry/geo.h"
// LArSoft Includes
#include "SimulationBase/MCFlux.h"
#include "SimulationBase/MCNeutrino.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "TDatabasePDG.h"
#include "T962/T962G4Ana/ChargedPion.h"
#include "MCCheater/BackTracker.h"

// ROOT includes
#include "TH1.h"

// C++ Includes
#include <iostream>
#include <cstring>
#include <sys/stat.h>

namespace cpion {

  //-----------------------------------------------------------------------
  // Constructor
  ChargedPion::ChargedPion(fhicl::ParameterSet const& pset) :
    fG4ModuleLabel(pset.get< std::string >("GeantModuleLabel")),
    fGenieModuleLabel(pset.get< std::string >("GenieModuleLabel"))
  {
    std::cout<<"constructor"<<std::endl;
  }

  //-----------------------------------------------------------------------
  // Destructor
  ChargedPion::~ChargedPion() 
  {
  }
  
  //-----------------------------------------------------------------------
  void ChargedPion::beginJob()
  {
    //Start art::ServiceHandle to save all trees
    art::ServiceHandle<art::TFileService> tfs;

    //Build your trees
    nTree=tfs->make<TTree>("MCNuTree","Nu Tree");
    gTree=tfs->make<TTree>("MCGenieTree", "Genie Tree");
    tTree=tfs->make<TTree>("MCGeantTree", "Geant Tree");
    // pTDPdg = new Int_t[pTNds];
    //pTNds4=pTNds*4;

    //Branch your tree
    nTree->Branch("MCEvt", &Evt, "MCEvt/I");
    nTree->Branch("MCNuPdg", &nPdg, "MCNuPdg/I");
    nTree->Branch("MCNuInteractionType", &nType, "MCNuInteractionType/I");
    nTree->Branch("MCNuVx", &nVx, "MCNuVx/D");
    nTree->Branch("MCNuVy", &nVy, "MCNuVy/D");
    nTree->Branch("MCNuVz", &nVz, "MCNuVz/D");
    nTree->Branch("MCNuEnergy", &nEnergy, "MCNuEnergy/D");
    nTree->Branch("MCNuQSquared", &nQSquared, "MCNuQSquared/D");

    gTree->Branch("MCEvt", &Evt, "MCEvt/I");
    gTree->Branch("MCGPdg", &gPdg, "MCGPdg/I");
    gTree->Branch("MCGTrackId", &gTrackId, "MCGTrackId/I");
    gTree->Branch("MCGMother", &gMother, "MCGMother/I");
    gTree->Branch("MCGStatus", &gStatus, "MCGStatus/I");
    gTree->Branch("MCGEnergy", &gEnergy, "MCGEnergy/D");
    gTree->Branch("MCGVx", &gVx, "MCGVx/D");
    gTree->Branch("MCGVy", &gVy, "MCGVy/D");
    gTree->Branch("MCGVz", &gVz, "MCGVz/D");
    //CCCCCC
    gTree->Branch("MCNuPdg", &nPdg, "MCNuPdg/I");
    gTree->Branch("MCNuPdg", &nPdg, "MCNuPdg/I");
    gTree->Branch("MCNuInteractionType", &nType, "MCNuInteractionType/I");
    gTree->Branch("MCNuVx", &nVx, "MCNuVx/D");
    gTree->Branch("MCNuVy", &nVy, "MCNuVy/D");
    gTree->Branch("MCNuVz", &nVz, "MCNuVz/D");
    gTree->Branch("MCNuEnergy", &nEnergy, "MCNuEnergy/D");
    gTree->Branch("MCNuQSquared", &nQSquared, "MCNuQSquared/D");

    tTree->Branch("MCEvt", &Evt, "MCEvt/I");
    tTree->Branch("MCTPdg", &tPdg, "MCTPdg/I");
    tTree->Branch("MCTTrackId", &tTrackId, "MCTTrackId/I");
    tTree->Branch("MCTMother", &tMother, "MCTMother/I");
    tTree->Branch("MCTNumDaughters", &tDaughters, "MCTNumDaughters/I");
    tTree->Branch("MCTDecay", &tDecay, "MCTDecay/I");
    tTree->Branch("MCTInelastic", &tInelastic, "MCTInelastic/I");
    tTree->Branch("MCTCapture", &tCapture, "MCTCapture/I");
    tTree->Branch("MCTEnergy", &tEnergy, "MCTEnergy/D");
    tTree->Branch("MCTVx", &tVx, "MCTVx/D");
    tTree->Branch("MCTVy", &tVy, "MCTVy/D");
    tTree->Branch("MCTVz", &tVz, "MCTVz/D");
    tTree->Branch("MCTEndx", &tEndx, "MCTEndx/D");
    tTree->Branch("MCTEndy", &tEndy, "MCTEndy/D");
    tTree->Branch("MCTEndz", &tEndz, "MCTEndz/D");
    tTree->Branch("MCTPx", &tPx, "MCTPx/D");
    tTree->Branch("MCTPy", &tPy, "MCTPy/D");
    tTree->Branch("MCTPz", &tPz, "MCTPz/D");
    tTree->Branch("MCTP", &tP, "MCTP/D");
    tTree->Branch("MCTPt", &tPt, "MCTPt/D");
    tTree->Branch("MCTDaughPdg", &tDaughPdg, "MCTDaughPdg[200]/D");
    tTree->Branch("MCTProtons", &tProtons, "MCTProtons/I");
    tTree->Branch("MCTNeutrons", &tNeutrons, "MCTNeutrons/I");
    tTree->Branch("MCTProtonE", &tProtonE, "MCTProtonE[10]/D");
    tTree->Branch("MCTNeutronE", &tNeutronE, "MCTNeutronE[10]/D");
    tTree->Branch("MCTProtonDist", &tProtonDist, "MCTProtonDist[10]/D");
    tTree->Branch("MCTMuPdg", &tMuPdg, "MCTMuPdg/I");
    tTree->Branch("MCTMuEnergy", &tMuEnergy, "MCTMuEnergy/D");
    tTree->Branch("MCTMuVx", &tMuVx, "MCTMuVx/D");
    tTree->Branch("MCTMuVy", &tMuVy, "MCTMuVy/D");
    tTree->Branch("MCTMuVz", &tMuVz, "MCTMuVz/D");
    tTree->Branch("MCTMuEndx", &tMuEndx, "MCTMuEndx/D");
    tTree->Branch("MCTMuEndy", &tMuEndy, "MCTMuEndy/D");
    tTree->Branch("MCTMuEndz", &tMuEndz, "MCTMuEndz/D");
    tTree->Branch("MCTMuPx", &tMuPx, "MCTMuPx/D");
    tTree->Branch("MCTMuPy", &tMuPy, "MCTMuPy/D");
    tTree->Branch("MCTMuPz", &tMuPz, "MCTMuPz/D");
    tTree->Branch("MCTTotalEnergy", &tTotalEnergy, "MCTTotalEnergy/D");
    tTree->Branch("MCTChargedEnergy", &tChargedEnergy, "MCTChargedEnergy/D");
    tTree->Branch("MCTTotalPPar", &tTotalPPar, "MCTTotalPPar/D");
    tTree->Branch("MCTTotalPPerp", &tTotalPPerp, "MCTTotalPPerp/D");
    tTree->Branch("MCTChargedPPar", &tChargedPPar, "MCTChargedPPar/D");
    tTree->Branch("MCTChargedPPerp", &tChargedPPerp, "MCTChargedPPerp/D");
    tTree->Branch("MCTNumPions", &tNumPions,"MCTNumPions/I");
    //CCCCCC
    tTree->Branch("MCNuPdg", &nPdg, "MCNuPdg/I");
    tTree->Branch("MCNuPdg", &nPdg, "MCNuPdg/I");
    tTree->Branch("MCNuInteractionType", &nType, "MCNuInteractionType/I");
    tTree->Branch("MCNuVx", &nVx, "MCNuVx/D");
    tTree->Branch("MCNuVy", &nVy, "MCNuVy/D");
    tTree->Branch("MCNuVz", &nVz, "MCNuVz/D");
    tTree->Branch("MCNuEnergy", &nEnergy, "MCNuEnergy/D");
    tTree->Branch("MCNuQSquared", &nQSquared, "MCNuQSquared/D");

    }

  //-----------------------------------------------------------------------
  void ChargedPion::analyze(const art::Event& evt) 
  {
    //TREES ARE ONLY FILLED FOR EVENTS WITH PI+

    //At the start of each event, zero all important variables
    PiLock=-999;
    nPdg=-999;
    nType=-999;
    nVx=-999;
    nVy=-999;
    nVz=-999;
    nEnergy=-999;
    nQSquared=-999;

    gPdg=-999;
    gTrackId=-999;
    gMother=-999;
    gStatus=-999;
    gEnergy=-999;
    gVx=-999;
    gVy=-999;
    gVz=-999;

    tPdg=-999;
    tTrackId=-999;
    tMother=-999;
    tDaughters=-999;
    tDecay=-999;
    tInelastic=-999;
    tEnergy=-999;
    tVx=-999;
    tVy=-999;
    tVz=-999;
    tEndx=-999;
    tEndy=-999;
    tEndz=-999; 
    tPx=-999;
    tPy=-999;
    tPz=-999;
    tP=-999;
    tPt=-999;  
    tDaughPdg[0]=-999;
    tProtons=-999;
    tNeutrons=-999;
    tMuPdg=-999;
    tMuEnergy=-999;
    tMuVx=-999;
    tMuVy=-999;
    tMuVz=-999;
    tMuEndx=-999;
    tMuEndy=-999;
    tMuEndz=-999;
    tMuPx=-999;
    tMuPy=-999;
    tMuPz=-999;
    tNumPions=0;

    tTotalEnergy=-999;
    tChargedEnergy=-999;
    tTotalPPar=-999;
    tTotalPPerp=-999;
    for (int pq=0; pq<10; pq++)
      {
	tProtonE[pq]=-999;
	tNeutronE[pq]=-999;
	tProtonDist[pq]=-999;
      }

    //Print out event information and set up handles
    std::cout<<"event    : "<<evt.id().event()<<std::endl;
    Evt=evt.id().event();
   
    //To get Genie information
    art::Handle< std::vector<simb::MCTruth> > mclist;
    evt.getByLabel(fGenieModuleLabel,mclist);
    const std::vector<simb::MCTruth>& vect(*mclist);
    art::Ptr<simb::MCTruth> n(mclist,0);
    if(vect.size()!=1)
      {
	std::cout<<"TROUBLE!! EVENT "<<Evt<<" HAS MORE THAN ONE NEUTRINO"<<std::endl;
      }
   
    //To get Geant information
    art::ServiceHandle<cheat::BackTracker> bt;
    std::vector<const simb::MCParticle*> pvec;

    for(unsigned int l=0; l < bt->ParticleList().size(); ++l)
      {
	pvec.push_back(bt->ParticleList().Particle(l));
      }

    if(mclist->size()==0) return;
   

    //Get && store neutrino information

    {
      nPdg=n->GetNeutrino().Nu().PdgCode();
      nType=n->GetNeutrino().InteractionType();
      nVx=n->GetNeutrino().Nu().Vx();
      nVy=n->GetNeutrino().Nu().Vy();
      nVz=n->GetNeutrino().Nu().Vz();
      nEnergy=n->GetNeutrino().Nu().E();
      nQSquared=n->GetNeutrino().QSqr();

      nTree->Fill();
    }
    //std::cout<<"NEUTRINO!!!!!!"<<n->GetNeutrino().Nu().PdgCode()<<" Vx "<<n->GetNeutrino().Nu().Vx()<<" VY "<<n->GetNeutrino().Nu().Vy()<<" VZ "<<n->GetNeutrino().Nu().Vz()<<std::endl;
     
     
    //Look at Genie information
    for(int r=0; r < n->NParticles(); r++)
      {	
      	simb::MCParticle part(n->GetParticle(r));

	//std::cout<<"PARTICLE: "<<part.PdgCode()<<" STATUS "<<part.StatusCode()<<" TRACK ID "<<part.TrackId()<<" MOTHER ID "<<(part.TrackId()*-1)-1<<" MOTHER "<<part.Mother()<<" ENERGY "<<part.E()<<" VX "<<part.Vx()<<" VY "<<part.Vy()<< " VZ "<<part.Vz()<<std::endl;
	gPdg=part.PdgCode();
	gTrackId=part.TrackId();
	gMother=part.Mother();
	gStatus=part.StatusCode();
	gEnergy=part.E();
	gVx=part.Vx();
	gVy=part.Vy();
	gVz=part.Vz();

	gTree->Fill();

	gPdg=-999;
	gTrackId=-999;
	gMother=-999;
	gStatus=-999;
	gEnergy=-999;
	gVx=-999;
	gVy=-999;
	gVz=-999;
      }
   
    //Now look at Geant info
    //Double_t length=0;
    tProtons=0;
    tNeutrons=0;
    tTotalEnergy=0;
    tChargedEnergy=0;
    tTotalPPar=0;
    tTotalPPerp=0;
    tnumPrimary=0;
    tChargedPPar=0;
    tChargedPPerp=0;

    //loop over all events and look at number of primary protons, neutrons, pions. Also sum the energy && momenta of all primary particles and all charged primary particles.
    for(unsigned int j=0;j<pvec.size();j++)
      {
	if(pvec[j]->Process()!="primary")
	  {
	    // break;
	  }
	if(pvec[j]->Process()=="primary"&&pvec[j]->PdgCode()<100000000)
	  {
	 //    length=(pow ( pow((pvec[j]->Vx()-pvec[j]->EndPosition()[0]),2) + pow((pvec[j]->Vy()-pvec[j]->EndPosition()[1]),2) + pow((pvec[j]->Vz()-pvec[j]->EndPosition()[2]),2), 0.5));

	    if(pvec[j]->PdgCode()==2212)
	      {
		tProtons++;
	      }
	    if(pvec[j]->PdgCode()==2112)
	      {
		tNeutrons++;
	      }
	    if(pvec[j]->PdgCode()==211)
	      {
		tNumPions++;
	      }

	    if(pvec[j]->PdgCode()!=2112 && pvec[j]->PdgCode()!=111 && pvec[j]->PdgCode()!=130 && pvec[j]->PdgCode()!=310 && pvec[j]->PdgCode()!=3212 && pvec[j]->PdgCode()!=3322 && pvec[j]->PdgCode()!=-3322 && pvec[j]->PdgCode()!=-421)
	      {
		//could also add a track length constraint using the length variable above...
		tChargedEnergy=tChargedEnergy+pvec[j]->E();
		tChargedPPar=tChargedPPar+pvec[j]->Pz();
		tChargedPPerp=tChargedPPerp+pow( (pow(pvec[j]->Px(),2)+pow(pvec[j]->Py(),2)) ,0.5);
		tnumPrimary++;
	//	std::cout<<" tChargedEnergy "<<tChargedEnergy<<" event "<<Evt<<" pdg "<<pvec[j]->PdgCode()<<std::endl;
	      }

	    if(pvec[j]->PdgCode()==2112 || pvec[j]->PdgCode()==111 || pvec[j]->PdgCode()==130 || pvec[j]->PdgCode()==310 || pvec[j]->PdgCode()==3212 || pvec[j]->PdgCode()==3322 || pvec[j]->PdgCode()==-3322 || pvec[j]->PdgCode()==-421)
	      {
		if(nType==1097)
		  {
		    std::cout<<"BAD! Neutral Coherent Particle! "<<pvec[j]->PdgCode()<<std::endl;
		  }
	      }

	    tTotalEnergy=tTotalEnergy+pvec[j]->E();
	    tTotalPPar=tTotalPPar+pvec[j]->Pz();
	    tTotalPPerp=tTotalPPerp+pow( pow(pvec[j]->Px(),2) + pow(pvec[j]->Py(),2) ,0.5);
	    if(nType==1097)
	      {
		//std::cout<<"COHERENT EVENT "<<Evt<<" PARTICLE "<<tnumPrimary<<" ENERGY "<<pvec[j]->E()<<" TOTAL ENERGY "<<tTotalEnergy<<" PDG "<<pvec[j]->PdgCode()<<std::endl;
	      }
	  }//if primary particle && NOT an excited nucleus
      }//loop over all particles in event

    //now we loop over all particles in an event again but look for pi+'s
    //probably a better way to do this than looping over all particles twice, if you think of it let me know!
    for(unsigned int i=0; i<pvec.size(); i++)
      {

	for(int kk=0;kk<10;kk++)
	  {
	    tProtonE[kk]=-999;
	    tProtonDist[kk]=-999;
	  }
	for(int jj=0;jj<200;jj++)
	  {
	    tDaughPdg[jj]=-999;
	  }

	//std::cout<<"PARTICLE: "<<pvec[i]->PdgCode()<<" EVENT "<<Evt<<" TRACK ID "<<pvec[i]->TrackId()<<" MOTHER "<<pvec[i]->Mother()<<" PROCESS "<<pvec[i]->Process()<<" ENERGY "<<pvec[i]->E()<<" PX "<<pvec[i]->Vx()<<" PY "<<pvec[i]->Vy()<< " PZ "<<pvec[i]->Vz()<<" NDAUGHTERS "<<pvec[i]->NumberDaughters()<<" POSITION "<<pvec[i]->Position()[0]<<" ENDPOINT "<<pvec[i]->EndPosition()[0]<<" "<<pvec[i]->EndPosition()[1]<<" "<<pvec[i]->EndPosition()[2]<<std::endl;

	if(pvec[i]->Process()=="primary" && (pvec[i]->PdgCode()==211 || pvec[i]->PdgCode()==13 || pvec[i]->PdgCode()==-13))
	  {
	    //want information on muons
	    if(pvec[i]->PdgCode()==13 && pvec[i]->Process()=="primary")
	      {
		tMuPdg=pvec[i]->PdgCode();
		tMuEnergy=pvec[i]->E();
		tMuVx=pvec[i]->Vx();
		tMuVy=pvec[i]->Vy();
		tMuVz=pvec[i]->Vz();
		tMuEndx=pvec[i]->EndPosition()[0];
		tMuEndy=pvec[i]->EndPosition()[1];
		tMuEndz=pvec[i]->EndPosition()[2];
		tMuPx=pvec[i]->Px();
		tMuPy=pvec[i]->Py();
		tMuPz=pvec[i]->Pz();
	      }
	    //or anti-muons
	    if(pvec[i]->PdgCode()==-13 && pvec[i]->Process()=="primary")
	      {
		tMuPdg=pvec[i]->PdgCode();
		tMuEnergy=pvec[i]->E();
		tMuVx=pvec[i]->Vx();
		tMuVy=pvec[i]->Vy();
		tMuVz=pvec[i]->Vz();
		tMuEndx=pvec[i]->EndPosition()[0];
		tMuEndy=pvec[i]->EndPosition()[1];
		tMuEndz=pvec[i]->EndPosition()[2];
		tMuPx=pvec[i]->Px();
		tMuPy=pvec[i]->Py();
		tMuPz=pvec[i]->Pz();
	      }
	    //collect everything you could possibly want to know about primary pi+'s
	    if(pvec[i]->PdgCode()==211&&pvec[i]->Process()=="primary")
	      {
		  {
		    PiLock=1;//to fill tree
		numPiPlus++;
		tPdg=pvec[i]->PdgCode();
		tTrackId=pvec[i]->TrackId();
		tMother=pvec[i]->Mother();
		tDaughters=pvec[i]->NumberDaughters();
		tEnergy=pvec[i]->E();
		tVx=pvec[i]->Vx();
		tVy=pvec[i]->Vy();
		tVz=pvec[i]->Vz();
		tEndx=pvec[i]->EndPosition()[0];
		tEndy=pvec[i]->EndPosition()[1];
		tEndz=pvec[i]->EndPosition()[2];
		tPx=pvec[i]->Px();
		tPy=pvec[i]->Py();
		tPz=pvec[i]->Pz();
		tP=pvec[i]->P();
		tPt=pvec[i]->Pt();
		
		tInelastic=0;
		tDecay=0;
		tCapture=0;
		//Now see if the pion has any daughters and how they were created.
		for (int s=0; s<pvec[i]->NumberDaughters(); s++)
	      {
		int daughTrack=pvec[i]->Daughter(s);
		//double check to make sure mother is pion...
		//note there is an off by one error (corrected for in variable offsets) between the track pulled via Daughter and which track I open...
		int momNum=pvec[daughTrack-1]->Mother();
		if(pvec[momNum-1]->PdgCode()!=211)
		  {
		    std::cout<<"MOTHER ERROR!!!!!! ERROR!!!! ERROR!!!!!"<<std::endl;
		    break;
		  }

		else if(pvec[momNum-1]->PdgCode()==211)
		  {
		  
		    //Check for a "New" Process
		    if(pvec[daughTrack-1]->Process()!="PionPlusInelastic"&&pvec[daughTrack-1]->Process()!="Decay"&&pvec[daughTrack-1]->Process()!="hIoni"&&pvec[daughTrack-1]->Process()!="hadElastic"&&pvec[daughTrack-1]->Process()!="phot"&&pvec[daughTrack-1]->Process()!="compt"&&pvec[daughTrack-1]->Process()!="eBrem"&&pvec[daughTrack-1]->Process()!="eIoni"&&pvec[daughTrack-1]->Process()!="annihil"&&pvec[daughTrack-1]->Process()!="hPairProd"&&pvec[daughTrack-1]->Process()!="hBrem")
		      {
			std::cout<<"KEY WORD!!!!!!!!!!! pdg "<<pvec[daughTrack-1]->PdgCode()<<" process "<<pvec[daughTrack-1]->Process()<<std::endl;
		      }


		    if(pvec[daughTrack-1]->Process()=="PionPlusInelastic")
		      {
			tDaughPdg[s]=pvec[daughTrack-1]->PdgCode();
			//std::cout<<"DAUGHTER PDG " <<tDaughPdg[s]<<std::endl;
			tInelastic++;
		      }
		    else if(pvec[daughTrack-1]->Process()=="Decay")
		      {
			tDecay++;
		      }
		  }
	      }
		  }
	      }
	if(PiLock==1)
	  {
		    tTree->Fill();
		    //std::cout<<"Filling: "<<Evt<<" tPdg "<<tPdg<<" tEnergy "<<tEnergy<<" tMuPdg "<<tMuPdg<<" tMuEnergy "<<tMuEnergy<<std::endl;
		    //reset everything after filling.
		    tPdg=-999;
		    tTrackId=-999;
		    tMother=-999;
		    tDaughters=-999;
		    tDecay=-999;
		    tInelastic=-999;
		    tCapture =-999;
		    tEnergy=-999;
		    tVx=-999;
		    tVy=-999;
		    tVz=-999;
		    tEndx=-999;
		    tEndy=-999;
		    tEndz=-999;    
		    tPx=-999;
		    tPy=-999;
		    tPz=-999;
		    tP=-999;
		    tPt=-999;  
		    PiLock=0;
	  }//if pilock==1
	  }//if primary && pi+

	//may put code to look for pi- here, most likely will just modify above to code to accept pdg==-211 & to look for daughter particles created via pi- capture
      
      }///loop over all particles (second time)
 
    if(tCapture==0 && tInelastic==0 && tDecay==0)
      {
	//std::cout<<" HM???? Event "<<Evt<<" Endx "<<tEndx<<" Endy "<<tEndy<<" Endz "<<tEndz<<std::endl;
      }
    std::cout<<"END OF EVENT!"<<std::endl;
    std::cout<<"NUMBER OF PIPLUS: "<<numPiPlus<<std::endl; 
    std::cout<<"NUMBER OF PIMINUS: "<<numPiMinus<<std::endl;
      
  }
}
// namespace ChargedPion


