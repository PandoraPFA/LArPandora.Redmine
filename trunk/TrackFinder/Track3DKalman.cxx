////////////////////////////////////////////////////////////////////////
//
// \file Track3DKalman.cxx
//
// \author echurch@fnal.gov
//
//  This algorithm is designed to reconstruct 3D tracks through  
//  GENFIT's Kalman filter.
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

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
#include "TrackFinder/Track3DKalman.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"
#include "Utilities/LArProperties.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"


// ROOT includes
#include "TVectorD.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"

// GENFIT includes
#include "Genfit/GFException.h"
#include "Genfit/GFAbsTrackRep.h"
#include "Genfit/GeaneTrackRep2.h"
#include "Genfit/RKTrackRep.h"
#include "Genfit/GFConstField.h"
#include "Genfit/GFFieldManager.h"
#include "Genfit/PointHit.h"
#include "Genfit/GFTrack.h"
#include "Genfit/GFKalman.h"
#include "Genfit/GFDaf.h"

static bool sp_sort_3dz(const recob::SpacePoint& h1, const recob::SpacePoint& h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[2] < xyz2[2];
}

//-------------------------------------------------
trkf::Track3DKalman::Track3DKalman(fhicl::ParameterSet const& pset) 
{

    this->reconfigure(pset);

    // get the random number seed, use a random default if not specified    
    // in the configuration file.  
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());

    createEngine(seed);

    produces< std::vector<recob::Track> >();
}

void trkf::Track3DKalman::reconfigure(fhicl::ParameterSet const& pset) 
  {

    fSpacePtsModuleLabel   = pset.get< std::string >("SpacePtsModuleLabel");
    fGenieGenModuleLabel   = pset.get< std::string >("GenieGenModuleLabel");
    fG4ModuleLabel         = pset.get< std::string >("G4ModuleLabel");

    fPosErr                = pset.get< std::vector < double >  >("PosErr3");   // resolution. cm
    fMomErr                = pset.get< std::vector < double >  >("MomErr3");   // GeV
    fMomStart              = pset.get< std::vector < double >  >("MomStart3"); // Will be unit norm'd.
    fGenfPRINT             = pset.get< bool >("GenfPRINT");

  }

//-------------------------------------------------
trkf::Track3DKalman::~Track3DKalman()
{

  /*
    delete stMCT ;
    delete covMCT;
    delete stREC;
    delete covREC;
    
    delete fpMCT;
    delete fpREC;
    delete fpRECL;
    delete fpRECt3D;
  */
}

//-------------------------------------------------
void trkf::Track3DKalman::beginJob()
{


  art::ServiceHandle<art::TFileService> tfs;
  
  stMCT  = new TMatrixT<Double_t>(5,1);
  covMCT = new TMatrixT<Double_t>(5,5);
  stREC  = new TMatrixT<Double_t>(5,1);
  covREC = new TMatrixT<Double_t>(5,5);
  
  fpMCT = new Float_t[4];
  fpREC = new Float_t[4];
  fpRECL = new Float_t[4];
  fpRECt3D = new Float_t[4];
  fDimSize = 400; // if necessary will get this from pset in constructor.

  fshx = new Float_t[fDimSize];
  fshy = new Float_t[fDimSize];
  fshz = new Float_t[fDimSize];

//   //TFile fileGENFIT("GENFITout.root","RECREATE");

  
  tree = tfs->make<TTree>("GENFITttree","GENFITttree");
  //tree->Branch("stMCT",&stMCT,"stMCT[5]/F"); // "TMatrixT<Double_t>"

  tree->Branch("stMCT","TMatrixD",&stMCT,64000,0);
  //tree->Branch("covMCT",&covMCT,"covMCT[25]/F");
  tree->Branch("covMCT","TMatrixD",&covMCT,64000,0);
  //tree->Branch("stREC",&stREC,"stREC[5]/F");
  tree->Branch("stREC","TMatrixD",&stREC,64000,0);
  //tree->Branch("covREC",&covREC,"covREC[25]/F");
  tree->Branch("covREC","TMatrixD",&covREC,64000,0);
  
  
  tree->Branch("chi2",&chi2,"chi2/F");
  tree->Branch("nfail",&nfail,"nfail/I");
  tree->Branch("ndf",&ndf,"ndf/I");
  tree->Branch("evtNo",&evtt,"evtNo/I");
  tree->Branch("chi2ndf",&chi2ndf,"chi2ndf/F");

  tree->Branch("trkNo",&nTrks,"trkNo/I");
  tree->Branch("ptsNo",&fptsNo,"ptsNo/I");
  tree->Branch("shx",fshx,"shx[ptsNo]/F");
  tree->Branch("shy",fshy,"shy[ptsNo]/F");
  tree->Branch("shz",fshz,"shz[ptsNo]/F");

  tree->Branch("pMCT",fpMCT,"pMCT[4]/F");
  tree->Branch("pRECKalF",fpREC,"pRECKalF[4]/F");
  tree->Branch("pRECKalL",fpRECL,"pRECKalL[4]/F");
  tree->Branch("pRECt3D",fpRECt3D,"pRECt3D[4]/F");
  

  //TGeoManager* geomGENFIT = new TGeoManager("Geometry", "Geane geometry");
  //TGeoManager::Import("config/genfitGeom.root");
  //  gROOT->Macro("config/Geane.C"); 
 
}

void trkf::Track3DKalman::endJob()
{
  if (!rep) delete rep;
  if (!repMC) delete repMC;
}


//------------------------------------------------------------------------------------//
void trkf::Track3DKalman::produce(art::Event& evt)
{ 

  rep=0;
  repMC=0;

  // get services
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  // get the random number generator service and make some CLHEP generators
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &engine = rng->getEngine();
  CLHEP::RandGaussQ gauss(engine);

  //////////////////////////////////////////////////////
  // Make a std::auto_ptr<> for the thing you want to put into the event
  // because that handles the memory management for you
  //////////////////////////////////////////////////////
  std::auto_ptr<std::vector<recob::Track> > tcol(new std::vector<recob::Track>);

  // define TPC parameters
  TString tpcName = geom->GetLArTPCVolumeName();


  // get input Hit object(s).
  art::Handle< std::vector<recob::Track> > trackListHandle;
  evt.getByLabel(fSpacePtsModuleLabel,trackListHandle);

  art::PtrVector<simb::MCTruth> mclist;
  std::vector<const sim::SimChannel*> simChannelHandle;

  if (!evt.isRealData())
    {

      //      std::cout << "Track3DKalman: This is MC." << std::endl;
      // std::cout<<"Run "<<evt.run()<<" Event "<<evt.id().event()<<std::endl;

      art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
      evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);

      for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
	{
	  art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
	  mclist.push_back(mctparticle);
	}
      evt.getView(fG4ModuleLabel, simChannelHandle);      

    }

  //create collection of spacepoints that will be used when creating the Track object
   std::vector<recob::SpacePoint> spacepoints;
  art::PtrVector<recob::Track> trackIn;
  // std::cout<<"Run "<<evt.run()<<" Event "<<evt.id().event()<<std::endl;
  mf::LogInfo("Track3DKalman: ") << "There are " <<  trackListHandle->size() << " Track3Dreco/SpacePt tracks/groups (whichever) in this event.";

  for(unsigned int ii = 0; ii < trackListHandle->size(); ++ii)
    {
      art::Ptr<recob::Track> track(trackListHandle, ii);
      trackIn.push_back(track);
      
    }
      
      TVector3 MCOrigin;
      TVector3 MCMomentum;
      // TVector3 posErr(.05,.05,.05); // resolution. 0.5mm
      // TVector3 momErr(.1,.1,0.2);   // GeV
      TVector3 posErr(fPosErr[0],fPosErr[1],fPosErr[2]); // resolution. 0.5mm
      TVector3 momErr(fMomErr[0],fMomErr[1],fMomErr[2]);   // GeV

      // This is strictly for MC
      if (!evt.isRealData())
	{
	  // Below breaks are stupid, I realize. But rather than keep all the MC
	  // particles I just take the first primary, e.g., muon and only keep its
	  // info in the branches of the Ttree. I could generalize later, ...
	  for( unsigned int ii = 0; ii < mclist.size(); ++ii )
	    {
	    //art::Ptr<const simb::MCTruth> mc(mctruthListHandle,i);
	    art::Ptr<simb::MCTruth> mc(mclist[ii]);
	    for(int jj = 0; jj < mc->NParticles(); ++jj)
	      {
		simb::MCParticle part(mc->GetParticle(jj));
		mf::LogInfo("Track3DKalman: ") << "FROM MC TRUTH, the particle's pdg code is: "<<part.PdgCode()<< " with energy = "<<part.E() <<", with energy = "<<part.E()<< " and vtx and momentum in Global (not volTPC) coords are " ;
		MCOrigin.SetXYZ(part.Vx(),part.Vy(),part.Vz()); // V for Vertex
		MCMomentum.SetXYZ(part.Px(),part.Py(),part.Pz());
		MCOrigin.Print();
		MCMomentum.Print();
		repMC = new genf::RKTrackRep(MCOrigin,
					     MCMomentum,
					     posErr,
					     momErr,
					     part.PdgCode());  
		break;
	      }
	    break;
	  }
	  //for saving of MC truth    
	  stMCT->ResizeTo(repMC->getState());
	  *stMCT = repMC->getState();
	  covMCT-> ResizeTo(repMC->getCov());
	  *covMCT = repMC->getCov();
	  mf::LogInfo("Track3DKalman: ") <<" repMC, covMC are ... " ;
	  repMC->getState().Print();
	  repMC->getCov().Print();

	} // !isRealData
      
      art::PtrVector<recob::Track>::const_iterator trackIter = trackIn.begin();
      
      nTrks=0;
    while(trackIter!=trackIn.end()) 
      {
	spacepoints.clear();
	spacepoints= (*trackIter)->SpacePoints();
	
	nTrks++;

	mf::LogInfo("Track3DKalman: ") << "found "<< spacepoints.size() <<" 3D spacepoint(s).";

// Add the 3D track to the vector of the reconstructed tracks
	if(spacepoints.size()>0)
	  {
      // Insert the GENFIT/Kalman stuff here then make the tracks. Units are cm, GeV.
	    const double resolution = 0.5; // dunno, 5 mm
	    const int numIT = 9; // 3->1, EC, 6-Jan-2011. Back, 7-Jan-2011.


	    //TVector3 mom(0.0,0.0,2.0);
	    TVector3 mom(fMomStart[0],fMomStart[1],fMomStart[2]);
	    //mom.SetMag(1.);
	    TVector3 momM(mom);
	    momM.SetX(gauss.fire(momM.X(),momErr.X()/* *momM.X() */));
	    momM.SetY(gauss.fire(momM.Y(),momErr.Y()/* *momM.Y() */));
	    momM.SetZ(gauss.fire(momM.Z(),momErr.Z()/* *momM.Z() */));
	    //std::cout << "Track3DKalman: sort spacepoints by x (volTPC coords) for tracker's sake." << std::endl;


	    std::sort(spacepoints.begin(), spacepoints.end(), sp_sort_3dz);
	    
	    //std::sort(spacepoints.begin(), spacepoints.end(), sp_sort_3dx); // Reverse sort!
	    
	    genf::GFFieldManager::getInstance()->init(new genf::GFConstField(0.,0.,0.0));
	    genf::GFDetPlane planeG((TVector3)(spacepoints[0].XYZ()),momM);
	    

      //      std::cout<<"Track3DKalman about to do GAbsTrackRep."<<std::endl;
      // Initialize with 1st spacepoint location and a guess at the momentum.
	    rep = new genf::RKTrackRep(//posM-.5/momM.Mag()*momM,
				       (TVector3)(spacepoints[0].XYZ()),
				       momM,
				       posErr,
				       momErr,
				       13);  // mu- hypothesis
      //      std::cout<<"Track3DKalman: about to do GFTrack. repDim is " << rep->getDim() <<std::endl;


	    genf::GFTrack fitTrack(rep);//initialized with smeared rep
      // Gonna sort in x cuz I want to essentially transform here to volTPC coords.
      // volTPC coords, cuz that's what the Geant3/Geane stepper wants, as that's its understanding
      // from the Geant4 geometry, which it'll use. EC, 7-Jan-2011.
	    int ihit = 0;
	    fptsNo = 0;
	    for (unsigned int point=0;point<spacepoints.size();++point)
	      {
		
		TVector3 spt3 = (TVector3)(spacepoints[point].XYZ());
		if (point%20) // Jump out of loop except on every 20th pt.
		  {
		    //continue;
		    // Icarus paper suggests we may wanna decimate our data in order to give
		    // trackfitter a better idea of multiple-scattering. EC, 7-Jan-2011.
		    //if (fabs(spt3[0]-spacepoints.at(point-1).XYZ()[0]) < 2) continue;
		  }
		if (fptsNo<fDimSize)
		  {
		    fshx[fptsNo] = spt3[0];
		    fshy[fptsNo] = spt3[1];
		    fshz[fptsNo] = spt3[2];
		  }
		fptsNo++;


		LOG_DEBUG("Track3DKalman: ") << "ihit xyz..." << spt3[0]<<","<< spt3[1]<<","<< spt3[2];
		fitTrack.addHit(new genf::PointHit(spt3,resolution),
				1,//dummy detector id
				ihit++
				);
	      }

	    //      std::cout<<"Track3DKalman about to do GFKalman."<<std::endl;
	    genf::GFKalman k;
      //k.setBlowUpFactor(500); // Instead of 500 out of box. EC, 6-Jan-2011.
      //k.setInitialDirection(+1); // Instead of 1 out of box. EC, 6-Jan-2011.
	    k.setNumIterations(numIT);
      //      std::cout<<"Track3DKalman back from setNumIterations."<<std::endl;
	    try{
	      //	std::cout<<"Track3DKalman about to processTrack."<<std::endl;
	      k.processTrack(&fitTrack);
	      //std::cout<<"Track3DKalman back from processTrack."<<std::endl;
	    }
	    catch(GFException& e){
	      mf::LogError("Track3DKalman: ") << "just caught a GFException."<<std::endl;
	      e.what();
	      mf::LogError("Track3DKalman: ") << "Exceptions won't be further handled ->exit(1) "<<__LINE__;

	      //	exit(1);
	    }

	    if(rep->getStatusFlag()==0) // 0 is successful completion
	      {
		LOG_DEBUG("Track3DKalman: ") << __FILE__ << " " << __LINE__ ;
		LOG_DEBUG("Track3DKalman: ") << "Track3DKalman.cxx: Original plane:";
	  
		if(fGenfPRINT) planeG.Print();
		LOG_DEBUG("Track3DKalman: ") << "Current (fit) reference Plane:";
		if(fGenfPRINT) rep->getReferencePlane().Print();
		
		LOG_DEBUG("Track3DKalman: ") << "Track3DKalman.cxx: Last reference Plane:";
		if(fGenfPRINT) rep->getLastPlane().Print();
		
		if(fGenfPRINT) 
		  {
		    if(planeG!=rep->getReferencePlane()) 
		      LOG_DEBUG("Track3DKalman: ")	<<"Track3DKalman: Original hit plane (not surprisingly) not current reference Plane!"<<std::endl;
		  }

		stREC->ResizeTo(rep->getState());
		*stREC = rep->getState();
		covREC->ResizeTo(rep->getCov());
		*covREC = rep->getCov();
		if(fGenfPRINT)
		  {
		    LOG_DEBUG("Track3DKalman: ") << " Final State and Cov:";
		    stREC->Print();
		    covREC->Print();
		  }
		chi2 = rep->getChiSqu();
		ndf = rep->getNDF();
		nfail = fitTrack.getFailedHits();
		chi2ndf = chi2/ndf;
		double dircoss[3],dircose[3];
		(*trackIter)->Direction(dircoss,dircose);      
		
		for (int ii=0;ii<3;++ii)
		  {
		    fpMCT[ii] = MCMomentum[ii]/MCMomentum.Mag();
		    fpREC[ii] = rep->getReferencePlane().getNormal()[ii];
		    fpRECL[ii] = rep->getLastPlane().getNormal()[ii];
		    fpRECt3D[ii] = dircoss[ii];
		  }
		fpMCT[3] = MCMomentum.Mag();
		fpREC[3] = -1.0/(*stREC)[0][0];
      
		evtt = (unsigned int) evt.id().event();
		mf::LogInfo("Track3DKalman: ") << "Track3DKalman about to do tree->Fill(). Chi2/ndf is " << chi2/ndf << ". All in volTPC coords .... pMCT[0-3] is " << fpMCT[0] << ", " << fpMCT[1] << ", " << fpMCT[2] << ", " << fpMCT[3] << ". pREC[0-3] is " << fpREC[0] << ", "<< fpREC[1] << ", " << fpREC[2] << ", " << fpREC[3] << ".";
	  
		tree->Fill();

	// Get the clusters associated to each track in induction and collection view
     		art::PtrVector<recob::Cluster> Icluster;
		art::PtrVector<recob::Cluster> Ccluster;
		
		art::ServiceHandle<geo::Geometry> geo;

		Icluster = (*trackIter)->Clusters(geo::kU); // induction
		Ccluster = (*trackIter)->Clusters(geo::kV); // collection
		
		art::PtrVector<recob::Cluster> clusters;
		clusters.clear();
		
		art::PtrVector<recob::Cluster>::const_iterator IclusterIter = Icluster.begin();
		while(IclusterIter!= Icluster.end() ){
		  clusters.push_back(*IclusterIter);
		  IclusterIter++;
		}
		
		art::PtrVector<recob::Cluster>::const_iterator CclusterIter = Ccluster.begin();
		while(CclusterIter!= Ccluster.end() ){
		  clusters.push_back(*CclusterIter);
		  CclusterIter++;
		}

	// Use new Track constructor... EC, 21-Apr-2011.	
		recob::Track  the3DTrack(clusters,spacepoints);
		double dircosF[3];
		double dircosL[3];
		for (int ii=0;ii<3;++ii)
		  {
		    dircosF[ii] = fpREC[ii];
		    dircosL[ii] = fpRECL[ii];
		  }
		the3DTrack.SetDirection(dircosF,dircosL);
		the3DTrack.SetID(tcol->size());
		tcol->push_back(the3DTrack);
	      } // getStatusFlag 
	  } // spacepoints.size()>0
	
  //
  //std::cout<<"Track3DKalman found "<< tcol->size() <<" 3D track(s)"<<std::endl;
	if(trackIter!=trackIn.end()) trackIter++;
	
      } // end loop over Track3Dreco/SpacePt tracks/groups (whichever) we brought into this event.
   
    evt.put(tcol);

   
}
