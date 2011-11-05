////////////////////////////////////////////////////////////////////////
/// \file  FileMuons.cxx
/// \brief Generator for muons from a file.
///
/// Module designed to produce a set list of particles for a MC event
///
/// \version $Id: FileMuons.cxx,v 1.4 2010/03/29 09:54:01 echurch Exp $
/// \author  echurch@fnal.gov
////////////////////////////////////////////////////////////////////////
// C++ includes.
#include <string>
#include <cmath>
#include <memory>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// nutools includes
#include "SimulationBase/simbase.h"

// lar includes
#include "EventGenerator/FileMuons.h"
#include "Geometry/geo.h"
#include "SummaryData/summary.h"

#include "TVector3.h"
#include "TDatabasePDG.h"

namespace evgen{

  //____________________________________________________________________________
  FileMuons::FileMuons(fhicl::ParameterSet const& pset)
     : fEventNumberOffset(pset.get<      int                 >("EventNumberOffset"))
     , fSeed(314159)					                          
     , fPDG              (pset.get< std::vector<int>	     >("PDG")              )	       
     , fFileName         (pset.get< std::string     	     >("FileName")         )      
     , fMuonsFileType    (pset.get< std::string              >("MuonsFileType")    )
     , fTreeName         (pset.get< std::string   	     >("TreeName")         )      
     , fBranchNames      (pset.get< std::vector<std::string> >("BranchNames")      )
  {

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();

  }

  //____________________________________________________________________________
  void FileMuons::beginJob()
  {
    countFile=fEventNumberOffset;
    LOG_DEBUG("FileMuons: starting at event ") << countFile <<std::endl;

    if (fMuonsFileType.compare("source")==0) 
      {
	std::cout << "FileMuons: Not yet equipped to walk through muons with TFS mojo."<< std::endl; 
      }
    else if (fMuonsFileType.compare("root")==0) 
      {
	std::cout << "FileMuons: You have chosen to read muons from Root File " << fFileName << std::endl; 
      }
    else if (fMuonsFileType.compare("text")==0) 
      {
	std::cout << "FileMuons: You have chosen to read muons from " << fFileName << "." << std::endl; 
      }
    else 
      {
	std::cout << "FileMuons: You must specify one of source/text/root file to read for muons."<< std::endl; 

      }

    if (fMuonsFileType.compare("text")==0) 
      {
	fMuonFile = new ifstream(fFileName.c_str());
	long begin = fMuonFile->tellg();
	fMuonFile->seekg (0, std::ios::end);
	long end = fMuonFile->tellg();
	std::cout << "FileMuons: " << fFileName << " size is: " << (end-begin) << " bytes.\n";
	fMuonFile->seekg (0, std::ios::beg);

	for (unsigned int header=0; header<3 && fMuonFile->good(); ++header)
	  {
	    std::string line;
	    getline(*fMuonFile,line);
	  }
	if (!fMuonFile->good())
	  {
	    std::cout << "FileMuons: Problem reading muon file header."<< std::endl; 
	  }
      }  // fMuonsFileType is a text file.
    else if (fMuonsFileType.compare("root")==0) 
      {
	fMuonFileR = new TFile(fFileName.c_str(),"READ");
	TNtuple = (TTree*)(fMuonFileR->Get(fTreeName.c_str()));

	TNtuple->SetBranchAddress("x", &xtmp, &b_x);
	TNtuple->SetBranchAddress("y", &ytmp, &b_y);
	TNtuple->SetBranchAddress("z", &ztmp, &b_z);
	TNtuple->SetBranchAddress("E", &E, &b_E);
	TNtuple->SetBranchAddress("costheta", &costheta, &b_costheta);
	TNtuple->SetBranchAddress("phi", &phi, &b_phi);
	TNtuple->SetBranchAddress("xdet", &xdet, &b_xdet);
	TNtuple->SetBranchAddress("ydet", &ydet, &b_ydet);
	TNtuple->SetBranchAddress("zdet", &zdet, &b_zdet);
	TNtuple->SetBranchAddress("px", &pxtmp, &b_px);
	TNtuple->SetBranchAddress("py", &pytmp, &b_py);
	TNtuple->SetBranchAddress("pz", &pztmp, &b_pz);
	TNtuple->SetBranchAddress("charge", &charge, &b_charge);

      }  // fMuonsFileType is a root file.
    
  }

  //____________________________________________________________________________
  void FileMuons::endJob()
  {
    if (fMuonsFileType.compare("text")==0) 
      fMuonFile->close();
    if (fMuonsFileType.compare("root")==0) 
      fMuonFileR->Close();
  }

  //____________________________________________________________________________
  void FileMuons::beginRun(art::Run& run)
  {

    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;

    geo::DetId_t detid = geo->DetId();

    std::auto_ptr<sumdata::RunData> runcol(new sumdata::RunData(detid));

    run.put(runcol);

    return;
  }

  //____________________________________________________________________________
  FileMuons::~FileMuons()
  {
  }

  //____________________________________________________________________________
  void FileMuons::produce(art::Event& evt)
  {

    ///auto_ptr allows ownership to be transferred to the art::Event after the put statement
    std::auto_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    ReadEvents(truth);

//     std::cout << "put mctruth into the vector" << std::endl;
    truthcol->push_back(truth);

//     std::cout << "add vector to the event " << truthcol->size() << std::endl;
    evt.put(truthcol);

    return;
  }

  //____________________________________________________________________________
  void FileMuons::ReadEvents(simb::MCTruth &mct) 
  {

//     std::cout << "size of particle vector is " << fPDG.size() << std::endl;

    ///every event will have one of each particle species in the fPDG array
    for (unsigned int i=0; i<fPDG.size(); ++i) {
      
      // Choose momentum
      //double p = 0.0;
      double m(0.108);

      TVector3 x;
      TVector3 p;
      Double_t q;
      Int_t pdgLocal;
      
      if (fMuonsFileType.compare("text")==0)
	{

	  std::string line;
	  getline(*fMuonFile,line);
	  if (!fMuonFile->good())
	    {
	      std::cout << "FileMuons: Problem reading muon file line ...."<< countFile << ". Perhaps you've exhausted the events in " << fFileName << std::endl; exit(0);
	    }
	  else
	    {
	      //	  std::cout << "FileMuons: getline() gives "<< line << " for event " << countFile << std::endl; 
	    }
	  countFile++;
	  
	  LOG_DEBUG("FileMuons: countFile is ") << countFile <<std::endl;
	  char * cstr, *ptok;
	  
      // Split this line into tokens
	  cstr = new char [line.size()+1];
	  strcpy (cstr, line.c_str());
	  // cstr now contains a c-string copy of str
	  ptok=strtok (cstr,"*");
	  unsigned int fieldCount = 0;
	  unsigned int posIndex = 0;
	  unsigned int pIndex = 0;
	  while (ptok!=NULL)
	    {
	      // std::cout << ptok << std::endl;
	      ptok=strtok(NULL,"*");
	      if (fieldCount==9 || fieldCount==10 || fieldCount==11)
		{
		  p[pIndex] = atof(ptok); pIndex++;
		  //   std::cout << ptok << std::endl;
		}
	      if (fieldCount==6 || fieldCount==7 || fieldCount==8)
		{
		  x[posIndex] = atof(ptok); 
		  // make the z axis point up for x, as with p
		  if (posIndex==2) {x[posIndex] = -1.0*x[posIndex];}
		  posIndex++;
		}
	      if (fieldCount==12)
		{
		  q = atof(ptok); 
		}
	      fieldCount++;
	    }
	  
	  delete[] cstr;  

	}
      else if (fMuonsFileType.compare("root")==0) // from root file
	{
	    /*
	      // Don't use this yet. Keep the specific branch-by-branch identification.
	      for (unsigned int ii=0;ii<fBranchNames.size();ii++)
	      {
	       TNtuple->SetBranchAddress(fBranchNames[ii], x+ii); 
	      }
	    */
	  //	  TNtuple->ResetBranchAddresses();
	  TNtuple->GetEntry(countFile);
	  //TNtuple->Show(countFile);

	  x.SetXYZ(xdet,ydet,-zdet); // as with txt file, make z point up.
	  p.SetXYZ(pxtmp,pytmp,pztmp);
	  q = charge;

	  countFile++;

	} // End read.
      
      static TDatabasePDG  pdgt;
      pdgLocal = -q*fPDG[i];
      
      TParticlePDG* pdgp = pdgt.GetParticle(pdgLocal);
      if (pdgp) m = pdgp->Mass();


      //       std::cout << "set the position "<<std::endl;
      // This gives coordinates at the center of the 300mx300m plate that is 3m above top of 
      // cavern. Got these by histogramming deJong's xdet,ydet,zdet.
      x[0] -= 875.;
      x[1] -= -1395.;
      x[2] += 1359.36+3+3+8.5; // 3 for plate height. 8.5 more to get to the center of TPC.
      // Now, must rotate to TPC coordinates. Let's orient TPC axis along z axis,
      // Cosmics, mostly going along deJong's +z axis must be going along TPC -y axis.
      x.RotateX(-M_PI/2);
      p.RotateX(-M_PI/2);
      //assume the position is relative to the center of the TPC
      art::ServiceHandle<geo::Geometry> geom;
      TVector3 off3 = geom->GetTPCFrontFaceCenter() * 0.01;
      x += off3;

      TLorentzVector pos(x[0]*100.0, x[1]*100.0, x[2]*100.0, 0.0);
      TLorentzVector pvec(p[0]*1000.0,p[1]*1000.0,p[2]*1000.0,sqrt(p.Mag2()*1000.0*1000.0+m*m));
      std::cout << "x[m] and p [TeV] are " << std::endl;
      x.Print();
      p.Print();
    
      int trackid = -1*(i+1); // set track id to -i as these are all primary particles and have id <= 0
      std::string primary("primary");
      simb::MCParticle part(trackid, pdgLocal, primary);
      part.AddTrajectoryPoint(pos, pvec);
      
      //       std::cout << "add the particle to the primary" << std::endl;
      
      mct.Add(part);

    }//end loop over particles
    
    return;
  }

}//end namespace evgen
////////////////////////////////////////////////////////////////////////
