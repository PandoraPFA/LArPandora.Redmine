////////////////////////////////////////////////////////////////////////
/// \file  GENIEHelper.h
/// \brief Wrapper for generating neutrino interactions with GENIE
///
/// \version $Id: GENIEHelper.cxx,v 1.13 2010/04/30 03:14:50 brebel Exp $
/// \author  brebel@fnal.gov
/// \update 2010/3/4 Sarah Budd added simple_flux
////////////////////////////////////////////////////////////////////////
#include <math.h>

//ROOT includes
#include "TH1.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCollection.h"
#include "TSystem.h"
#include "TString.h"
#include "TRandom3.h"

//GENIE includes
#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GMCJDriver.h"
#include "FluxDrivers/GCylindTH1Flux.h"
#include "FluxDrivers/GMonoEnergeticFlux.h"
#include "FluxDrivers/GNuMIFlux.h"
#ifndef  MISSING_GSIMPLENTPFLUX
#include "FluxDrivers/GSimpleNtpFlux.h"
#endif
#include "Geo/ROOTGeomAnalyzer.h"
#include "Interaction/InitialState.h"
#include "Interaction/Interaction.h"
#include "Interaction/Kinematics.h"
#include "Interaction/KPhaseSpace.h"
#include "Interaction/ProcessInfo.h"
#include "Interaction/XclsTag.h"
#include "GHEP/GHepParticle.h"

//NuTools includes
#include "EventGeneratorBase/inc/GENIEHelper.h"
#include "SimulationBase/inc/MCTruth.h"
#include "SimulationBase/inc/MCNeutrino.h"
#include "SimulationBase/inc/MCParticle.h"
#include "SimulationBase/inc/MCFlux.h"

//experiment includes - assumes every experiment has a Geometry package
#include "Geometry/inc/Geometry.h"

// Framework includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


namespace evgb{

  static const int kNue      = 0;
  static const int kNueBar   = 1;
  static const int kNuMu     = 2;
  static const int kNuMuBar  = 3;
  static const int kNuTau    = 4;
  static const int kNuTauBar = 5;

  //--------------------------------------------------
  GENIEHelper::GENIEHelper(edm::ParameterSet const& pset) :
    fFluxType          (pset.getParameter< std::string              >("FluxType")        ),
    fFluxFile          (pset.getParameter< std::string              >("FluxFile")        ),
    fBeamName          (pset.getParameter< std::string              >("BeamName")        ),
    fTopVolume         (pset.getParameter< std::string              >("TopVolume")       ),
    fTargetA           (pset.getParameter< double                   >("TargetA")         ),
    fEventsPerSpill    (pset.getParameter< double                   >("EventsPerSpill")  ),
    fPOTPerSpill       (pset.getParameter< double                   >("POTPerSpill")     ),
    fHistEventsPerSpill(0.),
    fMonoEnergy        (pset.getParameter< double                   >("MonoEnergy")      ),
    fPOTUsed           (0.),
    fBeamRadius        (pset.getParameter< double                   >("BeamRadius")      ),
    fSurroundingMass   (pset.getParameter< double                   >("SurroundingMass") ),
    fGlobalTimeOffset  (pset.getParameter< double                   >("GlobalTimeOffset")),
    fRandomTimeOffset  (pset.getParameter< double                   >("RandomTimeOffset")),
    fZCutOff           (pset.getParameter< double                   >("ZCutOff")         ),
    fEnvironment       (pset.getParameter< std::vector<std::string> >("Environment")     ),
    fDetLocation       (pset.getParameter< std::string              >("DetectorLocation"))
  {
    int ranseed(pset.getParameter< int >("RandomSeed"));

    std::vector<double> beamCenter   (pset.getParameter< std::vector<double> >("BeamCenter")   );
    std::vector<double> beamDirection(pset.getParameter< std::vector<double> >("BeamDirection"));
    fBeamCenter.SetXYZ(beamCenter[0], beamCenter[1], beamCenter[2]);
    fBeamDirection.SetXYZ(beamDirection[0], beamDirection[1], beamDirection[2]);

    std::vector<int>genFlavors(pset.getParameter< std::vector<int> >("GenFlavors"));

    for(unsigned int i = 0; i < genFlavors.size(); ++i) fGenFlavors.insert(genFlavors[i]);

    TString junk = "";
    TRandom3 *rand = new TRandom3();
    rand->SetSeed();
    gRandom = rand;
    if(ranseed == 0)
      junk += gRandom->Integer(100000000);
    else
      junk += ranseed;

    std::string seed(junk);
    fEnvironment.push_back("GSEED");
    fEnvironment.push_back(seed);
    ///set the environment, the vector should come in pairs of variable name, then value
    for(unsigned int i = 0; i < fEnvironment.size(); i += 2){
      gSystem->Setenv(fEnvironment[i].c_str(), fEnvironment[i+1].c_str());
      std::cout << "setting GENIE environment " << fEnvironment[i]
		<< " to " << fEnvironment[i+1] << std::endl; 
    }

    std::cout << "using " << fFluxType << " fluxes" << std::endl;

    ///make the histograms
    if(fFluxType.compare("histogram") == 0){
      std::cout << "setting beam direction and center at "
		<< fBeamDirection.X() << " " << fBeamDirection.Y() << " " << fBeamDirection.Z()
		<< " (" << fBeamCenter.X() << "," << fBeamCenter.Y() << "," << fBeamCenter.Z()
		<< ")" << std::endl;

      TDirectory *savedir = gDirectory;
    
      fFluxHistograms.clear();
      TFile tf(fFluxFile.c_str());
      tf.ls();

      for(std::set<int>::iterator flvitr = fGenFlavors.begin(); flvitr != fGenFlavors.end(); flvitr++){
	if(*flvitr ==  12) fFluxHistograms.push_back(dynamic_cast<TH1D *>(tf.Get("nue")));
	if(*flvitr == -12) fFluxHistograms.push_back(dynamic_cast<TH1D *>(tf.Get("nuebar")));
	if(*flvitr ==  14) fFluxHistograms.push_back(dynamic_cast<TH1D *>(tf.Get("numu")));
	if(*flvitr == -14) fFluxHistograms.push_back(dynamic_cast<TH1D *>(tf.Get("numubar")));
	if(*flvitr ==  16) fFluxHistograms.push_back(dynamic_cast<TH1D *>(tf.Get("nutau")));
	if(*flvitr == -16) fFluxHistograms.push_back(dynamic_cast<TH1D *>(tf.Get("nutaubar")));
      }

      for(unsigned int i = 0; i < fFluxHistograms.size(); ++i){
	fFluxHistograms[i]->SetDirectory(savedir);
	fTotalHistFlux += fFluxHistograms[i]->Integral();
      }

      std::cout << "total histogram flux over desired flavors = " << fTotalHistFlux << std::endl;

    }//end if getting fluxes from histograms

    std::cout << "Generating events for " << fFluxType.c_str() << " from "
	      << fFluxFile.c_str() << std::endl;  

    std::cout << "Generating the following flavors " << std::endl;
  
    for(std::set<int>::iterator itr = fGenFlavors.begin(); itr != fGenFlavors.end(); itr++)
      std::cout << "        " << *itr << std::endl;

    if(fFluxType.compare("mono")==0){
      fEventsPerSpill = 1;
      std::cout << "Generating monoenergetic (" << fMonoEnergy << " GeV) neutrinos " 
		<< "with " << fEventsPerSpill << " events per spill" << std::endl;
    }

    if(fEventsPerSpill != 0)
      std::cout << "Generating " << fEventsPerSpill << " events for each spill" << std::endl;
    else
      std::cout << "Using " << fPOTPerSpill << " pot for each spill" << std::endl;

    return;
  }

  //--------------------------------------------------
  GENIEHelper::~GENIEHelper()
  {
  }

  //--------------------------------------------------
  double GENIEHelper::TotalHistFlux() 
  {
    if(   fFluxType.compare("ntuple")       == 0
       || fFluxType.compare("mono")         == 0 
       || fFluxType.compare("simple_flux" ) == 0 ) return -999.;

    return fTotalHistFlux;
  }

  //--------------------------------------------------
  void GENIEHelper::Initialize()
  {

    ///initialize the Geometry and Flux drivers
    InitializeGeometry();
    InitializeFluxDriver();
    fDriver = new genie::GMCJDriver();
    fDriver->UseFluxDriver(fFluxD);
    fDriver->UseGeomAnalyzer(fGeomD);

    ///turn on following line to speed up driver initialization
    //driver->UseMaxPathLengths(***supply some xml file name***);

    fDriver->Configure();
    fDriver->UseSplines();
    fDriver->ForceSingleProbScale();

    if(fFluxType.compare("histogram") == 0){
      ///fluxes are assumed to be given in units of neutrinos/cm^2/1e20POT/energy in histogram bin width
      
      ///get the bin width to remove energy dependence
      ///bins are assumed to be in GeV while fluxes are in /MeV
      double binWidth = fFluxHistograms[0]->GetBinWidth(1)*1000.; 

      ///determine product of pot/spill, mass, and cross section
      ///events = flux * pot * 10^-38 cm^2 (xsec) * (mass detector (in kg) / carbon mass (in kg)) * energy bin size
      fXSecMassPOT  = 1.e-38*1.e-20*binWidth;
      fXSecMassPOT *= fPOTPerSpill*(fDetectorMass+fSurroundingMass)/(12*1.67262158e-27); 

      std::cout << "Number of events per spill will be base on poisson mean of "
		<< fXSecMassPOT*fTotalHistFlux << std::endl;

      fHistEventsPerSpill = gRandom->Poisson(fXSecMassPOT*fTotalHistFlux);
    }

    ///set the pot/event counters to zero
    fSpillTotal = 0.;
    fPOTUsed   = 0.;

    return;
  }

  //--------------------------------------------------
  void GENIEHelper::InitializeGeometry()
  {
    edm::Service<geo::Geometry> geo;
    genie::geometry::ROOTGeomAnalyzer *rgeom = new genie::geometry::ROOTGeomAnalyzer(geo->ROOTGeoManager());

    ///the detector geometry uses cgs units.
    rgeom->SetLengthUnits(genie::units::centimeter);
    rgeom->SetDensityUnits(genie::units::gram_centimeter3);
    rgeom->SetTopVolName(fTopVolume.c_str());
    rgeom->SetMixtureWeightsSum(1.);

    /// casting to the GENIE geometry driver interface
    fGeomD        = dynamic_cast<genie::GeomAnalyzerI *>(rgeom);
    fDetLength    = geo->DetLength()*0.01;  
    fDetectorMass = geo->TotalMass();

    return;
  }

  //--------------------------------------------------
  void GENIEHelper::InitializeFluxDriver()
  {

    if(fFluxType.compare("ntuple") == 0){

      genie::flux::GNuMIFlux* numiFlux = new genie::flux::GNuMIFlux();
      numiFlux->LoadBeamSimData(fFluxFile, fDetLocation);
    
      ///set the number of cycles to run
      ///+++++++++this is stupid - to really set it i have to get a 
      ///value from the MCJDriver and i am not even sure what i have 
      ///below is approximately correct.
      ///for now just run on a set number of events that is kept track of 
      ///in the sample method
      /// numiFlux->SetNumOfCycles(int(fPOT/fFluxNormalization));
    
      fFluxD =  dynamic_cast<genie::GFluxI *>(numiFlux);
    } //end if using ntuple flux files
    else if(fFluxType.compare("simple_flux")==0){

#ifdef MISSING_GSIMPLENTPFLUX
      std::cerr << "Not built with GSimpleNtpFlux enabled" << std::endl;
      assert(0);
#else
      genie::flux::GSimpleNtpFlux* simpleFlux = 
        new genie::flux::GSimpleNtpFlux();
      simpleFlux->LoadBeamSimData(fFluxFile, fDetLocation);

      fFluxD =  dynamic_cast<genie::GFluxI *>(simpleFlux);
#endif    
    
    } //end if using simple_flux flux files
    else if(fFluxType.compare("histogram") == 0){

      genie::flux::GCylindTH1Flux* histFlux = new genie::flux::GCylindTH1Flux();
    
      std::cout << "beam r = " << fBeamRadius << " centered at ("
		<< fBeamCenter.X() << "," << fBeamCenter.Y() << "," << fBeamCenter.Z() << ")" << std::endl; 

      ///now add the different fluxes - fluxes were added to the vector in the same 
      ///order that the flavors appear in fGenFlavors
      int ctr = 0;
      for(std::set<int>::iterator i = fGenFlavors.begin(); i != fGenFlavors.end(); i++){
	histFlux->AddEnergySpectrum(*i, fFluxHistograms[ctr]);
	++ctr;
      } //end loop to add flux histograms to driver

      histFlux->SetNuDirection(fBeamDirection);
      histFlux->SetBeamSpot(fBeamCenter);
      histFlux->SetTransverseRadius(fBeamRadius);
    
      fFluxD = dynamic_cast<genie::GFluxI *>(histFlux);
    } //end if using a histogram
    else if(fFluxType.compare("mono") == 0){

      ///weight each species equally in the generation
      double weight = 1./(1.*fGenFlavors.size());
      //make a map of pdg to weight codes
      std::map<int, double> pdgwmap;
      for(std::set<int>::iterator i = fGenFlavors.begin(); i != fGenFlavors.end(); i++)
	pdgwmap[*i] = weight;

      genie::flux::GMonoEnergeticFlux *monoflux = new genie::flux::GMonoEnergeticFlux(fMonoEnergy, pdgwmap);
      monoflux->SetDirectionCos(fBeamDirection.X(), fBeamDirection.Y(), fBeamDirection.Z());
      monoflux->SetRayOrigin(fBeamCenter.X(), fBeamCenter.Y(), fBeamCenter.Z());
      fFluxD = dynamic_cast<genie::GFluxI *>(monoflux);
    } //end if using monoenergetic beam

    return;
  }

  //--------------------------------------------------
  bool GENIEHelper::Stop()
  {
    //   std::cout << "in GENIEHelper::Stop(), fEventsPerSpill = " << fEventsPerSpill
    // 	    << " fPOTPerSpill = " << fPOTPerSpill << " fSpillTotal = " << fSpillTotal 
    // 	    << " fHistEventsPerSpill = " << fHistEventsPerSpill << std::endl;

    ///determine if we should keep throwing neutrinos for 
    ///this spill or move on
    if(fEventsPerSpill > 0){
      if(fSpillTotal < fEventsPerSpill) 
	return false;
    }
    else{
      if( ( fFluxType.compare("ntuple")      == 0 || 
            fFluxType.compare("simple_flux") == 0    ) && 
          fSpillTotal < fPOTPerSpill) return false;
      else if(fFluxType.compare("histogram") == 0      && 
              fSpillTotal < fHistEventsPerSpill) return false;
    }

    ///made it to here, means need to reset the counters
    fPOTUsed += fSpillTotal;
    fSpillTotal = 0.;
    fHistEventsPerSpill = gRandom->Poisson(fXSecMassPOT*fTotalHistFlux);
    return true;
  }

  //--------------------------------------------------
  bool GENIEHelper::Sample(simb::MCTruth &truth,
			   simb::MCFlux  &flux)
  {
    genie::EventRecord* record = fDriver->GenerateEvent();

    if ( !record ) return false;

    TLorentzVector *vertex = record->Vertex();
    if(vertex->Z() > fDetLength + fZCutOff) return false;

    PackMCTruth(record,truth); 

    ///pack the flux information
    if(fFluxType.compare("ntuple") == 0){
      fSpillTotal += dynamic_cast<genie::flux::GNuMIFlux *>(fFluxD)->UsedPOTs()/fDriver->GlobProbScale();
      ///check that the produced record is an allowed neutrino flavor as genie doesnt do that for you
      ///return false if it is not
      if(fGenFlavors.find(truth.GetNeutrino().Nu().PdgCode()) != fGenFlavors.end()){
	flux.fFluxType = simb::kNtuple;
	PackNuMIFlux(flux);
      }
      else return false;
    }
    else if ( fFluxType.compare("simple_flux")==0 ) {
 
#ifdef MISSING_GSIMPLENTPFLUX
      std::cerr << "Not built with GSimpleNtpFlux enabled" << std::endl;
      assert(0);
#else
      ///pack the flux information
      fSpillTotal += dynamic_cast<genie::flux::GSimpleNtpFlux *>(fFluxD)->UsedPOTs()/fDriver->GlobProbScale();
#endif
      ///check that the produced record is an allowed neutrino flavor as genie doesnt do that for you
      ///return false if it is not
      if ( fGenFlavors.find(truth.GetNeutrino().Nu().PdgCode()) != fGenFlavors.end()){
	flux.fFluxType = simb::kSimple_Flux;
	PackSimpleFlux(flux);
      }
      else return false;
    }
    else if(fFluxType.compare("histogram") == 0){
      ///set the flag in the parent object that says the 
      ///fluxes came from histograms and fill related values
      flux.fFluxType = simb::kHistPlusFocus;

      ///save the fluxes - fluxes were added to the vector in the same 
      ///order that the flavors appear in fGenFlavors
      int ctr = 0;
      int bin = fFluxHistograms[0]->FindBin(truth.GetNeutrino().Nu().E());
      std::vector<double> fluxes(6, 0.);
      for(std::set<int>::iterator i = fGenFlavors.begin(); i != fGenFlavors.end(); i++){
	if(*i ==  12) fluxes[kNue]      = fFluxHistograms[ctr]->GetBinContent(bin);
	if(*i == -12) fluxes[kNueBar]   = fFluxHistograms[ctr]->GetBinContent(bin);
	if(*i ==  14) fluxes[kNuMu]     = fFluxHistograms[ctr]->GetBinContent(bin);
	if(*i == -14) fluxes[kNuMuBar]  = fFluxHistograms[ctr]->GetBinContent(bin);
	if(*i ==  16) fluxes[kNuTau]    = fFluxHistograms[ctr]->GetBinContent(bin);
	if(*i == -16) fluxes[kNuTauBar] = fFluxHistograms[ctr]->GetBinContent(bin);
	++ctr;
      }

      ///get the flux for each neutrino flavor of this energy
      flux.SetFluxGen(fluxes[kNue],   fluxes[kNueBar],
		      fluxes[kNuMu],  fluxes[kNuMuBar],
		      fluxes[kNuTau], fluxes[kNuTauBar]);
    
      fSpillTotal += 1.;
    }
    else if(fFluxType.compare("mono") == 0){
      fSpillTotal += 1.;    
    }

    return true;
  }

  //--------------------------------------------------
  void GENIEHelper::PackNuMIFlux(simb::MCFlux &flux)
  {
    flux.Reset();

    ///cast the fFluxD pointer to be of the right type
    genie::flux::GNuMIFlux *gnf = dynamic_cast<genie::flux::GNuMIFlux *>(fFluxD);
    const genie::flux::GNuMIFluxPassThroughInfo& nflux = gnf->PassThroughInfo();

    ///check the particle codes and the units passed through
    /// nflux.pcodes: 0=original GEANT particle codes, 1=converted to PDG
    /// nflux.units:  0=original GEANT cm, 1=meters
    if(nflux.pcodes != 1 && nflux.units != 0)
      std::cerr << "either wrong particle codes or units from flux object - beware!!" 
		<< std::endl;

    // maintained variable names from gnumi ntuples
    // see http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/[/v19/output_gnumi.html]

    flux.frun      = nflux.run;
    flux.fevtno    = nflux.evtno;
    flux.fndxdz    = nflux.ndxdz;
    flux.fndydz    = nflux.ndydz;
    flux.fnpz      = nflux.npz;
    flux.fnenergy  = nflux.nenergy;
    flux.fndxdznea = nflux.ndxdznea;
    flux.fndydznea = nflux.ndydznea;
    flux.fnenergyn = nflux.nenergyn;
    flux.fnwtnear  = nflux.nwtnear;
    flux.fndxdzfar = nflux.ndxdzfar;
    flux.fndydzfar = nflux.ndydzfar;
    flux.fnenergyf = nflux.nenergyf;
    flux.fnwtfar   = nflux.nwtfar;
    flux.fnorig    = nflux.norig;
    flux.fndecay   = nflux.ndecay;
    flux.fntype    = nflux.ntype;
    flux.fvx       = nflux.vx;
    flux.fvy       = nflux.vy;
    flux.fvz       = nflux.vz;
    flux.fpdpx     = nflux.pdpx;
    flux.fpdpy     = nflux.pdpy;
    flux.fpdpz     = nflux.pdpz;
    flux.fppdxdz   = nflux.ppdxdz;
    flux.fppdydz   = nflux.ppdydz;
    flux.fpppz     = nflux.pppz;
    flux.fppenergy = nflux.ppenergy;
    flux.fppmedium = nflux.ppmedium;
    flux.fptype    = nflux.ptype;     // converted to PDG
    flux.fppvx     = nflux.ppvx;
    flux.fppvy     = nflux.ppvy;
    flux.fppvz     = nflux.ppvz;
    flux.fmuparpx  = nflux.muparpx;
    flux.fmuparpy  = nflux.muparpy;
    flux.fmuparpz  = nflux.muparpz;
    flux.fmupare   = nflux.mupare;
    flux.fnecm     = nflux.necm;
    flux.fnimpwt   = nflux.nimpwt;
    flux.fxpoint   = nflux.xpoint;
    flux.fypoint   = nflux.ypoint;
    flux.fzpoint   = nflux.zpoint;
    flux.ftvx      = nflux.tvx;
    flux.ftvy      = nflux.tvy;
    flux.ftvz      = nflux.tvz;
    flux.ftpx      = nflux.tpx;
    flux.ftpy      = nflux.tpy;
    flux.ftpz      = nflux.tpz;
    flux.ftptype   = nflux.tptype;   // converted to PDG
    flux.ftgen     = nflux.tgen;
    flux.ftgptype  = nflux.tgptype;  // converted to PDG
    flux.ftgppx    = nflux.tgppx;
    flux.ftgppy    = nflux.tgppy;
    flux.ftgppz    = nflux.tgppz;
    flux.ftprivx   = nflux.tprivx;
    flux.ftprivy   = nflux.tprivy;
    flux.ftprivz   = nflux.tprivz;
    flux.fbeamx    = nflux.beamx;
    flux.fbeamy    = nflux.beamy;
    flux.fbeamz    = nflux.beamz;
    flux.fbeampx   = nflux.beampx;
    flux.fbeampy   = nflux.beampy;
    flux.fbeampz   = nflux.beampz;    

    return;
  }

  //--------------------------------------------------
  void GENIEHelper::PackMCTruth(genie::EventRecord *record,
				simb::MCTruth &truth)
  {

    TLorentzVector *vertex = record->Vertex();

    //std::cout << "vertex loc " << vertex->X() << "," << vertex->Y() << "," << vertex->Z() << std::endl;  

    ///get the Interaction object from the record - this is the object
    ///that talks to the event information objects and is in m
    genie::Interaction *inter = record->Summary();
  
    ///get the different components making up the interaction
    const genie::InitialState &initState  = inter->InitState();
    const genie::ProcessInfo  &procInfo   = inter->ProcInfo();
    const genie::Kinematics   &kine       = inter->Kine();
    //const genie::XclsTag      &exclTag    = inter->ExclTag();
    //const genie::KPhaseSpace  &phaseSpace = inter->PhaseSpace();

    //choose a spill time (ns) to shift the vertex times by:

    double spillTime = fGlobalTimeOffset + gRandom->Uniform()*fRandomTimeOffset;

    ///add the particles from the interaction
    TIter partitr(record);
    genie::GHepParticle *part = 0;
    ///GHepParticles return units of GeV/c for p.  the V_i are all in fermis
    ///and are relative to the center of the struck nucleus.
    ///add the vertex X/Y/Z to the V_i for status codes 0 and 1
    int trackid = 0;
    std::string primary("primary");
    while( (part = dynamic_cast<genie::GHepParticle *>(partitr.Next())) ){
      --trackid;
      simb::MCParticle tpart(trackid, 
			     part->Pdg(), 
			     primary, 
			     part->FirstMother(), 
			     part->Mass(), 
			     part->Status());
      
      ///set the vertex location for the neutrino, nucleus and everything
      ///that is to be tracked.  vertex returns values in meters.
      if(part->Status() == 0 || part->Status() == 1){
	TLorentzVector pos(100.*(part->Vx()*1.e-15 + vertex->X()),
			   100.*(part->Vy()*1.e-15 + vertex->Y()),
			   100.*(part->Vz()*1.e-15 + vertex->Z()),
			   part->Vt() + spillTime);
	TLorentzVector mom(part->Px(), part->Py(), part->Pz(), part->E());
	tpart.AddTrajectoryPoint(pos,mom);
      }

      truth.Add(tpart);
        
    }///end loop to convert GHepParticles to MCParticles

    ///is the interaction NC or CC
    int CCNC = simb::kCC;
    if(procInfo.IsWeakNC()) CCNC = simb::kNC;

    ///what is the interaction type
    int mode = simb::kQE;
    if(procInfo.IsDeepInelastic()) mode = simb::kDIS;
    else if(procInfo.IsResonant()) mode = simb::kRes;
    else if(procInfo.IsCoherent()) mode = simb::kCoh;

    ///set the neutrino information in MCTruth
    truth.SetOrigin(simb::kBeamNeutrino);
    truth.SetNeutrino(CCNC, mode, 
		      initState.Tgt().Pdg(), 
		      initState.Tgt().HitNucPdg(), 
		      initState.Tgt().HitQrkPdg(),
		      kine.W(true), kine.x(true), kine.y(true), kine.Q2(true));


    return;
  }

  void GENIEHelper::PackSimpleFlux(simb::MCFlux &flux)
  {
#ifdef MISSING_GSIMPLENTPFLUX
    std::cerr << "Not built with GSimpleNtpFlux enabled" << std::endl;
    assert(0);
#else
    flux.Reset();

    ///cast the fFluxD pointer to be of the right type
    genie::flux::GSimpleNtpFlux *gsf = 
      dynamic_cast<genie::flux::GSimpleNtpFlux *>(fFluxD);
    
    // maintained variable names from gnumi ntuples
    // see http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/[/v19/output_gnumi.html]
    
    const genie::flux::GSimpleNtpEntry* nflux_entry = gsf->GetCurrentEntry();
    const genie::flux::GSimpleNtpNuMI*  nflux_numi  = gsf->GetCurrentNuMI();
    //const genie::flux::GSimpleNtpMeta*  nflux_meta  = gsf->GetCurrentMeta();
  
    flux.fntype  = nflux_entry->pdg;
    flux.fnimpwt = nflux_entry->wgt;

    if ( nflux_numi ) {
      flux.frun      = nflux_numi->run;
      flux.fevtno    = nflux_numi->evtno;
      flux.ftpx      = nflux_numi->tpx;
      flux.ftpy      = nflux_numi->tpy;
      flux.ftpz      = nflux_numi->tpz;
      flux.ftptype   = nflux_numi->tptype;   // converted to PDG
    }
    //   flux.fndxdz    = nflux.ndxdz;
    //   flux.fndydz    = nflux.ndydz;
    //   flux.fnpz      = nflux.npz;
    //   flux.fnenergy  = nflux.nenergy;
    //   flux.fndxdznea = nflux.ndxdznea;
    //   flux.fndydznea = nflux.ndydznea;
    //   flux.fnenergyn = nflux.nenergyn;
    //   flux.fnwtnear  = nflux.nwtnear;
    //   flux.fndxdzfar = nflux.ndxdzfar;
    //   flux.fndydzfar = nflux.ndydzfar;
    //   flux.fnenergyf = nflux.nenergyf;
    //   flux.fnwtfar   = nflux.nwtfar;
    //   flux.fnorig    = nflux.norig;
    //   flux.fndecay   = nflux.ndecay;
    //   flux.fntype    = nflux.ntype;
    //   flux.fvx       = nflux.vx;
    //   flux.fvy       = nflux.vy;
    //   flux.fvz       = nflux.vz;
    //   flux.fpdpx     = nflux.pdpx;
    //   flux.fpdpy     = nflux.pdpy;
    //   flux.fpdpz     = nflux.pdpz;
    //   flux.fppdxdz   = nflux.ppdxdz;
    //   flux.fppdydz   = nflux.ppdydz;
    //   flux.fpppz     = nflux.pppz;
    //   flux.fppenergy = nflux.ppenergy;
    //   flux.fppmedium = nflux.ppmedium;
    //   flux.fptype    = nflux.ptype;     // converted to PDG
    //   flux.fppvx     = nflux.ppvx;
    //   flux.fppvy     = nflux.ppvy;
    //   flux.fppvz     = nflux.ppvz;
    //   flux.fmuparpx  = nflux.muparpx;
    //   flux.fmuparpy  = nflux.muparpy;
    //   flux.fmuparpz  = nflux.muparpz;
    //   flux.fmupare   = nflux.mupare;
    //   flux.fnecm     = nflux.necm;
    //   flux.fnimpwt   = nflux.nimpwt;
    //   flux.fxpoint   = nflux.xpoint;
    //   flux.fypoint   = nflux.ypoint;
    //   flux.fzpoint   = nflux.zpoint;
    //   flux.ftvx      = nflux.tvx;
    //   flux.ftvy      = nflux.tvy;
    //   flux.ftvz      = nflux.tvz;
    //   flux.ftgen     = nflux.tgen;
    //   flux.ftgptype  = nflux.tgptype;  // converted to PDG
    //   flux.ftgppx    = nflux.tgppx;
    //   flux.ftgppy    = nflux.tgppy;
    //   flux.ftgppz    = nflux.tgppz;
    //   flux.ftprivx   = nflux.tprivx;
    //   flux.ftprivy   = nflux.tprivy;
    //   flux.ftprivz   = nflux.tprivz;
    //   flux.fbeamx    = nflux.beamx;
    //   flux.fbeamy    = nflux.beamy;
    //   flux.fbeamz    = nflux.beamz;
    //   flux.fbeampx   = nflux.beampx;
    //   flux.fbeampy   = nflux.beampy;
    //   flux.fbeampz   = nflux.beampz; 
#endif
    
    return;
  }

}
