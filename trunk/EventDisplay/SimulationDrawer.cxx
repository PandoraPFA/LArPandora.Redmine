///
/// \file    SimulationDrawer.cxx
/// \brief   Render the objects from the Simulation package
/// \author  messier@indiana.edu
/// \version $Id: SimulationDrawer.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
#include "TParticle.h"
#include "TLatex.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"

#include "EventDisplay/SimulationDrawer.h"
#include "EventDisplayBase/evdb.h"
#include "Geometry/geo.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "EventDisplay/Style.h"
#include "EventDisplay/SimulationDrawingOptions.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace {
  // Utility function to make uniform error messages.
  void writeErrMsg(const char* fcn,
		   cet::exception const& e)
  {
    mf::LogWarning("SimulationDrawer") << "SimulationDrawer::" << fcn
				       << " failed with message:\n"
				       << e;
  }
}

namespace evd{

  SimulationDrawer::SimulationDrawer() 
  {
  }

  //......................................................................

  SimulationDrawer::~SimulationDrawer()
  {
  }

  //......................................................................

  void SimulationDrawer::MCTruthShortText(const art::Event& evt,
					  evdb::View2D*     view) 
  {

    if( evt.isRealData() ) return;

    art::ServiceHandle<evd::SimulationDrawingOptions> drawopt;
    // Skip drawing if option is turned off
    if (!drawopt->fShowMCTruthText) return;

    art::PtrVector<simb::MCTruth> mctruth;
    this->GetMCTruth(evt, mctruth);
    
    for (unsigned int i=0; i<mctruth.size(); ++i) {
      std::string mctext;
      bool firstin  = true;
      bool firstout = true;
      std::string origin;
      std::string incoming;
      std::string outgoing;
      // Label cosmic rays -- others are pretty obvious
      if (mctruth[i]->Origin()==simb::kCosmicRay)  origin = "c-ray: "; 
      for (int j=0; j<mctruth[i]->NParticles(); ++j) {
	const simb::MCParticle& p = mctruth[i]->GetParticle(j);
	char buff[1024];
	if (p.P()>0.05) {
	  sprintf(buff,"#color[%d]{%s #scale[0.75]{[%.1f GeV/c]}}", 
		  Style::ColorFromPDG(p.PdgCode()),
		  Style::LatexName(p.PdgCode()), 
		  p.P());
	}
	else {
	  sprintf(buff,"#color[%d]{%s}", 
		  Style::ColorFromPDG(p.PdgCode()),
		  Style::LatexName(p.PdgCode()));
	}
	if (p.StatusCode()==0) {
	  if (firstin==false) incoming += " + ";
	  incoming += buff;
	  firstin = false;
	}
	if (p.StatusCode()==1) {
	  if (firstout==false) outgoing += " + ";
	  outgoing += buff;
	  firstout = false;
	}
      } // loop on j particles
      if (origin=="" && incoming=="") {
	mctext = outgoing;
      }
      else {
	mctext = origin+incoming+" #rightarrow "+outgoing;
      }
      TLatex& latex = view->AddLatex(0.03, 0.2, mctext.c_str());
      latex.SetTextSize(0.6);

    } // loop on i mctruth objects

  }

  //......................................................................

  void SimulationDrawer::MCTruthLongText(const art::Event& evt,
					 evdb::View2D* /*view*/) 
  {
    if( evt.isRealData() ) return;

    art::ServiceHandle<evd::SimulationDrawingOptions> drawopt;
    // Skip drawing if option is turned off
    if (!drawopt->fShowMCTruthText) return;

    art::PtrVector<simb::MCTruth> mctruth;
    this->GetMCTruth(evt, mctruth);
    
    for (unsigned int i=0; i<mctruth.size(); ++i) {
      std::cout << "\n";
      for (int j=0; j<mctruth[i]->NParticles(); ++j) {
	const simb::MCParticle& p = mctruth[i]->GetParticle(j);
	if(p.StatusCode() == 0 || p.StatusCode() == 1)
	  std::cout << Style::LatexName(p.PdgCode())
		    << "\t\t" << p.E() << " GeV"
		    << "\t"   << "(" << p.P() << " GeV/c)"
		    << std::endl;
      } // loop on j particles in list
    }
  }


  //......................................................................
  //this is the method you would use to color code hits by the MC truth pdg value
  void SimulationDrawer::MCTruthVectors2D(const art::Event& evt,
					  evdb::View2D*     view,
					  unsigned int      plane)
  {
    if( evt.isRealData() ) return;

    art::ServiceHandle<evd::SimulationDrawingOptions> drawopt;
    // If the option is turned off, there's nothing to do
    if (!drawopt->fShowMCTruthVectors) return;

    art::ServiceHandle<util::DetectorProperties> detprop;
    double sampleRate = detprop->SamplingRate();

    art::ServiceHandle<geo::Geometry>          geo;
    art::ServiceHandle<util::LArProperties>    larp;
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    // get the x position of the plane in question
    double xyz[3]  = {0.};
    double xyz2[3] = {0.};
    double loc[3]  = {0.};
    geo->Plane(plane).LocalToWorld(loc, xyz);
    double planex = xyz[0];

    // Unpack and draw the MC vectors
    art::PtrVector<simb::MCTruth> mctruth;
    this->GetMCTruth(evt, mctruth);
    
    for (unsigned int i=0; i<mctruth.size(); ++i) {
      for (int j=0; j<mctruth[i]->NParticles(); ++j) {
	const simb::MCParticle& p = mctruth[i]->GetParticle(j);
	
	// Skip all but incoming and out-going particles
	if (!(p.StatusCode()==0 || p.StatusCode()==1)) continue;
	
	double r  = p.P()*10.0;           // Scale length so 10 cm = 1 GeV/c
	if (r < 0.1) continue;            // Skip very short particles
	if (p.StatusCode() == 0) r = -r;  // Flip for incoming particles

  	xyz[0]  = planex;
	xyz[1]  = p.Vy();
	xyz[2]  = p.Vz();
	xyz2[0] = planex;
	xyz2[1] = xyz[1] + r * p.Py()/p.P();
	xyz2[2] = xyz[2] + r * p.Pz()/p.P();
	
	if(xyz2[2] < 0.)                          xyz2[2] = 0.;
	if(xyz2[2] > geo->DetLength() )           xyz2[2] = geo->DetLength();
	if(fabs(xyz2[1]) > geo->DetHalfHeight() ) xyz2[1] = geo->DetHalfHeight();

	unsigned int chan1 = geo->NearestChannel(xyz); 
	unsigned int chan2 = geo->NearestChannel(xyz2); 
	unsigned int pln   = 0;
	unsigned int tpc   = 0;
	unsigned int w1    = 0;
	unsigned int w2    = 0;
	geo->ChannelToWire(chan1, tpc, pln, w1);
	geo->ChannelToWire(chan2, tpc, pln, w2);
	
	// check that we have the correct plane and TPC
	if(pln != plane || tpc != rawopt->fTPC) continue;

	double time = p.Vx()*1e3/larp->DriftVelocity(larp->Efield(), larp->Temperature())/sampleRate;
	double time2 = (p.Vx() + r * p.Px()/p.P())*1e3/larp->DriftVelocity(larp->Efield(), larp->Temperature())/sampleRate;

	if(rawopt->fAxisOrientation < 1){
	  TLine& l = view->AddLine(w1, time, w2, time2);
	  evd::Style::FromPDG(l, p.PdgCode());
	}
	else{
	  TLine& l = view->AddLine(time, w1, time2, w2);
	  evd::Style::FromPDG(l, p.PdgCode());
	}

      } // loop on j particles in list
    } // loop on truths

  }

  //......................................................................
  //this method draws the true particle trajectories in 3D
  void SimulationDrawer::MCTruth3D(const art::Event& evt,
				   evdb::View3D*     view)
  {
    if( evt.isRealData() ) return;

    art::ServiceHandle<evd::SimulationDrawingOptions> drawopt;
    // If the option is turned off, there's nothing to do
    if (!drawopt->fShowMCTruthTrajectories) return;

    art::ServiceHandle<geo::Geometry> geom;

    // get the particles from the Geant4 step
    art::PtrVector<sim::Particle> plist;
    this->GetParticle(evt, plist);
    
    // loop over the LArVoxelList to get the true energy deposition
    // locations
    const sim::LArVoxelList voxels = sim::SimListUtils::GetLArVoxelList(evt,drawopt->fG4ModuleLabel);
								      
    // loop over all the particles
    for(size_t p = 0; p < plist.size(); ++p){

      // collect the points from this particle
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> z;

      // get the limits of the cryostat in the world coordinates
      double cryobounds[6] = {0.};
      geom->CryostatBoundaries(cryobounds);

      // this is slow, and there has to be a better way
      // loop over the voxel list and get the centers for each voxel
      // where the energy deposited by the current particle is above the 
      // threshold

      sim::LArVoxelList::const_iterator vxitr;
      for(vxitr = voxels.begin(); vxitr != voxels.end(); vxitr++){
	const sim::LArVoxelData &vxd = (*vxitr).second;
	if( vxd.find(plist[p]->TrackId()) != vxd.end() ){
	  if(vxd[plist[p]->TrackId()] > drawopt->fMinEnergyDeposition){                  	  
	    x.push_back(vxd.VoxelID().X());
	    y.push_back(vxd.VoxelID().Y());
	    z.push_back(vxd.VoxelID().Z());
	  }
	} // end if this track id is in the current voxel
      }// end loop over voxels
	  
      TPolyMarker3D& pm = view->AddPolyMarker3D(x.size(), 
						evd::Style::ColorFromPDG(plist[p]->PdgCode()), 
						1, 3);
      for(size_t i = 0; i < x.size(); ++i)
	pm.SetPoint(i, x[i], y[i], z[i]);
      
    }
    
    return;
  }

  //......................................................................
  int SimulationDrawer::GetParticle(const art::Event&              evt,
				    art::PtrVector<sim::Particle>& plist)
  {
    plist.clear();

    if( evt.isRealData() ) return 0;

    art::PtrVector<sim::Particle> temp;

    std::vector< art::Handle< std::vector<sim::Particle> > > plcol;
    // use get by Type because there should only be one collection of these in the event
    try{
      evt.getManyByType(plcol);
      for(size_t plc = 0; plc < plcol.size(); ++plc){
	for(unsigned int i = 0; i < plcol[plc]->size(); ++i){
	  art::Ptr<sim::Particle> mct(plcol[plc], i);
	  temp.push_back(mct);
	}
      }
      temp.swap(plist);
    }
    catch(cet::exception& e){
      writeErrMsg("GetRawDigits", e);
    }
  
    return plist.size();

  }

  //......................................................................

  int SimulationDrawer::GetMCTruth(const art::Event& evt,
				   art::PtrVector<simb::MCTruth>& mcvec) 
  {
    mcvec.clear();

    if( evt.isRealData() ) return 0;

    art::PtrVector<simb::MCTruth> temp;

    std::vector< art::Handle< std::vector<simb::MCTruth> > > mctcol;
    // use get by Type because there should only be one collection of these in the event
    try{
      evt.getManyByType(mctcol);
      for(size_t mctc = 0; mctc < mctcol.size(); ++mctc){
	for(unsigned int i = 0; i < mctcol[mctc]->size(); ++i){
	art::Ptr<simb::MCTruth> mct(mctcol[mctc], i);
	temp.push_back(mct);
      }
      }
      temp.swap(mcvec);
    }
    catch(cet::exception& e){
      writeErrMsg("GetMCTruth", e);
    }
  
    return mcvec.size();
  }


  //......................................................................

  void SimulationDrawer::HiLite(int trkId, bool dohilite)
  {
    fHighlite[trkId] = dohilite;
  }

}// namespace
////////////////////////////////////////////////////////////////////////
