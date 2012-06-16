////////////////////////////////////////////////////////////////////////
// $Id: EVD.cxx,v 1.2 2010/11/10 22:49:25 p-novaart Exp $
//
// EVD event display module for Argoneut
//
// brebel@fnal.gov
// msoderbe@syr.edu
//
//
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

//LArSoft includes
#include "EventDisplayBase/evdb.h"
#include "T962/ArgoneutEventDisplay/ArgoneutEVD.h"
#include "T962/ArgoneutEventDisplay/DisplayArgoneutView.h"
#include "EventDisplay/TWQProjectionView.h"
#include "EventDisplay/Display3DView.h"
#include "EventDisplay/Ortho3DView.h"


// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Builder for the TWQProjectionView canvas
static evdb::Canvas* mk_twqprojectionview_canvas(TGMainFrame* mf)
{
  return new evd::TWQProjectionView(mf);
}

// Builder for the Display3D view
static evdb::Canvas* mk_display3d_canvas(TGMainFrame* mf)
{
  return new evd::Display3DView(mf);
}

// Builder for the Ortho view
static evdb::Canvas* mk_ortho3d_canvas(TGMainFrame* mf)
{
  return new evd::Ortho3DView(mf);
}

// Builder for the ArgoNeuT 3D view
static evdb::Canvas* mk_displayargoneut_canvas(TGMainFrame* mf)
{
  return new argoevd::DisplayArgoneutView(mf);
}

// // Builder for the MCTruth view
// static evdb::ObjListCanvas* mk_mctrue_canvas(TGMainFrame* mf)
// {
//   return new evd::MCTrueView(mf);
// }


namespace argoevd{

  //----------------------------------------------------
  ArgoneutEVD::ArgoneutEVD(fhicl::ParameterSet const& pset)
  {

  }

  //----------------------------------------------------
  ArgoneutEVD::~ArgoneutEVD()
  {
  }

  //----------------------------------------------------
  void ArgoneutEVD::beginJob()
  {
    // Register the list of windows used by the event display
    evdb::DisplayWindow::Register("Time vs Wire, Charge View",
				  "Time vs Wire, Charge View",
				  700,
				  700,
				  mk_twqprojectionview_canvas);

    evdb::DisplayWindow::Register("Display3D",
				  "Display3D",
				  700,
				  700,
				  mk_display3d_canvas);

    evdb::DisplayWindow::Register("Ortho3D",
                                  "Ortho3D",
                                  700,
                                  700,
                                  mk_ortho3d_canvas);


  evdb::DisplayWindow::Register("DisplayArgoneut",
				  "DisplayArgoneut",
				  700,
				  700,
				  mk_displayargoneut_canvas);
    
//     evdb::ListWindow::Register("MC Particle List",
// 			       "MC Particle List",
// 			       400,
// 			       800,
// 			       mk_mctrue_canvas);
 
    // Open up the main display window and run
    evdb::DisplayWindow::OpenWindow(0);

  }

  //----------------------------------------------------
  void ArgoneutEVD::analyze(const art::Event& evt)
  {
  }

}//namespace
