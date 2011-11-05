////////////////////////////////////////////////////////////////////////
// $Id: RecoDrawingOption.cxx,v 1.16 2010/08/30 21:33:24 spitz7 Exp $
//
// Display parameters for the raw data 
//
// \author jpaley@indiana.edu
////////////////////////////////////////////////////////////////////////
#include "EventDisplay/RecoDrawingOptions.h"
#include <iostream>

namespace evd {

  //......................................................................
  RecoDrawingOptions::RecoDrawingOptions(fhicl::ParameterSet const& pset, 
					 art::ActivityRegistry& reg) 
  {
    this->reconfigure(pset);
  }
  
  //......................................................................
  RecoDrawingOptions::~RecoDrawingOptions() 
  {
  }

  //......................................................................
  void RecoDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
    fDrawHits        = pset.get< int >("DrawHits"       );
    fDrawClusters    = pset.get< int >("DrawClusters" 	);
    fDrawProngs      = pset.get< int >("DrawProngs"   	);
    fDrawTracks      = pset.get< int >("DrawTracks"   	);
    fDrawShowers     = pset.get< int >("DrawShowers"  	);
    fDrawVertices    = pset.get< int >("DrawVertices" 	);
    fDrawEvents      = pset.get< int >("DrawEvents"   	);
    fDraw2DEndPoints = pset.get< int >("Draw2DEndPoints");

    fHitLabels     = pset.get< std::vector<std::string> >("HitModuleLabels"    );
    fProngLabels   = pset.get< std::vector<std::string> >("ProngModuleLabels"  );
    fClusterLabels = pset.get< std::vector<std::string> >("ClusterModuleLabels");
    fTrackLabels   = pset.get< std::vector<std::string> >("TrackModuleLabels"  );
    fShowerLabels  = pset.get< std::vector<std::string> >("ShowerModuleLabels" );
    fVertexLabels  = pset.get< std::vector<std::string> >("VertexModuleLabels" );
    fEventLabels   = pset.get< std::vector<std::string> >("EventModuleLabels"  );
    fWireLabels    = pset.get< std::vector<std::string> >("WireModuleLabels"   );
  }
  
}
////////////////////////////////////////////////////////////////////////
