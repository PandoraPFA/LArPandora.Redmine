////////////////////////////////////////////////////////////////////////
//
// AggregateEvent module
//
// echurch@fnal.gov
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
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
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "EventFinder/AggregateEvent.h"
#include "Geometry/geo.h"

// ROOT includes
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TH2F.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"

namespace event {

  //-------------------------------------------------
  AggregateEvent::AggregateEvent(fhicl::ParameterSet const& pset) : 
    fVertexModuleLabel(pset.get< std::string >("VertexModuleLabel"))
  {

    produces< std::vector<recob::Event> >();

  }

  //-------------------------------------------------
  AggregateEvent::~AggregateEvent()
  {
  }

  //-------------------------------------------------
  void AggregateEvent::beginJob()
  {
  }

  //------------------------------------------------------------------------------------//
  void AggregateEvent::produce(art::Event& evt)
  { 

    std::auto_ptr<std::vector<recob::Event> > ecol(new std::vector<recob::Event>);
    
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    art::Handle< std::vector<recob::Vertex> > vertexListHandle;
    evt.getByLabel(fVertexModuleLabel,vertexListHandle);
    for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii){
      art::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
      fvertexlist.push_back(vertex); // class member
    }

  }

} // end namespace




