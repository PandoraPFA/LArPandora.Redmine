////////////////////////////////////////////////////////////////////////
/// \file    PhotonLibraryBuilder.h
/// \version $Id: CalHit.cxx,v 1.9 2011/01/24 23:18:19 p-novaart Exp $
/// \author  bpjones
////////////////////////////////////////////////////////////////////////

#ifndef __CINT__

#include "art/Framework/Core/EDAnalyzer.h"
#include "Simulation/PMTHit.h"
#include "PhotonPropagation/PhotonLibrary.h"
#include "Simulation/PhotonLibraryParameters.h"

// ROOT includes.
#include <Rtypes.h>
#ifndef PhotonLibraryBuilder_h
#define PhotonLibraryBuilder_h 1

///photon propagation tools
namespace phot {

  class PhotonLibraryBuilder : public art::EDAnalyzer{
    public:
      
      PhotonLibraryBuilder(const fhicl::ParameterSet&);
      virtual ~PhotonLibraryBuilder();
      
      void analyze(art::Event const&);
      
      void beginJob();
      
     
    private:
      
      // Parameters to read in

      std::string fLArG4InputModule;      // Input tag for PMT collection
      std::string fEvGenInputModule;      // Input tag for PMT collection

      PhotonLibrary * fLibrary;
    };
}

#endif
#endif
