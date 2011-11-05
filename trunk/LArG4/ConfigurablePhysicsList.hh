////////////////////////////////////////////////////////////////////////
/// \file CustomPhysicsList.h
//
/// \version $Id: PrimaryParticleInformation.cxx,v 1.3 2009/10/05 23:21:51 t962cvs Exp $
/// \author  bjpj@mit.edu
////////////////////////////////////////////////////////////////////////
// Configurable Physics List Class
//
// Ben Jones, MIT, 24/06/09
//
// Based on the QGSP_BERT physics list supplied with geant 4, but with
// options to switch on and off different physics processed from the config
// file. 
//
// To include a new physics process, create a new physics builder registering
// the relevant particles and interactions.  Then add the builder to the 
// GetPhysicsBuilders and GetDefaultSettings functions in this class,
// using the name by which the builder will be referenced in the config
// file.  Physics builders to be enabled are specified in the LArG4 config.
//
//


#ifndef TConfigurablePhysicsList_h
#define TConfigurablePhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"

namespace larg4 {

  template<class T>
  class TConfigurablePhysicsList: public T
  {
  public:
    TConfigurablePhysicsList(G4int ver = 1);
    virtual ~TConfigurablePhysicsList();
    bool UseCustomPhysics;
    std::vector<std::string> EnabledPhysics;
    std::vector<std::string> GetDefaultSettings();

  
  public:
    // SetCuts() 
    virtual void SetCuts();
    
  private:
    enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
  };

}
  
#include "LArG4/ConfigurablePhysicsList.icc"

namespace larg4 { 
  typedef TConfigurablePhysicsList<G4VModularPhysicsList> ConfigurablePhysicsList; 
}
  


#endif



// Sept 2009 - Ben Jones, MIT 

