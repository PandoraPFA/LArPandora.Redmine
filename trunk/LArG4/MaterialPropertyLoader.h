////////////////////////////////////////////////////////////////////////
/// \file MaterialPropertyLoader.h
//
/// \version $Id: PrimaryParticleInformation.cxx,v 1.3 2009/10/05 23:21:51 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
// Class to set material properties for different materials in
// the detector. Most likely use is via the child class
// MaterialPropertyLoaderFromXML to read material properties
// from an XML file.  But this base class is included in case
// another mechanism for loading material properties becomes
// necessary.

#ifndef LArG4_MaterialPropertyLoader_h
#define LArG4_MaterialPropertyLoader_h

#include <G4LogicalVolumeStore.hh>
#include <map>

namespace larg4 {
  
  class MaterialPropertyLoader
  {
  public:
    
    MaterialPropertyLoader() {}
    ~MaterialPropertyLoader() {}
    
    
    
  private:

    //         materials                properties            values
    std::map < std::string , std::map < std::string,double> > constPropertyList;

    //         materials                properties               energies  values
    std::map < std::string , std::map < std::string , std::map < double ,  double > > > propertyList;

    std::map<std::string, double> birksConstants;
    
    
  public:
    
    //Accessors - mostly trivial
    std::map<double,double> GetMaterialProperty(std::string Material,std::string Property) 
      {return propertyList[Material][Property];}
    
    double GetMaterialConstProperty(std::string Material, std::string Property)
    {return constPropertyList[Material][Property];}
    
    std::map<std::string,double> GetMaterialConstProperties(std::string Material)
      {return constPropertyList[Material];}
    
    std::map<std::string,std::map<double,double> >  GetMaterialProperties(std::string Material)
      {return propertyList[Material];}
    
     
    // For some reason birks constant is not set like other properties...
    void SetBirksConstant(std::string, double);
    
    
    // Setters. Include some checks on sizes of vectors etc
    void SetMaterialProperty(std::string,std::string,std::map<double,double>);
    void SetMaterialConstProperty(std::string,std::string,double);
    
    
    // Geometry updater.  Generate Geant4 material property tables and attach to detector geometry
    void UpdateGeometry(G4LogicalVolumeStore*);
    
    
  };  
  
}
#endif
