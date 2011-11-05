////////////////////////////////////////////////////////////////////////
/// \file MaterialPropertyLoader.cxx
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



#include "LArG4/MaterialPropertyLoader.h"
#include <G4Material.hh>
#include <G4MaterialPropertiesTable.hh>
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace larg4 {

  //----------------------------------------------
  void MaterialPropertyLoader::SetMaterialProperty(std::string Material,
						   std::string Property, 
						   std::map<double, double> PropertyVector)
  {
    propertyList[Material][Property]=PropertyVector;
    mf::LogInfo("MaterialPropertyLoader")<<"Added property " 
					 << Material<< "  " 
					 << Property;
  }
  
  //----------------------------------------------
  void MaterialPropertyLoader::SetMaterialConstProperty(std::string Material, 
							std::string Property, 
							double PropertyValue)
  {
    constPropertyList[Material][Property]=PropertyValue;
    mf::LogInfo("MaterialPropertyLoader") << "Added const property " 
					  << Material << "  " 
					  << Property;
  }
  
  //----------------------------------------------
  void MaterialPropertyLoader::SetBirksConstant(std::string Material, 
						double PropertyValue)
  {
    birksConstants[Material]=PropertyValue;	
    mf::LogInfo("MaterialPropertyLoader") << "Set Birks constant " 
					  << Material;
  }

  //----------------------------------------------  
  void MaterialPropertyLoader::UpdateGeometry(G4LogicalVolumeStore * lvs)
  {
    std::map<std::string,G4MaterialPropertiesTable*> MaterialTables;
    std::map<std::string,bool> MaterialsSet;
    
    mf::LogInfo("MaterialPropertyLoader") << "UPDATING GEOMETRY";
    
    // Loop over each material with a property vector and create a new material table for it
    for(std::map<std::string,std::map<std::string,std::map<double,double> > >::const_iterator i=propertyList.begin(); i!=propertyList.end(); i++){
      std::string Material=i->first;
      MaterialsSet[Material]=true;
      MaterialTables[Material]=new G4MaterialPropertiesTable;
    }
    
    // Loop over each material with a const property, 
    // if material table does not exist, create one
    for(std::map<std::string,std::map<std::string,double> >::const_iterator i=constPropertyList.begin(); i!=constPropertyList.end(); i++){
      std::string Material=i->first;
      if(!MaterialsSet[Material]){
	MaterialsSet[Material]=true;
	MaterialTables[Material]=new G4MaterialPropertiesTable;
      }
    }
    
    // For each property vector, convert to an array of g4doubles and 
    // feed to materials table Lots of firsts and seconds!  See annotation 
    // in MaterialPropertyLoader.h to follow what each element is
    
    for(std::map<std::string,std::map<std::string,std::map<double,double> > >::const_iterator i=propertyList.begin(); i!=propertyList.end(); i++){
      std::string Material=i->first;
      for(std::map<std::string,std::map<double,double> >::const_iterator j = i->second.begin(); j!=i->second.end(); j++){
	std::string Property=j->first;
	std::vector<G4double> g4MomentumVector;
	std::vector<G4double> g4PropertyVector;
	
	for(std::map<double,double>::const_iterator k=j->second.begin(); k!=j->second.end(); k++){
	  g4MomentumVector.push_back(k->first);
	  g4PropertyVector.push_back(k->second);
	}
	int NoOfElements=g4MomentumVector.size();
	MaterialTables[Material]->AddProperty(Property.c_str(),&g4MomentumVector[0], &g4PropertyVector[0],NoOfElements); 
	mf::LogInfo("MaterialPropertyLoader") << "Added property "
					      <<Property
					      <<" to material table " 
					      << Material;
      }
    }
    
    //Add each const property element
    for(std::map<std::string,std::map<std::string,double > >::const_iterator i = constPropertyList.begin(); i!=constPropertyList.end(); i++){
      std::string Material=i->first;
      for(std::map<std::string,double>::const_iterator j = i->second.begin(); j!=i->second.end(); j++){
	std::string Property=j->first;
	G4double PropertyValue=j->second;
	MaterialTables[Material]->AddConstProperty(Property.c_str(), PropertyValue); 
	mf::LogInfo("MaterialPropertyLoader") << "Added const property "
					      <<Property
					      <<" to material table " 
					      << Material;
      }
    }
    
    //Loop through geometry elements and apply relevant material table where materials match
    for ( G4LogicalVolumeStore::iterator i = lvs->begin(); i != lvs->end(); ++i ){
      G4LogicalVolume* volume = (*i);
      G4Material* TheMaterial = volume->GetMaterial();
      std::string Material = TheMaterial->GetName();
      for(std::map<std::string,G4MaterialPropertiesTable*>::const_iterator j=MaterialTables.begin(); j!=MaterialTables.end(); j++){
	if(Material==j->first){
	  TheMaterial->SetMaterialPropertiesTable(j->second);
	  //Birks Constant, for some reason, must be set separately
	  if(birksConstants[Material]!=0)
	    TheMaterial->GetIonisation()->SetBirksConstant(birksConstants[Material]);
	  volume->SetMaterial(TheMaterial);
	}
      }
    }
  }
  
}
