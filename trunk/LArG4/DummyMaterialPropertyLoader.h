// Sample material property loader - just hack in some random LAr Properties

#include "MaterialPropertyLoader.h"

namespace larg4 {

  class DummyMaterialPropertyLoader : public MaterialPropertyLoader
  {
  public:
    DummyMaterialPropertyLoader();
    ~DummyMaterialPropertyLoader(); 
  };
  
  
  DummyMaterialPropertyLoader::DummyMaterialPropertyLoader(){
    std::cout<<"MaterialPropertyLoader : Loading values" <<std::endl;
    
    std::vector<double> Energies;
    std::map<double,double> Scint;
    std::map<double,double> Rind;
    std::map<double,double> Absl;
    std::map<double,double> Rayleigh;
    std::map<double,double> Reflectance;
    std::map<double,double> DiffuseFraction;
    std::map<double,double> TPBAbsorptionLength;
    std::map<double,double> TPBEmissionSpectrum;
    
    
    Scint[9.5*eV]=0.5;
    Scint[9.7*eV]=1.0;
    Scint[9.9*eV]=0.5;
    this->SetMaterialProperty("LAr","FASTCOMPONENT",Scint);
    this->SetMaterialProperty("LAr","SLOWCOMPONENT",Scint);
    
    Rind[9.5*eV]=1.38;
    Rind[9.7*eV]=1.38;
    Rind[9.9*eV]=1.38;
    this->SetMaterialProperty("LAr","RINDEX",Rind);
    
    Absl[9.5*eV]=2000.*m; 
    Absl[9.7*eV]=2000.*m; 
    Absl[9.9*eV]=2000.*m; 
    this->SetMaterialProperty("LAr","ABSLENGTH", Absl);
    
    Absl[9.5*eV]=0.01*mm; 
    Absl[9.7*eV]=0.01*mm; 
    Absl[9.9*eV]=0.01*mm; 
    this->SetMaterialProperty("Glass","ABSLENGTH", Absl);
    
    Rayleigh[9.5*eV]=90*cm;
    Rayleigh[9.7*eV]=90*cm;
    Rayleigh[9.9*eV]=90*cm;
    this->SetMaterialProperty("LAr","RAYLEIGH", Rayleigh);
    
    Reflectance[9.5*eV]=0.25;
    Reflectance[9.7*eV]=0.25;
    Reflectance[9.9*eV]=0.25;
    this->SetMaterialProperty("LAr","REFLECTANCE_STEEL_STAINLESS_Fe7Cr2Ni", Reflectance);     
    
    DiffuseFraction[9.5*eV]=0.5;
    DiffuseFraction[9.7*eV]=0.5;
    DiffuseFraction[9.9*eV]=0.5;
    this->SetMaterialProperty("LAr","DIFFUSE_REFLECTANCE_FRACTION_STEEL_STAINLESS_Fe7Cr2Ni",DiffuseFraction);  
    
    TPBAbsorptionLength[9.5*eV]=1*nm;
    TPBAbsorptionLength[9.7*eV]=1*nm;
    TPBAbsorptionLength[9.9*eV]=1*nm;
    this->SetMaterialProperty("TPB","WLSABSLENGTH",TPBAbsorptionLength);
    
    TPBEmissionSpectrum[0.30*eV]=0; 
    TPBEmissionSpectrum[0.28*eV]=1;
    TPBEmissionSpectrum[0.30*eV]=0;   
    this->SetMaterialProperty("TPB","WLSCOMPONENT",TPBEmissionSpectrum);
    
    
    this->SetMaterialConstProperty("TPB","WLSTIMECONSTANT",0.5*ns);
    this->SetMaterialConstProperty("TPB","WLSMEANNUMBERPHOTONS",1.3);
    this->SetMaterialConstProperty("LAr","SCINTILLATIONYIELD", 24000.0*0.3/MeV); //for argon (BUT field dep)
    this->SetMaterialConstProperty("LAr","RESOLUTIONSCALE",0.005);
    this->SetMaterialConstProperty("LAr","FASTTIMECONSTANT",6.*ns);     //for argon
    this->SetMaterialConstProperty("LAr","SLOWTIMECONSTANT",1590.*ns);  //for argon
    this->SetMaterialConstProperty("LAr","YIELDRATIO",0.3);             //more work req
    
    SetBirksConstant("LAr",0.00322*cm/MeV);
  }
  
  DummyMaterialPropertyLoader::~DummyMaterialPropertyLoader(){
  }
}
