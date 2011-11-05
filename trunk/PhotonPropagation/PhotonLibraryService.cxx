#ifndef __CINT__

#include <sys/stat.h> 
#include "TFile.h"
#include "PhotonPropagation/PhotonLibrary.h"
#include "PhotonPropagation/PhotonLibraryService.h"
#include "Simulation/PhotonLibraryParameters.h"
#include "Geometry/geo.h"


namespace phot
{

  PhotonLibraryService::PhotonLibraryService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  {
    fFileName =  pset.get<std::string>("PhotonLibraryFile");
  }



  void PhotonLibraryService::InitializeBuildJob(sim::PhotonVoxelDef VoxelDef, bool ExtendExisting)
  {
    if(fBuildLibrary)
      if(!fExtendExisting)
	{
	  CreateLibraryFile(fFileName, VoxelDef);
	} 
    OpenLibraryFile(fFileName);
    CheckVoxelDef(VoxelDef);
    std::cout<<"PhotonLibraryService : file "<< fFileName<<" open for library building" <<std::cout;
  }




  void PhotonLibraryService::InitializeReadJob()
  {
    OpenLibraryFile(fFileName);
    std::cout<<"PhotonLibraryService : file "<< fFileName<<" open for reading" <<std::cout;
  }



  PhotonLibraryService::~PhotonLibraryService()
  {
    if(fTheLibraryFile) fTheLibraryFile->Close();
    std::cout<<"PhotonLibraryService : closing file " << fFileName<<std::endl;
  }

 


  void PhotonLibraryService::CreateLibraryFile(std::string FileName, sim::PhotonVoxelDef TheVoxelDef)
  {

    std::cout<<"PhotonLibraryService creating a new library file... " << FileName.c_str()<<std::endl;

    fVoxelDef = TheVoxelDef;
    fTheLibraryFile = new TFile(FileName.c_str(), "recreate");
    if(!fTheLibraryFile)
      {
	std::cout<<"PhotonLibraryService : Error opening library file." << std::endl;
	throw("");
      }
    else
      {
	std::cout<<"PhotonLibraryService : Successfully created new library file"<<std::endl;
      }

    // Create an empty photon library using the loaded voxel definition
    fTheLibrary = new PhotonLibrary(fVoxelDef);
    std::cout<<"Debug : making object"<<std::endl;
    
    fTheLibrary->SaveInFile(fTheLibraryFile);
    //fTheLibrary->Write("StoredPhotonLibrary");
    std::cout<<"Debug : object written" <<std::endl;

    // close the TFile so we can open it again, the standard way
    fTheLibraryFile->Close();
    std::cout<<"Debug : File closed" <<std::endl;
  }


  void PhotonLibraryService::OpenLibraryFile(std::string FileName)
  {
    std::cout<<"PhotonLibraryService attemping to open library file... " << FileName.c_str()<<std::endl;
    // Check our file exists
    struct stat stFileInfo;
    int intStat;
    intStat = stat(FileName.c_str(),&stFileInfo);
    if(intStat != 0) 
      {
	std::cout<<"PhotonLibraryService : Error, Library file does not exist! If you want to create a new one, set the CreateNewLibrary parameter to True" << std::endl;
	throw("");
      }
    else 
      {
	std::cout<<"PhotonLibraryService : File exists"<<std::endl;
      }
    fTheLibraryFile=TFile::Open(FileName.c_str());
    if(!fTheLibraryFile)
      {
	std::cout<<"PhotonLibraryService : Error opening library file." << std::endl;
	throw("");
      }

    // Extract library pointer from file
    fTheLibrary = new PhotonLibrary;
    fTheLibrary->LoadFromFile(fTheLibraryFile);
    std::cout<<"Library pointer extracted "<<std::endl;
    fVoxelDef = fTheLibrary->GetVoxelDef();
    std::cout<<"Voxel def extracted "<<std::endl;
  }


  void PhotonLibraryService::CheckVoxelDef(sim::PhotonVoxelDef TheDef)
  {
    // Check voxel def compatability
    if(TheDef != fVoxelDef)
      {
	std::cout<<"Photon voxel definition incompatible between stored library and config parameters. To use modify library, please modify job config accordingly"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Library voxel structure:                "<<std::endl<<std::endl;
	std::cout<<"           LIBRARY        THIS JOB      "<<std::endl;
	std::cout<<" XMin     "<< fVoxelDef.GetRegionLowerCorner()[0]  << "         " << TheDef.GetRegionLowerCorner()[0]  << std::endl;
	std::cout<<" XMax     "<< fVoxelDef.GetRegionUpperCorner()[0]  << "         " << TheDef.GetRegionUpperCorner()[0]  << std::endl;
	std::cout<<" YMin     "<< fVoxelDef.GetRegionLowerCorner()[1]  << "         " << TheDef.GetRegionLowerCorner()[1]  << std::endl;
	std::cout<<" YMax     "<< fVoxelDef.GetRegionUpperCorner()[1]  << "         " << TheDef.GetRegionUpperCorner()[1]  << std::endl;
	std::cout<<" ZMin     "<< fVoxelDef.GetRegionLowerCorner()[2]  << "         " << TheDef.GetRegionLowerCorner()[2]  << std::endl;
	std::cout<<" ZMax     "<< fVoxelDef.GetRegionUpperCorner()[2]  << "         " << TheDef.GetRegionUpperCorner()[2]  << std::endl;
	std::cout<<std::endl;
	std::cout<<" XSteps   "<< fVoxelDef.GetSteps()[0]  << "         " << TheDef.GetSteps()[0]  << std::endl;
	std::cout<<" YSteps   "<< fVoxelDef.GetSteps()[1]  << "         " << TheDef.GetSteps()[1]  << std::endl;
	std::cout<<" ZSteps   "<< fVoxelDef.GetSteps()[2]  << "         " << TheDef.GetSteps()[2]  << std::endl;

	throw("");
      }
  }
 

  PhotonLibrary * PhotonLibraryService::GetLibrary()
  {
    if(fTheLibrary) 
      return fTheLibrary;
    else
      {
	std::cout<<"PhotonLibraryService Error : No library is loaded"<<std::endl;
	throw("");
      }
  }

}

#endif
