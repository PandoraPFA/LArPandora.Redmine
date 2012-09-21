////////////////////////////////////////////////////////////////////////
//
//   MergeSimMINOS Class
//
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <dirent.h>
#include <time.h>
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>


#include "T962/MergeData/MergeSimMINOS.h"


namespace merge{
   static int event=1;
   //-------------------------------------------------
   MergeSimMINOS::MergeSimMINOS(fhicl::ParameterSet const& pset) : 
      fG4ModuleLabel(pset.get< std::string >("GeantModuleLabel"))
   {
      produces< std::vector<t962::MINOS> > ();
   }

   //-------------------------------------------------
   MergeSimMINOS::~MergeSimMINOS()
   {
   
   }

   //-------------------------------------------------
   void MergeSimMINOS::beginJob()
   {
      //get access to the TFile service
      art::ServiceHandle<art::TFileService> tfs;
   }

   //-------------------------------------------------
   void MergeSimMINOS::produce(art::Event& evt)
   {
      evt.getByLabel(fG4ModuleLabel, parHandle);
//       evt.getByLabel(fbeam_modulelabel,fbeam);
	  
      std::unique_ptr<std::vector<t962::MINOS> > MINOS_coll(new std::vector<t962::MINOS> );
      std::vector<t962::MINOS> vec_minos;

      if(MergeMINOS(vec_minos)){
 
         //std::cout << "No of MINOS objects saved is " << vec_minos.size() << std::endl;
         for(unsigned int i=0;i<vec_minos.size();i++)
         {
            MINOS_coll->push_back(vec_minos[i]);
         }
         evt.put(std::move(MINOS_coll));
      }
      else{
         std::cout <<" : No MINOS match." << std::endl;
      }

      vec_minos.clear();
      return;
   }
 
   //-------------------------------------------------
   bool MergeSimMINOS::MergeMINOS(std::vector<t962::MINOS> & vec_minos)
   {
      t962::MINOS minos;//new minos object;
  
      int flag=0;
      
      //loop through MINOS root files	  
      //std::string path = "/argoneut/data/outstage/spitz7/rashidmc/";
      std::string path(gSystem->Getenv("DIR"));//PWD
      
      path=( path + "/");
      

      DIR *pDIR;
      struct dirent *entry;
      if( (pDIR=opendir(path.c_str())) != NULL )
      {
     
         while((entry = readdir(pDIR)) != NULL)
         {;
            if( strcmp(entry->d_name, ".")==0 || strcmp(entry->d_name, "..")==00) continue;
            

            std::string filename(entry->d_name);    
            if((int)filename.find("fhc.root")==-1)
            continue;

            //grab initial/final timestamp info. from input MINOS file

            std::string file = (path + entry->d_name);
 
            TFile *f = new TFile(file.c_str());
            TTree *minitree = (TTree*)f->Get("minitree");
            minitree->SetBranchAddress("run",&minos.frun);  
            minitree->SetBranchAddress("subRun",&minos.fsubRun);
            minitree->SetBranchAddress("snarl",&minos.fsnarl);
            minitree->SetBranchAddress("utc",&minos.futc);
            minitree->SetBranchAddress("day",&minos.fday);
            minitree->SetBranchAddress("trkIndex",&minos.ftrkIndex);
            minitree->SetBranchAddress("trkChi2",&minos.ftrkChi2);
            minitree->SetBranchAddress("utc1",&minos.futc1);
            minitree->SetBranchAddress("nearns",&minos.fnearns);
            minitree->SetBranchAddress("nearsec",&minos.fnearsec);		  
            minitree->SetBranchAddress("trkE",&minos.ftrkE);
            minitree->SetBranchAddress("shwE",&minos.fshwE);
            minitree->SetBranchAddress("crateT0",&minos.fcrateT0);
            minitree->SetBranchAddress("tmframe",&minos.ftmframe);
            minitree->SetBranchAddress("year",&minos.fyear);		  
            minitree->SetBranchAddress("offset",&minos.foffset);
            minitree->SetBranchAddress("dtnear",&minos.fdtnear);		
            minitree->SetBranchAddress("trkErange",&minos.ftrkErange);
            minitree->SetBranchAddress("sgate53",&minos.fsgate53);
            minitree->SetBranchAddress("trkqp",&minos.ftrkqp);
            minitree->SetBranchAddress("trkeqp",&minos.ftrkeqp);
            minitree->SetBranchAddress("trkVtxX",&minos.ftrkVtxX);
            minitree->SetBranchAddress("trkVtxY",&minos.ftrkVtxY);
            minitree->SetBranchAddress("trkVtxZ",&minos.ftrkVtxZ);
            minitree->SetBranchAddress("trkVtxeX",&minos.ftrkVtxeX);
            minitree->SetBranchAddress("trkVtxeY",&minos.ftrkVtxeY);
            minitree->SetBranchAddress("charge",&minos.fcharge); 		  
            minitree->SetBranchAddress("trkmom",&minos.ftrkmom);	
            minitree->SetBranchAddress("trkVtxT",&minos.ftrkVtxT);	
            minitree->SetBranchAddress("trkTimeT0",&minos.ftrkTimeT0);  
            minitree->SetBranchAddress("trkdcosx",&minos.ftrkdcosx);
            minitree->SetBranchAddress("trkdcosy",&minos.ftrkdcosy);
            minitree->SetBranchAddress("trkdcosz",&minos.ftrkdcosz);
            minitree->SetBranchAddress("month",&minos.fmonth);
            minitree->SetBranchAddress("trtgtd",&minos.ftrtgtd);
            minitree->SetBranchAddress("tortgt",&minos.ftortgt);
            minitree->SetBranchAddress("tor101",&minos.ftor101);
            minitree->SetBranchAddress("tr101d",&minos.ftr101d);
            minitree->SetBranchAddress("goodspill",&minos.fgoodspill);
            minitree->SetBranchAddress("trkContained",&minos.ftrkcontained);
            minitree->SetBranchAddress("goodbeam",&minos.fgoodbeam);
            minitree->SetBranchAddress("ntrkstp",&minos.fntrkstp);
            // uncomment next 10 lines for MC 
	    // This is _only_ to be used for MC.
            minitree->SetBranchAddress("mcIndex",&minos.fmcIndex);
            minitree->SetBranchAddress("mcPDG",&minos.fmcPDG);
            minitree->SetBranchAddress("mcPx",&minos.fmcPx);
            minitree->SetBranchAddress("mcPy",&minos.fmcPy);
            minitree->SetBranchAddress("mcPz",&minos.fmcPz);
            minitree->SetBranchAddress("mcEne",&minos.fmcEne);
            minitree->SetBranchAddress("mcMass",&minos.fmcMass);
            minitree->SetBranchAddress("mcVtxX",&minos.fmcVtxX);
            minitree->SetBranchAddress("mcVtxY",&minos.fmcVtxY);
            minitree->SetBranchAddress("mcVtxZ",&minos.fmcVtxZ);
          
            double fmcPx, fmcPy, fmcPz;
             minitree->SetBranchAddress("mcPx",&fmcPx);
             minitree->SetBranchAddress("mcPy",&fmcPy);
             minitree->SetBranchAddress("mcPz",&fmcPz);
             minos.SetmcPx(fmcPx);
             minos.SetmcPy(fmcPy);
	     minos.SetmcPz(fmcPz);
            
            float trkstpX[1000];//don't like to hard-code the array size, but will suffice for now
            float trkstpY[1000];
            float trkstpZ[1000];
            float trkstpU[1000];
            float trkstpV[1000];
            minitree->SetBranchAddress("trkstpX",trkstpX);
            minitree->SetBranchAddress("trkstpY",trkstpY);
            minitree->SetBranchAddress("trkstpZ",trkstpZ);
            minitree->SetBranchAddress("trkstpU",trkstpU);
            minitree->SetBranchAddress("trkstpV",trkstpV);

            //------------------------------------------
            Long64_t nentries = minitree->GetEntries();
            Long64_t nbytes = 0;
		 		 
            for(int i=0;i<nentries;i++){

               nbytes+=minitree->GetEntry(i);
               //create correctly-sized vectors from arrays.
               std::vector<float> ftrkstpX(trkstpX,trkstpX+minos.fntrkstp);
               std::vector<float> ftrkstpY(trkstpY,trkstpY+minos.fntrkstp);
               std::vector<float> ftrkstpZ(trkstpZ,trkstpZ+minos.fntrkstp);
               std::vector<float> ftrkstpU(trkstpU,trkstpU+minos.fntrkstp);
               std::vector<float> ftrkstpV(trkstpV,trkstpV+minos.fntrkstp);
               
               //store in working MINOS object
               minos.ftrkstpX = ftrkstpX;
               minos.ftrkstpY = ftrkstpY;
               minos.ftrkstpZ = ftrkstpZ;
               minos.ftrkstpU = ftrkstpU;
               minos.ftrkstpV = ftrkstpV;

       if(minos.fsnarl==event)
       { 
          vec_minos.push_back(minos); 		  
          flag=1;          
       }
       else if(flag==1) break;

            }
  
            minitree->Delete(); 
            f->Close();
            f->Delete();
            
         }//while
         closedir(pDIR);
		  
      }
	  event++;
      if(vec_minos.size() > 0) return true;
      else return false;
    
   }


}//namespace merge
