////////////////////////////////////////////////////////////////////////
//
//   MergeDataMINOS Class
//
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Common/OrphanHandle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

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


#include "T962/MergeData/MergeDataMINOS.h"

namespace merge{

   //-------------------------------------------------
   MergeDataMINOS::MergeDataMINOS(fhicl::ParameterSet const& pset) : 
      fdaq_modulelabel(pset.get< std::string >("daq")),
      fbeam_modulelabel(pset.get< std::string >("beam"))
   {
      produces< std::vector<t962::MINOS> > ();
   }

   //-------------------------------------------------
   MergeDataMINOS::~MergeDataMINOS()
   {
   
   }

   //-------------------------------------------------
   void MergeDataMINOS::beginJob()
   {
      //get access to the TFile service
      art::ServiceHandle<art::TFileService> tfs;

      fproblemevent2d=tfs->make<TH2D>("fproblemevent2d","POT vs. Event number with no match", 400,0 ,40000,500,-1,45); 
      fPOTdiff_matched=tfs->make<TH1D>("fPOTdiff_matched","POT difference of match", 800,-40,40);
      fMINOSrun_event=tfs->make<TH2D>("MINOSrun_event","MINOS run # vs ArgoNeuT event # for match", 400,0 ,40000,500,15000 ,20000);
      futc1_tms_diff = tfs->make<TH1D>("futc1_tms_diff","fabs(utc1+500-(tms))", 200,0 ,2000);

   }

   //-------------------------------------------------
   void MergeDataMINOS::produce(art::Event& evt)
   {
      art::Handle< std::vector<raw::DAQHeader> > daqHandle;
      evt.getByLabel(fdaq_modulelabel,daqHandle);
      fdaq = art::Ptr<raw::DAQHeader>(daqHandle, daqHandle->size()-1);

      art::Handle< std::vector<raw::BeamInfo> > beamHandle;
      evt.getByLabel(fbeam_modulelabel,beamHandle);
      fbeam = art::Ptr<raw::BeamInfo>(beamHandle, beamHandle->size()-1);
	  
      std::auto_ptr<std::vector<t962::MINOS> > MINOS_coll(new std::vector<t962::MINOS> );
      std::vector<t962::MINOS> vec_minos;

      if(MergeMINOS(vec_minos)){
         //std::cout << "No of MINOS objects saved is " << vec_minos.size() << std::endl;
         for(int i=0;i<vec_minos.size();i++)
         {
            MINOS_coll->push_back(vec_minos[i]);
         }
         evt.put(MINOS_coll);
      }
      else{
         std::cout << "Run " << fdaq->GetRun() << " Event " << fdaq->GetEvent() << " : No MINOS time match." << std::endl;
      }

      vec_minos.clear();
 
      return;
   }
 
   //-------------------------------------------------
   bool MergeDataMINOS::MergeMINOS(std::vector<t962::MINOS> & vec_minos)
   {
      t962::MINOS minos;//new minos object;
  
      int flag=0;
	    
      double tms = (double)fbeam->get_t_ms();//millisecond timing information from BeamInfo object

      //loop through MINOS root files	  
      std::string path = "/argoneut/app/users/rmehdi/fhc/";
      DIR *pDIR;
      struct dirent *entry;
      if( (pDIR=opendir(path.c_str())) != NULL )
      {
         while((entry = readdir(pDIR)) != NULL)
         {
            if( strcmp(entry->d_name, ".")==0 || strcmp(entry->d_name, "..")==00) continue;
		 
            //grab initial/final timestamp info. from input MINOS file
            int firsttimestart = 0;
            int firsttimeend = 0;
            parse_filename(entry->d_name,firsttimestart,firsttimeend);
            //check that current ArgoNeuT event is within time window of current MINOS file
            if(!(firsttimestart<tms/1000.0  && tms/1000.<firsttimeend)) continue;
		 
            std::string file = (path + entry->d_name);
            //std::cout <<"***********FILE BEING SCANNED IS*********** "<< file << "\n";
 
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
		  		    
               //*************************************************************
               //Matching condition based on time info alone:
               double diff = fabs(minos.futc1 + 500 - tms);
          
               if(diff<1001){
                  if(minos.ftrkIndex==0){
                     // std::cout << "minos.futc1 = " << std::setprecision(12) << minos.futc1 
                     //           << " tms = " << std::setprecision(12) << tms << std::endl;
                     // std::cout << "diff = " << diff << std::endl;
                     fPOTdiff_matched->Fill(fbeam->get_tor101() - minos.ftor101);
                     fMINOSrun_event->Fill(fdaq->GetEvent(),minos.frun);
                     futc1_tms_diff->Fill(diff);
                  }
                  vec_minos.push_back(minos); 		  
                  flag=1;
               }
               else if(flag==1) break;//already past matching tracks, so stop looping over TTree

            }//loop over Minos root file		  

            minitree->Delete();
            f->Close();
            f->Delete();

            if(vec_minos.size()==0) fproblemevent2d->Fill(fdaq->GetEvent(),fbeam->get_tor101());
            else break;//already found our match, so why keep looping over files in directory?
            
         }//while
         closedir(pDIR);
		  
      }// if( pDIR=opendir("...
	 
      if(vec_minos.size() > 0) return true;
      else return false;
  
  
   }

   ////////////////////////////////////////////////////////////////////

   void  MergeDataMINOS::parse_filename(std::string filename, int & minos_t1, int & minos_t2)
   {
      std::string str("_");//delimiter
      size_t l1,l2,l3;
      
      l1 = filename.find(str);
      l2 = filename.find(str,l1+1);
      l3 = filename.find(str,l2+1);

      if(l2==std::string::npos || l3==std::string::npos){
         std::cout << "MergeData_MINOS::parse_filename Error: Couldn't identify start/end times for MINOS file " 
                   << filename << std::endl;
         return;
      }

      minos_t1 = atoi((filename.substr(l2+1,10)).c_str());//convert substring to int
      minos_t2 = atoi((filename.substr(l3+1,10)).c_str());//convert substring to int
      
      return;

   }

}//namespace merge
