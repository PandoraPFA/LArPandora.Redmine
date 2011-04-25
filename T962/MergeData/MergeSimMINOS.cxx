////////////////////////////////////////////////////////////////////////
//
//   MergeSimMINOS Class
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


#include "T962/MergeData/MergeSimMINOS.h"


namespace merge{

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

      fproblemevent2d=tfs->make<TH2D>("fproblemevent2d","POT vs. Event number with no match", 400,0 ,40000,500,-1,45); 

   }

   //-------------------------------------------------
   void MergeSimMINOS::produce(art::Event& evt)
   {
      evt.getByLabel(fG4ModuleLabel, parHandle);
//       evt.getByLabel(fbeam_modulelabel,fbeam);
	  
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
	    
//       double tms = (double)fbeam->get_t_ms();//millisecond timing information from BeamInfo object

      //loop through MINOS root files	  
      std::string path = "/argoneut/app/users/rmehdi/mcfhc/";
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
            
            
            double fmcPx, fmcPy, fmcPz;
            minitree->SetBranchAddress("mcPx",&fmcPx);
            minitree->SetBranchAddress("mcPy",&fmcPy);
            minitree->SetBranchAddress("mcPz",&fmcPz);
            minos.SetmcPx(fmcPx);
            minos.SetmcPy(fmcPy);
		    minos.SetmcPz(fmcPz); 
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
             // double diff = fabs(minos.futc1 + 500 - tms);
          
               //if(diff<1001){
               
               //if mcvertices are equal
               
 //////////////////////////////////////////////////////////////////////////////////////////                                             
//////////////////////////////////////////////////////////////////////////////////////////  

//        std::cout<<std::setprecision(9)<<first<<"  "
// 		<<run<<"  "<<event<<"  "<<tms<<" "
// 		<<tor101<<"  "<< tortgt <<"  "<< trtgtd <<std::endl;
   
       //if((fdaq->GetRun()==run)&&(fdaq->GetEvent()==event ) ){
 //       std::cout<<"MINOS "<<minos.fmcPx<<" "<<minos.fmcPy<<" "<<minos.fmcPz<<" "<<std::endl;
//        std::cout<<index<<" "<<stat<<" "<<par<<" "<<chil<<" "<<dummy<<" "<<px<<" "<<py<<" "<<pz<<std::endl;
//        

     // get the particles from the event handle


    art::PtrVector<sim::Particle> pvec;
    for(unsigned int i = 0; i < parHandle->size(); ++i){
      art::Ptr<sim::Particle> p(parHandle, i);      
      pvec.push_back(p);
    }    

    float minosenter_px=0.;
    float minosenter_py=0.;
    float minosenter_pz=0.;
    for(unsigned int i = 0; i < pvec.size(); ++i){
    
    if(pvec[i]->Process()!="primary")
    continue;
    
    simb::MCTrajectory trajectory = pvec[i]->Trajectory();
    TLorentzVector initposition = trajectory.Position(0);
    TLorentzVector initmomentum = trajectory.Momentum(0);
    int numberofpoints= pvec[i]->NumberTrajectoryPoints();   

    for(int j=1; j<numberofpoints; j++)
    {
      TLorentzVector prevposition = trajectory.Position(j-1);
      TLorentzVector position = trajectory.Position(j);
      TLorentzVector prevmomentum = trajectory.Momentum(j-1);
      TLorentzVector momentum = trajectory.Momentum(j);


    if(prevposition.Z()<154.22&&position.Z()>154.22&&abs(pvec[i]->PdgCode()!=14)&&abs(pvec[i]->PdgCode()!=12))
           {
       minosenter_px=prevmomentum.Px();
       minosenter_py=prevmomentum.Py();
       minosenter_pz=prevmomentum.Pz();
           }
    }   
  }

 std::cout<<"MINOS "<<fmcPx<<" "<<fmcPy<<" "<<fmcPz<<" "<<std::endl;
 std::cout<<"Event record "<<minosenter_px<<" "<<minosenter_py<<" "<<minosenter_pz<<std::endl;

       if(fmcPx==minosenter_px&&fmcPy==minosenter_py&&fmcPz==minosenter_pz
       &&minosenter_px&&minosenter_py&&minosenter_pz
       ){
  //             std::cout<<"MINOS "<<minos.fmcPx<<" "<<minos.fmcPy<<" "<<minos.fmcPz<<" "<<std::endl;
//        std::cout<<index<<" "<<stat<<" "<<par<<" "<<chil<<" "<<dummy<<" "<<px<<" "<<py<<" "<<pz<<std::endl;

//           beam.SetTOR101(tor101);
//           beam.SetTORTGT(tortgt);
//           beam.SetTRTGTD(trtgtd);
//           beam.SetT_MS(tms);
          vec_minos.push_back(minos); 		  
          flag=1;          
       }
        else if(flag==1) break;
            
//////////////////////////////////////////////////////////////////////////////////////////                                             
//////////////////////////////////////////////////////////////////////////////////////////                 
//                if(1)
//                {
//                   if(minos.ftrkIndex==0){
//                      // std::cout << "minos.futc1 = " << std::setprecision(12) << minos.futc1 
//                      //           << " tms = " << std::setprecision(12) << tms << std::endl;
//                      // std::cout << "diff = " << diff << std::endl;
// //                      fPOTdiff_matched->Fill(fbeam->get_tor101() - minos.ftor101);
// //                      fMINOSrun_event->Fill(fdaq->GetEvent(),minos.frun);
// //                      futc1_tms_diff->Fill(diff);
//                   }
//                   vec_minos.push_back(minos); 		  
//                   flag=1;
//                }
//                else if(flag==1) break;//already past matching tracks, so stop looping over TTree

            }
            //loop over Minos root file		  

            minitree->Delete();
            f->Close();
            f->Delete();

            if(vec_minos.size()==0) fproblemevent2d->Fill(1.,1.);
            else break;//already found our match, so why keep looping over files in directory?
            
         }//while
         closedir(pDIR);
		  
      }// if( pDIR=opendir("...
	 
      if(vec_minos.size() > 0) return true;
      else return false;
  
  
   }


}//namespace merge
