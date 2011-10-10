#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TParticle.h"
#include <iostream>
#include <iomanip>
#include <string>


#include "StandardNtuple/NtpStRecord.h"
#include "BeamDataNtuple/NtpBDLiteRecord.h"

#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "CandNtupleSR/NtpSRTimeStatus.h"
#include "CandNtupleSR/NtpSRDate.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREventSummary.h"

#include "RunQuality/RunQualityFinder.h"
#include "DcsUser/CoilTools.h"
#include "SpillTiming/SpillServerMonFinder.h"
#include "DataUtil/DataQualDB.h"

#include "TruthHelperNtuple/NtpTHTrack.h"
#include "Validity/VldContext.h"
#include "BeamDataUtil/BeamMonSpill.h"
#include "BeamDataUtil/BMSpillAna.h"
#include "BeamDataUtil/BDSpillAccessor.h"

#include "MCNtuple/NtpMCRecord.h"
#include "MCNtuple/NtpMCTruth.h"
#include "MCNtuple/NtpMCStdHep.h"


void microDst(TString inputfile, TString outputfile)
{	

  TFile fin(inputfile);

  TTree* tin=(TTree*)fin.Get("NtpSt");
  NtpStRecord* nsr=0;
  tin->SetBranchAddress("NtpStRecord", &nsr);

  /* 
  TTree* ntpBD=(TTree*)fin.Get("NtpBDLite");
  NtpBDLiteRecord* nbd=0;
  ntpBD->SetBranchAddress("NtpBDLiteRecord", &nbd);
  tin->AddFriend(ntpBD);


  int nspills = ntpBD->GetEntries();

  const int rec = ntpBD->GetEntry(0);
  if(rec<=0){
    cout << "*** Failed to get beam record " <<endl;
  }
  */
  
  TFile fout(outputfile, "RECREATE");
  TTree tout("minitree", "This is just the title to initialize");
 
  // Set up the variables you want here
  Int_t run, subRun, snarl, ntrkstp, goodSpillCnt, badSpillCnt, trkContained;
  Float_t trkE, shwE, trkErange, trkqp, trkeqp, trkIndex;
  Float_t trkVtxX, trkVtxY, trkVtxZ, trkdcosx, trkdcosy, trkdcosz;
  Float_t trkTimeT0, trkVtxT, trkChi2, trtgtd, tortgt, tor101, tr101d;
  Float_t trkstpX[100000],trkstpY[100000],trkstpZ[100000];
  Float_t trkstpU[100000],trkstpV[100000];
  Float_t charge, trkmom, dtnear, nearsec, mcIndex;
  Double_t crateT0, tmframe, sgate53, offset, trkVtxeX, trkVtxeY, utc, utc1;
  Double_t day, month, year, nearns;
  Double_t mcPDG, mcPx, mcPy, mcPz,  mcEne,  mcMass, mcVtxX, mcVtxY, mcVtxZ;
  Int_t ntracks, goodspill, goodbeam; 

#define BRANCH(branchname) tout.Branch(#branchname, &branchname)
  BRANCH(run); BRANCH(subRun); BRANCH(snarl); BRANCH(utc); BRANCH(day); 
  BRANCH(trkIndex); BRANCH(trkChi2);  BRANCH(utc1); BRANCH(nearns); BRANCH(nearsec); 
  BRANCH(trkE); BRANCH(shwE); BRANCH(crateT0); BRANCH(tmframe);BRANCH(year);
  BRANCH(trkErange);BRANCH(sgate53);
  BRANCH(trkqp);BRANCH(trkeqp);BRANCH(trkVtxX);BRANCH(trkVtxY);BRANCH(trkVtxZ);
  BRANCH(trkVtxeX);BRANCH(trkVtxeY); BRANCH(trkChi2); BRANCH(offset); BRANCH(dtnear);
  BRANCH(charge); BRANCH(trkmom); BRANCH(trkTimeT0); BRANCH(trkVtxT); 
  BRANCH(trkdcosx); BRANCH(trkdcosy); BRANCH(trkdcosz); BRANCH(month); 
  BRANCH(trtgtd); BRANCH(tortgt); BRANCH(tor101); BRANCH(tr101d); BRANCH(goodspill);
  BRANCH(ntrkstp); BRANCH(trkContained); BRANCH(goodbeam); BRANCH(ntracks);
  //truth
  BRANCH(mcIndex);
  BRANCH(mcPDG);BRANCH(mcPx);BRANCH(mcPy);BRANCH(mcPz); BRANCH(mcEne);
  BRANCH(mcMass);BRANCH(mcVtxX);BRANCH(mcVtxY);BRANCH(mcVtxZ);

  tout.Branch("trkstpX", &trkstpX, "trkstpX[ntrkstp]/F");
  tout.Branch("trkstpY", &trkstpY, "trkstpY[ntrkstp]/F");
  tout.Branch("trkstpZ", &trkstpZ, "trkstpZ[ntrkstp]/F");
  tout.Branch("trkstpU", &trkstpU, "trkstpU[ntrkstp]/F");
  tout.Branch("trkstpV", &trkstpV, "trkstpV[ntrkstp]/F");

  cout << tin->GetEntries() << " entries" << endl;


  tor101=0; tr101d=0; trtgtd=0; tortgt=0;
  utc=0; day=0; month=0; year=0; ntrkstp=0; trkContained=0;
  tmframe=0; crateT0=0; sgate53=0; nearsec=0; offset=0;
  nearns=0; dtnear=0; utc1=0; goodspill=0; goodSpillCnt=0; badSpillCnt=0;
 



  // Instantiate a BDSpillAccessor object, necessary to access the
  // BeamMonSpill database.


  //  BDSpillAccessor& sa = BDSpillAccessor::Get();

  Int_t index=0;
  Int_t snarlCount=0;
  Int_t utcS=0;
  Int_t utcL=0;
  goodbeam=0;

  // Loop over the snarls in the input file
  for(int iEntry=0; iEntry<tin->GetEntries(); ++iEntry){
    
    snarlCount++;
    tin->GetEntry(iEntry);
    Int_t nt = tin->GetEntry(iEntry);
    
    run=nsr->GetHeader().GetRun();
    subRun=nsr->GetHeader().GetSubRun();
    snarl=nsr->GetHeader().GetSnarl();
    
    const NtpSREventSummary& evthdr = nsr->evthdr;
    utc =evthdr.date.utc;
    day =evthdr.date.day;
    month =evthdr.date.month;
    year  =evthdr.date.year;
    
    // Get the time of the snarl, the value in the header is
    // accurate enough for this purpose
    VldTimeStamp snrl_vts = nsr->GetHeader().GetVldContext().GetTimeStamp();

    /*
    // Get the closest spill out of the database
    const BeamMonSpill* bms = sa.LoadSpill(snrl_vts);
    VldTimeStamp bms_vts;
    if (!bms) {
      cerr << "No BeamMonSpill for " << snrl_vts << endl;
      bms_vts = VldTimeStamp::GetEOT();
    }
    else {
      bms_vts = bms->SpillTime();
      //     cout << "bms_vts "<< bms_vts<<endl;
    }
    offset = bms_vts-snrl_vts;
    //    cout << "offset "<< offset<<endl;
    //    double streamSpillTimeDiff = offset.GetSeconds();
    
    //now determine if good beam
    Bool_t goodBeam=false;
    BMSpillAna bmana;
    bmana.UseDatabaseCuts();
    bmana.SetSpill(*bms);
    
    //set the time difference to determine if spill is in time
    VldTimeStamp delta_vts=snrl_vts-bms->SpillTime();
    bmana.SetTimeDiff(delta_vts.GetSeconds());
    
    //determine whether to select
    if(bmana.SelectSpill()) goodBeam=true;
    else goodBeam=false;
    */

    Bool_t goodBeam=true;
    //store the variable
    if (goodBeam) {
      
      ///////////////////////////////////////////
      goodbeam=1;
      const NtpSRTimeStatus& tstatus = nsr->timestatus;
      tmframe= tstatus.timeframe;
      crateT0= tstatus.crate_t0_ns;
      sgate53= tstatus.sgate_53mhz;
      
      cout<<setprecision(13);
      
      /*
      nearsec=nbd->nearest_sec;
      tor101=nbd->tor101;
      tr101d=nbd->tr101d;
      trtgtd=nbd->trtgtd;
      tortgt=nbd->tortgt;
      dtnear=nbd->dt_nearest;
      nearns=nbd->nearest_nsec;
      
      utc1=1000.*(utc+nearns/1000000000.);
      */

      // print out only relevant records for test
      // Double_t uu=(Int_t)utc%10000; 
      
      
      //    if (uu>=7730 && uu< 7757) {    
      //      cout<< "utc/crateT0/sgate53 "<< utc <<" "<<crateT0<<" "<<sgate53<<endl;
      //    cout<< "nearsec/dtnear/nearns "<< nearsec <<" "<< dtnear <<" "<<nearns<<" nearms "<< nearns/1000000<<endl;
      //   cout<< "utc/utc1 "<< utc <<" "<< utc1 <<" diff "<< utc1-utc <<endl;
      //    cout << "offset "<< offset <<endl;
      //    }
      
      if (snarlCount == 1) utcS=utc;
      
      // Commented out things that are only relevant for MC
      //vector<const NtpTHTrack*> thtrks=nsr->GetTHTracks();
      //vector<const NtpMCTruth*>  truths=nsr->GetMCTruth();
      //assert(trks.size()==thtrks.size());
      
      vector<const NtpSRTrack*>  trks=nsr->GetTracks();
      //       cout << "Retrieved " << trks.size() << " tracks. " <<endl;
      vector<const NtpSRShower*> shws=nsr->GetShowers();
      //    cout << "Retrieved " << shws.size() << " showers." <<endl;
      vector<const NtpSRStrip*>  strs=nsr->GetStrips();
      //       cout << "Retrieved " << strs.size() << " strips. " <<endl;
      
      
      // NB, done this way we loop over all the tracks in the snarl,
      // with no regard for whether they're part of events or     
      
      
      ntracks=(int)trks.size();
      
      int totnumtrkstp=0; 
      int ntra=0;
          vector<const NtpMCStdHep*>  truthstdhep=nsr->NtpStRecord::GetMCStdHeps();
    if ((int)truthstdhep.size() == 0) continue;     
    // go here only if stdhep truth info is in the file 

    int nmc = (int)truthstdhep.size();
    std::cout<<"nmc "<<nmc<<endl; 
    

    //    TClonesArray& mc_array=*(nsr->mc);  //mc not filled
    //    NtpMCTruth* mct = (NtpMCTruth*) mc_array[0];

    int icounter=0; 
      
      for(Short_t j=0; j<nmc; j++){
      
      icounter++;
      std::cout<<nmc<<" nmc j "<<j<<std::endl;
      
      TClonesArray& stdhepArray=*(nsr->stdhep);
      NtpMCStdHep* stdhep = (NtpMCStdHep*) stdhepArray[j];
      
      
      //if(stdhep->IstHEP==999) continue;
      if(stdhep->IstHEP!=1) continue;
      
      double p2 =0.0; for(int k=0; k<3; k++) p2+=stdhep->p4[k]*stdhep->p4[k];
      double e2 = stdhep->p4[3]*stdhep->p4[3];
      double diff=fabs(sqrt(e2)-sqrt(p2+stdhep->mass*stdhep->mass) );

      int child0=stdhep->child[0];
      
      int child1=stdhep->child[1];
      
      int parent0=stdhep->parent[0];
      
      int parent1=stdhep->parent[1];
      
      //if(stdhep->IstHEP!=1) continue;

      cout<<" parent1 "<<parent0<<" parent2 "<<parent1<<endl;
      cout<<" child1 "<<child0<<" child2 "<<child1<<endl;
      
      
      TParticle part(stdhep->IdHEP,stdhep->IstHEP,
		     parent0,parent1,child0,child1,0,0,0,0,0,0,0,0);
      double mass_term=stdhep->mass;
      cout<<" icount "<<icounter<<endl;
      cout<<" mass "<< mass_term<<endl;
      
      
      
      printf("%30s %+7.3f %+7.3f %+7.3f %7.3f %7.3f\n",
	     "",stdhep->vtx[0],stdhep->vtx[1],stdhep->vtx[2],stdhep->vtx[3]);
      
      mcEne = stdhep->p4[3];
      mcPDG = stdhep->IdHEP;
      mcMass = stdhep->mass;
      mcVtxX = stdhep->vtx[0];
      mcVtxY = stdhep->vtx[1];
      mcVtxZ = stdhep->vtx[2];
      mcPx = stdhep->p4[0];
      mcPy = stdhep->p4[1];
      mcPz = stdhep->p4[2];

      if(mcPDG == 13) { 
      printf("%30s %+7.3f %+7.3f %+7.3f %+7.3f %7.3f %7.3f\n",
	     "",mcEne, mcPDG, mcMass, mcPx, mcPy, mcPz); 
      
      cout<<" trkmom "<<trkmom<<" trkRange "<< trkErange<<endl;
      }
      break;
    }
      
      
      // loop over all tracks in the event
      for(int iTrk=0; iTrk < (int)trks.size(); iTrk++){
	
	const NtpSRTrack* trk=trks[iTrk];
	
	totnumtrkstp += trk->nstrip;
	ntra++;
      }
      
      //      cout << "Retrieved " << ntra << " # tracks. " <<endl; 
      //     if (totnumtrkstp!=0) {
      // loop over all tracks in the event
      for(int iTrk=0; iTrk < (int)trks.size(); ++iTrk){
	
	const NtpSRTrack* trk=trks[iTrk];
	
	trkIndex = trk->index;
	
	trkVtxX = trk->vtx.x;
	trkVtxY = trk->vtx.y;
	trkVtxZ = trk->vtx.z;
	trkVtxeX = trk->vtx.ex;
	trkVtxeY = trk->vtx.ey;
	trkVtxT = trk->vtx.t;
	
       trkdcosx = trk->vtx.dcosx;  // trk->vtx.dcosx or trk->lin.dcosx or trk->end.dcosx o
        trkdcosy = trk->vtx.dcosy;  // etc
        trkdcosz = trk->vtx.dcosz;
        trkContained = trk->contained;
	cout<<"old stopper "<<trkContained<<endl;
 

	//redefine the containment flag


	Int_t cont_flag=0;
        trkContained=0;
        double margin=0.12;
	// 1 = track is fully     contained in the upstream   region
	// 2 = track is partially contained in the upstream   region
	// 3 = track is fully     contained in the downstream region
	// 4 = track is partially contained in the downstream region

	double end_x = trk->end.x; //originally was evt.end.x - does it matter ?
	double end_y = trk->end.y; //
	double end_z = trk->end.z; //
	double end_u = trk->end.u; //
	double end_v = trk->end.v; //

	double rad   = TMath::Sqrt(end_x*end_x + end_y*end_y);

	//--- STOPPING
	// downstream
	bool down_stop = (end_y > -1.4 && end_y < 1.4 &&
			  end_x >-1.5 && end_x <  2.4 &&
			  end_u >-0.85 && end_u < 2.1 &&
			  end_v >-2.1 && end_v < 0.85 &&
			  rad   > 0.5 &&
			  end_z <  7);
	//ustream
	bool up_stop = (end_u >  0.3 && end_u < 1.8 &&
			end_v > -1.8 && end_v <-0.3 &&
			end_x <  2.4 && rad > 0.8 && end_z < 7);

	bool pass_uv = false;
        const double x = end_u/sqrt(2.0)-end_v/sqrt(2.0);

	// limits in meters
	const double ulow = -0.3;
        const double uhigh = 2.4;
        const double vlow = -2.4;
        const double vhigh = 0.3;
        const double vcoil = -0.2;
        const double ucoil = 0.2;
        const double xhigh = 2.8;

	// margin. cut into nominal fid area by this much
	const double mar = margin;
	
	if( (end_u > ulow+mar) && (end_u < uhigh-mar) &&
	    (end_v > vlow+mar) && (end_v < vhigh-mar) && (x < xhigh-mar) ){
	  if(end_u<ucoil+mar && end_v<vcoil-mar ){
	    pass_uv = true;
	  }
	  if(end_u>ucoil+mar && end_v>vcoil-mar){
	    pass_uv = true;
	  }
	  if(end_u>ucoil+mar  && end_v<vcoil-mar){
	    pass_uv = true;
	    cout<<"pass_uv is true"<<endl;
	  }
	} 
	//--- EXITING
	//downstream
// 	bool down_exit =  rad > 0.5 && (((end_y < -1.4 || end_y >  1.4  ||
// 					  end_x < -1.5 || end_x >  2.4 || end_u < -0.85 ||
// 					  end_u >  2.1 || end_v < -2.1 || end_v >  0.85) &&
// 					 end_z >  7 && end_z <16) || end_z >16);
// 	//upstream
// 	bool up_exit = rad > 0.5 && (end_u < 0.3 || end_u > 1.8 ||
// 				     end_v <-1.8 || end_v >-0.3 || end_x >2.4) && end_z <7;
	
	
	
	// set the track's containment flag
	
	if(up_stop) cont_flag = 1;
	else if(down_stop) cont_flag = 3;

// 	
// 	if (up_exit  ) cont_flag = 2;
// 	if (down_stop) cont_flag = 3;
// 	if (down_exit) cont_flag = 4;
	// if((cont_flag == 1 || cont_flag == 3) && pass_uv) {
	if((cont_flag == 1 || cont_flag == 3) && pass_uv) {
	  trkContained=1;
	  cout<<"stopper redefined"<<endl;
	}
	else if(cont_flag == 2 || cont_flag == 4) {
	  trkContained=0;
	  cout<<"exiter redefined"<<endl;
	}
	else {
	  cout<<"trkContainment flag has not been set"<<endl;
	  cout<<"end_u= "<<end_u<<" end_v= "<<end_v<<" rad= "<<rad<<endl;
	  cout<<" end_x= "<<end_x<<" end_y= "<<end_y<<" end_z= "<<end_z<<endl;   
	}  
	
	
	trkErange = trk->momentum.range;
        trkqp = trk->momentum.qp;
        trkeqp = trk->momentum.eqp;
        if (trkqp>0) charge=+1.;
        if (trkqp<0) charge=-1.;
        if (trkqp !=0) {
          trkmom = TMath::Abs(1./trkqp);
        }
        else {
          trkmom = 0;
        }
	
	trkChi2 = trk->fit.chi2;
	
	trkE = trk->ph.gev;
	
	trkTimeT0 = trk->time.t0;
	
	ntrkstp = trk->nstrip;
	// Copy the strip entries into the array that is linked with the branch
	for(int iTrkStp=0; iTrkStp<ntrkstp; ++iTrkStp) {
	  
	  trkstpX[iTrkStp] = trk->stpx[iTrkStp];
	  trkstpY[iTrkStp] = trk->stpy[iTrkStp];
	  trkstpZ[iTrkStp] = trk->stpz[iTrkStp];
	  trkstpU[iTrkStp] = trk->stpu[iTrkStp];
	  trkstpV[iTrkStp] = trk->stpv[iTrkStp];
	}
	
	// over events
	tout.Fill();
      } // over tracks
      //     }
    } // over good beam


	 //     }


  } // over snarls
  utcL=utc;
  
  tout.Write();
  cout << "Written output trees" << endl;
  cout << "utcS_"<<utcS<<"_"<<utcL<<endl;
}
