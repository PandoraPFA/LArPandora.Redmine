/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "GFKalman.h"

#include "assert.h"
#include <iostream>
#include <sstream>

#include "TMath.h"

#include "GFTrack.h"
#include "GFAbsRecoHit.h"
#include "GFAbsTrackRep.h"
#include "GFException.h"
  
#define COVEXC "cov_is_zero"


genf::GFKalman::GFKalman():fInitialDirection(1),fNumIt(3),fBlowUpFactor(500.){;}

genf::GFKalman::~GFKalman(){;}

void genf::GFKalman::processTrack(GFTrack* trk){
  int direction=fInitialDirection;
  assert(direction==1 || direction==-1);
  //  trk->clearGFBookkeeping();
  trk->clearRepAtHit();
  /*
  int nreps=trk->getNumReps();
  for(int i=0; i<nreps; ++i){ 
    trk->getBK(i)->setNhits(trk->getNumHits());
    trk->getBK(i)->bookMatrices("fPredSt");
    trk->getBK(i)->bookMatrices("fPredCov");
    trk->getBK(i)->bookMatrices("bPredSt");
    trk->getBK(i)->bookMatrices("bPredCov");
    trk->getBK(i)->bookGFDetPlanes("f");
    trk->getBK(i)->bookGFDetPlanes("b");
    trk->getBK(i)->bookNumbers("fChi2");
    trk->getBK(i)->bookNumbers("bChi2");
  }
  */

  /*why is there a factor of two here (in the for statement)?
    Because we consider one full iteration to be one back and
    one forth fitting pass */
  for(int ipass=0; ipass<2*fNumIt; ipass++){
    //    std::cout << "GFKalman:: numreps, numhits, ipass, direction"<< trk->getNumReps()<<", "<< trk->getNumHits()<<", "<< ipass<<", " <<direction << std::endl;
    if(ipass>0) blowUpCovs(trk);
    if(direction==1){
      trk->setNextHitToFit(0);
    }
    else {
      trk->setNextHitToFit(trk->getNumHits()-1);
    }
    fittingPass(trk,direction);
    
    //save first and last plane,state&cov after the fitting pass
    if(direction==1){//forward at last hit
      int nreps=trk->getNumReps();
      for(int i=0; i<nreps; ++i){
	trk->getTrackRep(i)->
	  setLastPlane( trk->getTrackRep(i)->getReferencePlane() );
	trk->getTrackRep(i)->
	  setLastState( trk->getTrackRep(i)->getState() );
	trk->getTrackRep(i)->
	  setLastCov( trk->getTrackRep(i)->getCov() );
      }
    }
    else{//backward at first hit
      int nreps=trk->getNumReps();
      for(int i=0; i<nreps; ++i){
	trk->getTrackRep(i)->
	  setFirstPlane( trk->getTrackRep(i)->getReferencePlane() );
	trk->getTrackRep(i)->
	  setFirstState( trk->getTrackRep(i)->getState() );
	trk->getTrackRep(i)->
	  setFirstCov( trk->getTrackRep(i)->getCov() );
      }
    }

    //switch direction of fitting and also inside all the reps
    if(direction==1) direction=-1;
    else direction=1;
    switchDirection(trk);

  }
  return;
}

void
genf::GFKalman::switchDirection(GFTrack* trk){
  int nreps=trk->getNumReps();
  for(int i=0; i<nreps; ++i){
    trk->getTrackRep(i)->switchDirection();
  }
}

void genf::GFKalman::blowUpCovs(GFTrack* trk){
  int nreps=trk->getNumReps();
  for(int irep=0; irep<nreps; ++irep){
    GFAbsTrackRep* arep=trk->getTrackRep(irep);
    //dont do it for already compromsied reps, since they wont be fitted anyway
    if(arep->getStatusFlag()==0) { 
      TMatrixT<Double_t> cov = arep->getCov();
      for(int i=0;i<cov.GetNrows();++i){
	for(int j=0;j<cov.GetNcols();++j){
	  //if(i!=j){//off diagonal
	  //  cov[i][j]=0.;
	  //}
	  //else{//diagonal
	  cov[i][j] = cov[i][j] * fBlowUpFactor;
	  //}
	}
      }
      arep->setCov(cov);
    }
  }  
}

void
genf::GFKalman::fittingPass(GFTrack* trk, int direction){
  //loop over hits

  unsigned int nhits=trk->getNumHits();
  unsigned int starthit=trk->getNextHitToFit();

  int nreps=trk->getNumReps();
  int ihit=(int)starthit;

  //clear chi2 sum and ndf sum in track reps
  for(int i=0;i<nreps;++i){
    GFAbsTrackRep* arep=trk->getTrackRep(i);
    arep->setChiSqu(0.);
    arep->setNDF(0);
  }

  //clear failedHits and outliers
  trk->getBK()->clearFailedHits();

  while((ihit<(int)nhits && direction==1) || (ihit>-1 && direction==-1)){
    //    GFAbsRecoHit* ahit=trk->getHit(ihit);
    // loop over reps
    for(int irep=0; irep<nreps; ++irep){
	  GFAbsTrackRep* arep=trk->getTrackRep(irep);
	  if(arep->getStatusFlag()==0) { 
	    try {
	      processHit(trk,ihit,irep,direction);
	    }
	    catch(GFException& e) {
	      trk->addFailedHit(irep,ihit);
	      std::cerr << e.what();
	      e.info();
	      if(e.isFatal()) {
		arep->setStatusFlag(1);
		continue; // go to next rep immediately
	      }
	    }	
	  }
    }// end loop over reps
    ihit+=direction;
  }// end loop over hits
  trk->setNextHitToFit(ihit-2*direction);
  //trk->printGFBookkeeping();
}

double genf::GFKalman::chi2Increment(const TMatrixT<Double_t>& r,const TMatrixT<Double_t>& H,
			     const TMatrixT<Double_t>& cov,const TMatrixT<Double_t>& V){

  // residuals covariances:R=(V - HCH^T)
  TMatrixT<Double_t> R(V);
  TMatrixT<Double_t> covsum1(cov,TMatrixT<Double_t>::kMultTranspose,H);
  TMatrixT<Double_t> covsum(H,TMatrixT<Double_t>::kMult,covsum1);

  R-=covsum;

  // chisq= r^TR^(-1)r
  double det=0.;
  TMatrixT<Double_t> Rsave(R);
  R.SetTol(1.0e-23); // to avoid inversion problem, EC, 8-Aug-2011.
  R.Invert(&det);
  TMatrixT<Double_t> residTranspose(r);
  residTranspose.T();
  TMatrixT<Double_t> chisq=residTranspose*(R*r);
  assert(chisq.GetNoElements()==1);

  if(TMath::IsNaN(chisq[0][0])){
	GFException exc("chi2 is nan",__LINE__,__FILE__);
	exc.setFatal();
	std::vector<double> numbers;
	numbers.push_back(det);
	exc.setNumbers("det",numbers);
	std::vector< TMatrixT<Double_t> > matrices;
	matrices.push_back(r);
	matrices.push_back(V);
	matrices.push_back(Rsave);
	matrices.push_back(R);
	matrices.push_back(cov);
	exc.setMatrices("r, V, Rsave, R, cov",matrices);
    throw exc;
  }

  return chisq[0][0];
}


double
genf::GFKalman::getChi2Hit(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{
  // get prototypes for matrices
  int repDim=rep->getDim();
  TMatrixT<Double_t> state(repDim,1);
  TMatrixT<Double_t> cov(repDim,repDim);;
  GFDetPlane pl=hit->getDetPlane(rep);
  rep->extrapolate(pl,state,cov);


  TMatrixT<Double_t> H = hit->getHMatrix(rep);
  // get hit covariances  
  TMatrixT<Double_t> V=hit->getHitCov(pl);

  TMatrixT<Double_t> r=hit->residualVector(rep,state,pl);
  assert(r.GetNrows()>0);

  //this is where chi2 is calculated
  double chi2 = chi2Increment(r,H,cov,V);

  return chi2/r.GetNrows();
}

  
void
genf::GFKalman::processHit(GFTrack* tr, int ihit, int irep,int direction){
  GFAbsRecoHit* hit = tr->getHit(ihit);
  GFAbsTrackRep* rep = tr->getTrackRep(irep);

  // get prototypes for matrices
  int repDim = rep->getDim();
  TMatrixT<Double_t> state(repDim,1);
  TMatrixT<Double_t> cov(repDim,repDim);;
  GFDetPlane pl;

  /* do an extrapolation, if the trackrep irep is not given
   * at this ihit position. This will usually be the case, but
   * not if the fit turnes around
   */
  //  std::cout << "GFKalman::ProcessHit(): ihit is " << ihit << std::endl;
  //  std::cout << "GFKalman::ProcessHit(): direction is "  << direction << std::endl;
  //rep->Print();

  // Take out this if/else condition. EC, 5-Jan-2011. Put it back. 6-Jan.

    if(ihit!=tr->getRepAtHit(irep)){
    //std::cout << "not same" << std::endl;
    // get the (virtual) detector plane
      pl=hit->getDetPlane(rep);

      //      std::cout << "GFKalman::ProcessHit(): hit is ... " << ihit << std::endl;
      //hit->Print();
      //std::cout << "GFKalman::ProcessHit(): plane is ... " <<  std::endl;
      //pl.Print();

    //do the extrapolation
      rep->extrapolate(pl,state,cov);
    }
      else{
    //std::cout << "same" << std::endl;
    //    return;

	pl = rep->getReferencePlane();
	state = rep->getState();
	cov = rep->getCov();
      }

  
  if(cov[0][0]<1.E-50){
        std::cout<<"GFKalman::processHit() 0. Calling Exception."<<std::endl;
    GFException exc(COVEXC,__LINE__,__FILE__);
        std::cout<<"GFKalman::processHit() 1. About to throw GFException."<<std::endl;
    throw exc;
        std::cout<<"GFKalman::processHit() 2. Back from GFException."<<std::endl;
  }
  
  /*
  if(direction==1){
    tr->getBK(irep)->setMatrix("fPredSt",ihit,state);
    tr->getBK(irep)->setMatrix("fPredCov",ihit,cov);
    tr->getBK(irep)->setGFDetPlane("f",ihit,pl);
  }
  else{
    tr->getBK(irep)->setMatrix("bPredSt",ihit,state);
    tr->getBK(irep)->setMatrix("bPredCov",ihit,cov);
    tr->getBK(irep)->setGFDetPlane("b",ihit,pl);
  }
  */

  //  TMatrixT<Double_t> origcov=rep->getCov();
  TMatrixT<Double_t> H=hit->getHMatrix(rep);
  // get hit covariances  
  TMatrixT<Double_t> V=hit->getHitCov(pl);
  // calculate kalman gain ------------------------------
  TMatrixT<Double_t> Gain(calcGain(cov,V,H));

	//	 std::cout << "GFKalman:: processHits(), state is  " << std::endl;
	//rep->getState().Print();

  TMatrixT<Double_t> res=hit->residualVector(rep,state,pl);
  // calculate update -----------------------------------

  TMatrixT<Double_t> update=Gain*res;
  state+=update; // prediction overwritten!
  cov-=Gain*(H*cov);

  // calculate filtered chisq
  // filtered residual
  TMatrixT<Double_t> r=hit->residualVector(rep,state,pl);
  double chi2 = chi2Increment(r,H,cov,V);
  int ndf = r.GetNrows();
  rep->addChiSqu( chi2 );
  rep->addNDF( ndf );

  /*
  if(direction==1){
    tr->getBK(irep)->setNumber("fChi2",ihit,chi2/ndf);
  }
  else{
    tr->getBK(irep)->setNumber("bChi2",ihit,chi2/ndf);
  }
  */

  // if we survive until here: update TrackRep
  //rep->setState(state);
  //rep->setCov(cov);
  //rep->setReferencePlane(pl);

  rep->setData(state,pl,&cov);
  tr->setRepAtHit(irep,ihit);

}


TMatrixT<Double_t>
genf::GFKalman::calcGain(const TMatrixT<Double_t>& cov, 
			 const TMatrixT<Double_t>& HitCov,
			 const TMatrixT<Double_t>& H){

  // calculate covsum (V + HCH^T)

  // Comment next 3 out for normal running.
  //cov.Print();
  //H.Print();
  //HitCov.Print();
  TMatrixT<Double_t> covsum1(cov,TMatrixT<Double_t>::kMultTranspose,H);
  //  covsum1.Print();
  TMatrixT<Double_t> covsum(H,TMatrixT<Double_t>::kMult,covsum1);

  //covsum.Print();

  covsum+=HitCov;
  //covsum = (covsum,TMatrixT<Double_t>::kPlus,HitCov);
  // covsum.Print();

  
  // invert
  double det=0;
  covsum.SetTol(1.0e-23); // to avoid inversion problem, EC, 8-Aug-2011.
  covsum.Invert(&det);
  //    std::cout << "GFKalman:: calGain(), det is  "<< det << std::endl;
  if(TMath::IsNaN(det)) {
    GFException e("Kalman Gain: det of covsum is nan",__LINE__,__FILE__);
    e.setFatal();
    throw e;
  }

  if(det==0){
	GFException exc("cannot invert covsum in Kalman Gain - det=0",
						__LINE__,__FILE__);
	exc.setFatal();
	std::vector< TMatrixT<Double_t> > matrices;
	matrices.push_back(cov);
	matrices.push_back(HitCov);
	matrices.push_back(covsum1);
	matrices.push_back(covsum);
	exc.setMatrices("cov, HitCov, covsum1 and covsum",matrices);
    throw exc;

  }

  // calculate gain
  TMatrixT<Double_t> gain1(H,TMatrixT<Double_t>::kTransposeMult,covsum);
  TMatrixT<Double_t> gain(cov,TMatrixT<Double_t>::kMult,gain1);

  return gain;
}







