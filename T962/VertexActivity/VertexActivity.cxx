////////////////////////////////////////////////////////////////////////
//
// VertexActivity class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to characterize the vertex activity associated with a neutrino event
////////////////////////////////////////////////////////////////////////


#include <iostream>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "T962/VertexActivity/VertexActivity.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "TMath.h"


#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
// #include "T962_MergeData/ScanInfo.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

//-----------------------------------------------------------------------------
vertex::VertexActivity::VertexActivity(fhicl::ParameterSet const& pset) :
  fDBScanModuleLabel            (pset.get< std::string >("DBScanModuleLabel")),
  fLArG4ModuleLabel             (pset.get< std::string >("LArG4ModuleLabel")),
  fHitsModuleLabel              (pset.get< std::string >("HitsModuleLabel")),  
  fGenieGenModuleLabel          (pset.get< std::string >("GenieGenModuleLabel")),
  fScanModuleLabel              (pset.get< std::string > ("ScanModuleLabel")),
  fCathodetimelocation          (pset.get< double >("Cathodetimelocation")),
  fDelta_Cathodetimelocation    (pset.get< double >("Delta_Cathodetimelocation")),
  fE_lifetime                   (pset.get< double >("E_lifetime")),
  fDelta_E_lifetime             (pset.get< double >("Delta_E_lifetime")),
  fRecombination_factor         (pset.get< double >("Recombination_factor")),
  fDelta_Recombination_factor   (pset.get< double >("Delta_Recombination_factor")),
  fWorkfunction_factor          (pset.get< double >("Workfunction_factor")),
  fDelta_Workfunction_factor    (pset.get< double >("Delta_Workfunction_factor")),
  fCalibration_factor           (pset.get< double >("Calibration_factor")),
  fDelta_Calibration_factor     (pset.get< double >("Delta_Calibration_factor")),
  fActivityRadius               (pset.get< double >("ActivityRadius"))
{
  produces< std::vector<recob::EndPoint2D> >();
}

//-----------------------------------------------------------------------------
vertex::VertexActivity::~VertexActivity()
{
}

//-----------------------------------------------------------------------------


void vertex::VertexActivity::beginJob()
{
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;
  fIndEfficiency    = tfs->make<TH2F>("Ind Efficiency",  "Ind Efficiency vs E_rec", 100, 0, 300, 100,0,2);
  fColEfficiency    = tfs->make<TH2F>("Col Efficiency",  "Col Efficiency vs E_rec", 100, 0, 300, 100,0,2);
  fEfficiency    = tfs->make<TH2F>("Efficiency",  "Efficiency vs E_rec", 100, 0, 300, 100,0,2);
  findcol    = tfs->make<TH2F>("IndCol_mips",  "IndCol_mips", 1000, 0, 1000, 1000,0,1000);
}

//-----------------------------------------------------------------------------
void vertex::VertexActivity::produce(art::Event& evt)
{
  double vertexcolwire=-1.;
  double vertexcoltime=-1.;
  double efficiency=1.;

  double elifetime_factor=0;
  double hitamplitude=0.;
  double vertex [3] = { 0, 0, 0 };
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larp;
  double electronlifetime=larp->ElectronLifetime();
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);
    
//   art::Handle< std::vector<recob::Vertex> > vertexListHandle;
//   evt.getByLabel(fVertexModuleLabel,vertexListHandle);
  
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);

  sim::LArVoxelList vxlistHandle = sim::SimListUtils::GetLArVoxelList(evt, fLArG4ModuleLabel);
  
  art::PtrVector<simb::MCTruth> mclist;
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    } 
   
  //  neutrinos  
//      for( unsigned int i = 0; i < mclist.size(); ++i ){
// 
//      art::Ptr<simb::MCTruth> mc(mclist[i]);
// 
// 	simb::MCNeutrino neut(mc->GetNeutrino());
// 
//     std::cout<<"vertex: "<<neut.Nu().Vx()<<" "<<neut.Nu().Vy()<<" "<<neut.Nu().Vz()<<std::endl;
//     vertex[0] =neut.Nu().Vx();
//     vertex[1] =neut.Nu().Vy();
//     vertex[2] =neut.Nu().Vz();
// 
//     }
//    
//    
   //muons
  for( unsigned int i = 0; i < mclist.size(); ++i ){

    art::Ptr<simb::MCTruth> mc(mclist[i]);

    simb::MCParticle neut(mc->GetParticle(i));

    // std::cout<<"vertex: "<<neut.Nu().Vx()<<" "<<neut.Nu().Vy()<<" "<<neut.Nu().Vz()<<std::endl;
    vertex[0] =neut.Vx();
    vertex[1] =neut.Vy();
    vertex[2] =neut.Vz();

  }
   
  // There's probably only one LArVoxelList per event, but FMWK
  // always reads a vector of pointers.  For each LArVoxelList:
  double Energy=0.;
  
  sim::LArVoxelList::const_iterator vxitr;
  for(vxitr = vxlistHandle.begin(); vxitr != vxlistHandle.end(); vxitr++){
    // Get the reference to the LArVoxelID in the LArVoxelList.
    const sim::LArVoxelData &voxel = (*vxitr).second;
    
    int numberParticles = voxel.NumberParticles();
      
    //  std::cout<<"numberParticles "<<numberParticles<<std::endl;
    for ( int i = 0; i != numberParticles; ++i )
      {
	// if(sqrt(pow(TMath::Abs(voxel.VoxelID().X()-vertex[0]),2)+pow(TMath::Abs(voxel.VoxelID().Y()-vertex[1]),2)+pow(TMath::Abs(voxel.VoxelID().Z()-vertex[2]),2))<fActivityRadius)
		
		
	if(TMath::Abs(voxel.VoxelID().X()-vertex[0])<=(fActivityRadius*.4)
	   &&((voxel.VoxelID().Z()-vertex[2])-((sqrt(3)/3)*(voxel.VoxelID().Y()-vertex[1]))+((8*sqrt(3)/3)*fActivityRadius*.1))>=0
	   &&((voxel.VoxelID().Z()-vertex[2])-((sqrt(3)/3)*(voxel.VoxelID().Y()-vertex[1]))-((8*sqrt(3)/3)*fActivityRadius*.1))<=0
	   &&((voxel.VoxelID().Z()-vertex[2])+((sqrt(3)/3)*(voxel.VoxelID().Y()-vertex[1]))-((8*sqrt(3)/3)*fActivityRadius*.1))<=0
	   &&((voxel.VoxelID().Z()-vertex[2])+((sqrt(3)/3)*(voxel.VoxelID().Y()-vertex[1]))+((8*sqrt(3)/3)*fActivityRadius*.1))>=0
	   ){
	  
	  // std::cout<<TMath::Abs(voxel.VoxelID().X()-vertex[0])<<" "<<TMath::Abs(voxel.VoxelID().Y()-vertex[1])<<" "<<TMath::Abs(voxel.VoxelID().Z()-vertex[2])<<" "<<voxel.Energy(i)*1000.<<std::endl;
	  
	  Energy+= (voxel.Energy(i)*1000.);
	  //   if(voxel.Energy(i))
	  //       std::cout<<(voxel.VoxelID().Z()-vertex[2])*10<<std::endl;
	}
      }
        
        
    //std::cout<<voxel.VoxelID().X()<<" "<<voxel.VoxelID().Y()<<" "<<voxel.VoxelID().Z()<<" "<<Energy<<std::endl;    
  }	      
  //std::cout<<"Energy: "<<Energy<<" vxlisthandlesize "<<vxlistHandle->size()<<std::endl;  


		      
  //hand scan info
  //   art::PtrVector<merge::ScanInfo> scanIn;
  //   scanIn.clear();
  //   art::Handle< std::vector<merge::ScanInfo> > scanHandle;
  //   evt.getByLabel(fScanModuleLabel,scanHandle);
  // 
  //  for(unsigned int i = 0; i < scanHandle->size(); ++i){     
  //      art::Ptr<merge::ScanInfo> scaninfo(scanHandle, i);
  //       scanIn.push_back(scaninfo);     
  //     }
  // 
  //     for(unsigned int i = 0; i < scanIn.size(); ++i){     
  //      vertexcolwire=scanIn[i]->Get_VertColWire();
  //      vertexcoltime=scanIn[i]->Get_VertColTime();   
  //     }
  //take into account automated vertex finding as well (not yet implemented)  
  //   art::PtrVector<recob::Vertex> vertIn;
  // 
  //   for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii)
  //     {
  //       art::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
  //       vertIn.push_back(vertex);
  //     }
  //   
  //    for(unsigned int i = 0; i < vertIn.size(); ++i){     
  //      std::cout<<vertIn[i]->WireNum()<<" "<<vertIn[i]->DriftTime()<<std::endl;
  //     }
  
  

  filter::ChannelFilter chanFilt;  
  art::PtrVector<recob::Hit> cHits;
  art::PtrVector<recob::Hit> hit;
   
//   art::PtrVector<recob::Cluster> clusIn;
//   for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
//     {
//       art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
//       clusIn.push_back(cluster);
//     }
//   
  art::Handle< std::vector<recob::Hit> > hitcol;
  evt.getByLabel(fHitsModuleLabel,hitcol);
  
  
  ///loop over all hits in the event and look for clusters (for each plane)
  
  
  art::PtrVector<recob::Hit> allhits;

  
  int numberwires;
  double numbertimesamples;

  unsigned int channel2,plane,plane2,wire,wire2,plane3;
  unsigned int p(0),w(0), t(0),channel(0);
  
      

  //std::cout<<"vertex0 "<<vertex[0]<<std::endl;
  double drifttick=(vertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./.198);
      

      

	  
	  
	  
  // std::cout<<"vertex transform collection "<<wire2<<" "<<drifttick<<std::endl;
  double energyofhits_ind=0.;
  double energyofhits_col=0.;
  std::vector<float> hit_ind;
  std::vector<float> hit_col;
  hit_ind.clear();
  hit_col.clear();
  for(int plane = 0; plane < geom->Nplanes(); plane++) 
    {
     allhits.clear();
        for(unsigned int i = 0; i< hitcol->size(); ++i){
  
	art::Ptr<recob::Hit> hit(hitcol, i);
  
  
	channel=hit->Wire()->RawDigit()->Channel();
	geom->ChannelToWire(channel,t,p,w);
    
	if(p == plane) allhits.push_back(hit);
   
      }

    
    
      if(plane==0)
	vertex[0]=.3;//force time coordinate to be closer to induction plane 
      else
	vertex[0]=-.3;//force time coordinate to be closer to collection plane
     
      channel2 = geom->NearestChannel(vertex);
      geom->ChannelToWire(channel2,t,plane2,wire2); 
      // std::cout<<"channel2 "<<channel2<<std::endl;
  //     art::PtrVector<recob::Hit> vHits;
//       art::PtrVector<recob::Cluster>::const_iterator clusterIter = clusIn.begin();
//       hit.clear();
//       cHits.clear();      
      
      
//       while(clusterIter!= clusIn.end() ) {
// 	cHits = (*clusterIter)->Hits(p);
// 	if(cHits.size() > 0)
// 	  for(int i = 0; i < cHits.size(); i++)
// 	    hit.push_back(cHits[i]);
//       
// 	clusterIter++;  
//       } 
// 
//    std::cout<<"allhits size "<<allhits.size()<<" "<<hitcol->size()<<std::endl;

      if(allhits.size() == 0) 
        continue;

      numberwires=geom->Nwires(0);
      numbertimesamples=allhits[0]->Wire()->fSignal.size();
      

	  
      for(unsigned int i=0;i < allhits.size(); i++)
	{
	
	  // std::cout<<"hit size: "<<hit.size()<<std::endl;
	
	   channel=allhits[i]->Wire()->RawDigit()->Channel();
	   geom->ChannelToWire(channel,t,plane3,wire);   
// // std::cout<<hit[i]->Charge()<<" "<<hit[i]->PeakTime()<<" "<<hit[i]->Channel()<<std::endl;
// 	  
// 	  if(plane!=p)
// 	    continue;

	  //units are cm:
	  // if(sqrt(pow(TMath::Abs(vertexcolwire-wire)*.0743,2)+pow(TMath::Abs(vertexcoltime-hit[i]->CrossingTime()),2)) < fActivityRadius) 
	  // std::cout<<"wires "<<wire<<" "<<wire2<<" "<<hit[i]->CrossingTime()<<" "<<drifttick<<std::endl;
	  // 	  
	  // 	  if(larp->DriftVelocity(larp->Efield(),larp->Temperature())*.198*sqrt(pow(TMath::Abs((int)(wire2-wire))*(1/(larp->DriftVelocity(larp->Efield(),larp->Temperature())*(1/.4)*.198)),2)+pow(TMath::Abs(drifttick-hit[i]->CrossingTime()),2)) < fActivityRadius) 
	  // 	  {
	  
	  // std::cout<<i<<" "<<wire<<" "<<wire2<<" "<<fActivityRadius<<" "<<drifttick<<" "<<allhits[i]->PeakTime()<<std::endl;
	  if(TMath::Abs((int)(wire2-wire))<=fActivityRadius&&TMath::Abs(drifttick-allhits[i]->PeakTime())*.0743<=fActivityRadius)
	    {
	       //std::cout<<"drifttick "<<drifttick<<" hit crossing time "<<allhits[i]->PeakTime()<<std::endl;
	  
	      hitamplitude=allhits[i]->Charge()/20.4;//20.4 is scale factor between charge and energy (empirical, from plot)
	      //elifetime_factor=TMath::Exp((-hit[i]->CrossingTime()*.198)/electronlifetime); 
	      elifetime_factor=1.;
	      //elifetime_factor=1;
	      //6241.5 electrons/fC
	  

	      if(plane==0)	  {energyofhits_ind+=(6241.5*fCalibration_factor*fWorkfunction_factor/fRecombination_factor)*(1/elifetime_factor)*hitamplitude*.000001*2.64;
		hit_ind.push_back(hitamplitude*2.64);//2.64=col/ind
	      }
	      if(plane==1){energyofhits_col+=((6241.5*fCalibration_factor*fWorkfunction_factor/fRecombination_factor)*(1/elifetime_factor)*hitamplitude*.000001);
		hit_col.push_back(hitamplitude);
	      }
 	      if(plane==0)	   
 		std::cout<<plane<<" "<<wire<<" "<<wire2<<" "<<larp->DriftVelocity(larp->Efield(),larp->Temperature())*.198*sqrt(pow(TMath::Abs((int)(wire2-wire))*(1/(larp->DriftVelocity(larp->Efield(),larp->Temperature())*(1/.4)*.198)),2))<<" "<<allhits[i]->PeakTime()<<" "<<drifttick<<" "<<hitamplitude*2.64<<" "<<energyofhits_ind<<std::endl;
//   	   
	      if(plane==1)
		std::cout<<plane<<" "<<wire<<" "<<wire2<<" "<<larp->DriftVelocity(larp->Efield(),larp->Temperature())*.198*sqrt(pow(TMath::Abs((int)(wire2-wire))*(1/(larp->DriftVelocity(larp->Efield(),larp->Temperature())*(1/.4)*.198)),2))<<" "<<allhits[i]->PeakTime()<<" "<<drifttick<<" "<<hitamplitude<<" "<<energyofhits_col<<std::endl;
  	   
	    }
	    else
	    {
	     if(plane==0)
		std::cout<<"not counted "<<plane<<" "<<wire<<" "<<wire2<<" "<<larp->DriftVelocity(larp->Efield(),larp->Temperature())*.198*sqrt(pow(TMath::Abs((int)(wire2-wire))*(1/(larp->DriftVelocity(larp->Efield(),larp->Temperature())*(1/.4)*.198)),2))<<" "<<allhits[i]->PeakTime()<<" "<<drifttick<<" "<<hitamplitude<<" "<<energyofhits_col<<std::endl;
	    
	    
	     if(plane==1)
		std::cout<<"not counted "<<plane<<" "<<wire<<" "<<wire2<<" "<<larp->DriftVelocity(larp->Efield(),larp->Temperature())*.198*sqrt(pow(TMath::Abs((int)(wire2-wire))*(1/(larp->DriftVelocity(larp->Efield(),larp->Temperature())*(1/.4)*.198)),2))<<" "<<allhits[i]->PeakTime()<<" "<<drifttick<<" "<<hitamplitude<<" "<<energyofhits_col<<std::endl;
	    
	    }
	  

	}
	
	
      // double A=TMath::Abs(energyofhits_ind-energyofhits_col)/sqrt(pow(energyofhits_ind,2)+pow(energyofhits_col,2));
	
      double A=TMath::Abs(energyofhits_ind-energyofhits_col)/sqrt(pow(energyofhits_ind,2)+pow(energyofhits_col,2)-(energyofhits_ind*energyofhits_col));
	
      double B=1-A;

      double averageen=(energyofhits_ind+energyofhits_col)/2.;
      double reduceden=(energyofhits_ind*energyofhits_col)/(energyofhits_ind+energyofhits_col);
	
	
      double finalenergy=((A*reduceden)+(B*averageen))/(A+B);
      if(plane==1)
	{
	  std::cout<<"Energy2 "<<Energy<<" MeV "<<finalenergy<<" "<<energyofhits_ind<<" "<<energyofhits_col<<std::endl;
	  fEfficiency->Fill(Energy,finalenergy/Energy);
	  fIndEfficiency->Fill(Energy,energyofhits_ind/Energy);
	  fColEfficiency->Fill(Energy,energyofhits_col/Energy);
	}
      //Energy=0.;
    
    
      //       for(int vertexnum=0;vertexnum<vertexListHandle->size();vertexnum++)
      // 	{
      // 
      // 
      // 	}

      allhits.clear();



    }
   
  int size=hit_ind.size();
  int size2=hit_col.size();
  int size3=0;
  if(size>size2)
    size3=size2;
  if(size<size2)
    size3=size;
  if(size==size2)
    size3=size;
    
 //  std::cout<<"size3 "<<size3<<" "<<size2<<" "<<size<<std::endl;
  for(int i=0;i<size3;i++)
    findcol->Fill(hit_ind[i],hit_col[i]);

 }



