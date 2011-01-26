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
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "VertexActivity.h"
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
#include "Simulation/LArVoxelCalculator.h"
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "T962_MergeData/ScanInfo.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

//-----------------------------------------------------------------------------
vertex::VertexActivity::VertexActivity(edm::ParameterSet const& pset) :
  fDBScanModuleLabel            (pset.getParameter< std::string >("DBScanModuleLabel")),
  fLArG4ModuleLabel         (pset.getParameter< std::string >("LArG4ModuleLabel")        ),
  fVertexModuleLabel            (pset.getParameter<std::string > ("VertexModuleLabel")),
  fGenieGenModuleLabel      (pset.getParameter< std::string >("GenieGenModuleLabel")     ),
  fScanModuleLabel              (pset.getParameter< std::string > ("ScanModuleLabel")),
  fCathodetimelocation          (pset.getParameter< double >("Cathodetimelocation")),
  fDelta_Cathodetimelocation    (pset.getParameter< double >("Delta_Cathodetimelocation")),
  fE_lifetime                   (pset.getParameter< double >("E_lifetime")),
  fDelta_E_lifetime             (pset.getParameter< double >("Delta_E_lifetime")),
  fRecombination_factor         (pset.getParameter< double >("Recombination_factor")),
  fDelta_Recombination_factor   (pset.getParameter< double >("Delta_Recombination_factor")),
  fWorkfunction_factor          (pset.getParameter< double >("Workfunction_factor")),
  fDelta_Workfunction_factor    (pset.getParameter< double >("Delta_Workfunction_factor")),
  fCalibration_factor           (pset.getParameter< double >("Calibration_factor")),
  fDelta_Calibration_factor     (pset.getParameter< double >("Delta_Calibration_factor")),
  fActivityRadius               (pset.getParameter< double >("ActivityRadius"))
{
  produces< std::vector<recob::Vertex> >();
}

//-----------------------------------------------------------------------------
vertex::VertexActivity::~VertexActivity()
{
}

//-----------------------------------------------------------------------------


void vertex::VertexActivity::beginJob(edm::EventSetup const&)
{
    // get access to the TFile service
    edm::Service<edm::TFileService> tfs;
    fIndEfficiency    = tfs->make<TH2F>("Ind Efficiency",  "Ind Efficiency vs E_rec", 100, 0, 100, 100,0,2);
    fColEfficiency    = tfs->make<TH2F>("Col Efficiency",  "Col Efficiency vs E_rec", 100, 0, 100, 100,0,2);
    fEfficiency    = tfs->make<TH2F>("Efficiency",  "Efficiency vs E_rec", 100, 0, 100, 100,0,2);
    findcol    = tfs->make<TH2F>("IndCol_mips",  "IndCol_mips", 1000, 0, 100, 1000,0,100);
}

//-----------------------------------------------------------------------------
void vertex::VertexActivity::produce(edm::Event& evt, edm::EventSetup const&)
{
double vertexcolwire=-1.;
double vertexcoltime=-1.;
double efficiency=1.;

double elifetime_factor=0;
double hitamplitude=0.;
double vertex [3] = { 0, 0, 0 };
  edm::Service<geo::Geometry> geom;
  edm::Service<util::LArProperties> larp;
  double electronlifetime=larp->ElectronLifetime();
  edm::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);
    
  edm::Handle< std::vector<recob::Vertex> > vertexListHandle;
  evt.getByLabel(fVertexModuleLabel,vertexListHandle);
  
  edm::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);
  
   // std::cout << "Vertex  list size = " << vertexListHandle->size() << " AND cluster hit size:" <<  clusterListHandle->size()<< std::endl;

    edm::Handle< std::vector<sim::LArVoxelData> > vxlistHandle;
    evt.getByLabel(fLArG4ModuleLabel,vxlistHandle);

    // std::cout<<"vxlistHandle->size() "<<vxlistHandle->size()<<std::endl;

   edm::PtrVector<simb::MCTruth> mclist;
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      edm::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    }
   
//  neutrinos  
//      for( unsigned int i = 0; i < mclist.size(); ++i ){
// 
//      edm::Ptr<simb::MCTruth> mc(mclist[i]);
// 
// 	simb::MCNeutrino neut(mc->GetNeutrino());
// 
//     std::cout<<"vertex: "<<neut.Nu().Vx()<<" "<<neut.Nu().Vy()<<" "<<neut.Nu().Vz()<<std::endl;
//     vertex[0] =neut.Nu().Vx();
//     vertex[1] =neut.Nu().Vy();
//     vertex[2] =neut.Nu().Vz();
// 
//     }
   
   
   
    for( unsigned int i = 0; i < mclist.size(); ++i ){

     edm::Ptr<simb::MCTruth> mc(mclist[i]);

	simb::MCParticle neut(mc->GetParticle(i));

    // std::cout<<"vertex: "<<neut.Nu().Vx()<<" "<<neut.Nu().Vy()<<" "<<neut.Nu().Vz()<<std::endl;
    vertex[0] =neut.Vx();
    vertex[1] =neut.Vy();
    vertex[2] =neut.Vz();

    }
   
    // There's probably only one LArVoxelList per event, but FMWK
    // always reads a vector of pointers.  For each LArVoxelList:
    double Energy=0.;
    for(unsigned int i = 0; i < vxlistHandle->size(); ++i){
      // Get the reference to the LArVoxelID in the LArVoxelList.
      edm::Ptr<sim::LArVoxelData> voxel(vxlistHandle, i);
      
        int numberParticles = voxel->NumberParticles();
	      
	      // std::cout<<"numberParticles "<<numberParticles<<std::endl;
	      for ( int i = 0; i != numberParticles; ++i )
		{
		// if(sqrt(pow(TMath::Abs(voxel->VoxelID().X()-vertex[0]),2)+pow(TMath::Abs(voxel->VoxelID().Y()-vertex[1]),2)+pow(TMath::Abs(voxel->VoxelID().Z()-vertex[2]),2))<fActivityRadius)
		
		
	  if(TMath::Abs(voxel->VoxelID().X()-vertex[0])<=(fActivityRadius*.4)
	  &&((voxel->VoxelID().Z()-vertex[2])-((sqrt(3)/3)*(voxel->VoxelID().Y()-vertex[1]))+((8*sqrt(3)/3)*fActivityRadius*.1))>=0
	  &&((voxel->VoxelID().Z()-vertex[2])-((sqrt(3)/3)*(voxel->VoxelID().Y()-vertex[1]))-((8*sqrt(3)/3)*fActivityRadius*.1))<=0
	  &&((voxel->VoxelID().Z()-vertex[2])+((sqrt(3)/3)*(voxel->VoxelID().Y()-vertex[1]))-((8*sqrt(3)/3)*fActivityRadius*.1))<=0
	  &&((voxel->VoxelID().Z()-vertex[2])+((sqrt(3)/3)*(voxel->VoxelID().Y()-vertex[1]))+((8*sqrt(3)/3)*fActivityRadius*.1))>=0
	  ){
	  
	  // std::cout<<TMath::Abs(voxel->VoxelID().X()-vertex[0])<<" "<<TMath::Abs(voxel->VoxelID().Y()-vertex[1])<<" "<<TMath::Abs(voxel->VoxelID().Z()-vertex[2])<<" "<<voxel->Energy(i)*1000.<<std::endl;
	  
      Energy+= (voxel->Energy(i)*1000.);
    //   if(voxel->Energy(i))
//       std::cout<<(voxel->VoxelID().Z()-vertex[2])*10<<std::endl;
      }
        }
        
        
       //std::cout<<voxel->VoxelID().X()<<" "<<voxel->VoxelID().Y()<<" "<<voxel->VoxelID().Z()<<" "<<Energy<<std::endl;    
	}	      
 //std::cout<<"Energy: "<<Energy<<" vxlisthandlesize "<<vxlistHandle->size()<<std::endl;  


		      
  //hand scan info
//   edm::PtrVector<merge::ScanInfo> scanIn;
//   scanIn.clear();
//   edm::Handle< std::vector<merge::ScanInfo> > scanHandle;
//   evt.getByLabel(fScanModuleLabel,scanHandle);
// 
//  for(unsigned int i = 0; i < scanHandle->size(); ++i){     
//      edm::Ptr<merge::ScanInfo> scaninfo(scanHandle, i);
//       scanIn.push_back(scaninfo);     
//     }
// 
//     for(unsigned int i = 0; i < scanIn.size(); ++i){     
//      vertexcolwire=scanIn[i]->Get_VertColWire();
//      vertexcoltime=scanIn[i]->Get_VertColTime();   
//     }
//take into account automated vertex finding as well (not yet implemented)  
//   edm::PtrVector<recob::Vertex> vertIn;
// 
//   for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii)
//     {
//       edm::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
//       vertIn.push_back(vertex);
//     }
//   
//    for(unsigned int i = 0; i < vertIn.size(); ++i){     
//      std::cout<<vertIn[i]->WireNum()<<" "<<vertIn[i]->DriftTime()<<std::endl;
//     }
  
  

  filter::ChannelFilter chanFilt;  
  edm::PtrVector<recob::Hit> cHits;
  edm::PtrVector<recob::Hit> hit;
   
  edm::PtrVector<recob::Cluster> clusIn;
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusIn.push_back(cluster);
    }


  int numberwires;
  double numbertimesamples;

  unsigned int channel,channel2,plane,plane2,wire,wire2;
  
      

      //std::cout<<"vertex0 "<<vertex[0]<<std::endl;
      double drifttick=(vertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./.198);
      

      

	  
	  
	  
	 // std::cout<<"vertex transform collection "<<wire2<<" "<<drifttick<<std::endl;
      double energyofhits_ind=0.;
      double energyofhits_col=0.;
      std::vector<float> hit_ind;
      std::vector<float> hit_col;
      hit_ind.clear();
      hit_col.clear();
  for(int p = 0; p < geom->Nplanes(); p++) 
    {
    
      if(p==0)
      vertex[0]=.3;//force time coordinate to be closer to induction plane 
      else
      vertex[0]=-.3;//force time coordinate to be closer to collection plane
     
	  channel2 = geom->NearestChannel(vertex);
	  geom->ChannelToWire(channel2,plane2,wire2); 
     std::cout<<"channel2 "<<channel2<<std::endl;
      edm::PtrVector<recob::Hit> vHits;
      edm::PtrVectorItr<recob::Cluster> clusterIter = clusIn.begin();
      hit.clear();
      cHits.clear();      
      
      
      while(clusterIter!= clusIn.end() ) {
	cHits = (*clusterIter)->Hits(p);
	if(cHits.size() > 0)
      for(int i = 0; i < cHits.size(); i++)
      hit.push_back(cHits[i]);
      
	clusterIter++;  
      } 
      if(hit.size() == 0) 
        continue;

      numberwires=geom->Nwires(0);
      numbertimesamples=hit[0]->Wire()->fSignal.size();
      

	  
      for(unsigned int i=0;i < hit.size(); i++)
	 {
	  channel=hit[i]->Wire()->RawDigit()->Channel();
	  geom->ChannelToWire(channel,plane,wire);   

	  
	  if(plane!=p)
	  continue;

	  //units are cm:
	  // if(sqrt(pow(TMath::Abs(vertexcolwire-wire)*.0743,2)+pow(TMath::Abs(vertexcoltime-hit[i]->CrossingTime()),2)) < fActivityRadius) 
	  // std::cout<<"wires "<<wire<<" "<<wire2<<" "<<hit[i]->CrossingTime()<<" "<<drifttick<<std::endl;
// 	  
// 	  if(larp->DriftVelocity(larp->Efield(),larp->Temperature())*.198*sqrt(pow(TMath::Abs((int)(wire2-wire))*(1/(larp->DriftVelocity(larp->Efield(),larp->Temperature())*(1/.4)*.198)),2)+pow(TMath::Abs(drifttick-hit[i]->CrossingTime()),2)) < fActivityRadius) 
// 	  {
	  
	  if(TMath::Abs((int)(wire2-wire))<=fActivityRadius&&TMath::Abs(drifttick-hit[i]->CrossingTime())*.0743<=fActivityRadius)
	  {
	  std::cout<<"drifttick "<<drifttick<<" hit crossing time "<<hit[i]->CrossingTime()<<std::endl;
	  
	  hitamplitude=hit[i]->MIPs();
	  //elifetime_factor=TMath::Exp((-hit[i]->CrossingTime()*.198)/electronlifetime); 
	  elifetime_factor=1.;
	  //elifetime_factor=1;
	  //6241.5 electrons/fC
	  

	  if(p==0)	  {energyofhits_ind+=(6241.5*fCalibration_factor*fWorkfunction_factor/fRecombination_factor)*(1/elifetime_factor)*hitamplitude*.000001*3.2;
	   hit_ind.push_back(hit[i]->MIPs());
	  }
	  if(p==1){energyofhits_col+=(6241.5*fCalibration_factor*fWorkfunction_factor/fRecombination_factor)*(1/elifetime_factor)*hitamplitude*.000001;
	  hit_col.push_back(hit[i]->MIPs());
	  }
	   if(p==0)
	   
  	   std::cout<<p<<" "<<wire<<" "<<wire2<<" "<<larp->DriftVelocity(larp->Efield(),larp->Temperature())*.198*sqrt(pow(TMath::Abs((int)(wire2-wire))*(1/(larp->DriftVelocity(larp->Efield(),larp->Temperature())*(1/.4)*.198)),2))<<" "<<hit[i]->CrossingTime()<<" "<<hitamplitude*3.2<<" "<<energyofhits_ind<<std::endl;
  	   
  	   if(p==1)
  	   std::cout<<p<<" "<<wire<<" "<<wire2<<" "<<larp->DriftVelocity(larp->Efield(),larp->Temperature())*.198*sqrt(pow(TMath::Abs((int)(wire2-wire))*(1/(larp->DriftVelocity(larp->Efield(),larp->Temperature())*(1/.4)*.198)),2))<<" "<<hit[i]->CrossingTime()<<" "<<hitamplitude<<" "<<energyofhits_col<<std::endl;
  	   
	  }
	  
	  
	}
	
	
	// double A=TMath::Abs(energyofhits_ind-energyofhits_col)/sqrt(pow(energyofhits_ind,2)+pow(energyofhits_col,2));
	
	double A=TMath::Abs(energyofhits_ind-energyofhits_col)/sqrt(pow(energyofhits_ind,2)+pow(energyofhits_col,2)-(energyofhits_ind*energyofhits_col));
	
	double B=1-A;

	double averageen=(energyofhits_ind+energyofhits_col)/2.;
	double reduceden=(energyofhits_ind*energyofhits_col)/(energyofhits_ind+energyofhits_col);
	
	
	double finalenergy=((A*reduceden)+(B*averageen))/(A+B);
	if(p==1)
	{
    std::cout<<"Energy2 "<<Energy<<" MeV "<<finalenergy<<" "<<std::endl;
    fEfficiency->Fill(finalenergy,finalenergy/Energy);
    fIndEfficiency->Fill(energyofhits_ind,energyofhits_ind/Energy);
    fColEfficiency->Fill(energyofhits_col,energyofhits_col/Energy);
    }
    //Energy=0.;
    
    
//       for(int vertexnum=0;vertexnum<vertexListHandle->size();vertexnum++)
// 	{
// 
// 
// 	}

      hit.clear();
      if(clusterIter!=clusIn.end()) clusterIter++;


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
    
    std::cout<<"size3 "<<size3<<" "<<size2<<" "<<size<<std::endl;
    for(int i=0;i<size3;i++)
    findcol->Fill(hit_ind[i],hit_col[i]);

}



