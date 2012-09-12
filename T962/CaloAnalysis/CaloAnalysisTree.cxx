////////////////////////////////////////////////////////////////////////
//
// Create a TTree for analysis
// This module is meant to analyze single protons. It compares how well our 
// reco is doing with reconstruction of track length as compared to MC truth. 
// It also evaluates reconstructed kinetic energy among other things.
//
// \author kinga.partyka@yale.edu
// 
////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>


#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/SubRun.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 


#include "T962/CaloAnalysis/CaloAnalysisTree.h"
#include "T962/T962_Objects/MINOS.h"
#include "T962/T962_Objects/ScanInfo.h"
#include "Geometry/geo.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "RecoBase/recobase.h"
#include "RawData/RawDigit.h"
#include "SummaryData/summary.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "MCCheater/BackTracker.h"

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

//-------------------------------------------------
t962::CaloAnalysisTree::CaloAnalysisTree(fhicl::ParameterSet const& pset) : 
  
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel")        ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fCalDataModuleLabel       (pset.get< std::string >("CalDataModuleLabel")      ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")        ),
  fEndPoint2DModuleLabel    (pset.get< std::string >("EndPoint2DModuleLabel")   ),
  fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel")       ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel")          ),
  fvertextrackWindow        (pset.get< double >("vertextrackWindow")            ),
  fvertexclusterWindow      (pset.get< double >("vertexclusterWindow")          ),
  fboundaryWindow           (pset.get< double >("boundaryWindow")               ),
  ntracks_reco(400),
  no_geant_particles(400)
 
  
 
  
{
}

//-------------------------------------------------
t962::CaloAnalysisTree::~CaloAnalysisTree()
{
delete twodvtx_w_reco;
delete twodvtx_t_reco;
delete twodvtx_w_truth;
delete twodvtx_t_truth;
delete pdg;
delete Eng;
delete Eng_at_endtrack;
delete Px;
 delete Py;
 delete Pz;
 delete StartPointx;
 delete StartPointy;
 delete StartPointz;
 delete EndPointx;
 delete EndPointy;
 delete EndPointz;
 delete NumberDaughters;
 delete trk_length_truth;
 delete trk_length_straight_line_truth;
 delete trk_length_reco;
 delete process_primary;
 delete Kin_Eng_reco;

 
}

void t962::CaloAnalysisTree::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("anatree","analysis tree");
  //make quick plots --------------------------------------------------
  diff_length_reco_truth = tfs->make<TH1F>("diff_length_reco_truth","trk_length_reco-trk_length_truth [cm]",200,-10.,10.);
  diff_KE_reco_truth = tfs->make<TH1F>("diff_KE_reco_truth","trk_KE_reco-trk_KE_truth [MeV]",100,-50.,50.);
  diff_length_truth1_vs_truth2 = tfs->make<TH1F>("diff_length_truth1_vs_truth2","trk_length_straight_line_truth-trk_length_truth [cm]",200,-10.,10.);
   
   diff_length_reco_truth_1geant = tfs->make<TH1F>("diff_length_reco_truth_1geant","trk_length_reco-trk_length_truth [cm] for 1 geant particle events",200,-10.,10.);
  diff_KE_reco_truth_1geant = tfs->make<TH1F>("diff_KE_reco_truth_1geant","trk_KE_reco-trk_KE_truth [MeV] for 1 geant particle events",100,-50.,50.);
  diff_length_truth1_vs_truth2_1geant = tfs->make<TH1F>("diff_length_truth1_vs_truth2_1geant","trk_length_straight_line_truth-trk_length_truth [cm] for 1 geant particle events",200,-10.,10.);
   
   
   
   
   
  //-------------------------------------------------------------------
  twodvtx_w_reco= new double[2];
  twodvtx_t_reco= new double[2];
  twodvtx_w_truth= new double[2];
  twodvtx_t_truth= new double[2];
  
  pdg= new int[no_geant_particles];
  Eng= new double[no_geant_particles];
  Eng_at_endtrack= new double[no_geant_particles];
  Px= new double[no_geant_particles];
  Py= new double[no_geant_particles];
  Pz= new double[no_geant_particles];
  StartPointx= new double[no_geant_particles];
  StartPointy= new double[no_geant_particles];
  StartPointz= new double[no_geant_particles];
  EndPointx= new double[no_geant_particles];
  EndPointy= new double[no_geant_particles];
  EndPointz= new double[no_geant_particles];
  NumberDaughters= new int[no_geant_particles];
  process_primary= new int[no_geant_particles];
 
  trk_length_truth = new double[no_geant_particles];
  trk_length_straight_line_truth = new double[no_geant_particles];
  trk_length_reco = new double[ntracks_reco];
  Kin_Eng_reco = new double[ntracks_reco];
 
  
 
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("pot",&pot,"pot/D");
  fTree->Branch("isdata",&isdata,"isdata/I");
  fTree->Branch("vtxx_reco",&vtxx_reco,"vtxx_reco/D");
  fTree->Branch("vtxy_reco",&vtxy_reco,"vtxy_reco/D");
  fTree->Branch("vtxz_reco",&vtxz_reco,"vtxz_reco/D");
  
  fTree->Branch("twodvtx_w_reco", twodvtx_w_reco, "twodvtx_w_reco[2]/D");
  fTree->Branch("twodvtx_t_reco", twodvtx_t_reco, "twodvtx_t_reco[2]/D");
  fTree->Branch("twodvtx_w_truth", twodvtx_w_truth, "twodvtx_w_truth[2]/D");
  fTree->Branch("twodvtx_t_truth", twodvtx_t_truth, "twodvtx_t_truth[2]/D");
  
  fTree->Branch("ntracks_reco",&ntracks_reco,"ntracks_reco/I");
  fTree->Branch("trackstart_dcosx_reco",&trackstart_dcosx_reco, "trackstart_dcosx_reco/D");
  fTree->Branch("trackstart_dcosy_reco",&trackstart_dcosy_reco, "trackstart_dcosy_reco/D");
  fTree->Branch("trackstart_dcosz_reco",&trackstart_dcosz_reco, "trackstart_dcosz_reco/D");
  fTree->Branch("trackexit_dcosx_reco",&trackexit_dcosx_reco, "trackexit_dcosx_reco/D");
  fTree->Branch("trackexit_dcosy_reco",&trackexit_dcosy_reco, "trackexit_dcosy_reco/D");
  fTree->Branch("trackexit_dcosz_reco",&trackexit_dcosz_reco, "trackexit_dcosz_reco/D");
  fTree->Branch("trackstart_x_reco",&trackstart_x_reco, "trackstart_x_reco/D");
  fTree->Branch("trackstart_y_reco",&trackstart_y_reco, "trackstart_y_reco/D");
  fTree->Branch("trackstart_z_reco",&trackstart_z_reco, "trackstart_z_reco/D");
  fTree->Branch("trackexit_x_reco",&trackexit_x_reco, "trackexit_x_reco/D");
  fTree->Branch("trackexit_y_reco",&trackexit_y_reco, "trackexit_y_reco/D");
  fTree->Branch("trackexit_z_reco",&trackexit_z_reco, "trackexit_z_reco/D");    
  
fTree->Branch("TrkPitchC", &fTrkPitchC, "TrkPitchC/F");
fTree->Branch("nhitsCOL",&fnhitsCOL,"nhitsCOL/I");



    //......................................................
// from geant4:

  fTree->Branch("no_geant_particles",&no_geant_particles,"no_geant_particles/I");
  fTree->Branch("pdg",pdg,"pdg[no_geant_particles]/I");
  fTree->Branch("Eng",Eng,"Eng[no_geant_particles]/D");
  fTree->Branch("Eng_at_endtrack",Eng_at_endtrack,"Eng_at_endtrack[no_geant_particles]/D");
  fTree->Branch("Px",Px,"Px[no_geant_particles]/D");
  fTree->Branch("Py",Py,"Py[no_geant_particles]/D");
  fTree->Branch("Pz",Pz,"Pz[no_geant_particles]/D");
  fTree->Branch("StartPointx",StartPointx,"StartPointx[no_geant_particles]/D");
  fTree->Branch("StartPointy",StartPointy,"StartPointy[no_geant_particles]/D");
  fTree->Branch("StartPointz",StartPointz,"StartPointz[no_geant_particles]/D");
  fTree->Branch("EndPointx",EndPointx,"EndPointx[no_geant_particles]/D");
  fTree->Branch("EndPointy",EndPointy,"EndPointy[no_geant_particles]/D");
  fTree->Branch("EndPointz",EndPointz,"EndPointz[no_geant_particles]/D");
  fTree->Branch("NumberDaughters",NumberDaughters,"NumberDaughters[no_geant_particles]/I");
  
  fTree->Branch("process_primary",process_primary,"process_primary[no_geant_particles]/I");
 fTree->Branch("trk_length_truth",trk_length_truth,"trk_length_truth[no_geant_particles]/D");
  fTree->Branch("trk_length_straight_line_truth",trk_length_straight_line_truth,"trk_length_straight_line_truth[no_geant_particles]/D");
  fTree->Branch("trk_length_reco",trk_length_reco,"trk_length_reco[ntracks_reco]/D");
  fTree->Branch("Kin_Eng_reco",Kin_Eng_reco,"Kin_Eng_reco[ntracks_reco]/D");
  fTree->Branch("Kin_Eng_truth",&Kin_Eng_truth,"Kin_Eng_truth/D");
  fTree->Branch("KE_sum_",&KE_sum_,"KE_sum_/D");
  
 //..................................

}


void t962::CaloAnalysisTree::beginSubRun(const art::SubRun& sr)
{

art::Handle< sumdata::POTSummary > potListHandle;
sr.getByLabel(fPOTModuleLabel,potListHandle);

if(sr.getByLabel(fPOTModuleLabel,potListHandle))
pot=potListHandle->totpot;
else
pot=0.;

}




void t962::CaloAnalysisTree::analyze(const art::Event& evt)
{
std::cout<<" IN *** MY *** CaloAnalysisTree *** ----------------"<<std::endl;
  ResetVars();

  run = evt.run();
  event = evt.id().event();

  if (evt.isRealData()){
    isdata = 1;
  }
  else isdata = 0;

  
  
 
 
  art::Handle< std::vector<recob::Track> > trackListHandle;
  evt.getByLabel(fTrackModuleLabel,trackListHandle);
  art::Handle< std::vector<recob::EndPoint2D> > endpointListHandle;
  evt.getByLabel(fEndPoint2DModuleLabel,endpointListHandle);
  art::Handle< std::vector<recob::Vertex> > vertexListHandle;
  evt.getByLabel(fVertexModuleLabel,vertexListHandle);
 

 


  art::PtrVector<recob::Track> tracklist;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle))
  for (unsigned int i = 0; i < trackListHandle->size(); ++i){
    art::Ptr<recob::Track> trackHolder(trackListHandle,i);
    tracklist.push_back(trackHolder);
  }

  art::PtrVector<recob::EndPoint2D> endpointlist;
  if(evt.getByLabel(fEndPoint2DModuleLabel,endpointListHandle))
    for (unsigned int i = 0; i < endpointListHandle->size(); ++i){
      art::Ptr<recob::EndPoint2D> endpointHolder(endpointListHandle,i);
      endpointlist.push_back(endpointHolder);
    }

  art::PtrVector<recob::Vertex> vertexlist;
  if(evt.getByLabel(fVertexModuleLabel,vertexListHandle))
  for (unsigned int i = 0; i < vertexListHandle->size(); ++i){
    art::Ptr<recob::Vertex> vertexHolder(vertexListHandle,i);
    vertexlist.push_back(vertexHolder);
  }

 
 
  
 
  art::ServiceHandle<geo::Geometry> geom;  
  art::ServiceHandle<util::LArProperties> larp;
   art::ServiceHandle<util::DetectorProperties> detprop;
  
 //--------------------------------------------------------------------
 //      FIGURE OUT THE TRACK LENGTH- SAME AS IN CALORIMETRY 
 //--------------------------------------------------------------------   
   fnhitsCOL = 0;
   double Trk_Length=0;
   int npC = 0;
   float time;
   double Kin_En = 0.;
   
    // Electronic calibration factor to convert from ADC to electrons
   double fElectronsToADC = detprop->ElectronsToADC();
  
   size_t trkCtr = 0;

   art::FindManyP<recob::Hit> fmh(trackListHandle, evt, fTrackModuleLabel);
   art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);

   for(art::PtrVector<recob::Track>::const_iterator trkIter = tracklist.begin(); trkIter != tracklist.end();  trkIter++){ 
     
     fTrkPitchC = (*trkIter)->ProjectedLength(geo::kV);

     // get the hits for this view
     std::vector< art::Ptr<recob::Hit> > hits = fmh.at(trkCtr);

     art::PtrVector<recob::Hit> hitlist;
     for(size_t h = 0; h < hits.size(); ++h) 
       if(hits[h]->View() == geo::kV) hitlist.push_back(hits[h]);
      
     //  fnhitsCOL = npC + hitlist.size();
//       
//       double dEdx_Coll_Tot = 0.;
//       double dEdx_Coll_Tot_5cm = 0.;

     Kin_En = 0.;
     double *xyz_previous = new double[3]; 
     
     std::cout<<" TrkPitchC="<<fTrkPitchC<<" nhitsCOL="<<fnhitsCOL<<std::endl;
  
     std::vector< art::Ptr<recob::SpacePoint> > sps = fmsp.at(trkCtr);
     //int nsp = sps.size();
     
     //loop over Collection hits
     for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
       //recover the Hit
       time = (*hitIter)->PeakTime() ;
       // 	stime = (*hitIter)->StartTime() ;
       // 	etime = (*hitIter)->EndTime();            
       // 	art::Ptr<recob::Wire> theWire = (*hitIter)->Wire();
       // 	channel = theWire->RawDigit()->Channel();
       // 	geom->ChannelToWire(channel,cstat,tpc,plane,wire);
       
       double MIPs   = (*hitIter)->Charge(true);   // in ADC
       double dQdx   = MIPs/fTrkPitchC;           // in ADC/cm
       // 	fdQdx_Coll->Fill(dQdx);
       
       double dQdx_e = dQdx/fElectronsToADC;  // Conversion from ADC/cm to e/cm
       
       dQdx_e *= LifetimeCorrection(time);   // Lifetime Correction (dQdx_e in e/cm)
       
       double dEdx = larp->BirksCorrection(dQdx_e);   // Correction for charge quenching (Recombination) dEdx in MeV/cm
       // 
       // 	fdEdx_Coll->Fill(dEdx);
       // 	fbirk->Fill(dQdx_e,dEdx);
       //fdEdx_Coll_vsXZangle->Fill(dEdx,TMath::ATan(trackCosStart[0]/trackCosStart[2]));
       
       //dEdx_Coll_Tot = dEdx_Coll_Tot + dEdx;
       Kin_En = Kin_En + dEdx * fTrkPitchC;
       
       
       //std::cout<<" Collection: npC, wire, time, ADC, dQdx (e/cm), dEdx (MeV/cm)"<<npC<<" "<<wire<<" "<<time<<" "<<MIPs<<" "<<dQdx_e<<" "<<dEdx<<std::endl;
       // 
       //  	fwireCOL[npC] = (int)wire;       
       // 	ftimeCOL[npC] = (double)time;
       // 	fstimeCOL[npC] = (double)stime;
       // 	fetimeCOL[npC] = (double)etime;
       // 	fMIPsCOL[npC] =  (double)MIPs;
       // 	fdEdxCOL[npC] = (double)dEdx;
       
       npC++;
       
       art::FindManyP<recob::Hit> fmsph(sps, evt, fTrackModuleLabel);
       
       // dE/dx vs Residual Range plot
       for(size_t s = 0; s < sps.size(); ++s) { 
	 // space point associations were made in fTrackModuleLabel
	 std::vector< art::Ptr<recob::Hit> > sphitsAll = fmsph.at(s);
	 std::vector< art::Ptr<recob::Hit> > sphits;
	 for(size_t sh = 0; sh < sphitsAll.size(); ++sh)
	   if(sphitsAll[sh]->View() == geo::kV) sphits.push_back(sphitsAll[sh]); 
	 
	 if(sphits.size()==0) continue;
	 
	 if((*hitIter)->Channel()!=sphits[0]->Channel()) continue;
	 //protect against multiple hits on the same wire in the cluster.
	 if((*hitIter)->PeakTime()!=sphits[0]->PeakTime()) continue;

	  const double *xyz = new double[3];
	  xyz = sps[s]->XYZ();
	  std::vector<double> larStart, larEnd;
	  (*trkIter)->Extent(larStart,larEnd);//put xyz coordinates at begin/end of track into vectors(?)

	//   double res_range = sqrt((larEnd[0]-xyz[0])*(larEnd[0]-xyz[0]) +
// 				  (larEnd[1]-xyz[1])*(larEnd[1]-xyz[1]) +
// 				  (larEnd[2]-xyz[2])*(larEnd[2]-xyz[2]));

	  // Track Length from I - last points
// 	  Trk_Length = sqrt((larStart[0]-larEnd[0])*(larStart[0]-larEnd[0]) +
// 	  		       (larStart[1]-larEnd[1])*(larStart[1]-larEnd[1]) +
// 	  		       (larStart[2]-larEnd[2])*(larStart[2]-larEnd[2]));

	  if(s == 0) 
	    Trk_Length = Trk_Length  + sqrt((larStart[0]-xyz[0])*(larStart[0]-xyz[0]) +
					    (larStart[1]-xyz[1])*(larStart[1]-xyz[1]) +
					    (larStart[2]-xyz[2])*(larStart[2]-xyz[2]));
	  else 
	    Trk_Length = Trk_Length  + sqrt((xyz_previous[0]-xyz[0])*(xyz_previous[0]-xyz[0]) +
					    (xyz_previous[1]-xyz[1])*(xyz_previous[1]-xyz[1]) +
					    (xyz_previous[2]-xyz[2])*(xyz_previous[2]-xyz[2]));

	 //  channel = theWire->RawDigit()->Channel();
	  // 	  geom->ChannelToWire(channel,cstat,tpc,plane,wire);
	  
	  // std::cout<<" Space Points: Coll. wire, time, charge, Ind. channel, time, charge "<<wire<<" "<<time<<" "<<MIPs<<" "<<sphitsI[0]->Channel()<<" "<<sphitsI[0]->PeakTime()<<" "<<sphitsI[0]->Charge(true)<<std::endl;
	  xyz_previous[0] = xyz[0];
	  xyz_previous[1] = xyz[1];
	  xyz_previous[2] = xyz[2];

           
	  /// First and last hits not included in dE/dx vs. residual range plots (because the actual extension of the track within the wire pitch is unknown)
	
	        
	  // npC++;
               
       }
     }// end of loop over the (geo::kV) cluster hits
      
      
      //fKinetic_En->Fill(Kin_En);

      

      std::cout<<"|-*                    Track Length="<<Trk_Length<<" cm"<<std::endl;
     // std::cout<<"  |-* Collection View Calorimetric Reco"<<std::endl;
      //std::cout<<"    |-* Hits="<<fnhitsCOL<<std::endl;
      //std::cout<<"    |-* <dE/dx>="<<dEdx_Coll_Mean<<" MeV/cm"<<std::endl;
      std::cout<<"    |-* Kinetic Energy deposited in LAr="<<Kin_En<<" MeV"<<std::endl;
     
 
    trk_length_reco[trkIter-tracklist.begin()]=Trk_Length;
    Kin_Eng_reco[trkIter-tracklist.begin()]=Kin_En;
    
    ++trkCtr;
   }//track list
 //--------------------------------------------------------------------
  //vertex information
  if(vertexlist.size())
  {
    double vtxxyz[3];
    vertexlist[0]->XYZ(vtxxyz);
    vtxx_reco = vtxxyz[0];
    vtxy_reco = vtxxyz[1];
    vtxz_reco = vtxxyz[2];
  }
  
  
  // 2d vertex information
 
  
  for (unsigned int j = 0; j<endpointlist.size();j++){
 std::cout<<"j="<<j<<" W_VERTEX_RECO= "<<endpointlist[j]->WireNum()<<" T_VERTEX_RECO= "<<endpointlist[j]->DriftTime()<<std::endl;
          twodvtx_w_reco[j]=endpointlist[j]->WireNum();
          twodvtx_t_reco[j]=endpointlist[j]->DriftTime();
          
    }
  
  //track information
     ntracks_reco=tracklist.size();
     double larStart[3];
     double larEnd[3];
 	 std::vector<double> trackStart;
  	 std::vector<double> trackEnd;
     trackStart.clear();
     trackEnd.clear();
    
   
     //grab information about where track started and ended and the dcos at those points
     trackStart.clear();
     trackEnd.clear();
     memset(larStart, 0, 3);
     memset(larEnd, 0, 3);
     for(unsigned int i=0; i<tracklist.size();++i){
     
    

       tracklist[i]->Direction(larStart,larEnd);
       tracklist[i]->Extent(trackStart,trackEnd);  
                     
       trackstart_dcosx_reco = larStart[0];
       trackstart_dcosy_reco = larStart[1];
       trackstart_dcosz_reco = larStart[2];       
       trackexit_dcosx_reco = larEnd[0];
       trackexit_dcosy_reco = larEnd[1];
       trackexit_dcosz_reco = larEnd[2];
           
       trackstart_x_reco=trackStart[0];
       trackstart_y_reco=trackStart[1];
       trackstart_z_reco=trackStart[2];
       trackexit_x_reco=trackEnd[0];
       trackexit_y_reco=trackEnd[1];
       trackexit_z_reco=trackEnd[2];  
      }

   
     
  
  
   //--------------------------------------------------------------//
 //        NOW I WILL GET INFO FROM GEANT4
 //       
 //--------------------------------------------------------------//
 double KE=0;
 double mass_proton=0.938;
 double mass_neutron=0.939;
 double KE_sum=0;
 
 
  if (!isdata){ 
 
 
  ////////////////////////////////////////////////////////////////////////
 
 //     GET DEPOSITED KINETIC ENERGY
 
 ///////////////////////////
double tot_energy=0.;
sim::ParticleList plist = sim::SimListUtils::GetParticleList(evt,fLArG4ModuleLabel);
art::ServiceHandle<geo::Geometry> geom;

    unsigned int cstat;    //hit cryostat number 
    unsigned int tpc;    //hit tpc number 
    unsigned int wire;   //hit wire number 
    unsigned int plane;  //hit plane number

   std::vector<const sim::SimChannel*> sccol;
   evt.getView(fLArG4ModuleLabel, sccol);

// loop over all sim::SimChannels in the event and make sure there are no
   // sim::IDEs with trackID values that are not in the sim::ParticleList
    for(size_t sc = 0; sc < sccol.size(); ++sc){
       //check whether the channel is collection - to make sure you sum only one wire plane basically.
       geom->ChannelToWire(sccol[sc]->Channel(),cstat,tpc,plane,wire);
       if(geom->Plane(plane).SignalType()!= geo::kCollection)
           continue;

       const std::map<unsigned short, std::vector<sim::IDE> >& tdcidemap = sccol[sc]->TDCIDEMap();
       std::map<unsigned short, std::vector<sim::IDE> >::const_iterator mapitr;
       for(mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
         const std::vector<sim::IDE> idevec = (*mapitr).second;
       for(size_t iv = 0; iv < idevec.size(); ++iv){
         if(plist.find( idevec[iv].trackID ) != plist.end()  //found the right track ID?
            && idevec[iv].trackID != sim::NoParticleId)  {
           // here you can put an if whether  idevec[iv].trackID== mother particle TrackID
            tot_energy+=idevec[iv].energy;
        }

           //  std::cout << "found track @? " << idevec[iv].trackID << " "<< idevec[iv].numElectrons << " sc: " << sc << std::endl;
         }
       }
     }

std::cout<<"tot_energy= "<<tot_energy<<std::endl;
Kin_Eng_truth=tot_energy;
////////////////////////////////////////
 
 
 art::Handle< std::vector<sim::Particle> > geant_list;
   if(evt.getByLabel (fLArG4ModuleLabel,geant_list));
 
  art::PtrVector<sim::Particle> geant_part;
   for (unsigned int ii = 0; ii <  geant_list->size(); ++ii)
    {
      art::Ptr<sim::Particle> p(geant_list,ii);
      geant_part.push_back(p);
    } 
    std::cout<<"No of geant part= "<<geant_list->size()<<std::endl;
 std::string pri ("primary");

  no_geant_particles=geant_part.size();
  
 std::cout<<"Geant4 list: "<<std::endl;
 
 for( unsigned int i = 0; i < geant_part.size(); ++i ){
 
  pdg[i]=geant_part[i]->PdgCode();
    
    Eng[i]=geant_part[i]->E();
    
    //energy at the last point of its trajectory:
    //for(unsigned int pt=0; pt<geant_part[i]->NumberTrajectoryPoints()-1; pt++){
    Eng_at_endtrack[i]=geant_part[i]->E(geant_part[i]->NumberTrajectoryPoints()-1);
    
    //}
    Px[i]=geant_part[i]->Px();
   
    Py[i]=geant_part[i]->Py();
    Pz[i]=geant_part[i]->Pz();
    
   StartPointx[i]=geant_part[i]->Vx();
   StartPointy[i]=geant_part[i]->Vy();
   StartPointz[i]=geant_part[i]->Vz();
   EndPointx[i]=geant_part[i]->EndPoint()[0];
   EndPointy[i]=geant_part[i]->EndPoint()[1];
   EndPointz[i]=geant_part[i]->EndPoint()[2];
   
   NumberDaughters[i]=geant_part[i]->NumberDaughters();
   
     if(geant_part[i]->Process()==pri){
       process_primary[i]=1;
       
       if(geant_part.size()==1) KE_sum=geant_part[i]->E()-mass_proton;
     } //if primary
    else { 
    process_primary[i]=0;
    
    if(geant_part[i]->PdgCode()==2212) KE=geant_part[i]->E()-mass_proton;
    if(geant_part[i]->PdgCode()==2112) KE=geant_part[i]->E()-mass_neutron;
    if(geant_part[i]->PdgCode()==22) KE=geant_part[i]->E();
    
    KE_sum+=KE;
    
    } 
    
    

     
     
     std::cout<<"KE_sum= "<<KE_sum*1000<<" MeV"<<std::endl;
     KE_sum_=KE_sum*1000;
     
    //-----------------------------------------------
    //  CALCULATE TRUE LENGTH OF A TRACK:
    double length=0.;
    //std::cout<<"NumberTrajectoryPoints= "<<geant_part[i]->NumberTrajectoryPoints()<<std::endl;
    for(unsigned int pt=0; pt<geant_part[i]->NumberTrajectoryPoints()-1; pt++){
    
    length+=sqrt((geant_part[i]->Position(pt+1).X()-geant_part[i]->Position(pt).X())*(geant_part[i]->Position(pt+1).X()-geant_part[i]->Position(pt).X()) + (geant_part[i]->Position(pt+1).Y()-geant_part[i]->Position(pt).Y())*(geant_part[i]->Position(pt+1).Y()-geant_part[i]->Position(pt).Y()) + (geant_part[i]->Position(pt+1).Z()-geant_part[i]->Position(pt).Z())*(geant_part[i]->Position(pt+1).Z()-geant_part[i]->Position(pt).Z()) );
    
    }
    trk_length_truth[i]=length;
     //-----------------------------------------------
     //  CALCULATE TRUE LENGTH OF A TRACK ASSUMING STRAIGHT LINE:
   trk_length_straight_line_truth[i]=sqrt((EndPointx[i]-StartPointx[i])*(EndPointx[i]-StartPointx[i]) + (EndPointy[i]-StartPointy[i])*(EndPointy[i]-StartPointy[i])+ (EndPointz[i]-StartPointz[i])*(EndPointz[i]-StartPointz[i]));
   
 // std::cout<<"pdg= "<<geant_part[i]->PdgCode()<<" trackId= "<<geant_part[i]->TrackId()<<" mother= "<<geant_part[i]->Mother()<<" NumberDaughters()= "<<geant_part[i]->NumberDaughters()<<" process= "<<geant_part[i]->Process()<<std::endl;
     
     
   std::cout<<"pdg= "<<geant_part[i]->PdgCode()<<" Process= "<<geant_part[i]->Process()<<" E= "<<geant_part[i]->E()<<" E last= "<<geant_part[i]->E(geant_part[i]->NumberTrajectoryPoints()-1)<<" P= "<<geant_part[i]->P()<<" "<<sqrt(geant_part[i]->Px()*geant_part[i]->Px() + geant_part[i]->Py()*geant_part[i]->Py()+ geant_part[i]->Pz()*geant_part[i]->Pz())<<" trk_length_truth= "<<trk_length_truth[i]<<" mother= "<<geant_part[i]->Mother()<<" trackID=" <<geant_part[i]->TrackId()<<std::endl;
   
   
   
   
    
    } //geant particles
    
    //Kin_Eng_truth=KE_sum*1000; //KE in MeV
    
    
      for(unsigned int k=0; k<tracklist.size();k++){
      std::cout<<std::endl;
      std::cout<<"________________________________"<<std::endl;
      std::cout<<"trk_length_truth = "<<trk_length_truth[k]<<" (assuming straight line track in truth= "<<trk_length_straight_line_truth[k]<<std::endl;
      std::cout<<"trk_length_reco = "<<trk_length_reco[k]<<std::endl;
      std::cout<<std::endl;
       std::cout<<"Kin_Eng (from summing geant p)= "<<KE_sum*1000<<" Kin_Eng_reco= "<<Kin_Eng_reco[k]<< "MeV"<<std::endl;
       std::cout<<"________________________________"<<std::endl;
       
       
       if(tracklist.size()==1 && trk_length_reco[k]!=0 && trk_length_truth[k]!=0){
        diff_length_reco_truth->Fill(trk_length_reco[k]-trk_length_truth[k]);
        diff_length_truth1_vs_truth2->Fill(trk_length_straight_line_truth[k]-trk_length_truth[k]);
        diff_KE_reco_truth->Fill(Kin_Eng_truth-Kin_Eng_reco[k]);
        
        if(geant_part.size()==1){
         diff_length_reco_truth_1geant->Fill(trk_length_reco[k]-trk_length_truth[k]);
        diff_length_truth1_vs_truth2_1geant->Fill(trk_length_straight_line_truth[k]-trk_length_truth[k]);
        diff_KE_reco_truth_1geant->Fill(Kin_Eng_truth-Kin_Eng_reco[k]);
        }
        }
      }

 
      //-------------------------------------------------------------------------
      //   FIGURE OUT WHAT PARTICLE IS IN EACH HIT, WANT TO FIND OUT PDG CODE OF
      //    SINGLE HITS AROUND PROTON TRACK
      //-------------------------------------------------------------------------
      art::ServiceHandle<cheat::BackTracker> bt;
      bt->SetEveIdCalculator(new sim::EmEveIdCalculator);
 
      art::Handle< std::vector<recob::Hit> > hitListHandle;
      evt.getByLabel(fHitsModuleLabel,hitListHandle);
      std::vector< art::Ptr<recob::Hit> > hits;
      art::fill_ptr_vector(hits, hitListHandle);
      
      std::vector< art::Ptr<recob::Hit> >::iterator itr = hits.begin();
      while(itr != hits.end()) {
	
	std::vector<cheat::TrackIDE> trackides = bt->HitToTrackID(*itr);
	std::vector<cheat::TrackIDE> eveides   = bt->HitToEveID(*itr);
		
	std::vector<cheat::TrackIDE>::iterator idesitr = trackides.begin();
	
	
	while( idesitr != trackides.end() ){
	  
	  

	  const sim::Particle* particle = bt->TrackIDToParticle((*idesitr).trackID);
		  
	  int pdg = particle->PdgCode();
	  
	  geom->ChannelToWire((*itr)->Channel(),cstat,tpc,plane,wire);
	  
	  std::cout<<"plane= "<<plane<<" w= "<<wire<<" time= "
		   <<(*itr)->PeakTime()<< " pdg= "<<pdg<<std::endl;
	  
	  //std::cout<<"pdg= "<<pdg<<std::endl;
   
	  idesitr++;
	  
	} //trackIDs  
	
	itr++;
      }
      
      
      //-------------------------------------------------------------------------
      //-------------------------------------------------------------------------
      
  } //if MC

 
  fTree->Fill();
}

  //---------------------------------------------------------------- 
void t962::CaloAnalysisTree::ResetVars(){

  run = -99999;
  event = -99999;
  isdata = -99999;
  vtxx_reco = -99999;
  vtxy_reco  = -99999;
  vtxz_reco  = -99999;
  trackstart_x_reco = -99999;
  trackstart_y_reco = -99999;
  trackstart_z_reco = -99999;
  trackexit_x_reco = -99999;
  trackexit_y_reco = -99999;
  trackexit_z_reco = -99999;  
  trackstart_dcosx_reco = -99999;
  trackstart_dcosy_reco = -99999;
  trackstart_dcosz_reco = -99999;       
  trackexit_dcosx_reco = -99999;
  trackexit_dcosy_reco = -99999;
  trackexit_dcosz_reco = -99999;
  ntracks_reco=-999;
  no_geant_particles=-999;
  Kin_Eng_truth=-999;
  KE_sum_=-999;
  

}

//------------------------------------------------------------------------//
double t962::CaloAnalysisTree::LifetimeCorrection(float time){

   float t = time;

   art::ServiceHandle<util::LArProperties> LArProp;
   art::ServiceHandle<util::DetectorProperties> detprop;

   double timetick = detprop->SamplingRate()*1.e-3;    //time sample in microsec
   double presamplings = detprop->TriggerOffset();

   t -= presamplings;
   time = t * timetick;  //  (in microsec)

   double tau = LArProp->ElectronLifetime();

   double correction = exp(time/tau);
   return correction;
}

//--------------------------------------------------  
  
  







