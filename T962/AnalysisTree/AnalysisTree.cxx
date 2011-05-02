////////////////////////////////////////////////////////////////////////
//
// Create a TTree for analysis
//
// \author tjyang@fnal.gov
// 
// 
////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>

#include "TTree.h"

#include "art/Framework/Core/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 


#include "T962/AnalysisTree/AnalysisTree.h"
#include "T962/T962_Objects/MINOS.h"
#include "T962/T962_Objects/MINOSTrackMatch.h"
#include "Geometry/geo.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "RecoBase/recobase.h"
#include "RawData/RawDigit.h"
#include "Utilities/LArProperties.h"

 
//-------------------------------------------------
t962::AnalysisTree::AnalysisTree(fhicl::ParameterSet const& pset) : 
  
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel")        ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fCalDataModuleLabel       (pset.get< std::string >("CalDataModuleLabel")      ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel")      ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")        ),
  fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel")       ),
  fMINOSModuleLabel         (pset.get< std::string >("MINOSModuleLabel")        ),
  fTrackMatchModuleLabel    (pset.get< std::string >("TrackMatchModuleLabel")   )
{
}

//-------------------------------------------------
t962::AnalysisTree::~AnalysisTree()
{
}

void t962::AnalysisTree::beginJob()
{

  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("anatree","analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("isdata",&isdata,"isdata/I");
  fTree->Branch("vtxx",&vtxx,"vtxx/D");
  fTree->Branch("vtxy",&vtxy,"vtxy/D");
  fTree->Branch("vtxz",&vtxz,"vtxz/D");
  fTree->Branch("ntrks",&ntrks,"ntrks/I");
  fTree->Branch("nclusu",&nclusu,"nclusu/I");
  fTree->Branch("nclusv",&nclusv,"nclusv/I");
  fTree->Branch("nclusw",&nclusw,"nclusw/I");
  fTree->Branch("matched",&matched,"matched/I");
  fTree->Branch("mtrk_mom",&mtrk_mom,"mtrk_mom/D");
  fTree->Branch("mtrk_charge",&mtrk_charge,"mtrk_charge/D");
  fTree->Branch("mtrk_dcosx",&mtrk_dcosx,"mtrk_dcosx/D");
  fTree->Branch("mtrk_dcosy",&mtrk_dcosy,"mtrk_dcosy/D");
  fTree->Branch("mtrk_dcosz",&mtrk_dcosz,"mtrk_dcosz/D");
  fTree->Branch("inu",&inu,"inu/I");
  fTree->Branch("ccnc",&ccnc,"ccnc/I");
  fTree->Branch("mode",&mode,"mode/I");
  fTree->Branch("enu",&enu,"enu/D");
  fTree->Branch("Q2",&Q2,"Q2/D");
  fTree->Branch("W",&W,"W/D");
  fTree->Branch("nuvtxx",&nuvtxx,"nuvtxx/D");
  fTree->Branch("nuvtxy",&nuvtxy,"nuvtxy/D");
  fTree->Branch("nuvtxz",&nuvtxz,"nuvtxz/D");
  fTree->Branch("lep_mom",&lep_mom,"lep_mom/D");
  fTree->Branch("lep_dcosx",&lep_dcosx,"lep_dcosx/D");
  fTree->Branch("lep_dcosy",&lep_dcosy,"lep_dcosy/D");
  fTree->Branch("lep_dcosz",&lep_dcosz,"lep_dcosz/D");

}

void t962::AnalysisTree::analyze(const art::Event& evt)
{

  ResetVars();

  run = evt.run();
  event = evt.id().event();

  if (evt.isRealData()){
    isdata = 1;
  }
  else isdata = 0;
  
  art::Handle< std::vector<raw::RawDigit> > rdListHandle;
  evt.getByLabel(fDigitModuleLabel,rdListHandle);
  art::Handle< std::vector<sim::SimChannel> > scListHandle;
  evt.getByLabel(fDigitModuleLabel,scListHandle);
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitsModuleLabel,hitListHandle);
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);
  art::Handle< std::vector<recob::Wire> > wireListHandle;
  evt.getByLabel(fCalDataModuleLabel,wireListHandle);
  art::Handle< std::vector<recob::Track> > trackListHandle;
  evt.getByLabel(fTrackModuleLabel,trackListHandle);
  art::Handle< std::vector<recob::Vertex> > vertexListHandle;
  evt.getByLabel(fVertexModuleLabel,vertexListHandle);
  art::Handle< std::vector<t962::MINOS> > minosListHandle;
  evt.getByLabel(fMINOSModuleLabel,minosListHandle);
  art::Handle< std::vector<t962::MINOSTrackMatch> > trackmatchListHandle;
  evt.getByLabel(fTrackMatchModuleLabel,trackmatchListHandle);

  art::PtrVector<simb::MCTruth> mclist;
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    }

  art::PtrVector<recob::Cluster> clusterlist;
  for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      clusterlist.push_back(clusterHolder);
    }


  art::PtrVector<recob::Track> tracklist;
  for (unsigned int i = 0; i < trackListHandle->size(); ++i){
    art::Ptr<recob::Track> trackHolder(trackListHandle,i);
    tracklist.push_back(trackHolder);
  }

  art::PtrVector<recob::Vertex> vertexlist;
  for (unsigned int i = 0; i < vertexListHandle->size(); ++i){
    art::Ptr<recob::Vertex> vertexHolder(vertexListHandle,i);
    vertexlist.push_back(vertexHolder);
  }

  art::PtrVector<t962::MINOS> minoslist;
  for (unsigned int i = 0; i < minosListHandle->size(); i++){
    art::Ptr<t962::MINOS> minosHolder(minosListHandle,i);
    minoslist.push_back(minosHolder);
  }

  art::PtrVector<t962::MINOSTrackMatch> trackmatchlist;
  for (unsigned int i = 0; i < trackmatchListHandle->size(); i++){
    art::Ptr<t962::MINOSTrackMatch> trackmatchHolder(trackmatchListHandle,i);
    trackmatchlist.push_back(trackmatchHolder);
  }


  art::ServiceHandle<geo::Geometry> geom;  

  //vertex information
  double vtxxyz[3];
  vertexlist[0]->XYZ(vtxxyz);
  vtxx = vtxxyz[0];
  vtxy = vtxxyz[1];
  vtxz = vtxxyz[2];

  //cluster information
  nclusu = 0;
  nclusv = 0;
  nclusw = 0;
  for(unsigned int i=0; i<clusterlist.size();++i){
    
    switch(clusterlist[i]->View()){
    case geo::kU :
      nclusu ++;
      break;
    case geo::kV :
      nclusv ++;
      break;
    case geo::kW :
      nclusw ++;
      break;
    default :
      break;
    }
  }

  //track information
  ntrks = trackmatchlist.size();

  //matching information
  matched = trackmatchlist.size();
  if (matched>0){
    mtrk_mom = minoslist[0]->ftrkmom;
    mtrk_charge = minoslist[0]->fcharge;
    mtrk_dcosx = minoslist[0]->ftrkdcosx;
    mtrk_dcosy = minoslist[0]->ftrkdcosy;
    mtrk_dcosz = minoslist[0]->ftrkdcosz;
  }

  if (!isdata){//mc truth information
    inu = mclist[0]->GetNeutrino().Nu().PdgCode();
    ccnc = mclist[0]->GetNeutrino().CCNC();
    mode = mclist[0]->GetNeutrino().Mode();
    Q2 = mclist[0]->GetNeutrino().QSqr();
    W = mclist[0]->GetNeutrino().W();
    enu = mclist[0]->GetNeutrino().Nu().E();
    nuvtxx = mclist[0]->GetNeutrino().Nu().Vx();
    nuvtxy = mclist[0]->GetNeutrino().Nu().Vy();
    nuvtxz = mclist[0]->GetNeutrino().Nu().Vz();
    lep_mom = mclist[0]->GetNeutrino().Lepton().P();
    if (mclist[0]->GetNeutrino().Lepton().P()){
      lep_dcosx = mclist[0]->GetNeutrino().Lepton().Px()/mclist[0]->GetNeutrino().Lepton().P();
      lep_dcosy = mclist[0]->GetNeutrino().Lepton().Py()/mclist[0]->GetNeutrino().Lepton().P();
      lep_dcosz = mclist[0]->GetNeutrino().Lepton().Pz()/mclist[0]->GetNeutrino().Lepton().P();
    }
  }

  fTree->Fill();

}

  //---------------------------------------------------------------- 
void t962::AnalysisTree::ResetVars(){

  run = -99999;
  event = -99999;
  isdata = -99999;
  vtxx = -99999;
  vtxy = -99999;
  vtxz = -99999;
  ntrks = -99999;
  nclusu = -99999;
  nclusv = -99999;
  nclusw = -99999;
  matched = -99999;
  mtrk_mom = -99999;
  mtrk_charge = -99999;
  mtrk_dcosx = -99999;
  mtrk_dcosy = -99999;
  mtrk_dcosz = -99999;
  inu = -99999;
  ccnc = -99999;
  mode = -99999;
  enu = -99999;
  Q2 = -99999;
  W = -99999;
  nuvtxx = -99999;
  nuvtxy = -99999;
  nuvtxz = -99999;
  lep_mom = -99999;
  lep_dcosx = -99999;
  lep_dcosy = -99999;
  lep_dcosz = -99999;
}
