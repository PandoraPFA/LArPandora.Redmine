////////////////////////////////////////////////////////////////////////
//
// Create a TTree for analysis
//
// \author tjyang@fnal.gov
// \author joshua.spitz@yale.edu
// 
////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>

#include "TTree.h"

#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/SubRun.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Principal/View.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 


#include "T962/AnalysisTree/AnalysisTree.h"
#include "T962/T962_Objects/MINOS.h"
#include "T962/T962_Objects/MINOSTrackMatch.h"
#include "T962/T962_Objects/ScanInfo.h"
#include "Geometry/geo.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "RecoBase/recobase.h"
#include "AnalysisBase/anabase.h"
#include "RawData/RawDigit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "SummaryData/summary.h"

 
//-------------------------------------------------
t962::AnalysisTree::AnalysisTree(fhicl::ParameterSet const& pset) : 
  
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel")        ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fCalDataModuleLabel       (pset.get< std::string >("CalDataModuleLabel")      ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel")      ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")        ),
  fEndPoint2DModuleLabel    (pset.get< std::string >("EndPoint2DModuleLabel")   ),
  fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel")       ),
  fMINOSModuleLabel         (pset.get< std::string >("MINOSModuleLabel")        ),
  fTrackMatchModuleLabel    (pset.get< std::string >("TrackMatchModuleLabel")   ),
  fScanModuleLabel          (pset.get< std::string >("ScanModuleLabel")         ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel")          ),
  fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel")  ),
  fvertextrackWindow        (pset.get< double >("vertextrackWindow")            ),
  fvertexclusterWindow      (pset.get< double >("vertexclusterWindow")          ),
  fboundaryWindow           (pset.get< double >("boundaryWindow")               )
  
{
}

//-------------------------------------------------
t962::AnalysisTree::~AnalysisTree()
{
}

void t962::AnalysisTree::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("anatree","analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("pot",&pot,"pot/D");
  fTree->Branch("isdata",&isdata,"isdata/I");
  fTree->Branch("vtxx_reco",&vtxx_reco,"vtxx_reco/D");
  fTree->Branch("vtxy_reco",&vtxy_reco,"vtxy_reco/D");
  fTree->Branch("vtxz_reco",&vtxz_reco,"vtxz_reco/D");
  fTree->Branch("nclusu_reco",&nclusu_reco,"nclusu_reco/I");
  fTree->Branch("nclusv_reco",&nclusv_reco,"nclusv_reco/I");
  fTree->Branch("nclusw_reco",&nclusw_reco,"nclusw_reco/I");
  fTree->Branch("ntracks_reco",&ntracks_reco,"ntracks_reco/I");
  fTree->Branch("nvertextracks_reco",&nvertextracks_reco,"nvertextracks_reco/I");
  fTree->Branch("nvertexclustersu_reco",&nvertexclustersu_reco,"nvertexclustersu_reco/I");
  fTree->Branch("nvertexclustersv_reco",&nvertexclustersv_reco,"nvertexclustersv_reco/I");
  fTree->Branch("nvertexclustersw_reco",&nvertexclustersw_reco,"nvertexclustersw_reco/I");
  fTree->Branch("ntrackendonboundary_reco",&ntrackendonboundary_reco,"ntrackendonboundary_reco/I");
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
  fTree->Branch("trackke",&trackke);
  fTree->Branch("trackrange",&trackrange);
  fTree->Branch("nmatched_reco",&nmatched_reco,"nmatched_reco/I");  
  fTree->Branch("trk_mom_minos",&trk_mom_minos,"trk_mom_minos/D");
  fTree->Branch("trk_charge_minos",&trk_charge_minos,"trk_charge_minos/D");
  fTree->Branch("trk_dcosx_minos",&trk_dcosx_minos,"trk_dcosx_minos/D");
  fTree->Branch("trk_dcosy_minos",&trk_dcosy_minos,"trk_dcosy_minos/D");
  fTree->Branch("trk_dcosz_minos",&trk_dcosz_minos,"trk_dcosz_minos/D");
  fTree->Branch("trk_vtxx_minos",&trk_vtxx_minos,"trk_vtxx_minos/D");
  fTree->Branch("trk_vtxy_minos",&trk_vtxy_minos,"trk_vtxy_minos/D");
  fTree->Branch("trk_vtxz_minos",&trk_vtxz_minos,"trk_vtxz_minos/D"); 
  fTree->Branch("mc_index_minos",&mc_index_minos,"mc_index_minos/F");
  fTree->Branch("mc_pdg_minos",&mc_pdg_minos,"mc_pdg_minos/D");
  fTree->Branch("mc_px_minos",&mc_px_minos,"mc_px_minos/D");
  fTree->Branch("mc_py_minos",&mc_py_minos,"mc_py_minos/D");
  fTree->Branch("mc_pz_minos",&mc_pz_minos,"mc_pz_minos/D");
  fTree->Branch("mc_ene_minos",&mc_ene_minos,"mc_ene_minos/D");
  fTree->Branch("mc_mass_minos",&mc_mass_minos,"mc_mass_minos/D");
  fTree->Branch("mc_vtxx_minos",&mc_vtxx_minos,"mc_vtxx_minos/D");
  fTree->Branch("mc_vtxy_minos",&mc_vtxy_minos,"mc_vtxy_minos/D");
  fTree->Branch("mc_vtxz_minos",&mc_vtxz_minos,"mc_vtxz_minos/D");   
  fTree->Branch("trkcontained_minos",&trkcontained_minos,"trkcontained_minos/I");  
  fTree->Branch("test_charge_minos",&test_charge_minos,"test_charge_minos/I"); 
  fTree->Branch("vtxx_scan", &vtxx_scan, "vtxx_scan/D");
  fTree->Branch("vtxy_scan", &vtxy_scan, "vtxy_scan/D");
  fTree->Branch("vtxz_scan", &vtxz_scan, "vtxz_scan/D");
  fTree->Branch("ntracks_scan", &ntracks_scan , "ntracks_scan/I");
  fTree->Branch("nshowers_scan", &nshowers_scan , "nshowers_scan/I");
  fTree->Branch("neutrino_scan", &neutrino_scan , "neutrino_scan/I");
  fTree->Branch("maybeneutrino_scan", &maybeneutrino_scan , "maybeneutrino_scan/I");    
  fTree->Branch("nuPDG_truth",&nuPDG_truth,"nuPDG_truth/I");
  fTree->Branch("ccnc_truth",&ccnc_truth,"ccnc_truth/I");
  fTree->Branch("mode_truth",&mode_truth,"mode_truth/I");
  fTree->Branch("enu_truth",&enu_truth,"enu_truth/D");
  fTree->Branch("Q2_truth",&Q2_truth,"Q2_truth/D");
  fTree->Branch("hitnuc_truth",&hitnuc_truth,"hitnuc_truth/I");
  fTree->Branch("W_truth",&W_truth,"W_truth/D");
  fTree->Branch("nuvtxx_truth",&nuvtxx_truth,"nuvtxx_truth/D");
  fTree->Branch("nuvtxy_truth",&nuvtxy_truth,"nuvtxy_truth/D");
  fTree->Branch("nuvtxz_truth",&nuvtxz_truth,"nuvtxz_truth/D");
  fTree->Branch("lep_mom_truth",&lep_mom_truth,"lep_mom_truth/D");
  fTree->Branch("lep_dcosx_truth",&lep_dcosx_truth,"lep_dcosx_truth/D");
  fTree->Branch("lep_dcosy_truth",&lep_dcosy_truth,"lep_dcosy_truth/D");
  fTree->Branch("lep_dcosz_truth",&lep_dcosz_truth,"lep_dcosz_truth/D");
}


void t962::AnalysisTree::beginSubRun(const art::SubRun& sr)
{

  art::Handle< sumdata::POTSummary > potListHandle;
  sr.getByLabel(fPOTModuleLabel,potListHandle);
  
  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
    pot=potListHandle->totpot;
  else
    pot=0.;
  
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
  art::Handle< std::vector<recob::EndPoint2D> > endpointListHandle;
  evt.getByLabel(fEndPoint2DModuleLabel,endpointListHandle);
  art::Handle< std::vector<recob::Vertex> > vertexListHandle;
  evt.getByLabel(fVertexModuleLabel,vertexListHandle);
  //art::Handle< std::vector<t962::MINOS> > minosListHandle;
  //evt.getByLabel(fMINOSModuleLabel,minosListHandle);
  art::View< t962::MINOS > minosListHandle;
  evt.getView(fMINOSModuleLabel,minosListHandle);
  art::Handle< std::vector<t962::MINOSTrackMatch> > trackmatchListHandle;
  evt.getByLabel(fTrackMatchModuleLabel,trackmatchListHandle);
  art::Handle< std::vector<t962::ScanInfo> > scanListHandle;
  evt.getByLabel(fScanModuleLabel,scanListHandle);

  art::PtrVector<simb::MCTruth> mclist;
  if(evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    }

  art::PtrVector<recob::Cluster> clusterlist;
  if(evt.getByLabel(fClusterModuleLabel,clusterListHandle))
  for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      clusterlist.push_back(clusterHolder);
    }

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

  //art::PtrVector<t962::MINOS> minoslist;
  std::vector<t962::MINOS*> minoslist;
  //if(evt.getByLabel(fMINOSModuleLabel,minosListHandle))
  //  for (unsigned int i = 0; i < minosListHandle->size(); i++){
  for (unsigned int i = 0; i < minosListHandle.vals().size(); i++){
    //    art::Ptr<t962::MINOS> minosHolder(minosListHandle,i);
    // minoslist.push_back(minosHolder);
    minoslist.push_back((t962::MINOS*)(minosListHandle.vals()[i]));
  }

  art::PtrVector<t962::MINOSTrackMatch> trackmatchlist;
  if(evt.getByLabel(fTrackMatchModuleLabel,trackmatchListHandle))
  for (unsigned int i = 0; i < trackmatchListHandle->size(); i++){
    art::Ptr<t962::MINOSTrackMatch> trackmatchHolder(trackmatchListHandle,i);
    trackmatchlist.push_back(trackmatchHolder);
  }
  
   art::PtrVector<t962::ScanInfo> scanlist;
  if(evt.getByLabel(fScanModuleLabel,scanListHandle))
  for (unsigned int i = 0; i < scanListHandle->size(); i++){
    art::Ptr<t962::ScanInfo> scanHolder(scanListHandle,i);
    scanlist.push_back(scanHolder);
  }
  art::ServiceHandle<geo::Geometry> geom;  

  //vertex information
  if(vertexlist.size())
  {
    double vtxxyz[3];
    vertexlist[0]->XYZ(vtxxyz);
    vtxx_reco = vtxxyz[0];
    vtxy_reco = vtxxyz[1];
    vtxz_reco = vtxxyz[2];
  }
  //cluster information
  nclusu_reco = 0;
  nclusv_reco = 0;
  nclusw_reco = 0;

  int nplanes = geom->Nplanes();
  std::vector<int> Cls[nplanes];

  for(unsigned int i=0; i<clusterlist.size();++i){
    
    switch(clusterlist[i]->View()){
    case geo::kU :
      nclusu_reco ++;
      Cls[0].push_back(i);
      break;
    case geo::kV :
      nclusv_reco ++;
      Cls[1].push_back(i);
      break;
    case geo::kW :
      nclusw_reco ++;
      Cls[2].push_back(i);
      break;
    default :
      break;
    }
  }


  for (int i = 0; i<nplanes; i++){
    int n_vertexclusters = 0;
    if (Cls[i].size()>0){
      int vtx2d_w = -99999;
      double vtx2d_t = -99999;
      bool found2dvtx = false;
      //find 2d vertex
      for (unsigned int j = 0; j<endpointlist.size();j++){
	if (endpointlist[j]->View() == clusterlist[Cls[i][0]]->View()){
	  vtx2d_w = endpointlist[j]->WireNum();
	  vtx2d_t = endpointlist[j]->DriftTime();
	  found2dvtx = true;
	  break;
	}
      }
      if (found2dvtx){
	for (unsigned j = 0; j<Cls[i].size(); j++){
	  double w = clusterlist[Cls[i][j]]->StartPos()[0];
	  double t = clusterlist[Cls[i][j]]->StartPos()[1];
	  double dtdw = clusterlist[Cls[i][j]]->dTdW();
	  double t_vtx = t+dtdw*(vtx2d_w-w);
	  double dis = TMath::Abs(vtx2d_t-t_vtx);
	  if (dis<fvertexclusterWindow){
	    n_vertexclusters++;
	  }
	}
      }
    }
    if (i==0){
      nvertexclustersu_reco = n_vertexclusters;
    }
    else if (i==1){
      nvertexclustersv_reco = n_vertexclusters;
    }
    else if (i==2){
      nvertexclustersw_reco = n_vertexclusters;
    }
  }

  //matching information  
  nmatched_reco = trackmatchlist.size();
  int ANTtrackID=-1;
  test_charge_minos=0.;
  
  for(unsigned int j = 0; j < minoslist.size(); j++)
    { 
      if (!isdata)
     	{
	  if(minoslist[j]->fcharge<0.)
	    test_charge_minos=-1;
        }
    }
  
  for(unsigned int i = 0; i < trackmatchlist.size(); i++)
    {
      
      for(unsigned int j = 0; j < minoslist.size(); j++)
	{
	  
	  if(minoslist[j]->ftrkIndex==trackmatchlist[i]->fMINOStrackid)
	    {
	      ANTtrackID=trackmatchlist[i]->fArgoNeuTtrackid;
	      
	      if(minoslist[j]->ftrkcontained)
		trk_mom_minos = minoslist[j]->ftrkErange;
	      else
		trk_mom_minos = minoslist[j]->ftrkmom;     
	      
	      trkcontained_minos = minoslist[j]->ftrkcontained; 
	      trk_charge_minos = minoslist[j]->fcharge;
	      trk_dcosx_minos = minoslist[j]->ftrkdcosx;
	      trk_dcosy_minos = minoslist[j]->ftrkdcosy;
	      trk_dcosz_minos = minoslist[j]->ftrkdcosz;
	      trk_vtxx_minos = minoslist[j]->ftrkVtxX;
	      trk_vtxy_minos = minoslist[j]->ftrkVtxY;
	      trk_vtxz_minos = minoslist[j]->ftrkVtxZ;        
	    }  
	  if (!isdata){       
	    mc_index_minos = minoslist[j]->fmcIndex;
	    mc_pdg_minos = minoslist[j]->fmcPDG;
	    mc_px_minos = minoslist[j]->fmcPx;
	    mc_py_minos = minoslist[j]->fmcPy;
	    mc_pz_minos = minoslist[j]->fmcPz;
	    mc_ene_minos = minoslist[j]->fmcEne;
	    mc_mass_minos = minoslist[j]->fmcMass;
	    mc_vtxx_minos = minoslist[j]->fmcVtxX;
	    mc_vtxy_minos = minoslist[j]->fmcVtxY;
	    mc_vtxz_minos = minoslist[j]->fmcVtxZ;
	  }
	}
    }

  //track information
  ntracks_reco=tracklist.size();
  double larStart[3];
  double larEnd[3];
  std::vector<double> trackStart;
  std::vector<double> trackEnd;
  trackStart.clear();
  trackEnd.clear();
  //grab information about whether a track is associated with a vertex and whether a track leaves the detector. these variables are powerful discriminators.
  int n_vertextracks=0;
  int n_endonboundarytracks=0;
  for(unsigned int i=0; i<tracklist.size();++i){
    tracklist[i]->Extent(trackStart,trackEnd);            
    trackstart_x_reco=trackStart[0];
    trackstart_y_reco=trackStart[1];
    trackstart_z_reco=trackStart[2];
    trackexit_x_reco=trackEnd[0];
    trackexit_y_reco=trackEnd[1];
    trackexit_z_reco=trackEnd[2];  
    if (!isdata){
      if(sqrt(pow(trackstart_x_reco-mclist[0]->GetNeutrino().Nu().Vx(),2)+pow(trackstart_y_reco-mclist[0]->GetNeutrino().Nu().Vy(),2)+pow(trackstart_z_reco-mclist[0]->GetNeutrino().Nu().Vz(),2))<fvertextrackWindow)
	n_vertextracks++; 
    }
    if(EndsOnBoundary(tracklist[i])) n_endonboundarytracks++;
    
    art::FindOneP<anab::Calorimetry> focal(trackListHandle, evt, fCalorimetryModuleLabel);
    trackke.push_back(focal.at(i)->KinematicEnergy());
    trackrange.push_back(focal.at(i)->Range());
  }
  nvertextracks_reco=n_vertextracks; 
  ntrackendonboundary_reco=n_endonboundarytracks;
  
  //grab information about where track started and ended and the dcos at those points
  trackStart.clear();
  trackEnd.clear();
  memset(larStart, 0, 3);
  memset(larEnd, 0, 3);
  for(unsigned int i=0; i<tracklist.size();++i){
    
    if(ANTtrackID!=tracklist[i]->ID())
      continue;
    
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
  
  //scan information (data only)
  double time= -99.;
  double y_vert=-99.;
  double z_vert=-99.;
  
  for(unsigned int i=0; i<scanlist.size();++i){ 
    time=(scanlist[i]->Get_VertIndTime()+scanlist[i]->Get_VertColTime())/2.;  
    int wire_I=scanlist[i]->Get_VertIndWire();    
    int wire_C=scanlist[i]->Get_VertColWire();   
    if(wire_I>-1 && wire_C>-1)
      {
	wire_I=geom->PlaneWireToChannel(0,wire_I);
	wire_C=geom->PlaneWireToChannel(1,wire_C);
	geom->ChannelsIntersect(wire_I,wire_C,y_vert,z_vert);
	
	vtxx_scan=(time-60)*.031;
	vtxy_scan=y_vert;
	vtxz_scan=z_vert;
	neutrino_scan=scanlist[i]->Get_IsNeutrino();
	maybeneutrino_scan=scanlist[i]->Get_IsMaybeNeutrino();
	ntracks_scan=scanlist[i]->Get_Track();
	nshowers_scan=scanlist[i]->Get_NumShower(); 
      }
  }
  //mc truth information
  if (!isdata){
    nuPDG_truth = mclist[0]->GetNeutrino().Nu().PdgCode();
    ccnc_truth = mclist[0]->GetNeutrino().CCNC();
    mode_truth = mclist[0]->GetNeutrino().Mode();
    Q2_truth = mclist[0]->GetNeutrino().QSqr();
    W_truth = mclist[0]->GetNeutrino().W();
    hitnuc_truth = mclist[0]->GetNeutrino().HitNuc();
    enu_truth = mclist[0]->GetNeutrino().Nu().E();
    nuvtxx_truth = mclist[0]->GetNeutrino().Nu().Vx();
    nuvtxy_truth = mclist[0]->GetNeutrino().Nu().Vy();
    nuvtxz_truth = mclist[0]->GetNeutrino().Nu().Vz();
    lep_mom_truth = mclist[0]->GetNeutrino().Lepton().P();
    if (mclist[0]->GetNeutrino().Lepton().P()){
      lep_dcosx_truth = mclist[0]->GetNeutrino().Lepton().Px()/mclist[0]->GetNeutrino().Lepton().P();
      lep_dcosy_truth = mclist[0]->GetNeutrino().Lepton().Py()/mclist[0]->GetNeutrino().Lepton().P();
      lep_dcosz_truth = mclist[0]->GetNeutrino().Lepton().Pz()/mclist[0]->GetNeutrino().Lepton().P();
    }
  }

  fTree->Fill();
}

  //---------------------------------------------------------------- 
void t962::AnalysisTree::ResetVars(){

  run = -99999;
  event = -99999;
  isdata = -99999;
  vtxx_reco = -99999;
  vtxy_reco  = -99999;
  vtxz_reco  = -99999;
  nclusu_reco  = -99999;
  nclusv_reco  = -99999;
  nclusw_reco  = -99999;
  ntracks_reco = -99999;
  nvertextracks_reco = -99999;
  nvertexclustersu_reco = -99999;
  nvertexclustersv_reco = -99999;
  nvertexclustersw_reco = -99999;
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
  nmatched_reco = -99999;
  trk_mom_minos = -99999;
  trk_charge_minos = -99999;
  trk_dcosx_minos = -99999;
  trk_dcosy_minos = -99999;
  trk_dcosz_minos = -99999;  
  trk_vtxx_minos = -99999;
  trk_vtxy_minos = -99999;
  trk_vtxz_minos = -99999;
  
  mc_index_minos=-99999;
  mc_pdg_minos=-99999;
  mc_px_minos=-99999;
  mc_py_minos=-99999;
  mc_pz_minos=-99999;
  mc_ene_minos=-99999;
  mc_mass_minos=-99999;
  mc_vtxx_minos=-99999;
  mc_vtxy_minos=-99999;
  mc_vtxz_minos=-99999;
  test_charge_minos=-99999;
  trkcontained_minos=-99999;
  vtxx_scan=-99999;
  vtxy_scan=-99999;
  vtxz_scan=-99999;
  neutrino_scan=-99999;
  maybeneutrino_scan=-99999;
  ntracks_scan=-99999;
  nshowers_scan=-99999;         
  nuPDG_truth = -99999;
  ccnc_truth = -99999;
  mode_truth = -99999;
  enu_truth = -99999;
  Q2_truth = -99999;
  W_truth = -99999;
  hitnuc_truth = -99999;
  nuvtxx_truth = -99999;
  nuvtxy_truth = -99999;
  nuvtxz_truth = -99999;
  lep_mom_truth = -99999;
  lep_dcosx_truth = -99999;
  lep_dcosy_truth = -99999;
  lep_dcosz_truth = -99999;
}

bool t962::AnalysisTree::EndsOnBoundary(art::Ptr<recob::Track> lar_track)
{
      std::vector<double> larStart, larEnd;
      lar_track->Extent(larStart,larEnd);//put xyz coordinates at begin/end of track into vectors(?)

 	if(fabs(larEnd[0])<fboundaryWindow 
	|| fabs(47.-larEnd[0])<fboundaryWindow 
	|| fabs(larEnd[1]+20.)<fboundaryWindow
	|| fabs(20.-larEnd[1])<fboundaryWindow 
	|| fabs(larEnd[2])<fboundaryWindow 
	|| fabs(90.-larEnd[2])<fboundaryWindow  )   
	return true;  
      else return false;
}








