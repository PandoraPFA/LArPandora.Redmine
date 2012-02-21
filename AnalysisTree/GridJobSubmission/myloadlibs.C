void myloadlibs() {
  gSystem->Load("libBeamDataNtuple.so");
  gSystem->Load("libBeamDataNtuple");
  gSystem->Load("libBeamDataUtil");

  gSystem->Load("libCandNtupleSR.so");
  gSystem->Load("libMCNtuple.so");
  gSystem->Load("libTruthHelperNtuple.so");
  gSystem->Load("libStandardNtuple.so");


  gSystem->Load("libCandNtupleEM");
  gSystem->Load("libNeugenInterface");
  gSystem->Load("libMCReweight");
  gSystem->Load("libAnalysisNtuples");
  gSystem->Load("libAnalysisNtuplesModule");
  gSystem->Load("libDcsUser");
  gSystem ->Load("libMuonRemoval.so");

  gSystem->Load("libMinuit2.so");
  gSystem -> Load("libPhysicsNtuple.so");
  gSystem -> Load("libPhysicsNtupleFill.so");
  gSystem -> Load("libPhysicsNtupleSelect.so");
  gSystem -> Load("libPhysicsNtupleStore.so");
  gSystem -> Load("libPhysicsNtuplekNNAlg.so");
  
  gSystem->Load("libMad.so");  
  gSystem->Load("libBeamDataUtil.so");
  gSystem->Load("libSpillTiming.so");
  gSystem->Load("libNtpFitSA.so");
  gSystem->Load("libNuBarPID.so");
  gSystem->Load("libNtupleUtils.so");
  gSystem->Load("libNuMuBar.so");
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libAstroUtil.so");
  gSystem->Load("libAstroUtiltest.so");

}
