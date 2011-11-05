/////////////////////////////////////////////////////////////////
//  DBSCANfinder.h
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include <iostream>
#include "art/Framework/Core/EDProducer.h"

class TH1F;

namespace cluster{
   
  //--------------------------------------------------------------- 
  class DBcluster : public art::EDProducer
  {
  public:
    explicit DBcluster(fhicl::ParameterSet const& pset); 
    ~DBcluster();
    void produce(art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:
       
    TH1F *fhitwidth;
    TH1F *fhitwidth_ind_test;  
    TH1F *fhitwidth_coll_test;  
    
    std::string fhitsModuleLabel;
    
  };

}
