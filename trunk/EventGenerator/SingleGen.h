////////////////////////////////////////////////////////////////////////
/// \file SingleGen.h
///
/// Module to produce a set list of particles for a MC event
///
/// \version $Id: SingleGen.h,v 1.2 2010/02/15 19:10:40 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
// Framework includes.
#include <vector>
#include <string>
#include "art/Framework/Core/EDProducer.h"

namespace simb { class MCTruth; }

namespace evgen {

  /// module to produce single or multiple specified particles in the detector
  class SingleGen : public art::EDProducer {

  public:
    explicit SingleGen(fhicl::ParameterSet const& pset);
    virtual ~SingleGen();

    // This is called for each event.
    void produce(art::Event& evt);
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    void SampleOne(unsigned int i, simb::MCTruth &mct);        
    void Sample(simb::MCTruth &mct);        

    static const int kUNIF = 0;    
    static const int kGAUS = 1;    

    unsigned int        fSeed;           ///< random number seed    
    int                 fMode;           ///< Particle Selection Mode 
                                         ///< 0--generate a list of all particles, 
                                         ///< 1--generate a single particle selected randomly from the list
    std::vector<int>    fPDG;            ///< PDG code of particles to generate    
    std::vector<double> fP0;             ///< Central momentum (GeV/c) to generate    
    std::vector<double> fSigmaP;         ///< Variation in momenta (GeV/c)    
    int                 fPDist;          ///< How to distribute momenta (gaus or uniform)    
    std::vector<double> fX0;             ///< Central x position (cm) in world coordinates 
    std::vector<double> fY0;             ///< Central y position (cm) in world coordinates
    std::vector<double> fZ0;             ///< Central z position (cm) in world coordinates
    std::vector<double> fSigmaX;         ///< Variation in x position (cm)    
    std::vector<double> fSigmaY;         ///< Variation in y position (cm)    
    std::vector<double> fSigmaZ;         ///< Variation in z position (cm)    
    int                 fPosDist;        ///< How to distribute xyz (gaus, or uniform)        
    std::vector<double> fTheta0XZ;       ///< Angle in XZ plane (degrees)    
    std::vector<double> fTheta0YZ;       ///< Angle in YZ plane (degrees)    
    std::vector<double> fSigmaThetaXZ;   ///< Variation in angle in XZ plane    
    std::vector<double> fSigmaThetaYZ;   ///< Variation in angle in YZ plane    
    int                 fAngleDist;      ///< How to distribute angles (gaus, uniform)

  };
};
////////////////////////////////////////////////////////////////////////
