////////////////////////////////////////////////////////////////////////
/// $Id: SimChannel.h,v 1.3 2010/03/26 20:08:36 brebel Exp $
///
/// \file  SimChannel.h
///
/// \brief object containing MC truth information necessary for making RawDigits 
/// and doing back tracking
///
/// \author  seligman@nevis.columbia.edu
///
////////////////////////////////////////////////////////////////////////

#ifndef SIM_SIMCHANNEL_H
#define SIM_SIMCHANNEL_H

#include <vector>
#include <map>

namespace sim {


  class IDE{
  public:
    
    IDE();
    virtual ~IDE();

    int    trackID;      ///< Geant4 supplied track ID
    double numElectrons; ///< total number of electrons for this track ID and time
    double x;            ///< x position of ionization
    double y;            ///< y position of ionization
    double z;            ///< z position of ionization
  };
  
  class SimChannel
  {
  public:

    // Default constructor
    SimChannel();
    explicit SimChannel(unsigned short channel);
    
    // Destructor.
    virtual ~SimChannel();

    // method to add ionization electrons to this channel
    void AddIonizationElectrons(int trackID,
				unsigned int tdc,
				double numberElectrons,
				double *xyz);

    unsigned short Channel() const { return fChannel; }

    // method to return a collection of IDE structs for all geant4
    // track ids represented between startTDC and endTDC
    std::vector<sim::IDE> TrackIDsAndEnergies(unsigned int startTDC,
					      unsigned int endTDC) const;

    const std::map<unsigned short, std::vector<sim::IDE> >& TDCIDEMap() const { return fTDCIDEs; }

    // The number of ionization electrons associated with this channel for the 
    // specified TDC.
    double Charge(unsigned int tdc) const;

    bool operator< (const SimChannel& other)     const {return fChannel < other.Channel();}

  private:
    
    unsigned short                                      fChannel; ///< electronics channel associated with these sim::Electrons
    std::map< unsigned short, std::vector< sim::IDE > > fTDCIDEs; ///< vector of IDE structs for each TDC with signal

  };

} // namespace sim

#endif // SIM_SIMCHANNEL_H

////////////////////////////////////////////////////////////////////////
