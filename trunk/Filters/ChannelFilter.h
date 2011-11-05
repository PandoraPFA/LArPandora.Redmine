////////////////////////////////////////////////////////////////////////
//
// ChannelFilter class:
//
// This class provides methods for returning the condition of
// a wire, as to whether it is bad and is to be ignored or perhaps
// if it has some known problem.  This allows the removal of detector
// specific code from a few places in LArSoft.  Right now is is only 
// implemented for Argoneut.
//  
//
// pagebri3@msu.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef CHANNELFILTER_H
#define CHANNELFILTER_H

#include "Geometry/geo.h"

///filters for channels, events, etc
namespace filter {

  class ChannelFilter {

  public:

    ChannelFilter();
    ~ChannelFilter();

    bool BadChannel(unsigned int channel);
    bool NoisyChannel(unsigned int channel);
    std::set<unsigned int> SetOfBadChannels()   const { return fBadChannels;   }
    std::set<unsigned int> SetOfNoisyChannels() const { return fNoisyChannels; }
  private:

    std::set<unsigned int>            fBadChannels;   ///< list of bad channels
    std::set<unsigned int>            fNoisyChannels; ///< list of bad channels

  }; //class ChannelFilter
}
#endif // CHANNELFILTER_H


