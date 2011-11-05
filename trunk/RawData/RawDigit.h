////////////////////////////////////////////////////////////////////////
// $Id: RawDigit.h,v 1.15 2010/03/26 20:06:04 brebel Exp $
//
// Definition of basic raw digits
//
// brebel@fnal.gov
//
// -modified RawDigit class slightly to be compatible with binary output of DAQ software,
//  and to include #samples/channel explicity, instead of via sizeof() methods. 
//  Mitch Soderberg  2/19/2008
//
// - Made the destructor virtual to allow for a Simulation class to inherit from this one.
//   <seligman@nevis.columbia.edu> 05Jun2009
//
////////////////////////////////////////////////////////////////////////

#ifndef RAWDATA_RAWDIGIT_H
#define RAWDATA_RAWDIGIT_H

#include <vector>
#include <iosfwd>

///Raw data description
namespace raw {
  
  typedef enum _compress {
    kNone,       ///< no compression 
    kHuffman,    ///< Huffman Encoding
    kDynamicDec  ///< Dynamic decimation
  } Compress_t;

  class RawDigit {

  public:
    RawDigit(); // Default constructor
    RawDigit(unsigned short channel, 
	     unsigned short samples,
	     std::vector<short> adclist,
	     raw::Compress_t compression=raw::kNone);
    RawDigit(unsigned int channel,
	     std::vector<short> adclist,
	     raw::Compress_t compression=raw::kNone);
    
    
    virtual ~RawDigit();
    
    // Set Methods
    void             SetADC(int i, short iADC);
    void             SetChannel(unsigned short iChan)    { fChannel = iChan; }
    void             SetSamples(unsigned short iSamples) { fSamples = iSamples; }
    
    void             Set(unsigned int channel, 
			 int i, 
			 short adc);
    
    void             SetPedestal(double ped);
    
    // Get Methods
    unsigned int    NADC()        const { return fADC.size();  }
    short           ADC(int i)    const;
    unsigned short  Channel()     const { return fChannel;     }
    unsigned short  Samples()     const { return fSamples;     }
    double          GetPedestal() const { return fPedestal;    } 
    double          GetSigma()    const { return fSigma;       } 
    raw::Compress_t Compression() const { return fCompression; }

    std::vector<short> fADC;

  private:

    unsigned short  fChannel;     ///< channel in the readout
    unsigned short  fSamples;     ///< number of ticks of the clock
    
    double          fPedestal;    ///< pedestal for this channel
    double          fSigma;       ///< sigma of the pedestal counts for this channel

    raw::Compress_t fCompression; ///< compression scheme used for the ADC vector
    
  };
}

#endif // RAWDATA_RAWDIGIT_H

////////////////////////////////////////////////////////////////////////
