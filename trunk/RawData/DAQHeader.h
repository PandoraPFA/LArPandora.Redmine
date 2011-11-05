////////////////////////////////////////////////////////////////////////
// $Id: DAQHeader.h,v 1.7 2009/02/19 23:07:30 soderber Exp $
//
// Definition of basic DAQ header information
//
// brebel@fnal.gov
//
// -modified DAQHeader class to save all information available in 
//  binary version of DAQ480 software.  - Mitch Soderberg 2/19/09  
//
////////////////////////////////////////////////////////////////////////

#ifndef DAQHEADER_H
#define DAQHEADER_H

#include <vector>
#include <iosfwd>
#include <time.h>

namespace raw {

  class DAQHeader {
    public:
      DAQHeader(); // Default constructor
      DAQHeader(unsigned int status); // Ascii DAQ constructor
      DAQHeader(unsigned int status,  // Binary DAQ constructor
		int fixed,
		unsigned short format,
		unsigned short software,
		unsigned short run,
		unsigned short event,
		time_t time,
		short spare,
		unsigned short nchan);
      
      ~DAQHeader();

      // Set Methods
      void             SetStatus(unsigned int i){fStatus = i;}
      void             SetFixedWord(int i){fFixed = i;}
      void             SetFileFormat(unsigned short i){fFormat = i;}
      void             SetSoftwareVersion(unsigned short i){fSoftware = i;}
      void             SetRun(unsigned short i){fRun = i;}
      void             SetEvent(unsigned short i){fEvent = i;}
      void             SetTimeStamp(time_t t){fTime = t;}
      void             SetSpareWord(short s){fSpare = s;}
      void             SetNChannels(unsigned short i){fNchan = i;}

      

      // Get Methods
      unsigned int     GetStatus()          const { return fStatus; }
      int              GetFixedWord()       const { return fFixed; }
      unsigned short   GetFileFormat()      const { return fFormat; }
      unsigned short   GetSoftwareVersion() const { return fSoftware; }
      unsigned short   GetRun()             const { return fRun; }
      unsigned short   GetEvent()           const { return fEvent; }
      time_t           GetTimeStamp()       const { return fTime; }
      short            GetSpareWord()       const { return fSpare; }
      unsigned short   GetNChannels()       const { return fNchan; }
     

    private:

      unsigned int   fStatus;
      int            fFixed;
      unsigned short fFormat;
      unsigned short fSoftware;
      unsigned short fRun;
      unsigned short fEvent;
      time_t         fTime;
      short          fSpare;
      unsigned short fNchan;

    };
}

#endif // DAQHEADER_H

////////////////////////////////////////////////////////////////////////
