///
/// \file    HeaderDrawer.cxx
/// \brief   Render the objects from the Simulation package
/// \author  messier@indiana.edu
/// \version $Id: HeaderDrawer.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
#include "EventDisplay/HeaderDrawer.h"
#include "EventDisplayBase/evdb.h"
#include "TText.h"
#include "TTimeStamp.h"
#include "EventDisplayBase/EventHolder.h"
#include "art/Framework/Principal/Event.h"

namespace evd{

  HeaderDrawer::HeaderDrawer() { }

  //......................................................................

  HeaderDrawer::~HeaderDrawer() { }


  //......................................................................

  void HeaderDrawer::Header(evdb::View2D* view)
  {
    TText& titlet = view->AddText(0.03,0.80,"LArSoft");
    titlet.SetTextSize(0.13);
    titlet.SetTextFont(72);

    // get the event
    const art::Event* evt = evdb::EventHolder::Instance()->GetEvent();
    if(!evt) return;

    int run   = evt->run();
    int srun  = evt->subRun();
    int event = evt->id().event();

    unsigned int year, month, day, dayofweek;
    unsigned int hour, minute, second;
    int          nano;

    // get the time stamp.  art::Timestamp::value() returns a TimeValue_t which is a typedef to unsigned long long.
    // The conventional use is for the upper 32 bits to have the seconds since 1970 epoch and the lower 32 bits to be
    // the number of microseconds with the current second.
    unsigned long long int tsval = evt->time().value();
  
    // taking it apart
    // the masking isn't strictly necessary *if* "long" is truly 32bits
    // but this can vary w/ compiler/platform
    const unsigned long int mask32 = 0xFFFFFFFFUL;
    unsigned long int lup = ( tsval >> 32 ) & mask32;
    unsigned long int llo = tsval & mask32;
    TTimeStamp ts(lup, (int)llo);
    
    ts.GetDate(kTRUE,0,&year,&month,&day);
    ts.GetTime(kTRUE,0,&hour,&minute,&second);
    nano = ts.GetNanoSec();
    dayofweek = ts.GetDayOfWeek();
    char eventbuff[256];
    char runbuff[256];
    char datebuff[256];
    char timebuff[256];
  
    // Skip first one since ROOT returns these numbers starting from 1 not 0
    static const char* days[] = {"",
				 "Mon","Tue","Wed","Thu","Fri","Sat","Sun"
    };
    static const char* months[] = {"",
				   "Jan","Feb","Mar","Apr","May","Jun",
				   "Jul","Aug","Sep","Oct","Nov","Dec"
    };
  
    sprintf(runbuff,  "Run:   %d/%d",run,srun);
    sprintf(eventbuff,"Event: %d",event);
    sprintf(datebuff, "UTC %s %s %d, %d", 
	    days[dayofweek],
	    months[month],
	    day,
	    year);
    sprintf(timebuff, "%.2d:%.2d:%2.9f", 
	    hour, 
	    minute, 
	    (float)second+(float)nano/1.0E9);

    TText& runt   = view->AddText(0.04,0.60, runbuff);
    TText& eventt = view->AddText(0.04,0.45, eventbuff);
    TText& datet  = view->AddText(0.04,0.25, datebuff);
    TText& timet  = view->AddText(0.04,0.10, timebuff);
	
    runt.SetTextSize(0.13);
    runt.SetTextFont(42);

    eventt.SetTextSize(0.13);
    eventt.SetTextFont(42);

    datet.SetTextSize(0.12);
    datet.SetTextFont(42);

    timet.SetTextSize(0.12);
    timet.SetTextFont(42);
  }

}// namespace
////////////////////////////////////////////////////////////////////////
