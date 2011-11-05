////////////////////////////////////////////////////////////////////////////
// \version $Id: Prong.h,v 1.4 2010/06/10 16:21:31 antonm Exp $
//
// \brief Definition of event object for LArSoft
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////

#ifndef RB_EVENT_H
#define RB_EVENT_H

#include <vector>
#include <iosfwd>
#include <string>
#include <iostream>
#include <iomanip>

#include "RecoBase/Vertex.h"
#include "art/Persistency/Common/PtrVector.h"


namespace recob {
  
  class Event  {

  public:
    
    Event();  // Default constructor
    Event(art::PtrVector<recob::Vertex> &vertices,
	  int id=-999);
    ~Event();

    art::PtrVector<recob::Hit>           Hits()          const;
    const art::PtrVector<recob::Vertex>& Vertices()      const { return fVertices; }
    double                               Energy()        const;			
    double                            	 SigmaEnergy()   const;			
    art::Ptr<recob::Vertex>           	 PrimaryVertex() const { return fPrimary;  }
    const int                         	 ID()            const { return fID;       }

    friend bool          operator <   (const Event & a, const Event & b);

  protected:
     virtual std::ostream&          Print(std::ostream& stream) const;
     friend  std::ostream& operator <<   (std::ostream& o, const Event & a);

    
  private:
    
    art::PtrVector<recob::Vertex> fVertices; ///< collection of Vertex objects
    art::Ptr<recob::Vertex>       fPrimary;  ///< primary vertex in the event
    int                           fID;       ///< id for this event
  };
}

#endif // RB_EVENT_H
