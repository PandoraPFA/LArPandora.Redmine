// PMTHit, PMTPhoton and PMTHitCollection implementation.
//
// These objects are primarily storage containers, so not much
// function code is actually required.
//
// See comments at head of PMTHit.h for more details.
//
// Ben Jones, MIT, 06/04/2010
//
/// \version $Id: ParticleHistory.cxx,v 1.1 2010/04/29 15:38:01 seligman Exp $


#include "Simulation/sim.h"
#include <iostream>

namespace sim
{
   PMTHitCollection::PMTHitCollection()
  {
  }

  PMTHitCollection::~PMTHitCollection()
  {
    for(PMTHitCollection::const_iterator it=this->begin(); it!=this->end(); it++)
      delete ((*it).second);
  }

  PMTHit * PMTHitCollection::GetHit(int key)
  {
    if( !((*this)[key]) )
      (*this)[key] = new PMTHit();
    return (*this)[key];
      
  }
  
  PMTHitCollection & PMTHitCollection::operator+=(const PMTHitCollection &rhs)
  {
    for(PMTHitCollection::const_iterator it = rhs.begin(); it!=rhs.end(); it++)
      {
	GetHit(it->first)->operator+=(*(it->second));
      }
    SetSDName("CompositeHitCollection");
    return *this;
  }
  
  const PMTHitCollection PMTHitCollection::operator+(const PMTHitCollection &rhs) const
  {
    return PMTHitCollection(*this)+=rhs;
  }


  PMTPhoton::PMTPhoton()
  {
  }

  PMTPhoton::~PMTPhoton()
  {
  }


 
  PMTHit::PMTHit()
  {
  }

  PMTHit::~PMTHit()
  {
  }

 
  PMTHit & PMTHit::operator+=(const PMTHit &rhs)
  {
    for(PMTHit::const_iterator it = rhs.begin(); it!=rhs.end(); it++)
      {
	push_back(*it);
      }
    return *this;
  }

  const PMTHit PMTHit::operator+(const PMTHit &rhs) const
  {
    return PMTHit(*this)+=rhs;
  }


}
