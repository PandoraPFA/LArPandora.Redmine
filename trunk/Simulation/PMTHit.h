////////////////////////////////////////////////////////////////////////
/// \file  PMTHit.h
/// \brief contains objects relating to PMT hits
///
/// \version $Id: ParticleList.h,v 1.13 2010/05/13 16:12:20 seligman Exp $
/// \author  Ben Jones
////////////////////////////////////////////////////////////////////////
// This file contains the definitions of the classes which
// are stored in the event representing PMT hits.
//
// A PMT Hit stores data for each photon which steps inside the PMT 
// volume.  Currently the quantities stored are 4 potition, 4 momentum
// and TrackID.  A PMTHitCollection is a set of PMTHits, one per PMT
// in the collection. 
// 
// The PMTHit is filled in by the PMTSensitiveDetector class in LArG4 
// and will be used to generate the PMT response later in the simulation
// chain.
//
// PMTPhoton, PMTHit and PMTHitCollection are all persistent under
// ROOT I/O.  For compilation to succeed, the relevant pragma lines
// must be present in LinkDef.h.
//
// The current implementation resembles that of an STL container in 
// some respects but needs more work before it is polished product.
//
// Ben Jones, MIT, 06/04/2010
//

#include "TLorentzVector.h"
#include <map>

#ifndef PMTHit_h
#define PMTHit_h 1

namespace sim
{

  // This structure contains all the information per photon
  // which entered the sensitive PMT volume.

  class PMTPhoton 
  {
  public:
    PMTPhoton();
    virtual ~PMTPhoton();

    bool           SetInSD;
    TLorentzVector Position;
    TLorentzVector Momentum;
  };
  

  // Define a PMT Hit as a list of PMT photons which were
  // recorded in the PMT volume.

  class PMTHit : public std::vector<PMTPhoton> 
    {
    public:
      PMTHit();
      virtual ~PMTHit();
      typedef std::vector<PMTPhoton>             list_type;
      typedef list_type::value_type              value_type;
      typedef list_type::iterator                iterator;
      typedef list_type::const_iterator          const_iterator;
      typedef list_type::reverse_iterator        reverse_iterator;
      typedef list_type::const_reverse_iterator  const_reverse_iterator;
      typedef list_type::size_type               size_type;
      typedef list_type::difference_type         difference_type;

      // define addition operators for combining hits
      //   (add all photons to one vector)
      PMTHit& operator+=(const PMTHit &rhs);
      const PMTHit operator+(const PMTHit &rhs) const;

      const int PMTID()       const { return fPMTID; }
      void      SetID(int id)       { fPMTID = id; }

      int  fPMTID;  /// volume number for the PMT
      
    };
 


  // The PMT Hit collection is the set of all PMT Hits indexed
  // by PMT ID 

  class PMTHitCollection : public std::map<int, PMTHit *>
    {
    public:
      typedef std::map<int,PMTHit *>             list_type;
      typedef list_type::key_type                key_type;
      typedef list_type::mapped_type             mapped_type;
      typedef list_type::value_type              value_type;
      typedef list_type::iterator                iterator;
      typedef list_type::const_iterator          const_iterator;
      typedef list_type::reverse_iterator        reverse_iterator;
      typedef list_type::const_reverse_iterator  const_reverse_iterator;
      typedef list_type::size_type               size_type;
      typedef list_type::difference_type         difference_type;
      typedef list_type::key_compare             key_compare;
      typedef list_type::allocator_type          allocator_type;
 
      PMTHit * GetHit(int);

      // define addition operators for combining hit collections
      //   (add each hit in the collection)
      PMTHitCollection& operator+=(const PMTHitCollection &rhs);
      const PMTHitCollection operator+(const PMTHitCollection &rhs) const; 

    public:
      PMTHitCollection();
      virtual ~PMTHitCollection();
      
    public:
      void SetSDName(std::string TheSDName)  { fTheSDName = TheSDName; }
      std::string GetSDName()                { return fTheSDName; }
     

    private:
      std::string fTheSDName;


    };

}

#endif
