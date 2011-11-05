////////////////////////////////////////////////////////////////////////////
// \version $Id: Cluster.cxx,v 1.7 2010/06/12 21:46:34 spitz7 Exp $
//
// \brief Definition of cluster object for LArSoft
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////

#include "RecoBase/Cluster.h"
#include <string>
#include <iostream>

namespace recob{

  //----------------------------------------------------------------------
  Cluster::Cluster() : 
    fdTdW(0.),
    fdQdW(0.),
    fID(-1)
  {
    std::cerr << "WARNING: using default Cluster ctor - bad!" << std::endl;
  }

  //----------------------------------------------------------------------
  Cluster::Cluster(art::PtrVector<recob::Hit> &hits,
		   double startWire, double sigmaStartWire,
		   double startTime, double sigmaStartTime,
		   double endWire, double sigmaEndWire, 
		   double endTime, double sigmaEndTime,
		   double dTdW, double sigmadTdW,
		   double dQdW, double sigmadQdW,
		   int id) :
    fdTdW(dTdW),
    fdQdW(dQdW),
    fSigmadTdW(sigmadTdW),
    fSigmadQdW(sigmadQdW),
    fID(id),
    fHits(hits)
  {
    fStartPos.push_back(startWire);
    fStartPos.push_back(startTime);
    fSigmaStartPos.push_back(sigmaStartWire);
    fSigmaStartPos.push_back(sigmaStartTime);

    fEndPos.push_back(endWire);
    fEndPos.push_back(endTime);
    fSigmaEndPos.push_back(sigmaEndWire);
    fSigmaEndPos.push_back(sigmaEndTime);
  }

  //----------------------------------------------------------------------
  Cluster::~Cluster()
  {
  }

  //----------------------------------------------------------------------
  ///if allwires is true then return all hits for this view and wire if specified
  ///just return those hits
  art::PtrVector<recob::Hit> Cluster::Hits(unsigned int wire, bool allwires) const
  {
  
    unsigned int t = 0;
    unsigned int p = 0;
    unsigned int w = 0;

    art::PtrVector<recob::Hit> hits;
    art::ServiceHandle<geo::Geometry> geo;

    for(size_t i = 0; i < fHits.size(); ++i){
      geo->ChannelToWire(fHits[i]->Wire()->RawDigit()->Channel(), t, p, w);
      if(w == wire || allwires) hits.push_back(fHits[i]);
    }
    
    return hits;
  }

  //----------------------------------------------------------------------
  double Cluster::Charge() const
  {
    double charge = 0.;

    for(size_t i = 0; i < fHits.size(); ++i) charge += fHits[i]->Charge();
  
    return charge;
  }
  
  //----------------------------------------------------------------------
  // Print information for all Hits in a Cluster.  
  //
  void Cluster::PrintHits() 
  {
    for(size_t i = 0; i < fHits.size(); ++i) 
      std::cout << "Hit #" << i << " : " << *fHits[i] << std::endl;

    return;
  }
  
  

  //----------------------------------------------------------------------
  //  Addition operator.  
  //

  Cluster Cluster::operator +(Cluster a)
  {

    // check the start and end positions - for now the 
    // smallest wire number means start position, largest means end position
    std::vector<double> astart(a.StartPos());
    std::vector<double> aend  (a.EndPos()  );
    std::vector<double> start(StartPos());
    std::vector<double> end  (EndPos()  );
    std::vector<double> sigstart(SigmaStartPos());
    std::vector<double> sigend  (SigmaEndPos()  );

    if(astart[0] < fStartPos[0]){
      start = astart;
      sigstart = a.SigmaStartPos();
    }
    
    if(aend[0] > fEndPos[0]){
      end = aend;
      sigend = a.SigmaEndPos();
    }
  
    //inelegant merging of Hits since there is no "+" operator for PtrVector;
    //also seems a bit kludgey that to grab the Hits requires knowing the View.
    art::PtrVector<recob::Hit> hits = Hits();
    art::PtrVector<recob::Hit> ihits = a.Hits();
    
    //take weighted mean in obtaining average slope and differential charge, based on number of hits in each cluster
    double dtdw = ((hits.size()*dTdW()) + (ihits.size()*a.dTdW()))/(hits.size()+ihits.size());
    double dqdw = ((hits.size()*dQdW()) + (ihits.size()*a.dQdW()))/(hits.size()+ihits.size());
    
    double sigdtdw = TMath::Max(SigmadTdW(), a.SigmadTdW());
    double sigdqdw = TMath::Max(SigmadQdW(), a.SigmadQdW());
   
    art::PtrVector<recob::Hit>::const_iterator hitIter = ihits.begin();
    while (hitIter!=ihits.end()){	
      hits.push_back(*hitIter);
      hitIter++;
    }
    
    hits.sort();//sort the PtrVector to organize Hits of new Cluster
  
    Cluster sum(hits, 
		start[0], sigstart[0], 
		start[1], sigstart[1],
		end[0],   sigend[0],
		end[1],   sigend[1],
		dtdw, sigdtdw,
		dqdw, sigdqdw,
		ID());
    
    return sum;
    
  }

  //----------------------------------------------------------------------
  // ostream operator.  
  //

  std::ostream& operator<< (std::ostream& o, const Cluster& c)
  {
    o << std::setiosflags(std::ios::fixed) << std::setprecision(2);
    o << "Cluster ID "    << std::setw(5)  << std::right << c.ID() 
      << " : View = "     << std::setw(3)  << std::right << c.View()  
      << " StartWire = "  << std::setw(7)  << std::right << c.StartPos()[0] 
      << " EndWire = "    << std::setw(7)  << std::right << c.EndPos()[0] 
      << " StartTime = "  << std::setw(9)  << std::right << c.StartPos()[1]
      << " EndTime = "    << std::setw(9)  << std::right << c.EndPos()[1] 
      << " dTdQ = "       << std::setw(9)  << std::right << c.dTdW()
      << " Charge = "     << std::setw(10) << std::right << c.Charge()
      << " #Hits = "      << std::setw(5)  << std::right << c.Hits().size();

    return o;
  }


  //----------------------------------------------------------------------
  // < operator.  
  //
  bool operator < (const Cluster & a, const Cluster & b)
  {
    if(a.View() != b.View())
      return a.View()<b.View();
    if(a.ID() != b. ID())
      return a.ID()<b.ID();
    if(a.StartPos()[0] != b.StartPos()[0])
      return a.StartPos()[0] < b.StartPos()[0];
    if(a.EndPos()[0] != b.EndPos()[0])
      return a.EndPos()[0] < b.EndPos()[0];

    return false; //They are equal
  }

}

