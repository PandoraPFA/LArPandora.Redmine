///
/// \file    SimulationDrawer.h
/// \brief   Render the objects from the Simulation package
/// \author  messier@indiana.edu
/// \version $Id: SimulationDrawer.h,v 1.2 2010/11/10 22:38:34 p-novaart Exp $
///
#ifndef EVD_SIMULATIONDRAWER_H
#define EVD_SIMULATIONDRAWER_H
#include <string>
#include <vector>
#include <map>
#ifdef __CINT__
namespace art {
  class Ptr;
  class PtrVector;
}
#else 
#include "art/Persistency/Common/PtrVector.h"
#endif

namespace art  { class Event;       }
namespace evdb { class View2D;      }
namespace evdb { class View3D;      }
namespace geo  { class Geometry;    }
namespace simb { class MCTruth;     }
namespace sim  { class Particle;    }

namespace evd {
  class SimulationDrawer {
  public:
    SimulationDrawer();
    ~SimulationDrawer();

  public:
    // Drawing functions
    void MCTruthShortText(const art::Event& evt,
			  evdb::View2D*     view);
    void MCTruthLongText(const art::Event& evt,
			 evdb::View2D*     view);
    void MCTruthVectors2D(const art::Event& evt,
			  evdb::View2D*     view,
			  unsigned int      plane);
    void MCTruth3D(const art::Event& evt,
		   evdb::View3D*     view);

#ifndef __CINT__

#endif

    void HiLite(int trkId, bool hlt=true);

  private:
    int GetMCTruth(const art::Event&              evt,
		   art::PtrVector<simb::MCTruth>& mctruth);

    int GetParticle(const art::Event&              evt,
		    art::PtrVector<sim::Particle>& plist);

  private:
    std::map<int,bool>       fHighlite;

  };
}

#endif
////////////////////////////////////////////////////////////////////////
