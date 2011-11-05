#ifndef RecoCheckAna_h
#define RecoCheckAna_h
////////////////////////////////////////////////////////////////////////
// Class:       RecoCheckAna
// Module Type: analyzer
// File:        RecoCheckAna.h
//
// Generated at Fri Jul 15 09:54:26 2011 by Brian Rebel using artmod
// from art v0_07_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "RecoBase/recobase.h"
#include "Simulation/sim.h"

namespace cheat {
  class RecoCheckAna;
}

class TH1D;

class cheat::RecoCheckAna : public art::EDAnalyzer {
public:
  explicit RecoCheckAna(fhicl::ParameterSet const &p);
  virtual ~RecoCheckAna();

  virtual void analyze(art::Event const &e);

  virtual void reconfigure(fhicl::ParameterSet const & p);
  virtual void beginRun(art::Run const &r);

  

private:

  void CheckRecoClusters(std::vector< art::Ptr<recob::Cluster> >  const& clscol,
			 std::vector<const sim::SimChannel*>           & scs,
			 sim::ParticleList                        const& plist,
			 art::PtrVector<recob::Hit>              const& allhits);
  void CheckRecoProngs(std::vector< art::Ptr<recob::Prong> >      const& prgcol,
		       std::vector<const sim::SimChannel*>             & scs,
		       sim::ParticleList                          const& plist,
		       art::PtrVector<recob::Hit>                const& allhits,
		       TH1D*                                             purity,
		       TH1D*                                             efficiency);
  void CheckRecoVertices(std::vector< art::Ptr<recob::Vertex> >   const& vtxcol,
			 std::vector<const sim::SimChannel*>           & scs,
			 sim::ParticleList                        const& plist,
			 art::PtrVector<recob::Hit>              const& allhits);
  void CheckRecoEvents(std::vector< art::Ptr<recob::Event> >      const& evtcol,
		       std::vector<const sim::SimChannel*>             & scs,
		       sim::ParticleList                          const& plist,
		       art::PtrVector<recob::Hit>                const& allhits);

  std::string fHitModuleLabel;     ///< lable for module making the hits
  std::string fClusterModuleLabel; ///< label for module making the clusters
  std::string fShowerModuleLabel;  ///< label for module making the showers
  std::string fTrackModuleLabel;   ///< label for module making the tracks
  std::string fVertexModuleLabel;  ///< label for module making the vertices
  std::string fEventModuleLabel;   ///< label for module making the events
  std::string fG4ModuleLabel;      ///< label for module making the geant4 simulation

  bool        fCheckClusters;      ///< should we check the reconstruction of clusters?
  bool        fCheckShowers;       ///< should we check the reconstruction of showers?
  bool        fCheckTracks;        ///< should we check the reconstruction of tracks?
  bool        fCheckVertices;      ///< should we check the reconstruction of vertices?
  bool        fCheckEvents;        ///< should we check the reconstruction of events?

  TH1D*       fClusterPurity;      ///< histogram of cluster purity
  TH1D*       fClusterEfficiency;  ///< histogram of cluster efficiency
  TH1D*       fShowerPurity;       ///< histogram of shower purity
  TH1D*       fShowerEfficiency;   ///< histogram of shower efficiency
  TH1D*       fTrackPurity;        ///< histogram of track purity
  TH1D*       fTrackEfficiency;    ///< histogram of track efficiency
  TH1D*       fVertexPurity;       ///< histogram of vertex purity
  TH1D*       fVertexEfficiency;   ///< histogram of vertex efficiency
  TH1D*       fEventPurity;        ///< histogram of event purity
  TH1D*       fEventEfficiency;    ///< histogram of event efficiency

};
#endif /* RecoCheckAna_h */
