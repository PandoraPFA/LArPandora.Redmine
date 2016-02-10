/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraInput.h
 *
 *  @brief  Helper functions for providing inputs to pandora
 */

#ifndef LAR_PANDORA_INPUT_H
#define LAR_PANDORA_INPUT_H 1

#include "LArPandoraInterface/ILArPandora.h"
#include "LArPandoraInterface/LArPandoraHelper.h"

namespace lar_pandora
{

/**
 *  @brief  LArPandoraInput class
 */
class LArPandoraInput
{
public:
    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        /**
         *  @brief  Default constructor
         */
        Settings();

        const pandora::Pandora *m_pPrimaryPandora;          ///< 
        const ILArPandora      *m_pILArPandora;             ///<
        bool                    m_useHitWidths;             ///<
        int                     m_uidOffset;                ///<
        double                  m_dx_cm;                    ///<
        double                  m_int_cm;                   ///<
        double                  m_rad_cm;                   ///<
        double                  m_dEdX_max;                 ///<
        double                  m_dEdX_mip;                 ///<
        double                  m_mips_to_gev;              ///<
        double                  m_recombination_factor;     ///<
    };

    /**
     *  @brief  Create the Pandora 2D hits from the ART hits
     *
     *  @param  settings the settings
     *  @param  hits the input list of ART hits for this event
     *  @param  idToHitMap to receive the mapping from Pandora hit ID to ART hit
     */
    static void CreatePandoraHits2D(const Settings &settings, const HitVector &hitVector, IdToHitMap &idToHitMap);

    /**
     *  @brief  Create the Pandora 3D hits from the ART space points
     *
     *  @param  settings the settings
     *  @param  spacepoints the input list of ART spacepoints for this event
     *  @param  spacePointsToHits the mapping between ART space points and ART hits
     *  @param  spacePointMap mapping from Pandora hit addresses to ART space points
     */
    static void CreatePandoraHits3D(const Settings &settings, const SpacePointVector &spacePointVector,
        const SpacePointsToHits &spacePointsToHits, SpacePointMap &spacePointMap);

    /**
     *  @brief  Create the Pandora MC particles from the MC particles
     *
     *  @param  settings the settings
     *  @param  truthToParticles  mapping from MC truth to MC particles
     *  @param  particlesToTruth  mapping from MC particles to MC truth
     */
    static void CreatePandoraMCParticles(const Settings &settings, const MCTruthToMCParticles &truthToParticles,
        const MCParticlesToMCTruth &particlesToTruth);

    /**
     *  @brief  Create 2D projections of the Pandora MC particles
     *
     *  @param  settings the settings
     *  @param  particleVector the input vector of MC particles
     */
    static void CreatePandoraMCParticles2D(const Settings &settings, const MCParticleVector &particleVector);

    /**
     *  @brief  Create links between the 2D hits and Pandora MC particles
     *
     *  @param  settings the settings
     *  @param  hitMap mapping from Pandora hit addresses to ART hits
     *  @param  hitToParticleMap mapping from each ART hit to its underlying G4 track ID
     */
    static void CreatePandoraMCLinks2D(const Settings &settings, const HitMap &hitMap, const HitsToTrackIDEs &hitToParticleMap);

private:
    /**
     *  @brief  Loop over MC trajectory points and identify start and end points within the detector
     *
     *  @param  settings the settings
     *  @param  volumeID the drift volume
     *  @param  particle the true particle
     *  @param  startT the first trajectory point in the detector
     *  @param  endT the last trajectory point in the detector
     */
    static void GetTrueStartAndEndPoints(const Settings &settings, const int volumeID, const art::Ptr<simb::MCParticle> &particle,
        int &startT, int &endT);

    /**
     *  @brief  Loop over MC trajectory points and identify start and end points within a given cryostat and TPC
     *
     *  @param  cstat the cryostat
     *  @param  tpc the TPC
     *  @param  particle the true particle
     *  @param  startT the first trajectory point in the detector
     *  @param  endT the last trajectory point in the detector
     */
    static void GetTrueStartAndEndPoints(const unsigned int cstat, const unsigned int tpc, const art::Ptr<simb::MCParticle> &particle,
        int &startT, int &endT);

    /**
     *  @brief  Use detector and time services to get a true X offset for a given trajectory point
     *
     *  @param  particle the true particle
     *  @param  nT the trajectory point
     */
    static float GetTrueX0(const art::Ptr<simb::MCParticle> &particle, const int nT);

    /**
     *  @brief  Convert charge in ADCs to approximate MIPs
     *
     *  @param  settings the settings
     *  @param  hit_Charge the input charge
     *  @param  hit_View the input view number
     */
    static double GetMips(const Settings &settings, const double hit_Charge, const geo::View_t hit_View);
};

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_INPUT_H
