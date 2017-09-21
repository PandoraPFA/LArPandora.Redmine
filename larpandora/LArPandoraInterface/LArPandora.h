/**
 *  @file   larpandora/LArPandoraInterface/LArPandora.h
 *
 *  @brief  Base producer module for reconstructing recob::PFParticles from recob::Hits
 *
 */

#ifndef LAR_PANDORA_H
#define LAR_PANDORA_H 1

#include "larpandora/LArPandoraInterface/ILArPandora.h"
#include "larpandora/LArPandoraInterface/LArPandoraInput.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "larreco/Calorimetry/LinearEnergyAlg.h"

#include <string>
#include <memory> // std::unique_ptr<>

namespace lar_pandora
{

/**
 *  @brief  LArPandora class
 */
class LArPandora : public ILArPandora
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset the parameter set
     */
    LArPandora(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~LArPandora();

    void beginJob();
    void produce(art::Event &evt);

protected:
    void DeletePandoraInstances();
    void CreatePandoraInput(art::Event &evt, IdToHitMap &idToHitMap);
    void ProcessPandoraOutput(art::Event &evt, const IdToHitMap &idToHitMap);
    void RunPandoraInstances();
    void ResetPandoraInstances();

    /**
     *  @brief  Create a new Pandora instance and register lar content algs and plugins
     *
     *  @return the address of the new Pandora instance
     */
    const pandora::Pandora *CreateNewPandora() const;

protected:
    std::string                     m_configFile;               ///< The config file
    std::string                     m_stitchingConfigFile;      ///< The stitching config file (multi drift volumes only)

    LArDriftVolumeList              m_driftVolumeList;          ///< The drift volume list
    LArDriftVolumeMap               m_driftVolumeMap;           ///< The map from volume id to drift volume

private:
    LArPandoraInput::Settings       m_inputSettings;            ///< The lar pandora input settings
    LArPandoraOutput::Settings      m_outputSettings;           ///< The lar pandora output settings
    LArPandoraGeometry::Settings    m_geometrySettings;         ///< The lar pandora geometry settings

    std::unique_ptr<calo::LinearEnergyAlg> m_showerEnergyAlg;   ///< Optional cluster energy algorithm.
    bool                            m_runStitchingInstance;     ///< Whether to run the pandora stitching instance (multi drift volumes only)
    bool                            m_enableProduction;         ///< Whether to persist output products
    bool                            m_enableDetectorGaps;       ///< Whether to pass detector gap information to Pandora instances
    bool                            m_lineGapsCreated;          ///< Book-keeping: whether line gap creation has been called
    bool                            m_enableMCParticles;        ///< Whether to pass mc information to Pandora instances to aid development

    std::string                     m_geantModuleLabel;         ///< The geant module label
    std::string                     m_hitfinderModuleLabel;     ///< The hit finder module label
    std::string                     m_mvaModuleLabel;           ///< The mva module label
};

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_H
