
# where should the scripts/..xml file be installed?  Perhaps in bin?

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{PANDORA_INC} )

art_make( BASENAME_ONLY
          LIBRARY_NAME  LArPandoraInterface
          LIB_LIBRARIES DetectorProperties_service
                        LArProperties_service
                        Geometry
                        Geometry_service
                        TimeService_service
                        Simulation
                        RawData
                        RecoBase
                        Utilities
                        ClusterRecoUtil
                        RecoAlg
                        ClusterFinder
                        ${PANDORASDK}
                        ${PANDORAMONITORING}
                        ${LARCONTENT}
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES LArPandoraInterface 
          )

install_headers()
install_fhicl()
install_source()

add_subdirectory(scripts)

