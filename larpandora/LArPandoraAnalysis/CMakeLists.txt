
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( 
          
          LIB_LIBRARIES lardata_Utilities_DetectorProperties_service
                        lardata_Utilities_LArProperties_service
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities_TimeService_service
                        larsim_Simulation
                        lardata_RawData
                        lardata_RecoBase
                        lardata_Utilities
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
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
           MODULE_LIBRARIES larpandora_LArPandoraInterface
          )

install_headers()
install_fhicl()
install_source()

