
art_make(	
      
      LIB_LIBRARIES
            ${ART_FRAMEWORK_CORE}
            art_Utilities canvas
            ${MF_MESSAGELOGGER}
            cetlib cetlib_except
            ${Boost_SYSTEM_LIBRARY}
            ${ROOT_BASIC_LIB_LIST}
        
        MODULE_LIBRARIES
            lardataobj_RecoBase
            nusimdata_SimulationBase
            ${LARDATA_LIB}
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
            art_Persistency_Common canvas
            art_Persistency_Provenance canvas
            art_Utilities canvas
            ${MF_MESSAGELOGGER}
            ${FHICLCPP}
            cetlib cetlib_except
            ${CLHEP}
            ${ROOT_GEOM}
            ${ROOT_XMLIO}
            ${ROOT_GDML}
            ${ROOT_BASIC_LIB_LIST}
            ${LARCOREOBJ_LIB}
            ${LARDATAOBJ_LIB}
            ${LARSIM_LIB}
            ${DK2NUDATA_LIB}
            lardataobj_Simulation
            lardataobj_MCBase
            nusimdata_SimulationBase
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
            larpandora_LArPandoraInterface
            lardata_RecoObjects
            lardata_Utilities
        )

install_headers()
install_fhicl()
install_source()