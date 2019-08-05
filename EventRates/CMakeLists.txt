art_make( BASENAME_ONLY
	  LIB_LIBRARIES

		${ART_FRAMEWORK_CORE}
                ${Boost_SYSTEM_LIBRARY}
		${MF_MESSAGELOGGER}
                ${ROOT_BASIC_LIB_LIST}

		art_Utilities canvas
		cetlib cetlib_except
                lardata_RecoObjects 
		lardataobj_RawData
		lardataobj_Simulation
		larreco_Deprecated
		larreco_RecoAlg
		larcore_Geometry_Geometry_service
		larcorealg_Geometry
		larsim_MCCheater_BackTrackerService_service
		larsim_MCCheater_ParticleInventoryService_service
		larsim_Simulation
                larpandora_LArPandoraInterface
		ubana_AnalysisTree_MCTruth_BackTrackerTruth_tool
                ubana_AnalysisTree_MCTruth_AssociationsTruth_tool
		ubana_AnalysisTree_MCTruth_MCTruthBase
		ubana_BasicTool_GeoAlgo
		#ubana_UBXSec
		ubobj_UBXSec
		ubraw_RawData
	
	MODULE_LIBRARIES 

		art_Persistency_Common canvas
		art_Persistency_Provenance canvas
		art_Utilities canvas
		cetlib cetlib_except
                #larcore_Geometry                                                
		larcore_Geometry_Geometry_service
		larcorealg_Geometry
		larcoreobj_SummaryData 
		lardata_Utilities
		lardataobj_RecoBase
		lardataobj_AnalysisBase 
		lardataobj_RawData
		lardata_RecoObjects 
		lardataobj_Simulation
		larevt_Filters
		larreco_RecoAlg
                larpandora_LArPandoraInterface
		larsim_MCCheater_BackTrackerService_service
		larsim_MCCheater_ParticleInventoryService_service
		larsim_Simulation 
		nutools_ParticleNavigation
		nusimdata_SimulationBase
		ubana_AnalysisTree_MCTruth_BackTrackerTruth_tool
                ubana_AnalysisTree_MCTruth_AssociationsTruth_tool
                ubana_AnalysisTree_MCTruth_MCTruthBase
		ubana_BasicTool_GeoAlgo
		ubana_MyClasses
                #ubana_UBXSec
                ubobj_UBXSec
		ubobj_Trigger
		ubraw_RawData

		${LARRECO_LIB}
		${LARDATA_LIB}
		${ART_FRAMEWORK_CORE}
		${ART_FRAMEWORK_PRINCIPAL}
		${ART_FRAMEWORK_SERVICES_REGISTRY}
		${ART_FRAMEWORK_SERVICES_OPTIONAL}
		${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
		${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
		${MF_MESSAGELOGGER}
		${FHICLCPP}
		${CLHEP}
		${ROOT_GEOM}
		${ROOT_XMLIO}
		${ROOT_GDML}
		${ROOT_BASIC_LIB_LIST}

)

install_headers()
install_fhicl()
install_source()
