////////////////////////////////////////////////////////////////////////
// Class:       NuMIEventRates
// Plugin Type: analyzer (art v3_01_02)
// File:        NuMIEventRates_module.cc
// A module to get the event rate distribution of NuMI events in MicroBooNE
// to validate the input flux prediction
// Generated at Sun Aug  4 21:35:12 2019 by Krishan Mistry using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

// Default Art Includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

//cetlib
#include "cetlib/cpu_timer.h"
#include "cetlib_except/exception.h"

// LArSoft Includes

// larcore
#include "larcore/Geometry/Geometry.h"

// larcorealg
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

// larcoreobj
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::Compress_t, raw::Channel_t
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/ParticleID_VariableTypeEnums.h"

// lardata
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/GeometryUtilities.h"

// lardataobj
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"

// larsim
#include "larsim/MCCheater/ParticleInventoryService.h"

// nusim
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

// nutools
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nutools/ParticleNavigation/EveIdCalculator.h"

// C++ Includes
#include <cmath>
#include <map>
#include <iostream>
#include <algorithm>

// ROOT includes
#include "TH1.h" 
#include "TH2.h" 
#include "TTree.h" 
#include "TDatabasePDG.h" 
#include "TParticlePDG.h" 
#include "TCanvas.h" 
#include "TF1.h" 
#include "TMath.h"

#define PI 3.14159265

class NuMIEventRates;
//______________________________________________________________________________
class NuMIEventRates : public art::EDAnalyzer {
public:
    explicit NuMIEventRates(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    NuMIEventRates(NuMIEventRates const&) = delete;
    NuMIEventRates(NuMIEventRates&&) = delete;
    NuMIEventRates& operator=(NuMIEventRates const&) = delete;
    NuMIEventRates& operator=(NuMIEventRates&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;
    void endSubRun(art::SubRun const &sr) override;
    void endJob() override;
    
    double Calc_Theta(double Pz, double P);
    double Calc_Phi  (double Px, double Py, double P);

private:

    // Declare member data here.

    // MC Truth/Particle product labels
    std::string fMC_Particle_Module_Label;
    std::string fMC_Truth_Module_Label;

    // SW Trigger product label
    std::string fSoftware_Trigger_Algo;
    std::string fSoftware_Trigger_Module_Label;

    // POT product labels
    std::string _potsum_producer_mc;
    std::string _potsum_producer_data;
    std::string _potsum_instance;
    
    // Histograms
    // Nue all
    TH1D* 	hNue_Energy_All;	 
    TH1D* 	hNue_Theta_All;	  
    TH1D* 	hNue_Phi_All;	

    TH2D* 	hNue_E_vs_Theta_All;	 
    TH2D* 	hNue_E_vs_Phi_All;   

    // Nue
    TH1D* 	hNue_Energy;	
    TH1D*	hNue_Theta;		 
    TH1D* 	hNue_Phi;			  

    // Nue bar
    TH1D* 	hNue_bar_Energy;	
    TH1D*	hNue_bar_Theta;		 
    TH1D* 	hNue_bar_Phi;			 

    // NuMu all
    TH1D* 	hNuMu_Energy_All;	 
    TH1D* 	hNuMu_Theta_All;	  
    TH1D* 	hNuMu_Phi_All;	    

    TH2D* 	hNuMu_E_vs_Theta_All;	 
    TH2D* 	hNuMu_E_vs_Phi_All;  

    // NuMu
    TH1D* 	hNuMu_Energy;	
    TH1D*	hNuMu_Theta;		 
    TH1D* 	hNuMu_Phi;	

    // NuMu bar
    TH1D* 	hNuMu_bar_Energy;	
    TH1D*	hNuMu_bar_Theta;		 
    TH1D* 	hNuMu_bar_Phi;	

    int NC_event{0};

    // Add Ttree variables 
    TTree *MCTree;
    int run, subrun, evt;
    double NueEnergy,  NueTheta,  NuePhi;
    double NuMuEnergy, NuMuTheta, NuMuPhi;
    int PDG;

    TTree * pot_tree;
    double pot = 0.0;

    TTree* _sr_tree;
    int _sr_run, _sr_subrun;
    double _sr_begintime, _sr_endtime;
    double _sr_pot;
    std::ofstream _run_subrun_list_file;

};
//______________________________________________________________________________
NuMIEventRates::NuMIEventRates(fhicl::ParameterSet const& p) : EDAnalyzer{p} {
    
    // MC Truth/Particle product labels
    fMC_Particle_Module_Label = p.get<std::string>("MCParticleModuleLabel");
    fMC_Truth_Module_Label    = p.get<std::string>("MCTruthModuleLabel");

    // SW Trigger product label
    // fSoftware_Trigger_Algo         = p.get<std::string>("SoftwareTriggerAlgo");
    // fSoftware_Trigger_Module_Label = p.get<std::string>("SoftwareTriggerModuleLabel");

    // POT product labels
    _potsum_producer_mc             = p.get<std::string>("POTSummaryProducerMC");
	_potsum_instance                = p.get<std::string>("POTSummaryInstance");

}
//______________________________________________________________________________
double NuMIEventRates::Calc_Theta(double Pz, double P){

  double Z_dir{ Pz / P }; // Z direction
  double Theta{ std::acos(Z_dir) * (180. / PI) }; // Theta

  return Theta; 
}
//______________________________________________________________________________
double NuMIEventRates::Calc_Phi(double Px, double Py, double P){

  double Y_dir {Py / P}; // Y direction 
  double X_dir {Px / P}; // X direction
  
  double Phi{std::atan2(Y_dir, X_dir) * (180. / PI)}; // Phi

  return Phi; 
}
//______________________________________________________________________________
void NuMIEventRates::beginJob() {
   
    // Art TFile Service
    art::ServiceHandle<art::TFileService> tfs;

    // POT Trees
     pot_tree = tfs->make<TTree>("pot_tree", "pot_per_subrun");
    _sr_tree = tfs->make<TTree>("pottree","");
    pot_tree->Branch("pot", &pot, "pot/D");

	_sr_tree->Branch("run",                &_sr_run,                "run/I");
	_sr_tree->Branch("subrun",             &_sr_subrun,             "subrun/I");
	_sr_tree->Branch("begintime",          &_sr_begintime,          "begintime/D");
	_sr_tree->Branch("endtime",            &_sr_endtime,            "endtime/D");
	_sr_tree->Branch("pot",                &_sr_pot,                "pot/D");
    _run_subrun_list_file.open ("run_subrun_list.txt", std::ofstream::out | std::ofstream::trunc);

    // Create the TTree and add relavent branches
    MCTree = tfs->make<TTree>("EventTree","EventTree");
        
    // Make directories
    art::TFileDirectory Nue_All_dir = tfs->mkdir( "Nue_All_dir" );
    art::TFileDirectory Nue_dir     = tfs->mkdir( "Nue_dir"     );
    art::TFileDirectory Nue_bar_dir = tfs->mkdir( "Nue_bar_dir" );

    art::TFileDirectory NuMu_All_dir  = tfs->mkdir( "NuMu_All_dir" );
    art::TFileDirectory NuMu_dir      = tfs->mkdir( "NuMu_dir"  );
    art::TFileDirectory NuMu_bar_dir  = tfs->mkdir( "NuMu_bar_dir"   );

    // Histograms 

    // Nue All
    hNue_E_vs_Theta_All = Nue_All_dir.make<TH2D>("Nue_E_vs_Theta_All","Nue_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",15., 0., 7. , 10., 0., 180);
    hNue_E_vs_Phi_All   = Nue_All_dir.make<TH2D>("Nue_E_vs_Phi_All","Nue_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",10., 0., 10. , 10., -180., 180);
    hNue_E_vs_Theta_All ->SetOption("COLZ,TEXT");
    hNue_E_vs_Phi_All   ->SetOption("COLZ,TEXT");
    
    hNue_Energy_All = Nue_All_dir.make<TH1D>("Nue_Energy_All","Nue_Energy_All; E [GeV]; Events",400., 0., 20.);
    hNue_Theta_All  = Nue_All_dir.make<TH1D>("Nue_Theta_All","Nue_Theta_All; Theta [Degrees]; Events", 100., 0., 180);
    hNue_Phi_All    = Nue_All_dir.make<TH1D>("Nue_Phi_All","Nue_Phi_All; Phi [Degrees]; Events", 10., -180., 180);
    hNue_Energy_All ->SetOption("HIS");
    hNue_Theta_All  ->SetOption("HIST,TEXT00");
    hNue_Phi_All    ->SetOption("HIST,TEXT00");

    // Nue
    hNue_Energy = Nue_dir.make<TH1D>("Nue_Energy","Nue_Energy; E [GeV]; Events",400., 0., 20.);
    hNue_Theta  = Nue_dir.make<TH1D>("Nue_Theta","Nue_Theta; Theta [Degrees]; Events", 200., 0., 180);
    hNue_Phi    = Nue_dir.make<TH1D>("Nue_Phi","Nue_Phi; Phi [Degrees]; Events", 10., -180., 180);
    hNue_Energy ->SetOption("HIST,TEXT00");
    hNue_Theta  ->SetOption("HIST,TEXT00");
    hNue_Phi    ->SetOption("HIST,TEXT00");

    // Nue bar
    hNue_bar_Energy = Nue_bar_dir.make<TH1D>("Nue_bar_Energy","Nue_bar_Energy; E [GeV]; Events",400., 0., 20.);
    hNue_bar_Theta  = Nue_bar_dir.make<TH1D>("Nue_bar_Theta","Nue_bar_Theta; Theta [Degrees]; Events", 10., 0., 180);
    hNue_bar_Phi    = Nue_bar_dir.make<TH1D>("Nue_bar_Phi","Nue_bar_Phi; Phi [Degrees]; Events", 10., -180., 180);
    hNue_bar_Energy ->SetOption("HIS");
    hNue_bar_Theta  ->SetOption("HIST,TEXT00");
    hNue_bar_Phi    ->SetOption("HIST,TEXT00");  

    // NuMu all
    hNuMu_E_vs_Theta_All = NuMu_All_dir.make<TH2D>("NuMu_E_vs_Theta_All","NuMu_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",20., 0., 10. , 10., 0., 180);
    hNuMu_E_vs_Phi_All   = NuMu_All_dir.make<TH2D>("NuMu_E_vs_Phi_All","NuMu_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",20., 0., 10. , 10., -180., 180);
    hNuMu_E_vs_Theta_All->SetOption("COLZ,TEXT");
    hNuMu_E_vs_Phi_All->SetOption("COLZ,TEXT");
    
    hNuMu_Energy_All = NuMu_All_dir.make<TH1D>("NuMu_Energy_All","NuMu_Energy_All; E [GeV]; Events",400., 0., 20.);
    hNuMu_Theta_All  = NuMu_All_dir.make<TH1D>("NuMu_Theta_All","NuMu_Theta_All; Theta [Degrees]; Events", 10., 0., 180);
    hNuMu_Phi_All    = NuMu_All_dir.make<TH1D>("NuMu_Phi_All","NuMu_Phi_All; Phi [Degrees]; Events", 10., -180., 180);
    hNuMu_Energy_All  ->SetOption("HIS");
    hNuMu_Theta_All   ->SetOption("HIST,TEXT00");
    hNuMu_Phi_All     ->SetOption("HIST,TEXT00");

    // NuMu
    hNuMu_Energy = NuMu_dir.make<TH1D>("NuMu_Energy","NuMu_Energy; E [GeV]; Events",400., 0., 20.);
    hNuMu_Theta  = NuMu_dir.make<TH1D>("NuMu_Theta","NuMu_Theta; Theta [Degrees]; Events", 10., 0., 180);
    hNuMu_Phi    = NuMu_dir.make<TH1D>("NuMu_Phi","NuMu_Phi; Phi [Degrees]; Events", 10., -180., 180);
    hNuMu_Energy ->SetOption("HIS");
    hNuMu_Theta ->SetOption("HIST,TEXT00");
    hNuMu_Phi   ->SetOption("HIST,TEXT00");

    // NuMu_bar
    hNuMu_bar_Energy = NuMu_bar_dir.make<TH1D>("NuMu_bar_Energy","NuMu_bar_Energy; E [GeV]; Events",400., 0., 20.);
    hNuMu_bar_Theta  = NuMu_bar_dir.make<TH1D>("NuMu_bar_Theta","NuMu_bar_Theta; Theta [Degrees]; Events", 10., 0., 180);
    hNuMu_bar_Phi    = NuMu_bar_dir.make<TH1D>("NuMu_bar_Phi","NuMu_bar_Phi; Phi [Degrees]; Events", 10., -180., 180);
    hNuMu_bar_Energy->SetOption("HIS");
    hNuMu_bar_Theta->SetOption("HIST,TEXT00");
    hNuMu_bar_Phi->SetOption("HIST,TEXT00");

    // Add Tree Information
    MCTree->Branch("run",   &run);
    MCTree->Branch("subrun",&subrun);
    MCTree->Branch("event", &evt);
    MCTree->Branch("PDG",   &PDG);

    MCTree->Branch("NueEnergy", &NueEnergy);
    MCTree->Branch("NueTheta",  &NueTheta);
    MCTree->Branch("NuePhi",    &NuePhi);
    
    MCTree->Branch("NuMuEnergy", &NuMuEnergy);
    MCTree->Branch("NuMuTheta",  &NuMuTheta);
    MCTree->Branch("NuMuPhi",    &NuMuPhi);

}
//______________________________________________________________________________
void NuMIEventRates::analyze(art::Event const& e) {

    std::cout << "Running over Event:\t" << e.event() << std::endl;

    // Determine run, subrun, event IDs
    run    = e.id().run();
    subrun = e.id().subRun();
    evt    = e.id().event();

    double ParticleE{ 0. };
    double Theta{ 0. }; 
    double Phi{ 0. }; 
    NC_event = 0; // reset NC counter

    // Get the MC Truth vector
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
  
    if (e.getByLabel("generator",mctruthListHandle)) 
        art::fill_ptr_vector(mclist, mctruthListHandle);

    int iList{ 0 }; // 1 nu int per spill

    for (int p = 0; p < mclist[iList]->NParticles(); p++) { // Loop over MCTruth Particles

        //if (mclist[iList]->GetNeutrino().CCNC() == 1) NC_event ++; // The event was NC and we dont want to include the additional scattered neutrino
        if (mclist[iList]->GetNeutrino().CCNC() == 1) continue; // The event was NC then skip the event

        simb::MCParticle particle{mclist[iList]->GetParticle(p)}; // Get a MC Particle

        ParticleE = particle.E() ;                               // Energy

        if ( NC_event >= 2 ) continue;// Skip fill for the second neutrino if NC

        if (mclist[iList]->Origin() == simb::kBeamNeutrino){ // Require a beam neutrino
            
            if (particle.PdgCode() < 2212) std::cout << "PDG:\t"<< particle.PdgCode()<< std::endl;
            
            // BEGIN SELECTING PARTICLE BLOCK
            if (particle.PdgCode() == 12){ // nue in the event
                
                // Fill a histogram with the pdg code, energy, theta, phi
                hNue_Energy_All             ->Fill( ParticleE );
                hNue_Theta_All              ->Fill(Theta); 
                hNue_Phi_All                ->Fill(Phi); 
                hNue_E_vs_Theta_All         ->Fill(ParticleE, Theta);
                hNue_E_vs_Phi_All           ->Fill(ParticleE, Phi);
                
                hNue_Energy ->Fill( ParticleE );
                hNue_Theta  ->Fill(Theta); 
                hNue_Phi    ->Fill(Phi);

                // Set laues for filling TTree
                NueEnergy =  ParticleE;
                NueTheta  =  Theta; 
                NuePhi    =  Phi;
                PDG       =  particle.PdgCode();

                
            } 
            else if (particle.PdgCode() == -12){ // nue bar in the event

                // Fill a histogram with the pdg code, energy, theta, phi
                hNue_Energy_All       ->Fill( ParticleE );
                hNue_Theta_All        ->Fill( Theta ); 
                hNue_Phi_All          ->Fill( Phi ); 
                hNue_E_vs_Theta_All   ->Fill( ParticleE, Theta );
                hNue_E_vs_Phi_All     ->Fill( ParticleE, Phi );

                hNue_bar_Energy ->Fill( ParticleE );
                hNue_bar_Theta  ->Fill( Theta ); 
                hNue_bar_Phi    ->Fill( Phi );

                // Set laues for filling TTree
                NueEnergy =  ParticleE;
                NueTheta  =  Theta; 
                NuePhi    =  Phi;
                PDG       =  particle.PdgCode();

            } 
            else if (particle.PdgCode() == 14){ // numu in the event

                // Fill a histogram with the pdg code, energy, theta, phi
                hNuMu_Energy_All     ->Fill( ParticleE );
                hNuMu_Theta_All      ->Fill( Theta );
                hNuMu_Phi_All        ->Fill( Phi );
                hNuMu_E_vs_Theta_All ->Fill( ParticleE, Theta );
                hNuMu_E_vs_Phi_All   ->Fill( ParticleE, Phi );

                hNuMu_Energy   ->Fill( ParticleE );
                hNuMu_Theta    ->Fill( Theta );
                hNuMu_Phi      ->Fill( Phi );

                // Set laues for filling TTree
                NuMuEnergy =  ParticleE;
                NuMuTheta  =  Theta; 
                NuMuPhi    =  Phi;
                PDG        =  particle.PdgCode();

            } 
            else if (particle.PdgCode() == -14){ // numu bar in the event

                // Fill a histogram with the pdg code, energy, theta, phi
                hNuMu_Energy_All     ->Fill( ParticleE );
                hNuMu_Theta_All      ->Fill( Theta );
                hNuMu_Phi_All        ->Fill( Phi );
                hNuMu_E_vs_Theta_All ->Fill( ParticleE, Theta );
                hNuMu_E_vs_Phi_All   ->Fill( ParticleE, Phi );

                hNuMu_bar_Energy   ->Fill( ParticleE );
                hNuMu_bar_Theta    ->Fill( Theta );
                hNuMu_bar_Phi      ->Fill( Phi );

                // Set laues for filling TTree
                NuMuEnergy =  ParticleE;
                NuMuTheta  =  Theta; 
                NuMuPhi    =  Phi;
                PDG        =  particle.PdgCode();

            } // END IF CONDITION BLOCK  
    
        } // END if a beam neutrino

    } // END loop over mclist

    MCTree->Fill();

}

void NuMIEventRates::endSubRun(art::SubRun const & sr) {

    
    auto const & POTSummaryHandle = sr.getValidHandle < sumdata::POTSummary >("generator");
    auto const & POTSummary(*POTSummaryHandle);
    const double total_pot = POTSummary.totpot;
    std::cout << "----------------------------" << std::endl;
    std::cout << "Total POT / subRun: " << total_pot << std::endl;
    std::cout << "----------------------------" << std::endl;

    pot = total_pot;
    pot_tree->Fill();

    // Saving run and subrun number on file so that we can run Zarko's script easily
	_run_subrun_list_file << sr.run() << " " << sr.subRun() << std::endl;
	
    _sr_run       = sr.run();
    _sr_subrun    = sr.subRun();
    _sr_begintime = sr.beginTime().value();
    _sr_endtime   = sr.endTime().value();

    art::Handle<sumdata::POTSummary> potsum_h;

    // MC
    if(sr.getByLabel(_potsum_producer_mc, potsum_h)) {
        _sr_pot = potsum_h->totpot;
    }
    else
        _sr_pot = 0.;

    _sr_tree->Fill();

}

//______________________________________________________________________________
void NuMIEventRates::endJob() {

}
//______________________________________________________________________________
DEFINE_ART_MODULE(NuMIEventRates)
