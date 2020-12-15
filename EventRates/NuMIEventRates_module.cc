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


#include "larcoreobj/SummaryData/POTSummary.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"


// ROOT includes
#include "TH1.h" 
#include "TH2.h" 
#include "TTree.h" 


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
    
private:

    // Declare member data here.

    // MC Truth/Particle product labels
    std::string fMC_Particle_Module_Label;
    std::string fMC_Truth_Module_Label;

    // POT product labels
    std::string _potsum_producer_mc;
    std::string _potsum_instance;

    // Add Ttree variables 
    TTree *MCTree;
    int run, subrun, evt;
    double vx, vy, vz;
    double energy;
    int Interaction;
    int PDG;

    TTree* _sr_tree;
    double _sr_pot;

    int iteration{0};

};
//______________________________________________________________________________
NuMIEventRates::NuMIEventRates(fhicl::ParameterSet const& p) : EDAnalyzer{p} {
    
    // MC Truth/Particle product labels
    fMC_Particle_Module_Label = p.get<std::string>("MCParticleModuleLabel");
    fMC_Truth_Module_Label    = p.get<std::string>("MCTruthModuleLabel");

    // POT product labels
    _potsum_producer_mc       = p.get<std::string>("POTSummaryProducerMC");
	_potsum_instance          = p.get<std::string>("POTSummaryInstance");

}
//______________________________________________________________________________
void NuMIEventRates::beginJob() {
 
    std::cout << "Starting begin Job" << std::endl;
    
    // Art TFile Service
    art::ServiceHandle<art::TFileService> tfs;

    // POT Trees
    _sr_tree = tfs->make<TTree>("pottree","");
	_sr_tree->Branch("pot",                &_sr_pot,                "pot/D");

    // Create the TTree and add relavent branches
    MCTree = tfs->make<TTree>("EventTree","EventTree");
        
    // Add Tree Information
    MCTree->Branch("run",   &run);
    MCTree->Branch("subrun",&subrun);
    MCTree->Branch("event", &evt);
    MCTree->Branch("PDG",   &PDG);
    MCTree->Branch("Interaction",   &Interaction);
    MCTree->Branch("energy",   &energy);

    
    std::cout << "Ending begin Job" << std::endl;

}
//______________________________________________________________________________
void NuMIEventRates::analyze(art::Event const& e) {

    iteration+=1;
    std::cout << "Running over Event:\t" << iteration << std::endl;

    // Determine run, subrun, event IDs
    run    = e.id().run();
    subrun = e.id().subRun();
    evt    = e.id().event();

    // Get the MC Truth vector
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
  
    if (e.getByLabel("generator",mctruthListHandle)) 
        art::fill_ptr_vector(mclist, mctruthListHandle);

    int iList{ 0 }; // 1 nu int per spill

    for (int p = 0; p < mclist[iList]->NParticles(); p++) { // Loop over MCTruth Particles



        // Get the interation type
        Interaction = mclist[iList]->GetNeutrino().Mode();

        simb::MCParticle particle{mclist[iList]->GetParticle(p)}; // Get a MC Particle

        if (mclist[iList]->Origin() == simb::kBeamNeutrino){ // Require a beam neutrino

            if (particle.PdgCode() == 12 || particle.PdgCode() == -12){ // nue in the event
            
                // Set values for filling TTree
                PDG =  particle.PdgCode();
                vx  =  particle.Vx();
                vy  =  particle.Vy();
                vz  =  particle.Vz();
                energy = particle.E();
            }

        } // END if a beam neutrino

    } // END loop over mclist

    // Check whether the interaction was in the FV and if the energy is > 250 MeV
     if ( vx   >= 20.0   && vx <= 236.35 && vy >= -96.5 && vy <= 96.5 && vz   >= 20.0   && vz <= 1016.8  ){
        // Was in FV

        // And the energy was greater than the threshold
        if (energy > 0.250){
             MCTree->Fill();
        }

    }  
}

void NuMIEventRates::endSubRun(art::SubRun const & sr) {
	
    art::Handle<sumdata::POTSummary> potsum_h;

    // MC
    if(sr.getByLabel(_potsum_producer_mc, potsum_h)) {
        _sr_pot = potsum_h->totpot;
    }
    else
        _sr_pot = 0.;

    std::cout << "----------------------------" << std::endl;
    std::cout << "Total POT / subRun: " << _sr_pot << std::endl;
    std::cout << "----------------------------" << std::endl;


    _sr_tree->Fill();

}

//______________________________________________________________________________
void NuMIEventRates::endJob() {

}
//______________________________________________________________________________
DEFINE_ART_MODULE(NuMIEventRates)
