////////////////////////////////////////////////////////////////////////
// Class:       EventWeightReader
// Plugin Type: analyzer (art v2_05_01)
// File:        EventWeightReader_module.cc
//
// Generated at Wed Aug 22 05:28:06 2018 by Krishan Mistry using cetskelgen
// from cetlib version v1_21_00.
//
// This LArsoft module will read in the MC weight object from a root file 
// and put it into a ttree format for easier reading into code
////////////////////////////////////////////////////////////////////////

// Default art includes
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

#include "larsim/EventWeight/Base/MCEventWeight.h"

// ROOT includes
#include "TTree.h" 
#include "TString.h"
#include <iostream>

// -----------------------------------------------------------------------------------------------------
class EventWeightReader;
class EventWeightReader : public art::EDAnalyzer {
public:
    explicit EventWeightReader(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    EventWeightReader(EventWeightReader const &) = delete;
    EventWeightReader(EventWeightReader &&) = delete;
    EventWeightReader & operator = (EventWeightReader const &) = delete;
    EventWeightReader & operator = (EventWeightReader &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;
    void AddWeights(std::vector<Event_List> &N, art::Event const & e, bool &resize_flag);

private:

    // Declare member data here.
    int run, subrun, evt;
    std::vector<std::map<std::string, std::vector<double>>> weights; // Map (Model Name, Weight Vector in each universe)

    int Iterations{1};          // Number of times looped event to see how quickly module is running

    // TTree
    TTree *tree;

};
// -----------------------------------------------------------------------------------------------------
// A function that loops over all the parameter weights and universes and re-weights the desired events. 
void EventWeightReader::AddWeights(std::vector<Event_List> &N, art::Event const & e, bool &resize_flag){

    auto GenieEW_Handle = e.getValidHandle<std::vector<evwgh::MCEventWeight>>("mcweight"); // Request the mcweight data product

    if (GenieEW_Handle.isValid())  std::cout << "[Analyze] GenieEW_Handle is valid" << std::endl; 
    
    std::vector<evwgh::MCEventWeight> const& GenieEWvec(*GenieEW_Handle);

    for (evwgh::MCEventWeight const& GenieEW: GenieEWvec) weights = GenieEW.fWeight; // Grab the weights 

    // Loop over the labels e.g qevec, qema, all etc.
    for (unsigned int j=0; j < N.size(); j++){
        
        // Loop over the labels in the EventWeight object
        for (auto const& it : weights) {
                
            if (it.first.find(N.at(j).label.c_str()) != std::string::npos) { // Match up the labels from input list and EventWeight

                std::cout << "N.at(j).label.c_str():\t" << N.at(j).label.c_str() << std::endl;

                // unisim
                if (it.second.size() == 2 && N.at(j).mode == "unisim"){
                    // std::cout << "unisim input!" << std::endl;

                    // Do this once per file
                    if (resize_flag == false){ 
                        N_size = N.at(j).N_reweight.size();                   // The current size (the number of universes already ran over)
                        N.at(j).N_reweight.resize(N_size + it.second.size()); // resize to the number of universes if we havent already done so
                    }

                    // Loop over each universe for a parameter and weight the event
                    for (unsigned int i = N_size; i < N_size + it.second.size(); i++)
                        N.at(j).N_reweight.at(i) += it.second.at(i - N_size) ; // Weight the event

                }
                // mutisim
                if (it.second.size() != 2 && N.at(j).mode == "multisim"){
                    // std::cout << "multisim input!" << std::endl;
                    // Do this once per file
                    if (resize_flag == false){ 
                        N_size = N.at(j).N_reweight.size();                   // The current size (the number of universes already ran over)
                        N.at(j).N_reweight.resize(N_size + it.second.size()); // resize to the number of universes if we havent already done so
                    }
                
                    // Loop over each universe for a parameter and weight the event
                    for (unsigned int i = N_size; i < N_size + it.second.size(); i++)
                        N.at(j).N_reweight.at(i) += it.second.at(i - N_size) ; // Weight the event

                }
            
            }
        } // END loop over parameters
    }
    resize_flag = true;
}
// -----------------------------------------------------------------------------------------------------
EventWeightReader::EventWeightReader(fhicl::ParameterSet const & p) : EDAnalyzer(p) {}
// -----------------------------------------------------------------------------------------------------
void EventWeightReader::beginJob() {
    
    // Access ART's TFileService, which will handle histograms/trees/etc.
    art::ServiceHandle<art::TFileService> tfs;

    // Create the TTree and add relavent branches
    tree     = tfs->make<TTree>("weight_tree","weight_ree");
    tree    ->Branch("weights", "std::vector<std::map<std::string, std::vector<double>>>", &weights);
    tree    ->Branch("run",    &run); 
    tree    ->Branch("subrun", &subrun); 
    tree    ->Branch("evt",    &evt); 

    

}
void EventWeightReader::analyze(art::Event const & e) {
    
    // Determine event ID, run and subrun 
    run =     e.id().run();
    subrun =  e.id().subRun();
    evt =     e.id().event();


    // Get the MC weight dataproduct
    auto GenieEW_Handle = e.getValidHandle<std::vector<evwgh::MCEventWeight>>("mcweight"); // Request the mcweight data product
    if (GenieEW_Handle.isValid())  std::cout << "[Analyze] GenieEW_Handle is valid" << std::endl; 
    std::vector<evwgh::MCEventWeight> const& GenieEWvec(*GenieEW_Handle);
    
    weight = *GenieEWvec;
    

    tree->Fill();

    
}
// -----------------------------------------------------------------------------------------------------
void EventWeightReader::endJob() {

}
// -----------------------------------------------------------------------------------------------------
DEFINE_ART_MODULE(EventWeightReader)
