////////////////////////////////////////////////////////////////////////
// Class:       GeneratorReader
// Plugin Type: analyzer (art v3_01_02)
// File:        GeneratorReader_module.cc
// A module read in truth information from a generator for cross section comparisons
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

#include "nusimdata/SimulationBase/MCFlux.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/calcLocationWeights.h"




// ROOT includes
#include "TH1.h" 
#include "TH2.h" 
#include "TTree.h" 
#include "TFile.h" 

#define PI 3.14159265

class GeneratorReader;
//______________________________________________________________________________
class GeneratorReader : public art::EDAnalyzer {
public:
    explicit GeneratorReader(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    GeneratorReader(GeneratorReader const&) = delete;
    GeneratorReader(GeneratorReader&&) = delete;
    GeneratorReader& operator=(GeneratorReader const&) = delete;
    GeneratorReader& operator=(GeneratorReader&&) = delete;

    // Required functions.
    // Selected optional functions.
    void beginJob() override;
    void endSubRun(art::SubRun const &sr) override;
    void endJob() override;

    void analyze(art::Event const& e) override;

    bool in_fv(double x, double y, double z);
    

private:

    TFile* output;
    TTree* Tree;
    TTree* POTTree;

    double pot;

    int ievent = 0;

    double elec_e;
    double beta;
    double cosbeta;
    double ppfx_cv;

    bool infv;

    int ccnc;
    int interaction;
    double nu_e;
    int nu_pdg;

    double true_nu_vtx_t;
    double true_nu_vtx_x;
    double true_nu_vtx_y;
    double true_nu_vtx_z;
    double true_nu_px;
    double true_nu_py;
    double true_nu_pz;

    double elec_px;
    double elec_py;
    double elec_pz;

    std::string fEventWeightLabel00, fEventWeightLabel01;

    art::InputTag fMCTproducer;     // MCTruth from neutrino generator
    art::InputTag fMCPproducer;     // MCParticle from Geant4 stage

};
//______________________________________________________________________________
GeneratorReader::GeneratorReader(fhicl::ParameterSet const& p) : EDAnalyzer{p} {
    
    // MC Truth/Particle product labels
    fMCTproducer = p.get<art::InputTag>("MCTproducer");
    fMCPproducer = p.get<art::InputTag>("MCPproducer");

}
//______________________________________________________________________________
void GeneratorReader::beginJob() {
 
    std::cout << "Starting begin Job" << std::endl;

    output = new TFile("events.root", "RECREATE");

    POTTree = new TTree("pottree","pottree");
    POTTree->SetDirectory(output);
    POTTree -> Branch("pot", &pot);

    // Tree for POT counting
    Tree = new TTree("tree","tree");
    Tree->SetDirectory(output);
    Tree -> Branch("nu_e", &nu_e);
    Tree -> Branch("nu_pdg", &nu_pdg);
    Tree -> Branch("ppfx_cv", &ppfx_cv);
    Tree -> Branch("ccnc", &ccnc);
    Tree -> Branch("interaction", &interaction);
    Tree -> Branch("true_nu_vtx_t", &true_nu_vtx_t);
    Tree -> Branch("true_nu_vtx_x", &true_nu_vtx_x);
    Tree -> Branch("true_nu_vtx_y", &true_nu_vtx_y);
    Tree -> Branch("true_nu_vtx_z", &true_nu_vtx_z);
    Tree -> Branch("true_nu_px", &true_nu_px);
    Tree -> Branch("true_nu_py", &true_nu_py);
    Tree -> Branch("true_nu_pz", &true_nu_pz);
    Tree -> Branch("infv", &infv);
    Tree -> Branch("elec_e", &elec_e);
    Tree -> Branch("elec_px", &elec_px);
    Tree -> Branch("elec_py", &elec_py);
    Tree -> Branch("elec_pz", &elec_pz);
    Tree -> Branch("beta", &beta);
    Tree -> Branch("cosbeta", &cosbeta);

}
//______________________________________________________________________________
void GeneratorReader::analyze(art::Event const& e) {

    elec_e = 0.0;
    beta = 0.0;
    cosbeta = 0.0;
    ppfx_cv = 0.0;
    infv = false;
    ccnc = 0;
    interaction = 0;
    nu_e = 0.0;
    nu_pdg = 0;
    true_nu_vtx_x = 0.0;
    true_nu_vtx_y = 0.0;
    true_nu_vtx_t = 0.0;
    true_nu_vtx_z = 0.0;
    true_nu_px = 0.0;
    true_nu_py = 0.0;
    true_nu_pz = 0.0;

    std::vector<art::InputTag> vecTag;
    art::InputTag eventweight_tag_00("eventweight","","EventWeightSept24");
    //vecTag.push_back(eventweight_tag_00);


    for(auto& thisTag : vecTag){

        art::Handle<std::vector<evwgh::MCEventWeight>> eventweights_handle;
        e.getByLabel(thisTag, eventweights_handle);

        if(eventweights_handle.isValid()){

            // std::cout << " [ EventWeightTree ]" << " isValid! " << std::endl;

            std::vector<art::Ptr<evwgh::MCEventWeight>> eventweights;
            art::fill_ptr_vector(eventweights, eventweights_handle);
            std::map<std::string, std::vector<double>> evtwgt_map = eventweights.at(0)->fWeight;

            // Get the PPFX Central Value
            if(evtwgt_map.find("ppfx_cv_UBPPFXCV") != evtwgt_map.end()){
                ppfx_cv = evtwgt_map.find("ppfx_cv_UBPPFXCV")->second[0];
                std::cout << "ppfx cv weight: "<< ppfx_cv <<  std::endl;
            }
        }
    }


    // load MCTruth
    auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);

    // load MCTruth [from geant]
    // auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);


    auto mct = mct_h->at(0);
    auto neutrino = mct.GetNeutrino();
    auto nu = neutrino.Nu();

    ccnc = neutrino.CCNC();
    interaction = neutrino.Mode();
    nu_pdg = nu.PdgCode();
    nu_e = nu.Trajectory().E(0);

    true_nu_vtx_t = nu.T();
    true_nu_vtx_x = nu.Vx();
    true_nu_vtx_y = nu.Vy();
    true_nu_vtx_z = nu.Vz();
    true_nu_px = nu.Px();
    true_nu_py = nu.Py();
    true_nu_pz = nu.Pz();


    infv = in_fv(true_nu_vtx_x, true_nu_vtx_y, true_nu_vtx_z);


    elec_px     = (neutrino.Lepton().Px() / neutrino.Lepton().P());
    elec_py     = (neutrino.Lepton().Py() / neutrino.Lepton().P());
    elec_pz     = (neutrino.Lepton().Pz() / neutrino.Lepton().P());
    elec_e   = neutrino.Lepton().E();

    TVector3 elec_dir(elec_px, elec_py, elec_pz); // elec direction
    elec_dir.Unit();

    TVector3 nu_dir(true_nu_px, true_nu_py, true_nu_pz); // nu direction
    nu_dir.Unit();

    beta = elec_dir.Angle(nu_dir) * 180 / 3.14159;
    cosbeta = std::cos(elec_dir.Angle(nu_dir));

    std::cout << 
    "Nu PDG: " << nu_pdg << "\n" <<
    "Nu E: " << nu_e*1000 << "\n" <<
    "Elec E: " << elec_e*1000 << "\n" <<
    "Beta: " << beta << "\n" <<
    "cosbeta: " << cosbeta << "\n" <<
    "FV: " << infv << "\n" <<
    std::endl;
    

    Tree->Fill();

}

void GeneratorReader::endSubRun(art::SubRun const & sr) {

    auto const & POTSummaryHandle = sr.getValidHandle < sumdata::POTSummary >("generator");
    auto const & POTSummary(*POTSummaryHandle);
    const double total_pot = POTSummary.totpot;
    std::cout << "----------------------------" << std::endl;
    std::cout << "Total POT / subRun: " << total_pot << std::endl;
    std::cout << "----------------------------" << std::endl;

    pot = total_pot;
    POTTree->Fill();

}

//______________________________________________________________________________
void GeneratorReader::endJob() {

    output->cd();
    
    POTTree->Write();
    Tree->Write();

}
//___________________________________________________________________________

bool GeneratorReader::in_fv(double x, double y, double z){
    
    // Require the vertex is in the boundary
    if ( x   >= 8.45   && x <= 244.8
      && y   >= -106.5 && y <= 106.5
      && z   >= 5      && z <= 1031.8  ){
        return true;   // pass
    }  
    else return false; // fail
}

//______________________________________________________________________________
DEFINE_ART_MODULE(GeneratorReader)
