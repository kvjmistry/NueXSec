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
    double GetNuMIAngle(double px, double py, double pz);
    double SetPPFXCVWeight(double _nu_e, double _nu_ang, int _nu_pdg);
    

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
    double nu_e, nu_ang;
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

    std::string fPPFX_origin;

    enum flav {
    k_numu,
    k_numubar,
    k_nue,
    k_nuebar,
    k_FLAV_MAX
    };

    std::vector<std::string> flav_str = {
        "numu",
        "numubar",
        "nue",
        "nuebar"
    };

    // Hiatograms to store the flux histograms ratios
    std::vector<TH2D*> hist_ratio;
    std::vector<TH2D*> hist_ratio_uw;

    TFile *f_flux;

};
//______________________________________________________________________________
GeneratorReader::GeneratorReader(fhicl::ParameterSet const& p) : EDAnalyzer{p} {
    
    // MC Truth/Particle product labels
    fMCTproducer = p.get<art::InputTag>("MCTproducer");
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
    fPPFX_origin = p.get<std::string>("PPFX_origin");

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
    Tree -> Branch("nu_ang", &nu_ang);
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

    TFile *f_flux = TFile::Open("/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/FHC/output_uboone_fhc_run0_merged.root", "READ");

    hist_ratio.resize(k_FLAV_MAX);
    hist_ratio_uw.resize(k_FLAV_MAX);

    // Get the flux histograms
    for (unsigned int f = 0; f < flav_str.size(); f++){
        hist_ratio.at(f)      = (TH2D*) f_flux->Get(Form("%s/Detsmear/%s_CV_AV_TPC_2D", flav_str.at(f).c_str(), flav_str.at(f).c_str()));
        hist_ratio_uw.at(f)   = (TH2D*) f_flux->Get(Form("%s/Detsmear/%s_unweighted_AV_TPC_2D", flav_str.at(f).c_str(), flav_str.at(f).c_str()));

        hist_ratio.at(f)->SetDirectory(0);
        hist_ratio_uw.at(f)->SetDirectory(0);

        // Take the ratio
        hist_ratio.at(f)->Divide(hist_ratio_uw.at(f));

    }

    f_flux->Close();

    output->cd();

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
    nu_ang = 0.0;
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
    if (fPPFX_origin == "artroot")
        vecTag.push_back(eventweight_tag_00);


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

    nu_ang = GetNuMIAngle(true_nu_px, true_nu_py, true_nu_pz);

    infv = in_fv(true_nu_vtx_x, true_nu_vtx_y, true_nu_vtx_z);

    if (fPPFX_origin == "calculate")
        ppfx_cv = SetPPFXCVWeight(nu_e, nu_ang, nu_pdg);


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

    // Skip NC interactions
    if (ccnc == 1)
        continue;

    std::cout << 
    "Nu PDG: " << nu_pdg << "\n" <<
    "Nu E: " << nu_e*1000 << "\n" <<
    "Nu Ang: " << nu_ang << "\n" <<
    "Elec E: " << elec_e*1000 << "\n" <<
    "Beta: " << beta << "\n" <<
    "cosbeta: " << cosbeta << "\n" <<
    "FV: " << infv << "\n" <<
    "CCNC: " << ccnc << "\n" <<
    "PPFX CV: " << ppfx_cv << "\n" <<
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

//___________________________________________________________________________
double GeneratorReader::GetNuMIAngle(double px, double py, double pz){

    // Variables
    TRotation RotDet2Beam;             // Rotations
    TVector3  detxyz, BeamCoords;      // Translations
    std::vector<double> rotmatrix;     // Inputs

    // input detector coordinates to translate
    detxyz = {px, py, pz};     

    // From beam to detector rotation matrix
    rotmatrix = {
        0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021,
        4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359,
        -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291 };

    // Return the TRotation
    TVector3 newX, newY, newZ;
    newX = TVector3(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
    newY = TVector3(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
    newZ = TVector3(rotmatrix[6], rotmatrix[7], rotmatrix[8]);

    RotDet2Beam.RotateAxes(newX, newY, newZ); // Return the TRotation now det to beam
    // RotDet2Beam.Invert(); // Invert back to the beam to det

    // Rotate to beam coords
    BeamCoords = RotDet2Beam * detxyz;

    TVector3 beamdir = {0 , 0 , 1};;
    
    // Get the angle wrt to the beam
    beamdir = {0 , 0 , 1};

    double angle = BeamCoords.Angle(beamdir) * 180 / 3.1415926;

    return angle;
}

double GeneratorReader::SetPPFXCVWeight(double _nu_e, double _nu_ang, int _nu_pdg){
    
    float weight = 1.0;

    double xbin{1.0},ybin{1.0};

    if (_nu_pdg == 14) {
        xbin =  hist_ratio.at(k_numu)->GetXaxis()->FindBin(_nu_e);
        ybin =  hist_ratio.at(k_numu)->GetYaxis()->FindBin(_nu_ang);
        weight =  hist_ratio.at(k_numu)->GetBinContent(xbin, ybin);
    }
    if (_nu_pdg == -14) {
        xbin =  hist_ratio.at(k_numubar)->GetXaxis()->FindBin(_nu_e);
        ybin = hist_ratio.at(k_numubar)->GetYaxis()->FindBin(_nu_ang);
        weight = hist_ratio.at(k_numubar)->GetBinContent(xbin, ybin);
    }
    if (_nu_pdg == 12) {
        xbin = hist_ratio.at(k_nue)->GetXaxis()->FindBin(_nu_e);
        ybin = hist_ratio.at(k_nue)->GetYaxis()->FindBin(_nu_ang);
        weight = hist_ratio.at(k_nue)->GetBinContent(xbin, ybin);
    }
    if (_nu_pdg == -12) {
        xbin = hist_ratio.at(k_nuebar)->GetXaxis()->FindBin(_nu_e);
        ybin = hist_ratio.at(k_nuebar)->GetYaxis()->FindBin(_nu_ang);
        weight = hist_ratio.at(k_nuebar)->GetBinContent(xbin, ybin);
    }

    // Add some catches to remove unphysical weights
    if (std::isinf(weight))      weight = 1.0; 
    if (std::isnan(weight) == 1) weight = 1.0;
    if (weight > 100)            weight = 1.0;

    return weight;

}

//______________________________________________________________________________
DEFINE_ART_MODULE(GeneratorReader)