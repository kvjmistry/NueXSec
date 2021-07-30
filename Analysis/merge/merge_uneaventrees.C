// Script to merge the mc, data, ext and dirt ttrees to one file

void merge_uneaventrees(std::string run_type, bool intrinsic_mode, bool fake_intrinsic_mode, std::string mc, std::string data, std::string ext, std::string dirt, std::string detvar) {

    enum types {
        k_mc,
        k_data,
        k_ext,
        k_dirt,
        k_types_MAX
    };

    int resize;

    // If intrinsic mode, we have another tree to loop over
    if (intrinsic_mode){
        resize = k_types_MAX + 1;
    }
    else
        resize = k_types_MAX;

    std::vector<TFile*> files(resize);
    std::vector<TTree*> trees(resize);
    std::vector<std::string> filenames;
    
    // Using intrinsic nue in MC only
    if (intrinsic_mode && !fake_intrinsic_mode){
        std::cout << "MC Intrinsic Only"<< std::endl;
        filenames = {mc, mc, data, ext, dirt};
    }
    // Using intrinsic nue in fake data only
    else if (!intrinsic_mode && fake_intrinsic_mode){
        std::cout << "Fake Data Intrinsic Only"<< std::endl;
        filenames = {mc, data, data, ext, dirt};
    }
    // Using intrinsic nue in MC and fake data
    else if (intrinsic_mode && fake_intrinsic_mode){
        std::cout << "MC and Data Intrinsic"<< std::endl;
        filenames = {mc, mc, data, data, ext, dirt};
    }
    else
        filenames = {mc, data, ext, dirt};
    

    // Tree variables
    int run{0}, subrun{0}, event{0}, _run{0}, _subrun{0}, _event{0};;
    std::string classification, *_classification = NULL;
    bool gen{false},                      _gen{false};
    bool passed_selection{false},         _passed_selection{false};           
    double weight{0.0},                   _weight{0.0}; 
    double true_energy,                   _true_energy;
    double reco_energy,                   _reco_energy;
    float shr_tkfit_dedx_Y{0.0},          _shr_tkfit_dedx_Y{0.0};
    float n_showers{0},                   _n_showers{0};
    float n_tracks{0},                    _n_tracks{0};
    float shr_theta{0.0},                 _shr_theta{0.0};
    float shr_phi{0.0},                   _shr_phi{0.0};
    float shr_energy_cali{0.0},           _shr_energy_cali{0.0};
    float shrmoliereavg{0.0},             _shrmoliereavg{0.0};
    float shr_hits_max{0.0},              _shr_hits_max{0.0};
    float elec_e{0.0},                    _elec_e{0.0};
    float ppfx_cv{1.0},                   _ppfx_cv{0.0};
    float weightSplineTimesTune{1.0},     _weightSplineTimesTune{1.0};
    float numi_ang{0.0},                  _numi_ang{0.0};
    int nu_pdg{0},                        _nu_pdg{0};
    float shr_bkt_purity{0.0},            _shr_bkt_purity{0.0};
    float shr_bkt_completeness{0.0},      _shr_bkt_completeness{0.0};
    float shr_bkt_E{0.0},                 _shr_bkt_E{0.0};
    int shr_bkt_pdg{0},                   _shr_bkt_pdg{0};
    int npi0{0},                          _npi0{0};
    double pi0_e{0.0},                    _pi0_e{0};
    int interaction{0},                   _interaction{0};
    double effective_angle{0.0},          _effective_angle{0.0};
    double cos_effective_angle{0.0},      _cos_effective_angle{0.0};
    double true_effective_angle{0.0},     _true_effective_angle{0.0};
    double cos_true_effective_angle{0.0}, _cos_true_effective_angle{0.0};
    std::vector<float> all_shr_hits;
    std::vector<float> all_shr_energies;
    std::vector<float> *_all_shr_hits = NULL;
    std::vector<float> *_all_shr_energies = NULL;
    int ccnc{0}, _ccnc{0};

    std::vector<unsigned short> weightsGenie;
    std::vector<unsigned short> weightsReint;
    std::vector<unsigned short> weightsFlux ;
    std::vector<unsigned short> *_weightsGenie  = NULL;
    std::vector<unsigned short> *_weightsReint  = NULL;
    std::vector<unsigned short> *_weightsFlux   = NULL;
    
    double knobRPAup{0.0},              _knobRPAup{0.0};
    double knobCCMECup{0.0},            _knobCCMECup{0.0};
    double knobAxFFCCQEup{0.0},         _knobAxFFCCQEup{0.0};
    double knobVecFFCCQEup{0.0},        _knobVecFFCCQEup{0.0};
    double knobDecayAngMECup{0.0},      _knobDecayAngMECup{0.0};
    double knobThetaDelta2Npiup{0.0},   _knobThetaDelta2Npiup{0.0};
    double knobThetaDelta2NRadup{0.0},  _knobThetaDelta2NRadup{0.0};
    double knobRPA_CCQE_Reducedup{0.0}, _knobRPA_CCQE_Reducedup{0.0};
    double knobNormCCCOHup{0.0},        _knobNormCCCOHup{0.0};
    double knobNormNCCOHup{0.0},        _knobNormNCCOHup{0.0};
    double knobxsr_scc_Fv3up{0.0},      _knobxsr_scc_Fv3up{0.0};
    double knobxsr_scc_Fa3up{0.0},      _knobxsr_scc_Fa3up{0.0};
    double knobRPAdn{0.0},              _knobRPAdn{0.0};
    double knobCCMECdn{0.0},            _knobCCMECdn{0.0};
    double knobAxFFCCQEdn{0.0},         _knobAxFFCCQEdn{0.0};
    double knobVecFFCCQEdn{0.0},        _knobVecFFCCQEdn{0.0};
    double knobDecayAngMECdn{0.0},      _knobDecayAngMECdn{0.0};
    double knobThetaDelta2Npidn{0.0},   _knobThetaDelta2Npidn{0.0};
    double knobThetaDelta2NRaddn{0.0},  _knobThetaDelta2NRaddn{0.0};
    double knobRPA_CCQE_Reduceddn{0.0}, _knobRPA_CCQE_Reduceddn{0.0};
    double knobNormCCCOHdn{0.0},        _knobNormCCCOHdn{0.0};
    double knobNormNCCOHdn{0.0},        _knobNormNCCOHdn{0.0};
    double knobxsr_scc_Fv3dn{0.0},      _knobxsr_scc_Fv3dn{0.0};
    double knobxsr_scc_Fa3dn{0.0},      _knobxsr_scc_Fa3dn{0.0};



    std::vector<std::string> treenames;
    
    if (intrinsic_mode && !fake_intrinsic_mode){
        treenames = {"mc_nue_tree", "mc_tree", "data_tree", "ext_tree", "dirt_tree"};
    }
    else if (!intrinsic_mode && fake_intrinsic_mode){
        treenames = {"mc_tree", "data_nue_tree", "data_tree", "ext_tree", "dirt_tree"};
    }
    else if (intrinsic_mode && fake_intrinsic_mode){
        treenames = {"mc_nue_tree", "mc_tree", "data_nue_tree", "data_tree", "ext_tree", "dirt_tree"};
    }
    else
        treenames = {"mc_tree", "data_tree", "ext_tree", "dirt_tree"};
    
    

    TFile *outfile;
    
    if (detvar == "") outfile = new TFile(Form("./files/trees/nuexsec_tree_merged_run%s.root", run_type.c_str()), "UPDATE");
    else outfile = new TFile(Form("./files/trees/nuexsec_tree_merged_run%s_%s.root", run_type.c_str(), detvar.c_str()), "UPDATE");
    
    TTree* outtree = new TTree("tree", "UPDATE");

    outtree->Branch("run",              &run,    "run/I");
    outtree->Branch("subrun",           &subrun, "subrun/I");
    outtree->Branch("event",            &event,  "event/I");
    outtree->Branch("gen",              &gen,    "gen/O");
    outtree->Branch("passed_selection", &passed_selection,    "passed_selection/O");
    outtree->Branch("weight",           &weight, "weight/D");
    outtree->Branch("true_energy",      &true_energy, "true_energy/D");
    outtree->Branch("reco_energy",      &reco_energy, "reco_energy/D");
    outtree->Branch("classification",   &classification);
    outtree->Branch("shr_tkfit_dedx_Y", &shr_tkfit_dedx_Y);
    outtree->Branch("n_showers",        &n_showers);
    outtree->Branch("n_tracks",         &n_tracks);
    outtree->Branch("shr_theta",        &shr_theta);
    outtree->Branch("shr_phi",          &shr_phi);
    outtree->Branch("shr_energy_cali",  &shr_energy_cali);
    outtree->Branch("shrmoliereavg",    &shrmoliereavg);
    outtree->Branch("shr_hits_max",     &shr_hits_max);
    outtree->Branch("elec_e",           &elec_e,  "elec_e/F");
    outtree->Branch("ppfx_cv",          &ppfx_cv,  "ppfx_cv/F");
    outtree->Branch("weightSplineTimesTune",  &weightSplineTimesTune,  "weightSplineTimesTune/F");
    outtree->Branch("numi_ang",             &numi_ang,  "numi_ang/F");
    outtree->Branch("nu_pdg",               &nu_pdg,  "nu_pdg/I");
    outtree->Branch("shr_bkt_purity",       &shr_bkt_purity);
    outtree->Branch("shr_bkt_completeness", &shr_bkt_completeness);
    outtree->Branch("shr_bkt_E",            &shr_bkt_E);
    outtree->Branch("shr_bkt_pdg",          &shr_bkt_pdg);
    outtree->Branch("npi0",                 &npi0,  "npi0/I");
    outtree->Branch("pi0_e",                &pi0_e,  "pi0_e/D");
    outtree->Branch("all_shr_hits", "std::vector<float>",     &all_shr_hits);
    outtree->Branch("all_shr_energies", "std::vector<float>", &all_shr_energies);
    outtree->Branch("interaction",                            &interaction,  "interaction/I");
    outtree->Branch("effective_angle",          &effective_angle);
    outtree->Branch("cos_effective_angle",      &cos_effective_angle);
    outtree->Branch("true_effective_angle",     &true_effective_angle);
    outtree->Branch("cos_true_effective_angle", &cos_true_effective_angle);
    outtree->Branch("ccnc", &ccnc);
    
    outtree->Branch("weightsGenie", "std::vector<unsigned short>", &weightsGenie);
    outtree->Branch("weightsReint", "std::vector<unsigned short>", &weightsReint);
    outtree->Branch("weightsFlux", "std::vector<unsigned short>",  &weightsFlux);
    outtree->Branch("knobRPAup",             &knobRPAup);
    outtree->Branch("knobRPAdn",             &knobRPAdn);
    outtree->Branch("knobCCMECup",           &knobCCMECup);
    outtree->Branch("knobCCMECdn",           &knobCCMECdn);
    outtree->Branch("knobAxFFCCQEup",        &knobAxFFCCQEup);
    outtree->Branch("knobAxFFCCQEdn",        &knobAxFFCCQEdn);
    outtree->Branch("knobVecFFCCQEup",       &knobVecFFCCQEup);
    outtree->Branch("knobVecFFCCQEdn",       &knobVecFFCCQEdn);
    outtree->Branch("knobDecayAngMECup",     &knobDecayAngMECup);
    outtree->Branch("knobDecayAngMECdn",     &knobDecayAngMECdn);
    outtree->Branch("knobThetaDelta2Npiup",  &knobThetaDelta2Npiup);
    outtree->Branch("knobThetaDelta2Npidn",  &knobThetaDelta2Npidn);
    outtree->Branch("knobThetaDelta2NRadup", &knobThetaDelta2NRadup);
    outtree->Branch("knobThetaDelta2NRaddn", &knobThetaDelta2NRaddn);
    outtree->Branch("knobRPA_CCQE_Reducedup",&knobRPA_CCQE_Reducedup);
    outtree->Branch("knobRPA_CCQE_Reduceddn",&knobRPA_CCQE_Reduceddn);
    outtree->Branch("knobNormCCCOHup",       &knobNormCCCOHup);
    outtree->Branch("knobNormCCCOHdn",       &knobNormCCCOHdn);
    outtree->Branch("knobNormNCCOHup",       &knobNormNCCOHup);
    outtree->Branch("knobNormNCCOHdn",       &knobNormNCCOHdn);
    outtree->Branch("knobxsr_scc_Fv3up",     &knobxsr_scc_Fv3up);
    outtree->Branch("knobxsr_scc_Fv3dn",     &knobxsr_scc_Fv3dn);
    outtree->Branch("knobxsr_scc_Fa3up",     &knobxsr_scc_Fa3up);
    outtree->Branch("knobxsr_scc_Fa3dn",     &knobxsr_scc_Fa3dn);

    // Loop over the flles
    for (unsigned int k = 0; k < files.size(); k++){

        files.at(k) = TFile::Open(filenames.at(k).c_str());
        
        if (files.at(k) == NULL){
            std::cout << "Couldn't get the file: " << filenames.at(k) << std::endl;
        }

        // Get the tree
        trees.at(k) = (TTree*)files.at(k)->Get(treenames.at(k).c_str());
    
        trees.at(k)->SetBranchAddress("run",              &_run);
        trees.at(k)->SetBranchAddress("subrun",           &_subrun);
        trees.at(k)->SetBranchAddress("event",            &_event);
        trees.at(k)->SetBranchAddress("gen",              &_gen);
        trees.at(k)->SetBranchAddress("passed_selection", &_passed_selection);
        trees.at(k)->SetBranchAddress("weight",           &_weight);
        trees.at(k)->SetBranchAddress("true_energy",      &_true_energy);
        trees.at(k)->SetBranchAddress("reco_energy",      &_reco_energy);
        trees.at(k)->SetBranchAddress("classification",   &_classification);
        trees.at(k)->SetBranchAddress("shr_tkfit_dedx_Y", &_shr_tkfit_dedx_Y);
        trees.at(k)->SetBranchAddress("n_showers",        &_n_showers);
        trees.at(k)->SetBranchAddress("n_tracks",         &_n_tracks);
        trees.at(k)->SetBranchAddress("shr_theta",        &_shr_theta);
        trees.at(k)->SetBranchAddress("shr_phi",          &_shr_phi);
        trees.at(k)->SetBranchAddress("shr_energy_cali",  &_shr_energy_cali);
        trees.at(k)->SetBranchAddress("shrmoliereavg",    &_shrmoliereavg);
        trees.at(k)->SetBranchAddress("shr_hits_max",     &_shr_hits_max);
        trees.at(k)->SetBranchAddress("elec_e",           &_elec_e);
        trees.at(k)->SetBranchAddress("ppfx_cv",          &_ppfx_cv);
        trees.at(k)->SetBranchAddress("weightSplineTimesTune",  &_weightSplineTimesTune);
        trees.at(k)->SetBranchAddress("numi_ang",         &_numi_ang);
        trees.at(k)->SetBranchAddress("nu_pdg",           &_nu_pdg);
        trees.at(k)->SetBranchAddress("shr_bkt_purity",   &_shr_bkt_purity);
        trees.at(k)->SetBranchAddress("shr_bkt_completeness", &_shr_bkt_completeness);
        trees.at(k)->SetBranchAddress("shr_bkt_E",        &_shr_bkt_E);
        trees.at(k)->SetBranchAddress("shr_bkt_pdg",      &_shr_bkt_pdg);
        trees.at(k)->SetBranchAddress("all_shr_hits",     &_all_shr_hits);
        trees.at(k)->SetBranchAddress("all_shr_energies", &_all_shr_energies);
        trees.at(k)->SetBranchAddress("npi0",             &_npi0);
        trees.at(k)->SetBranchAddress("pi0_e",            &_pi0_e);
        trees.at(k)->SetBranchAddress("interaction",      &_interaction);
        trees.at(k)->SetBranchAddress("effective_angle",            &_effective_angle);
        trees.at(k)->SetBranchAddress("cos_effective_angle",        &_cos_effective_angle);
        trees.at(k)->SetBranchAddress("true_effective_angle",       &_true_effective_angle);
        trees.at(k)->SetBranchAddress("cos_true_effective_angle",   &_cos_true_effective_angle);
        trees.at(k)->SetBranchAddress("ccnc",   &_ccnc);
        
        trees.at(k)->SetBranchAddress("weightsGenie",          &_weightsGenie);
        trees.at(k)->SetBranchAddress("weightsReint",          &_weightsReint);
        trees.at(k)->SetBranchAddress("weightsFlux",           &_weightsFlux);
        trees.at(k)->SetBranchAddress("knobRPAup",             &_knobRPAup);
        trees.at(k)->SetBranchAddress("knobRPAdn",             &_knobRPAdn);
        trees.at(k)->SetBranchAddress("knobCCMECup",           &_knobCCMECup);
        trees.at(k)->SetBranchAddress("knobCCMECdn",           &_knobCCMECdn);
        trees.at(k)->SetBranchAddress("knobAxFFCCQEup",        &_knobAxFFCCQEup);
        trees.at(k)->SetBranchAddress("knobAxFFCCQEdn",        &_knobAxFFCCQEdn);
        trees.at(k)->SetBranchAddress("knobVecFFCCQEup",       &_knobVecFFCCQEup);
        trees.at(k)->SetBranchAddress("knobVecFFCCQEdn",       &_knobVecFFCCQEdn);
        trees.at(k)->SetBranchAddress("knobDecayAngMECup",     &_knobDecayAngMECup);
        trees.at(k)->SetBranchAddress("knobDecayAngMECdn",     &_knobDecayAngMECdn);
        trees.at(k)->SetBranchAddress("knobThetaDelta2Npiup",  &_knobThetaDelta2Npiup);
        trees.at(k)->SetBranchAddress("knobThetaDelta2Npidn",  &_knobThetaDelta2Npidn);
        trees.at(k)->SetBranchAddress("knobThetaDelta2NRadup", &_knobThetaDelta2NRadup);
        trees.at(k)->SetBranchAddress("knobThetaDelta2NRaddn", &_knobThetaDelta2NRaddn);
        trees.at(k)->SetBranchAddress("knobRPA_CCQE_Reducedup",&_knobRPA_CCQE_Reducedup);
        trees.at(k)->SetBranchAddress("knobRPA_CCQE_Reduceddn",&_knobRPA_CCQE_Reduceddn);
        trees.at(k)->SetBranchAddress("knobNormCCCOHup",       &_knobNormCCCOHup);
        trees.at(k)->SetBranchAddress("knobNormCCCOHdn",       &_knobNormCCCOHdn);
        trees.at(k)->SetBranchAddress("knobNormNCCOHup",       &_knobNormNCCOHup);
        trees.at(k)->SetBranchAddress("knobNormNCCOHdn",       &_knobNormNCCOHdn);
        trees.at(k)->SetBranchAddress("knobxsr_scc_Fv3up",     &_knobxsr_scc_Fv3up);
        trees.at(k)->SetBranchAddress("knobxsr_scc_Fv3dn",     &_knobxsr_scc_Fv3dn);
        trees.at(k)->SetBranchAddress("knobxsr_scc_Fa3up",     &_knobxsr_scc_Fa3up);
        trees.at(k)->SetBranchAddress("knobxsr_scc_Fa3dn",     &_knobxsr_scc_Fa3dn);

        int tree_entries = trees.at(k)->GetEntries();

        // Event loop
        for (int ievent = 0; ievent < tree_entries; ievent++){
                trees.at(k)->GetEntry(ievent); 

                // If intrinsic mode is turned on and we are using the default mc, then skip the fill for generated events which are already covered
                if (intrinsic_mode && !fake_intrinsic_mode){
                     if (k == 1 && _gen)
                        continue;
                }
                else if (!intrinsic_mode && fake_intrinsic_mode){
                     if (k == 2 && _gen)
                        continue;
                }
                else if (intrinsic_mode && fake_intrinsic_mode){
                     if (( k == 1 || k == 3) && _gen)
                        continue;
                }

                run              = _run;
                subrun           = _subrun;
                event            = _event;
                gen              = _gen;
                passed_selection = _passed_selection;
                classification   = *_classification;
                weight           = _weight;
                true_energy      = _true_energy;
                reco_energy      = _reco_energy;
                shr_tkfit_dedx_Y = _shr_tkfit_dedx_Y;
                n_showers        = _n_showers;
                n_tracks         = _n_tracks;
                shr_theta        = _shr_theta;
                shr_phi          = _shr_phi;
                shr_energy_cali  = _shr_energy_cali;
                shrmoliereavg    = _shrmoliereavg;
                shr_hits_max     = _shr_hits_max;
                elec_e           = _elec_e;
                ppfx_cv          = _ppfx_cv;
                weightSplineTimesTune = _weightSplineTimesTune;
                numi_ang         = _numi_ang;
                nu_pdg           = _nu_pdg;
                shr_bkt_pdg      = _shr_bkt_pdg;
                shr_bkt_purity   = _shr_bkt_purity;
                shr_bkt_completeness = _shr_bkt_completeness;
                shr_bkt_E        = _shr_bkt_E;
                npi0             = _npi0;
                pi0_e            = _pi0_e;
                interaction      = _interaction;
                effective_angle          = _effective_angle;
                cos_effective_angle      = _cos_effective_angle;
                true_effective_angle     = _true_effective_angle;
                cos_true_effective_angle = _cos_true_effective_angle;
                ccnc = _ccnc;

                if (_all_shr_hits != NULL)    all_shr_hits = *_all_shr_hits;
                if (_all_shr_energies != NULL)all_shr_energies = *_all_shr_energies;

                if (_weightsGenie != NULL) weightsGenie           = *_weightsGenie; // If these aren't set by default then bad things happen in memory land
                if (_weightsReint != NULL) weightsReint           = *_weightsReint;
                if (_weightsFlux  != NULL) weightsFlux            = *_weightsFlux;
                knobRPAup =              _knobRPAup;
                knobCCMECup =            _knobCCMECup;
                knobAxFFCCQEup =         _knobAxFFCCQEup;
                knobVecFFCCQEup =        _knobVecFFCCQEup;
                knobDecayAngMECup =      _knobDecayAngMECup;
                knobThetaDelta2Npiup =   _knobThetaDelta2Npiup;
                knobThetaDelta2NRadup =  _knobThetaDelta2NRadup;
                knobRPA_CCQE_Reducedup = _knobRPA_CCQE_Reducedup;
                knobNormCCCOHup =        _knobNormCCCOHup;
                knobNormNCCOHup =        _knobNormNCCOHup;
                knobRPAdn =              _knobRPAdn;
                knobCCMECdn =            _knobCCMECdn;
                knobAxFFCCQEdn =         _knobAxFFCCQEdn;
                knobVecFFCCQEdn =        _knobVecFFCCQEdn;
                knobDecayAngMECdn =      _knobDecayAngMECdn;
                knobThetaDelta2Npidn =   _knobThetaDelta2Npidn;
                knobThetaDelta2NRaddn =  _knobThetaDelta2NRaddn;
                knobRPA_CCQE_Reduceddn = _knobRPA_CCQE_Reduceddn;
                knobNormCCCOHdn =        _knobNormCCCOHdn;
                knobNormNCCOHdn =        _knobNormNCCOHdn;
                knobxsr_scc_Fv3up =      _knobxsr_scc_Fv3up;
                knobxsr_scc_Fv3dn =      _knobxsr_scc_Fv3dn;
                knobxsr_scc_Fa3up =      _knobxsr_scc_Fa3up;
                knobxsr_scc_Fa3dn =      _knobxsr_scc_Fa3dn;

                outtree->Fill();
        }

    }

    outfile->cd();
    outtree->Write("", TObject::kOverwrite);
    outfile->Close();

    std::cout << "\nMerged output ttree woo!!\n" << std::endl;

}
