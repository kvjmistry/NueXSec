// Script to merge the mc, data, ext and dirt ttrees to one file

void merge_uneaventrees(std::string run_type, std::string mc, std::string data, std::string ext, std::string dirt) {

    enum types {
        k_mc,
        k_data,
        k_ext,
        k_dirt,
        k_types_MAX
    };

    std::vector<TFile*> files(k_types_MAX);
    std::vector<TTree*> trees(k_types_MAX);
    std::vector<std::string> filenames{mc, data, ext, dirt};

    // Tree variables
    int run{0}, subrun{0}, event{0}, _run{0}, _subrun{0}, _event{0};;
    std::string classifcation, *_classifcation = NULL;
    bool gen{false}, _gen{false};           
    double weight{0.0}, _weight{0.0}; 
    double true_energy, _true_energy;
    double reco_energy, _reco_energy;

    std::vector<std::string> treenames{"mc_tree", "data_tree", "ext_tree", "dirt_tree"};

    TFile *outfile = new TFile(Form("./files/trees/nuexsec_tree_merged_run%s.root", run_type.c_str()), "UPDATE");
    TTree* outtree = new TTree("tree", "UPDATE");

    outtree->Branch("run",    &run,    "run/I");
    outtree->Branch("subrun", &subrun, "subrun/I");
    outtree->Branch("event",  &event,  "event/I");
    outtree->Branch("gen",    &gen,    "gen/O");
    outtree->Branch("weight", &weight, "weight/D");
    outtree->Branch("true_energy", &true_energy, "true_energy/D");
    outtree->Branch("reco_energy", &reco_energy, "reco_energy/D");
    outtree->Branch("classifcation",   &classifcation);

    // Loop over the flles
    for (unsigned int k = 0; k < files.size(); k++){
        files.at(k) = TFile::Open(filenames.at(k).c_str());
        if (files.at(k) == NULL){
            std::cout << "Couldn't get the file: " << filenames.at(k) << std::endl;
        }

        // Get the tree
        trees.at(k) = (TTree*)files.at(k)->Get(treenames.at(k).c_str());

        trees.at(k)->SetBranchAddress("run",    &_run);
        trees.at(k)->SetBranchAddress("subrun", &_subrun);
        trees.at(k)->SetBranchAddress("event",  &_event);
        trees.at(k)->SetBranchAddress("gen",    &_gen);
        trees.at(k)->SetBranchAddress("weight", &_weight);
        trees.at(k)->SetBranchAddress("true_energy", &_true_energy);
        trees.at(k)->SetBranchAddress("reco_energy", &_reco_energy);
        trees.at(k)->SetBranchAddress("classifcation",   &_classifcation);

        int tree_entries = trees.at(k)->GetEntries();

        // Event loop
        for (int ievent = 0; ievent < tree_entries; ievent++){
                trees.at(k)->GetEntry(ievent); 

                run = _run;
                subrun = _subrun;
                event = _event;
                classifcation = *_classifcation;
                weight = _weight;
                true_energy = _true_energy;
                reco_energy = _reco_energy;

                outtree->Fill();
        }

    }

    outfile->cd();
    outtree->Write("", TObject::kOverwrite);
    outfile->Close();

    std::cout << "Merged output ttree woo!!" << std::endl;

}