////////////////////////////////////////////////////////////////////////
// Class:       FlashValidation
// Plugin Type: analyzer (art v3_01_02)
// File:        FlashValidation_module.cc
//
// Generated at Fri Dec  6 09:53:59 2019 by Krishan Mistry using cetskelgen
// from cetlib version v3_05_01.
// A module to get the flash information from a sample to plot
////////////////////////////////////////////////////////////////////////

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

#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"


class FlashValidation;


class FlashValidation : public art::EDAnalyzer {
public:
  explicit FlashValidation(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FlashValidation(FlashValidation const&) = delete;
  FlashValidation(FlashValidation&&) = delete;
  FlashValidation& operator=(FlashValidation const&) = delete;
  FlashValidation& operator=(FlashValidation&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endSubRun(art::SubRun const &sr) override;

private:

  // Optical 
  TTree * optical_tree;
  int fEvent;
  int fRun;
  int fSubRun;
  int fOpFlashPE;
  double fOpFlashTime;
  double fOpFlashWidthY;
  double fOpFlashWidthZ;
  double fOpFlashCenterY;
  double fOpFlashCenterZ;
  int fHasRecoFlash;


  // POT 
  TTree * pot_tree;
  double pot = 0.0;
  bool run_pot_counting = false;

  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot;
  std::ofstream _run_subrun_list_file;

  std::string _potsum_producer_mc;
  std::string _potsum_producer_data;
  std::string _potsum_instance;

  // Other
  int iteration{0}; // index to count number of entries run over
  std::string _mode;
  bool _is_mc{false};
  bool _is_data{false};
  bool _is_overlay{false};
  bool _cosmic_only{false};

  // Histograms
  TH1D* h_flash_time;
  TH1D* h_flash_pe;


  bool shift_times{true};

};


FlashValidation::FlashValidation(fhicl::ParameterSet const& p) : EDAnalyzer{p} {
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  _potsum_producer_mc    = p.get<std::string>("POTSummaryProducerMC");
  _potsum_producer_data  = p.get<std::string>("POTSummaryProducerData");
  _potsum_instance       = p.get<std::string>("POTSummaryInstance");
  _mode                  = p.get<std::string>("Mode");
  shift_times            = p.get<bool>("Shift_Times"); 

}

void FlashValidation::beginJob() {
  // Implementation of optional member function here.

  art::ServiceHandle<art::TFileService> fs;

  // ----------------------- OPTICAL TREE ------------------------------------
  optical_tree = fs->make<TTree>("optical_tree", "optical_objects");
  optical_tree->Branch("Event",            &fEvent,               "fEvent/I"         );
  optical_tree->Branch("Run",              &fRun,                 "fRun/I"           );
  optical_tree->Branch("SubRun",           &fSubRun,              "fSubRun/I"        );
  optical_tree->Branch("OpFlashPE",        &fOpFlashPE,           "fOpFlashPE/I"     );
  optical_tree->Branch("OpFlashTime",      &fOpFlashTime,         "fOpFlashTime/D"   );
  optical_tree->Branch("OpFlashWidhtY",    &fOpFlashWidthY,       "fOpFlashWidthY/D" );
  optical_tree->Branch("OpFlashWidthZ",    &fOpFlashWidthZ,       "fOpFlashWidthZ/D" );
  optical_tree->Branch("OpFlashCenterY",   &fOpFlashCenterY,      "fOpFlashCenterY/D");
  optical_tree->Branch("OpFlashCenterZ",   &fOpFlashCenterZ,      "fOpFlashCenterZ/D");
  optical_tree->Branch("HasRecoFlash",     &fHasRecoFlash,        "fHasRecoFlash/I"  );


  pot_tree = fs->make<TTree>("pot_tree", "pot_per_subrun");
  pot_tree->Branch("pot", &pot, "pot/D");

  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run",       &_sr_run,       "run/I"      );
  _sr_tree->Branch("subrun",    &_sr_subrun,    "subrun/I"   );
  _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree->Branch("endtime",   &_sr_endtime,   "endtime/D"  );
  _sr_tree->Branch("pot",       &_sr_pot,       "pot/D"      );
  _run_subrun_list_file.open ("run_subrun_list.txt", std::ofstream::out | std::ofstream::trunc);

  h_flash_time = fs->make<TH1D>("h_flash_time", "Flash Time; Flash Time [us]; Entries",100, 0, 25);
  h_flash_pe = fs->make<TH1D>("h_flash_pe", "Flash PE; Flash PE; Entries",100, 0, 10000);

}

void FlashValidation::analyze(art::Event const& e) {
  
  // Implementation of required member function here.
  int run      = e.id().run();
  int event    = e.id().event();
  int subrun   = e.id().subRun();
    
  if      (_mode == "EXT")     _cosmic_only = true;
  else if (_mode == "Data")    _is_data     = true;
  else if (_mode == "Overlay"){
    _is_mc       = true;
    _is_overlay  = true;
  } 
  else _is_mc = true;

  std::cout << "[Analyze] ------------------------------------------- [Analyze]" << std::endl;
  std::cout << "[Analyze] Running over entry: " << iteration << std::endl;
  iteration++;

  if(_cosmic_only == true) {std::cout << "[Analyze] Running in Cosmic Only Configuration! " << std::endl;}
  
  if(_is_mc == true)       {std::cout << "[Analyze] Running with Monte Carlo " << std::endl;}
  
  if(_is_data == true)     {std::cout << "[Analyze] Running with Data " << std::endl;}
  
  if(_is_data == true)     {run_pot_counting = false; std::cout << "[Analyze] Do Not Count MC POT" << std::endl;}
  
  if(_is_mc == true)       {run_pot_counting = true;  std::cout << "[Analyze] Count MC POT" << std::endl;}

  // Get the optical information --  Maybe make them filled at the same place as the other - so it's a per event
  std::string beam_flash_tag = "simpleFlashBeam";
  auto const & beam_opf = e.getValidHandle<std::vector < recob::OpFlash> >(beam_flash_tag);
  auto const & beam_opflashes(*beam_opf);
  std::cout << "[Analyze] [OPTICAL] " << beam_flash_tag << " in this event: " << beam_opflashes.size() << std::endl;

  // If there is no optical activity in this event then I have to check where the true nu vtx is
  // If it's a true nue without reco optical event, but interacts in the FV then this is a loss in efficiency!
  if(beam_opflashes.size() == 0) {
      std::cout << "[Analyze] [Optical] No Optical Activity in this Event!" << std::endl;
      fHasRecoFlash   = 0;
      fEvent          = event;
      fRun            = run;
      fSubRun         = subrun;
      fOpFlashPE      = -999;
      fOpFlashTime    = -999;
      fOpFlashWidthY  = -999;
      fOpFlashWidthZ  = -999;
      fOpFlashCenterY = -999;
      fOpFlashCenterZ = -999;
      optical_tree->Fill();
  }

  for (auto const & opflsh : beam_opflashes) {
      fHasRecoFlash   = 1;
      fEvent          = event;
      fRun            = run;
      fSubRun         = subrun;
      fOpFlashPE      = opflsh.TotalPE();
      fOpFlashTime    = opflsh.Time();
      fOpFlashWidthY  = opflsh.YWidth();
      fOpFlashWidthZ  = opflsh.ZWidth();
      fOpFlashCenterY = opflsh.YCenter();
      fOpFlashCenterZ = opflsh.ZCenter();
      optical_tree->Fill();

      // Shift the ext by -359 ns to match beam on
      if (_cosmic_only && shift_times)  fOpFlashTime+=-0.359;

      // Shift the MC by +55ns to ma
      if (_is_mc && shift_times)        fOpFlashTime+=0.055;

      // Fill histogram if > 50 PE
      if (fOpFlashPE > 50) h_flash_time->Fill(fOpFlashTime);
      //h_flash_time->Fill(fOpFlashTime);
      //

      // Now fill flash PE if > 8 us (cuts out the ext)
      if (fOpFlashTime > 8) h_flash_pe->Fill(fOpFlashPE);
      if (fOpFlashTime > 8) std::cout <<"Flash PE: " << fOpFlashPE << std::endl;
  }

}

void FlashValidation::endSubRun(art::SubRun const & sr) {

  bool _debug = false;

    if(run_pot_counting == true) {
        auto const & POTSummaryHandle = sr.getValidHandle < sumdata::POTSummary >("generator");
        auto const & POTSummary(*POTSummaryHandle);
        const double total_pot = POTSummary.totpot;
        std::cout << "----------------------------" << std::endl;
        std::cout << "Total POT / subRun: " << total_pot << std::endl;
        std::cout << "----------------------------" << std::endl;

        pot = total_pot;
        pot_tree->Fill();
    }
    
    if (_debug) std::cout << "[Analysis::endSubRun] Starts" << std::endl;

    // Saving run and subrun number on file so that we can run Zarko's script easily
    _run_subrun_list_file << sr.run() << " " << sr.subRun() << std::endl;

    _sr_run       = sr.run();
    _sr_subrun    = sr.subRun();
    _sr_begintime = sr.beginTime().value();
    _sr_endtime   = sr.endTime().value();

    art::Handle<sumdata::POTSummary> potsum_h;

    // MC
    if (_is_mc) {
        
        if (_debug) std::cout << "[Analysis::endSubRun] Getting POT for MC" << std::endl;
        
        if(sr.getByLabel(_potsum_producer_mc, potsum_h)) {
            if (_debug) std::cout << "[Analysis::endSubRun] POT are valid" << std::endl;
            _sr_pot = potsum_h->totpot;
        }
        else
            _sr_pot = 0.;
    }

    // Data
    if (_is_data) {
        
        if (_debug) std::cout << "[Analysis::endSubRun] Getting POT for DATA, producer " << _potsum_producer_data << ", instance " << _potsum_instance << std::endl;
        
        if (sr.getByLabel(_potsum_producer_data, _potsum_instance, potsum_h)) {
            if (_debug) std::cout << "[Analysis::endSubRun] POT are valid" << std::endl;
            _sr_pot = potsum_h->totpot;
        }
        else
            _sr_pot = 0;
    
    }

    _sr_tree->Fill();

    if (_debug) std::cout << "[Analysis::endSubRun] Ends" << std::endl;
}

DEFINE_ART_MODULE(FlashValidation)
