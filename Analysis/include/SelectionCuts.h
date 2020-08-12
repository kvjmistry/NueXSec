#ifndef SELECTIONCUTS_H
#define SELECTIONCUTS_H

#include "Utility.h"
#include "SliceContainer.h"

// Class for applying the selection cuts. Cut classes will inherit
// from this one and override function
class SelectionCuts{

    private:

    public:
    SelectionCuts(){}; // Default constructor

    Utility _util;
    // -------------------------------------------------------------------------
    void Initalise(Utility _utility);
    // -------------------------------------------------------------------------
    // Software Trigger -- MC Only
    bool swtrig(SliceContainer &SC, int type);
    // -------------------------------------------------------------------------
    // Common Optical Filter PE -- MC Only
    bool opfilt_pe(SliceContainer &SC, int type);
    // -------------------------------------------------------------------------
    // Common Optical Filter Michel Veto -- MC Only
    bool opfilt_veto(SliceContainer &SC, int type);
    // -------------------------------------------------------------------------
    // Require pass for the slice ID
    bool slice_id(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Require the space charge corrected reco nu vertex in the FV
    bool in_fv(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Electron Candidate Cut
    bool e_candidate(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Cut on the topological score
    bool topo_score(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Cut on the shower score
    bool shower_score(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Cut on the total shower energy to reject michels
    bool michel_rej(SliceContainer &SC);   
    // -------------------------------------------------------------------------
    // Cut on the shower to nu vertex distance
    bool shr_distance(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Cut on the ratio of shower hits to slice hits
    bool shr_hitratio(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Impact parameter of object thats not a pfp to the slice
    bool shr_cosmic_IP(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Slice Contained Fraction
    bool contained_frac(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Shower Moliere Average
    bool shr_moliere_avg(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // 2D cut in dEdx and shower vertex distance
    bool shr_dist_dEdx_max(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // dEdx cut in the case of no tracks
    bool dEdx_max_no_tracks(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Pi0 Selection Cuts
    bool pi_zero_cuts(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // NuMu Selection Cuts
    bool numu_cuts(SliceContainer &SC);
    // -------------------------------------------------------------------------


}; // End Class Selection Cuts

#endif
