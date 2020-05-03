#ifndef SELECTION_CUTS_h
#define SELECTION_CUTS_h

#include "utility.h"
#include "SliceContainer.h"

// Class for applying the selection cuts. Cut classes will inherit
// from this one and override function
class selection_cuts{

    private:

    public:
    selection_cuts(){}; // Default constructor

    utility _util;
    // -------------------------------------------------------------------------
    void Initalise(utility _utility);
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
    // Cut on the cluster fraction in the slice
    bool cluster_frac(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Cut on the shower score
    bool shower_score(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Cut on the total shower energy to reject michels
    bool michel_rej(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Cut on the dEdx y plane
    bool dEdx_y(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Cut on the dEdx y plane
    bool dEdx_v(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Cut on the dEdx y plane
    bool dEdx_u(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Cut on the total hits for the leading shower
    bool shr_hits(SliceContainer &SC);
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
    // Is the leading shower contained
    bool shr_contained(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // Shower Moliere Average
    bool shr_moliere_avg(SliceContainer &SC);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------


}; // End Class Selection Cuts

#endif
