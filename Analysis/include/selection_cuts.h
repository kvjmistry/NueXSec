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
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

}; // End Class Selection Cuts

#endif
