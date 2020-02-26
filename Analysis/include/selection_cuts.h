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
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

}; // End Class Selection Cuts

#endif
