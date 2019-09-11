#include "RecoTrueHelper.h"

namespace nue_xsec
{

void recotruehelper::Configure(art::Event const & e,
                               std::string _pfp_producer,
                               std::string _spacepointLabel,
                               std::string _hitfinderLabel,
                               std::string _geantModuleLabel) {

    // Collect hits
    lar_pandora::HitVector hitVector;
    lar_pandora::LArPandoraHelper::CollectHits(e, _hitfinderLabel, hitVector);

    // Collect PFParticles and match Reco Particles to Hits
    lar_pandora::PFParticleVector recoParticleVector;
    lar_pandora::PFParticleVector recoNeutrinoVector;
    lar_pandora::PFParticlesToHits recoParticlesToHits;
    lar_pandora::HitsToPFParticles recoHitsToParticles;

    lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
    lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel, recoParticlesToHits,
                                                          recoHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);

    if (_verbose) {
        std::cout << "[RecoTrueHelper] [Configure] RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
        std::cout << "[RecoTrueHelper] [Configure] RecoParticles: " << recoParticleVector.size() << std::endl;
    }

    // Collect MCParticles and match True Particles to Hits
    lar_pandora::MCParticleVector trueParticleVector;
    lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;
    lar_pandora::MCParticlesToHits trueParticlesToHits;
    lar_pandora::HitsToMCParticles trueHitsToParticles;

    if (!e.isRealData()) {
        lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, trueParticleVector);
        lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, truthToParticles, particlesToTruth);
        lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(e, _geantModuleLabel, hitVector, trueParticlesToHits,
                                                              trueHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);
    }

    if (_verbose) {
        std::cout << "[RecoTrueHelper] [Configure] TrueParticles: " << particlesToTruth.size() << std::endl;
        std::cout << "[RecoTrueHelper] [Configure] TrueEvents: " << truthToParticles.size() << std::endl;
    }

    // Now set the things we need for the future
    _trueHitsToParticles = trueHitsToParticles;
    _recoParticlesToHits = recoParticlesToHits;

    if (_debug) { // yes, don't do it
        std::cout << "[McPfpMatch] This is event " << e.id().run() << std::endl;
        //art::ServiceHandle<cheat::BackTracker> bt;
        std::cout << "[McPfpMatch] Number of MCParticles matched to hits: " << trueParticlesToHits.size() << std::endl;
        for (const auto & iter : trueParticlesToHits) {
            const art::Ptr<simb::MCTruth> mc_truth = TrackIDToMCTruth(e, _geantModuleLabel, (iter.first)->TrackId());
            //bt->TrackIDToMCTruth((iter.first)->TrackId());
            std::cout << "[McPfpMatch] MCParticle with pdg " << (iter.first)->PdgCode()
                      << " and origin " << (mc_truth->Origin() == 1 ? "neutrino" : "cosmic")
                      << " has " << (iter.second).size() << " hits ass." << std::endl;
            if (mc_truth->Origin() == 1) {
                lar_pandora::HitVector hits = (iter.second);
                for (const auto & hit : hits) {
                    std::cout << "[McPfpMatch]   > Hit on plane " << hit->View()
                              << " on wire " << hit->WireID()
                              << " with time " << hit->PeakTime() << std::endl;
                }
            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------
void recotruehelper::Configure(art::Event const & e,
                               std::string _pfp_producer,
                               std::string _spacepointLabel,
                               std::string _hitfinderLabel,
                               std::string _geantModuleLabel,
                               std::string _hit_mcp_producer,
                               lar_pandora::LArPandoraHelper::DaughterMode daughterMode) {

    // Collect hits
    lar_pandora::HitVector hitVector;
    //lar_pandora::LArPandoraHelper::CollectHits(e, _hitfinder_producer, hitVector);
    
    // Get the hits from the event
    art::Handle<std::vector<recob::Hit> > hit_h;
    e.getByLabel(_hitfinderLabel, hit_h);
    
    if (!hit_h.isValid()) {
        std::cout << "[RecoTrueHelper] [Configure] Hit Handle is not valid." << std::endl;
        throw std::exception();
    }
    
    art::fill_ptr_vector(hitVector, hit_h);

    // Collect PFParticles and match Reco Particles to Hits
    lar_pandora::PFParticleVector  recoParticleVector;
    lar_pandora::PFParticleVector  recoNeutrinoVector;
    lar_pandora::PFParticlesToHits recoParticlesToHits;
    lar_pandora::HitsToPFParticles recoHitsToParticles;

    // Get the reco PFP from the event
    lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
    
    // Select reconstructed neutrino particles from a list of all reconstructed particles
    lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
    
    // Build mapping between PFP and Hits using PFParticle/SpacePoint/Hit maps.
    lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel,
                                                          recoParticlesToHits, recoHitsToParticles, daughterMode, true); // Use clusters to go from pfp to hits

    if (_verbose) {
        std::cout << "[RecoTrueHelper] [Configure] RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
        std::cout << "[RecoTrueHelper] [Configure] RecoParticles: " << recoParticleVector.size() << std::endl;
    }

    // Collect MCParticles and match True Particles to Hits
    lar_pandora::MCParticleVector     trueParticleVector;
    lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;
    lar_pandora::MCParticlesToHits    trueParticlesToHits;
    lar_pandora::HitsToMCParticles    trueHitsToParticles; 

    if (!e.isRealData()) {
        
        // Create map from MCTruth to MCParticle
        lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, trueParticleVector);

        // Create map from MCParticle to MCTruth
        lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, truthToParticles, particlesToTruth); 

        // Create Map track ID to MCParticle
        lar_pandora::MCParticleMap particleMap;

        // Loop over the map of MC Truth to MC Particles
        for (lar_pandora::MCTruthToMCParticles::const_iterator iter1 = truthToParticles.begin(), iterEnd1 = truthToParticles.end(); iter1 != iterEnd1; ++iter1) {
            
            // Get the MC Particle vector
            const lar_pandora::MCParticleVector &particleVector = iter1->second;
            
            // Loop over the MC Particle Vector and map to the track ID
            for (lar_pandora::MCParticleVector::const_iterator iter2 = particleVector.begin(), iterEnd2 = particleVector.end(); iter2 != iterEnd2; ++iter2) {
                const art::Ptr<simb::MCParticle> particle = *iter2;
                particleMap[particle->TrackId()] = particle;
            }
        }

        // anab::BackTrackerHitMatchingData: Store information from a reco-truth matching module (which is based on the BackTracker service)
        // Stores the cleanliness and completeness of a match
        // Cleanliness  = charge in reco object from true object / total charge in reco object
        // Completeness = charge in reco object from true object / total charge deposited by true object
        art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> mcps_from_hit(hit_h, e, _hit_mcp_producer);
        std::vector<art::Ptr<simb::MCParticle> >             mcp_v;   // MC Particle vector
        std::vector<anab::BackTrackerHitMatchingData const*> match_v; // Data reference

        // Loop over the hits, get the ass MCP, and then true to link
        for (auto hit : hitVector) {

            mcp_v.clear(); match_v.clear();
            mcps_from_hit.get(hit.key(), mcp_v, match_v);

            double max_energy = -1;
            int best_match_id = -1;

            // Loop over hit matched and get hit with largest energy of track ID
            for (size_t m = 0; m < match_v.size(); m++) {
                double this_energy = match_v[m]->energy;
                
                if (this_energy > max_energy) {
                    best_match_id = m;
                    max_energy = this_energy;
                }
            }

            if (best_match_id > -1) {

                try {

                    const art::Ptr<simb::MCParticle> thisParticle = mcp_v.at(best_match_id);
                    const art::Ptr<simb::MCParticle> primaryParticle(lar_pandora::LArPandoraHelper::GetFinalStateMCParticle(particleMap, thisParticle));

                    // If using the add daughter mode (daughters absorbed into parent), then return the primary particle, else return this particle
                    const art::Ptr<simb::MCParticle> selectedParticle((lar_pandora::LArPandoraHelper::kAddDaughters == daughterMode) ? primaryParticle : thisParticle); 

                    if ((lar_pandora::LArPandoraHelper::kIgnoreDaughters == daughterMode) && (selectedParticle != primaryParticle))
                        continue;

                    // Skip the ones which are not visible
                    if (!(lar_pandora::LArPandoraHelper::IsVisible(selectedParticle)))
                        continue;

                    trueHitsToParticles[hit] = selectedParticle; // Assign the mapping of MC Particle to the hit

                } 
                
                catch (cet::exception &e) {}
                // const auto mct = UBXSecHelper::TrackIDToMCTruth(e, "largeant", selectedParticle->TrackId());
                // if (mct->Origin() == simb::kBeamNeutrino && selectedParticle->PdgCode() == 13 && selectedParticle->Mother() == 0) {
                //     std::cout << "Muon from neutrino ass to hit " << hit->PeakTime() << ", "<< hit->WireID () << std::endl;
                // }
            }

        } // End loop over hits
    
    }

    // Now set the things we need for the future
    _trueHitsToParticles = trueHitsToParticles; // True Hits to True MC Particles
    _recoParticlesToHits = recoParticlesToHits; // Reco Particles to True Hits

    std::cout << "[RecoTrueHelper] [Configure] trueHitsToParticles size: " << trueHitsToParticles.size() << std::endl;
    //_configured = true;
}
//------------------------------------------------------------------------------------------------------------------------------------------
void recotruehelper::BuildTrueNeutrinoHitMaps(const lar_pandora::MCTruthToMCParticles &truthToParticles, const lar_pandora::MCParticlesToHits &trueParticlesToHits,
                                              lar_pandora::MCTruthToHits &trueNeutrinosToHits, lar_pandora::HitsToMCTruth &trueHitsToNeutrinos) const
{
    for (lar_pandora::MCTruthToMCParticles::const_iterator iter1 = truthToParticles.begin(), iterEnd1 = truthToParticles.end();
         iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<simb::MCTruth> trueNeutrino = iter1->first;
        const lar_pandora::MCParticleVector &trueParticleVector = iter1->second;

        for (lar_pandora::MCParticleVector::const_iterator iter2 = trueParticleVector.begin(), iterEnd2 = trueParticleVector.end(); iter2 != iterEnd2; ++iter2)
        {
            const lar_pandora::MCParticlesToHits::const_iterator iter3 = trueParticlesToHits.find(*iter2);
            if (trueParticlesToHits.end() == iter3)
                continue;

            const lar_pandora::HitVector &hitVector = iter3->second;

            for (lar_pandora::HitVector::const_iterator iter4 = hitVector.begin(), iterEnd4 = hitVector.end(); iter4 != iterEnd4; ++iter4)
            {
                const art::Ptr<recob::Hit> hit = *iter4;
                trueHitsToNeutrinos[hit] = trueNeutrino;
                trueNeutrinosToHits[trueNeutrino].push_back(hit);
            }
        }
    }
}


//--------------------------------------------------------------------------------------------------------------------------------------------
void recotruehelper::BuildRecoNeutrinoHitMaps(const lar_pandora::PFParticleMap &recoParticleMap, const lar_pandora::PFParticlesToHits &recoParticlesToHits,
                                              lar_pandora::PFParticlesToHits &recoNeutrinosToHits, lar_pandora::HitsToPFParticles &recoHitsToNeutrinos) const
{
    for (lar_pandora::PFParticleMap::const_iterator iter1 = recoParticleMap.begin(), iterEnd1 = recoParticleMap.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> recoParticle = iter1->second;
        const art::Ptr<recob::PFParticle> recoNeutrino = lar_pandora::LArPandoraHelper::GetParentPFParticle(recoParticleMap, recoParticle);

        if (!lar_pandora::LArPandoraHelper::IsNeutrino(recoNeutrino))
            continue;

        const lar_pandora::PFParticlesToHits::const_iterator iter2 = recoParticlesToHits.find(recoParticle);
        if (recoParticlesToHits.end() == iter2)
            continue;

        const lar_pandora::HitVector &hitVector = iter2->second;

        for (lar_pandora::HitVector::const_iterator iter3 = hitVector.begin(), iterEnd3 = hitVector.end(); iter3 != iterEnd3; ++iter3)
        {
            const art::Ptr<recob::Hit> hit = *iter3;
            recoHitsToNeutrinos[hit] = recoNeutrino;
            recoNeutrinosToHits[recoNeutrino].push_back(hit);
        }
    }
}

//___________________________________________________________________________________________________
void recotruehelper::GetRecoToTrueMatches(lar_pandora::MCParticlesToPFParticles &matchedParticles,
                                          lar_pandora::MCParticlesToHits &matchedParticleHits) {


    GetRecoToTrueMatches(_recoParticlesToHits,
                         _trueHitsToParticles,
                         matchedParticles,
                         matchedParticleHits);
}


//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoNeutrinosToHits, const lar_pandora::HitsToMCTruth &trueHitsToNeutrinos,
                                          lar_pandora::MCTruthToPFParticles &matchedNeutrinos, lar_pandora::MCTruthToHits &matchedNeutrinoHits) const
{
    PFParticleSet recoVeto; MCTruthSet trueVeto;

    this->GetRecoToTrueMatches(recoNeutrinosToHits, trueHitsToNeutrinos, matchedNeutrinos, matchedNeutrinoHits, recoVeto, trueVeto);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoNeutrinosToHits, const lar_pandora::HitsToMCTruth &trueHitsToNeutrinos,
                                          lar_pandora::MCTruthToPFParticles &matchedNeutrinos, lar_pandora::MCTruthToHits &matchedNeutrinoHits,
                                          PFParticleSet &vetoReco, MCTruthSet &vetoTrue) const
{
    bool foundMatches(false);

    // Loop over map of PFP to True Hits
    for (lar_pandora::PFParticlesToHits::const_iterator iter1 = recoNeutrinosToHits.begin(), iterEnd1 = recoNeutrinosToHits.end(); iter1 != iterEnd1; ++iter1) {
        
        const art::Ptr<recob::PFParticle> recoNeutrino = iter1->first;
        if (vetoReco.count(recoNeutrino) > 0)
            continue;

        const lar_pandora::HitVector &hitVector = iter1->second;

        lar_pandora::MCTruthToHits truthContributionMap;

        // Loop over the hit vector
        for (lar_pandora::HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2) {
            const art::Ptr<recob::Hit> hit = *iter2;

            lar_pandora::HitsToMCTruth::const_iterator iter3 = trueHitsToNeutrinos.find(hit);
            if (trueHitsToNeutrinos.end() == iter3)
                continue;

            const art::Ptr<simb::MCTruth> trueNeutrino = iter3->second;
            if (vetoTrue.count(trueNeutrino) > 0)
                continue;

            truthContributionMap[trueNeutrino].push_back(hit);
        }

        lar_pandora::MCTruthToHits::const_iterator mIter = truthContributionMap.end();

        for (lar_pandora::MCTruthToHits::const_iterator iter4 = truthContributionMap.begin(), iterEnd4 = truthContributionMap.end();
             iter4 != iterEnd4; ++iter4)
        {
            if ((truthContributionMap.end() == mIter) || (iter4->second.size() > mIter->second.size()))
            {
                mIter = iter4;
            }
        }

        if (truthContributionMap.end() != mIter)
        {
            const art::Ptr<simb::MCTruth> trueNeutrino = mIter->first;

            lar_pandora::MCTruthToHits::const_iterator iter5 = matchedNeutrinoHits.find(trueNeutrino);

            if ((matchedNeutrinoHits.end() == iter5) || (mIter->second.size() > iter5->second.size()))
            {
                matchedNeutrinos[trueNeutrino] = recoNeutrino;
                matchedNeutrinoHits[trueNeutrino] = mIter->second;
                foundMatches = true;
            }
        }
    }

    if (!foundMatches)
        return;

    for (lar_pandora::MCTruthToPFParticles::const_iterator pIter = matchedNeutrinos.begin(), pIterEnd = matchedNeutrinos.end();
         pIter != pIterEnd; ++pIter)
    {
        vetoTrue.insert(pIter->first);
        vetoReco.insert(pIter->second);
    }

    if (_recursiveMatching)
        this->GetRecoToTrueMatches(recoNeutrinosToHits, trueHitsToNeutrinos, matchedNeutrinos, matchedNeutrinoHits, vetoReco, vetoTrue);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                          lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits) const
{
    PFParticleSet recoVeto; MCParticleSet trueVeto;

    this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, recoVeto, trueVeto);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                          lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits,
                                          PFParticleSet &vetoReco, MCParticleSet &vetoTrue) const
{
    bool foundMatches(false);

    for (lar_pandora::PFParticlesToHits::const_iterator iter1 = recoParticlesToHits.begin(), iterEnd1 = recoParticlesToHits.end();
         iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> recoParticle = iter1->first;
        if (vetoReco.count(recoParticle) > 0)
            continue;

        const lar_pandora::HitVector &hitVector = iter1->second;

        lar_pandora::MCParticlesToHits truthContributionMap;

        for (lar_pandora::HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::Hit> hit = *iter2;

            lar_pandora::HitsToMCParticles::const_iterator iter3 = trueHitsToParticles.find(hit);
            if (trueHitsToParticles.end() == iter3)
                continue;

            const art::Ptr<simb::MCParticle> trueParticle = iter3->second;
            if (vetoTrue.count(trueParticle) > 0)
                continue;

            truthContributionMap[trueParticle].push_back(hit);
        }

        lar_pandora::MCParticlesToHits::const_iterator mIter = truthContributionMap.end();

        for (lar_pandora::MCParticlesToHits::const_iterator iter4 = truthContributionMap.begin(), iterEnd4 = truthContributionMap.end();
             iter4 != iterEnd4; ++iter4)
        {
            if ((truthContributionMap.end() == mIter) || (iter4->second.size() > mIter->second.size()))
            {
                mIter = iter4;
            }
        }

        if (truthContributionMap.end() != mIter)
        {
            const art::Ptr<simb::MCParticle> trueParticle = mIter->first;

            lar_pandora::MCParticlesToHits::const_iterator iter5 = matchedHits.find(trueParticle);

            if ((matchedHits.end() == iter5) || (mIter->second.size() > iter5->second.size()))
            {
                matchedParticles[trueParticle] = recoParticle;
                matchedHits[trueParticle] = mIter->second;
                foundMatches = true;
            }
        }
    }

    if (!foundMatches)
        return;

    for (lar_pandora::MCParticlesToPFParticles::const_iterator pIter = matchedParticles.begin(), pIterEnd = matchedParticles.end();
         pIter != pIterEnd; ++pIter)
    {
        vetoTrue.insert(pIter->first);
        vetoReco.insert(pIter->second);
    }

    if (_recursiveMatching)
        this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, vetoReco, vetoTrue);
}



//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::BuildRecoParticleMap(const lar_pandora::PFParticleVector &particleVector, lar_pandora::PFParticleMap &particleMap) const
{
    for (lar_pandora::PFParticleVector::const_iterator iter = particleVector.begin(), iterEnd = particleVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = *iter;
        particleMap[particle->Self()] = particle;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::BuildTrueParticleMap(const lar_pandora::MCParticleVector &particleVector, lar_pandora::MCParticleMap &particleMap) const
{
    for (lar_pandora::MCParticleVector::const_iterator iter = particleVector.begin(), iterEnd = particleVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<simb::MCParticle> particle = *iter;
        particleMap[particle->TrackId()] = particle;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int recotruehelper::CountHitsByType(const int view, const lar_pandora::HitVector &hitVector) const
{
    int nHits(0);

    for (lar_pandora::HitVector::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::Hit> hit = *iter;
        if (hit->View() == view)
            ++nHits;
    }

    return nHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void recotruehelper::GetStartAndEndPoints(const art::Ptr<simb::MCParticle> particle, int &startT, int &endT) const
{
    art::ServiceHandle<geo::Geometry> theGeometry;

    bool foundStartPosition(false);

    const int numTrajectoryPoints(static_cast<int>(particle->NumberTrajectoryPoints()));

    for (int nt = 0; nt < numTrajectoryPoints; ++nt)
    {
        try
        {
            double pos[3] = {particle->Vx(nt), particle->Vy(nt), particle->Vz(nt)};
            unsigned int which_tpc(std::numeric_limits<unsigned int>::max());
            unsigned int which_cstat(std::numeric_limits<unsigned int>::max());
            theGeometry->PositionToTPC(pos, which_tpc, which_cstat);

            // TODO: Apply fiducial cut due to readout window

            endT = nt;
            if (!foundStartPosition)
            {
                startT = endT;
                foundStartPosition = true;
            }
        }
        catch (cet::exception &e) {
            continue;
        }
    }

    if (!foundStartPosition)
        throw cet::exception("LArPandora");
}
//------------------------------------------------------------------------------------------------------------------------------------------
double recotruehelper::GetLength(const art::Ptr<simb::MCParticle> particle, const int startT, const int endT) const
{
    if (endT <= startT)
        return 0.0;

    double length(0.0);

    for (int nt = startT; nt < endT; ++nt)
    {
        const double dx(particle->Vx(nt+1) - particle->Vx(nt));
        const double dy(particle->Vy(nt+1) - particle->Vy(nt));
        const double dz(particle->Vz(nt+1) - particle->Vz(nt));
        length += sqrt(dx * dx + dy * dy + dz * dz);
    }

    return length;
}
//------------------------------------------------------------------------------------------------------------------------------------------
art::Ptr<simb::MCTruth> recotruehelper::TrackIDToMCTruth(art::Event const & e, std::string _geant_producer, int geant_track_id)
{
    lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;

    lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geant_producer, truthToParticles, particlesToTruth);

    for (auto iter : particlesToTruth) {
        if (iter.first->TrackId() == geant_track_id) {
            return iter.second;
        }
    }

    art::Ptr<simb::MCTruth> null_ptr;
    return null_ptr;
}
//------------------------------------------------------------------------------------------------------------------------------------------
}//end namespace lar_pandora