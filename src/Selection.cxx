#include "../include/Selection.h"

#include <iostream>

// Constructor
Selection::Selection(const Utility &util): _utility{ util } {
	
	std::cout << "Initialising Selection Class" << std::endl;

}

// ------------------------------------------------------------------------------
bool Selection::ApplyLanternSelection(EventContainer &_evt, Utility::FileTypeEnums type, Utility::RunPeriodEnums runPeriod) {
	// populate derived variables
    _evt.populateDerivedVariables(type, runPeriod);

	// determine event classification
    _evt.EventClassifier(type);

	// determine interaction classification
	_evt.InteractionClassifier(type);

	// determine particle classification
	_evt.ParticleClassifier(type);
	_evt.ParticleClassifierWC(type);
	_evt.ParticleClassifierLantern(type);

	// determine event weight
    _evt.calculateCVEventWeight(type, runPeriod); 

	// neutrino vertex found
	if (_evt.lantern_foundVertex == 0) return false;

	// vertex in fiducial volume
	if(!ApplyVertexFVCut(_evt.lantern_vtxX, _evt.lantern_vtxY, _evt.lantern_vtxZ)) return false;

	// vertex score
	if (_evt.lantern_vtxScore < 0.8) return false;

	// containment 
	if (_evt.lantern_vtxContainment != 2) return false;

	// cosmic rejection
	if (_evt.lantern_vtxFracHitsOnCosmic > 0.2) return false;

	// signal topology
	// count primary tracks and showers
	// identify longest track
	int n_primary_tracks = 0;
	int n_primary_showers = 0;
	float longest_track_length = 0.0;
	int longest_track_id = -1;
	for (int idx_track = 0; idx_track < _evt.lantern_nTracks; idx_track++) {
		if (_evt.lantern_trackIsSecondary[idx_track] == 0) { // primary
			n_primary_tracks++; // track
			// calculate length
			float track_length = std::sqrt(
				std::pow(_evt.lantern_trackEndPosX[idx_track] - _evt.lantern_trackStartPosX[idx_track], 2) +
				std::pow(_evt.lantern_trackEndPosY[idx_track] - _evt.lantern_trackStartPosY[idx_track], 2) +
				std::pow(_evt.lantern_trackEndPosZ[idx_track] - _evt.lantern_trackStartPosZ[idx_track], 2) );
			if (track_length > longest_track_length) {
				longest_track_length = track_length;
				longest_track_id = idx_track;
			}
		}
	}
	_evt.lantern_primaryTrackIndex = longest_track_id;
	_evt.lantern_primaryTrackLength = longest_track_length;
	_evt.lantern_primaryTrackVertexSeparation = std::sqrt(
		std::pow(_evt.lantern_trackStartPosX[longest_track_id] - _evt.lantern_vtxX, 2) +
		std::pow(_evt.lantern_trackStartPosY[longest_track_id] - _evt.lantern_vtxY, 2) +
		std::pow(_evt.lantern_trackStartPosZ[longest_track_id] - _evt.lantern_vtxZ, 2) );

	for (int idx_shower = 0; idx_shower < _evt.lantern_nShowers; idx_shower++) {
		if (_evt.lantern_showerIsSecondary[idx_shower] == 0) { // primary
			n_primary_showers++; // track
		}
	}

	// at least one primary track
	if (n_primary_tracks < 1) return false; 

	// zero primary showers
	if (n_primary_showers != 0) return false;

	// WC generic neutrino selection, cosmic rejection and containment
	if (!ApplyVertexFVCut(_evt.wc_reco_nuvtxX, _evt.wc_reco_nuvtxY, _evt.wc_reco_nuvtxZ)) return false; // WC vertex in fiducial volume
	if (_evt.wc_match_isFC == 0) return false; // WC fully contained
	if (_evt.wc_numu_cc_flag < 0) return false;		// WC generic neutrino selection
	
	// WC at least one primary particle identified 
	int n_particle_WC = 0;
	for (unsigned int wc_idx = 0; wc_idx < _evt.wc_reco_Ntrack; wc_idx++) {

		// check primary
		if (_evt.wc_reco_mother[wc_idx] != 0) continue; // require WC primary track

		// skip WC pseudo-particles
		if (_evt.wc_reco_pdg[wc_idx] == 22 || _evt.wc_reco_pdg[wc_idx] == 2112) continue; // skip pseudo-particles

		n_particle_WC++;
	}

	if (n_particle_WC < 1) return false;

	// classify primary tracks
	int n_lantern_larpid_mu = 0;
	int n_lantern_larpid_pi = 0;
	int n_lantern_larpid_pr = 0;
	int n_lantern_larpid_other = 0;
	int n_lantern_larpid_unclassified = 0;

	std::vector<int> lantern_pion_candidate_indices;

	for (int idx_track = 0; idx_track < _evt.lantern_nTracks; idx_track++) {
		if (_evt.lantern_trackIsSecondary[idx_track] == 0) { // primary

			// skip short tracks
			float track_length = std::sqrt(
				std::pow(_evt.lantern_trackEndPosX[idx_track] - _evt.lantern_trackStartPosX[idx_track], 2) +
				std::pow(_evt.lantern_trackEndPosY[idx_track] - _evt.lantern_trackStartPosY[idx_track], 2) +
				std::pow(_evt.lantern_trackEndPosZ[idx_track] - _evt.lantern_trackStartPosZ[idx_track], 2) );
			
			if (track_length < 5) continue;
		
			// check if classified
			if (_evt.lantern_trackClassified[idx_track] == 0) {
				n_lantern_larpid_unclassified++;
				continue;
			}

			// check PID
			if (std::abs(_evt.lantern_trackPID[idx_track]) == 13) {
				n_lantern_larpid_mu++;
				lantern_pion_candidate_indices.push_back(idx_track);
			}
			else if (std::abs(_evt.lantern_trackPID[idx_track]) == 211) {
				n_lantern_larpid_pi++;
				lantern_pion_candidate_indices.push_back(idx_track);
			}
			else if (std::abs(_evt.lantern_trackPID[idx_track]) == 2212) n_lantern_larpid_pr++;
			else {
				n_lantern_larpid_other++;
			}
		
		}
	}
	
	// require exactly one muon or pion track
	if (n_lantern_larpid_mu + n_lantern_larpid_pi != 1) return false;

	// no unclassified particles
	if (n_lantern_larpid_unclassified != 0) return false;

	// sanity check, should only be one candidate at this point
	if (lantern_pion_candidate_indices.size() != 1) {
		std::cout << "Warning: more than one pion candidate found! Logic issue." << std::endl;
		return false;
	}
	// save pion candidate index
	_evt.lantern_pion_candidate_index = lantern_pion_candidate_indices[0];

	// WC cross-check 
	// count primary particles as classified by WC
	int n_wc_reco_el = 0;
	int n_wc_reco_mu = 0;
	int n_wc_reco_pi = 0;
	int n_wc_reco_pr = 0;
	int n_wc_reco_other = 0;

	for (unsigned int wc_idx = 0; wc_idx < _evt.wc_reco_Ntrack; wc_idx++) {

		// skip WC pseudo-particles
		if (_evt.wc_reco_pdg[wc_idx] == 22 || _evt.wc_reco_pdg[wc_idx] == 2112) continue; // skip pseudo-particles

		// check primary
		if (_evt.wc_reco_mother[wc_idx] != 0) continue; // require WC primary particle

		// showers (electrons pdg)
		if (_evt.wc_reco_pdg[wc_idx] == 11) {
			// apply energy threshold
			if (_evt.wc_reco_startMomentum[wc_idx][3] > 0.01) n_wc_reco_el++;			
			continue;
		}

		// impose length requirement
		// calculate track length
		float wc_track_length = std::sqrt(
		std::pow(_evt.wc_reco_startXYZT[wc_idx][0] - _evt.wc_reco_endXYZT[wc_idx][0], 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][1] - _evt.wc_reco_endXYZT[wc_idx][1], 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][2] - _evt.wc_reco_endXYZT[wc_idx][2], 2) );

		if (wc_track_length < 5.0) continue; // 5 cm

		if (_evt.wc_reco_pdg[wc_idx] == 13) {
			n_wc_reco_mu++;
			continue;
		} else if (_evt.wc_reco_pdg[wc_idx] == 211) {
			n_wc_reco_pi++;
			continue;
		} else if (_evt.wc_reco_pdg[wc_idx] == 2212) {
			n_wc_reco_pr++;
			continue;
		}

		// if here, unclassified track
		n_wc_reco_other++;		

	}

	// Classical WC Selection
	// reject events with WC identified electron
	if (n_wc_reco_el > 0) return false; // reject events with primary shower

	// require at least one WC identified muon or pion
	if (n_wc_reco_mu + n_wc_reco_pi == 0) return false; // require at least one muon or pion

	// max one WC identified muon or pion
	if (n_wc_reco_mu + n_wc_reco_pi > 1) return false; // max one muon or pion

	// LArPID LLR cuts
	float mu_pi_llr_norm = std::tanh(0.5 * (_evt.lantern_trackPiScore[_evt.lantern_pion_candidate_index] - _evt.lantern_trackMuScore[_evt.lantern_pion_candidate_index]));
	float pr_pi_llr_norm = std::tanh(0.5 * (_evt.lantern_trackPiScore[_evt.lantern_pion_candidate_index] - _evt.lantern_trackPrScore[_evt.lantern_pion_candidate_index]));
	float pr_mu_llr_norm = std::tanh(0.5 * (_evt.lantern_trackMuScore[_evt.lantern_pion_candidate_index] - _evt.lantern_trackPrScore[_evt.lantern_pion_candidate_index]));

	_evt.sel_passInitialSelection_ = true;
	_evt.sel_LanternPID_llr_mu_pi_ = mu_pi_llr_norm;
	_evt.sel_LanternPID_llr_pr_pi_ = pr_pi_llr_norm;

	// require pion candidate to have mu/pi LLR < 0.4 (muon-like rejected)
	if (mu_pi_llr_norm < 0.6) return false;

	_evt.sel_passMuPiLLR_= true;

	// require pion candidate to have pr/pi LLR < 0.4 (proton-like rejected)
	if (pr_pi_llr_norm < 0.6) return false;
	
	// LArPID cuts	
	//if (_evt.lantern_trackPiScore[_evt.lantern_pion_candidate_index] < -0.5) return false;
	//if (_evt.lantern_trackMuScore[_evt.lantern_pion_candidate_index] > -1.5) return false;
	//if (_evt.lantern_trackPrScore[_evt.lantern_pion_candidate_index] > -2.0) return false;

	// second muon/pion veto
	int nSecondMuonPion = 0;
	int nUnclassified = 0;
	// loop through primary tracks
	for (int idx_track = 0; idx_track < _evt.lantern_nTracks; idx_track++) {
		if (_evt.lantern_trackIsSecondary[idx_track] == 0) { // primary

			// skip candidate track
			if (idx_track == _evt.lantern_pion_candidate_index) continue;	
			
			// skip short tracks
			float track_length = std::sqrt(
				std::pow(_evt.lantern_trackEndPosX[idx_track] - _evt.lantern_trackStartPosX[idx_track], 2) +
				std::pow(_evt.lantern_trackEndPosY[idx_track] - _evt.lantern_trackStartPosY[idx_track], 2) +
				std::pow(_evt.lantern_trackEndPosZ[idx_track] - _evt.lantern_trackStartPosZ[idx_track], 2) );
			
			if (track_length < 2) continue;

			// look for non proton-like	tracks

			if ( std::exp(_evt.lantern_trackPrScore[idx_track]) < 0.2) { // not proton-like
				nSecondMuonPion++;
				continue;
			}
		}
	}

	if (nSecondMuonPion != 0) return false;

	_evt.sel_NC1pi_ = true;

	return true;
}

// ------------------------------------------------------------------------------

bool Selection::ApplyWCSelection(EventContainer &_evt, Utility::FileTypeEnums type, Utility::RunPeriodEnums runPeriod) {

	// populate derived variables
    _evt.populateDerivedVariables(type, runPeriod);

	// determine event classification
    _evt.EventClassifier(type);

	// determine interaction classification
	_evt.InteractionClassifier(type);

	// determine particle classification
	_evt.ParticleClassifier(type);
	_evt.ParticleClassifierWC(type);

    // determine event weight
    _evt.calculateCVEventWeight(type, runPeriod); 

	// generic neutrino selection
	if (_evt.wc_numu_cc_flag < 0) return false;		// WC generic neutrino selection

	// vertex in fiducial volume
	if(!ApplyVertexFVCut(_evt.wc_reco_nuvtxX, _evt.wc_reco_nuvtxY, _evt.wc_reco_nuvtxZ)) return false;

	// require fully contained
	if (_evt.wc_match_isFC == 0) return false; // WC fully contained

	// WC
	// loop through particles to find longest track and count number of tracks/showers
	int wc_primaryTrackIndex_loop = -1;
	float wc_primaryTrackLength_loop = 0.0;	
	int wc_n_primary_particles = 0;	

	for (unsigned int wc_idx = 0; wc_idx < _evt.wc_reco_Ntrack; wc_idx++) {

		// check primary
		if (_evt.wc_reco_mother[wc_idx] != 0) continue; // require WC primary track

		// calculate track length
		float wc_track_length = std::sqrt(
		std::pow(_evt.wc_reco_startXYZT[wc_idx][0] - _evt.wc_reco_endXYZT[wc_idx][0], 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][1] - _evt.wc_reco_endXYZT[wc_idx][1], 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][2] - _evt.wc_reco_endXYZT[wc_idx][2], 2) );	
		
		if (wc_track_length > wc_primaryTrackLength_loop) {
			wc_primaryTrackIndex_loop = wc_idx;
			wc_primaryTrackLength_loop = wc_track_length;
		}

		wc_n_primary_particles++;		
	}

	_evt.wc_primaryTrackIndex = wc_primaryTrackIndex_loop;
	_evt.wc_primaryTrackLength = wc_primaryTrackLength_loop;

	_evt.wc_n_primary_particles = wc_n_primary_particles;
	
	// calculate primary track vertex separation
	_evt.wc_primaryTrackVertexSeparation = std::sqrt(
		std::pow(_evt.wc_reco_startXYZT[_evt.wc_primaryTrackIndex][0] - _evt.wc_reco_nuvtxX, 2) +
		std::pow(_evt.wc_reco_startXYZT[_evt.wc_primaryTrackIndex][1] - _evt.wc_reco_nuvtxY, 2) +
		std::pow(_evt.wc_reco_startXYZT[_evt.wc_primaryTrackIndex][2] - _evt.wc_reco_nuvtxZ, 2) );

	// at least one WC track
	if(_evt.wc_primaryTrackIndex == -1) return false;

	// require primary track within 5cm of vertex
	if(_evt.wc_primaryTrackVertexSeparation > 5.0) return false; // primary track within 5cm of vertex

	// track length
	if (_evt.wc_primaryTrackLength < 5.0) return false; // WC primary track length > 10 cm

	// count primary particles classified by LArPID
	int n_wc_reco_larpid_el = 0;
	int n_wc_reco_larpid_ph = 0;
	int n_wc_reco_larpid_mu = 0;
	int n_wc_reco_larpid_pi = 0;
	int n_wc_reco_larpid_pr = 0;
	int n_wc_reco_larpid_unclassified = 0;

	// pion candidate index
	std::vector<unsigned int> wc_pion_candidate_indices;
	
	for (unsigned int wc_idx = 0; wc_idx < _evt.wc_reco_Ntrack; wc_idx++) {

		// skip WC pseudo-particles
		if (_evt.wc_reco_pdg[wc_idx] == 22 || _evt.wc_reco_pdg[wc_idx] == 2112) continue; // skip pseudo-particles

		// check primary
		if (_evt.wc_reco_mother[wc_idx] != 0) continue; // require WC primary particle

		if (_evt.wc_reco_larpid_pdg[wc_idx] == 11) {
			n_wc_reco_larpid_el++;
			continue;
		} else if (_evt.wc_reco_larpid_pdg[wc_idx] == 22) {
			n_wc_reco_larpid_ph++;
			continue;
		}

		// dealing with tracks from this point
		// not interested in tracks disconected from vertex, potentially miss-IDed pion daughters
		float wc_TrackVertexSeparation = std::sqrt(
		std::pow(_evt.wc_reco_startXYZT[wc_idx][0] - _evt.wc_reco_nuvtxX, 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][1] - _evt.wc_reco_nuvtxY, 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][2] - _evt.wc_reco_nuvtxZ, 2) );

		if (wc_TrackVertexSeparation > 5.0) continue;
		
		// impose length requirement
		// calculate track length
		float wc_track_length = std::sqrt(
		std::pow(_evt.wc_reco_startXYZT[wc_idx][0] - _evt.wc_reco_endXYZT[wc_idx][0], 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][1] - _evt.wc_reco_endXYZT[wc_idx][1], 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][2] - _evt.wc_reco_endXYZT[wc_idx][2], 2) );

		if (wc_track_length < 5.0) continue; // 5 cm

		if (_evt.wc_reco_larpid_pdg[wc_idx] == 13) {
			wc_pion_candidate_indices.push_back(wc_idx);
			n_wc_reco_larpid_mu++;
			continue;
		} else if (_evt.wc_reco_larpid_pdg[wc_idx] == 211) {
			wc_pion_candidate_indices.push_back(wc_idx);
			n_wc_reco_larpid_pi++;
			continue;
		} else if (_evt.wc_reco_larpid_pdg[wc_idx] == 2212) {
			n_wc_reco_larpid_pr++;
			continue;
		}

		// if here, unclassified track
		n_wc_reco_larpid_unclassified++;		

	}

	// count primary particles as classified by WC
	int n_wc_reco_el = 0;
	int n_wc_reco_mu = 0;
	int n_wc_reco_pi = 0;
	int n_wc_reco_pr = 0;
	int n_wc_reco_other = 0;

	// pion candidate index
	std::vector<unsigned int> wc_pion_candidate_indices_classical;
	
	for (unsigned int wc_idx = 0; wc_idx < _evt.wc_reco_Ntrack; wc_idx++) {

		// skip WC pseudo-particles
		if (_evt.wc_reco_pdg[wc_idx] == 22 || _evt.wc_reco_pdg[wc_idx] == 2112) continue; // skip pseudo-particles

		// check primary
		if (_evt.wc_reco_mother[wc_idx] != 0) continue; // require WC primary particle

		// showers (electrons pdg)
		if (_evt.wc_reco_pdg[wc_idx] == 11) {
			// apply energy threshold
			if (_evt.wc_reco_startMomentum[wc_idx][3] > 0.01) n_wc_reco_el++;			
			continue;
		}

		// dealing with tracks from this point
		// not interested in tracks disconected from vertex, potentially miss-IDed pion daughters
		float wc_TrackVertexSeparation = std::sqrt(
		std::pow(_evt.wc_reco_startXYZT[wc_idx][0] - _evt.wc_reco_nuvtxX, 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][1] - _evt.wc_reco_nuvtxY, 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][2] - _evt.wc_reco_nuvtxZ, 2) );

		if (wc_TrackVertexSeparation > 5.0) continue;
		
		// impose length requirement
		// calculate track length
		float wc_track_length = std::sqrt(
		std::pow(_evt.wc_reco_startXYZT[wc_idx][0] - _evt.wc_reco_endXYZT[wc_idx][0], 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][1] - _evt.wc_reco_endXYZT[wc_idx][1], 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][2] - _evt.wc_reco_endXYZT[wc_idx][2], 2) );

		if (wc_track_length < 5.0) continue; // 5 cm

		if (_evt.wc_reco_pdg[wc_idx] == 13) {
			wc_pion_candidate_indices_classical.push_back(wc_idx);
			n_wc_reco_mu++;
			continue;
		} else if (_evt.wc_reco_pdg[wc_idx] == 211) {
			wc_pion_candidate_indices_classical.push_back(wc_idx);
			n_wc_reco_pi++;
			continue;
		} else if (_evt.wc_reco_pdg[wc_idx] == 2212) {
			n_wc_reco_pr++;
			continue;
		}

		// if here, unclassified track
		n_wc_reco_other++;		

	}

	// store number of classified particles
	_evt.n_wc_reco_larpid_unclassified = n_wc_reco_larpid_unclassified;	
	
	// Classical WC Selection
	// reject events with WC identified electron
	if (n_wc_reco_el > 0) return false; // reject events with primary shower

	// require at least one WC identified muon or pion
	if (n_wc_reco_mu + n_wc_reco_pi == 0) return false; // require at least one muon or pion

	// max one WC identified muon or pion
	if (n_wc_reco_mu + n_wc_reco_pi > 1) return false; // max one muon or pion

	// sanity check
	if (wc_pion_candidate_indices_classical.size() != 1) {
		std::cout << "Error: Expected exactly one WC pion candidate, found " << wc_pion_candidate_indices_classical.size() << std::endl;
		return false;
	}
	// store candidate index
	_evt.wc_pion_candidate_index = wc_pion_candidate_indices_classical[0];

	// LArPID selection
	// check candidate is LARPID classified
	if (_evt.wc_reco_larpid_classified[_evt.wc_pion_candidate_index] == 0) return false;

	// reject events with unclassified non-stub tracks
	if (n_wc_reco_larpid_unclassified > 0) return false;

	// reject events with LArPID identified electron or photon
	if (n_wc_reco_larpid_el > 0 || n_wc_reco_larpid_ph > 0) return false; // reject events with electron or photon

	// construct LLR variables for pion candidate
	// LArPID LLR cuts
	float mu_pi_llr_norm = std::tanh(0.5 * (_evt.wc_reco_larpid_pidScore_pi[_evt.wc_pion_candidate_index] - _evt.wc_reco_larpid_pidScore_mu[_evt.wc_pion_candidate_index]));
	float pr_pi_llr_norm = std::tanh(0.5 * (_evt.wc_reco_larpid_pidScore_pi[_evt.wc_pion_candidate_index] - _evt.wc_reco_larpid_pidScore_pr[_evt.wc_pion_candidate_index]));
	float pr_mu_llr_norm = std::tanh(0.5 * (_evt.wc_reco_larpid_pidScore_mu[_evt.wc_pion_candidate_index] - _evt.wc_reco_larpid_pidScore_pr[_evt.wc_pion_candidate_index]));

	_evt.sel_passInitialSelection_ = true;
	_evt.sel_LanternPID_llr_mu_pi_ = mu_pi_llr_norm;
	_evt.sel_LanternPID_llr_pr_pi_ = pr_pi_llr_norm;

	// require pion candidate to have mu/pi LLR < 0.6 (muon-like rejected)
	if (mu_pi_llr_norm < 0.825) return false;

	_evt.sel_passMuPiLLR_= true;

	// require pion candidate to have pr/pi LLR < 0.5 (proton-like rejected)
	if (pr_pi_llr_norm < 0.6) return false;

	// veto second MIP-like track
	int n_mip_like_tracks = 0;
	for (unsigned int wc_idx = 0; wc_idx < _evt.wc_reco_Ntrack; wc_idx++) {

		// skip WC pseudo-particles
		if (_evt.wc_reco_pdg[wc_idx] == 22 || _evt.wc_reco_pdg[wc_idx] == 2112) continue; // skip pseudo-particles

		// check primary
		if (_evt.wc_reco_mother[wc_idx] != 0) continue; // require WC primary particle

		// skip candidate track
		if (wc_idx == _evt.wc_pion_candidate_index) continue;

		// dealing with tracks from this point
		// not interested in tracks disconected from vertex, potentially miss-IDed pion daughters
		float wc_TrackVertexSeparation = std::sqrt(
		std::pow(_evt.wc_reco_startXYZT[wc_idx][0] - _evt.wc_reco_nuvtxX, 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][1] - _evt.wc_reco_nuvtxY, 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][2] - _evt.wc_reco_nuvtxZ, 2) );

		if (wc_TrackVertexSeparation > 5.0) continue;
		
		// impose length requirement
		// calculate track length
		float wc_track_length = std::sqrt(
		std::pow(_evt.wc_reco_startXYZT[wc_idx][0] - _evt.wc_reco_endXYZT[wc_idx][0], 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][1] - _evt.wc_reco_endXYZT[wc_idx][1], 2) +
		std::pow(_evt.wc_reco_startXYZT[wc_idx][2] - _evt.wc_reco_endXYZT[wc_idx][2], 2) );

		if (wc_track_length < 3.0) continue; // 5 cm

		// check classified by LArPID
		if (_evt.wc_reco_larpid_classified[wc_idx] == 0) continue; // only consider classified tracks

		// count number muon-like tracks
		//if ( std::exp(_evt.wc_reco_larpid_pidScore_mu[wc_idx]) > 0.5 || std::exp(_evt.wc_reco_larpid_pidScore_pi[wc_idx]) > 0.5) {
		if ( std::exp(_evt.wc_reco_larpid_pidScore_pr[wc_idx]) < 0.25) { // not proton-like
			n_mip_like_tracks++;
			continue;
		}		
	}

	if (n_mip_like_tracks > 0) return false; // veto second MIP-like track

	// calculate pion momentum for candidate track -- range
	float wc_candidate_track_length = std::sqrt(
		std::pow(_evt.wc_reco_startXYZT[_evt.wc_pion_candidate_index][0] - _evt.wc_reco_endXYZT[_evt.wc_pion_candidate_index][0], 2) +
		std::pow(_evt.wc_reco_startXYZT[_evt.wc_pion_candidate_index][1] - _evt.wc_reco_endXYZT[_evt.wc_pion_candidate_index][1], 2) +
		std::pow(_evt.wc_reco_startXYZT[_evt.wc_pion_candidate_index][2] - _evt.wc_reco_endXYZT[_evt.wc_pion_candidate_index][2], 2) );

	float pionMomentumRange = _evt.CalculatePionMomentumRange(wc_candidate_track_length);

	// calculate pion momentum for candidate track -- hypfit
	float pionMomentumHypfit = _evt.CalculatePionMomentumHypfit(_evt.wc_pion_candidate_index);

	_evt.sel_NC1pi_ = true;

	// event passes
    return true;
}

// ------------------------------------------------------------------------------

bool Selection::ApplyGoodShowerSelection(const EventContainer &_evt) {

	if(_evt.shr_pfpgeneration != 2) return false; 	// require pfp generation 2

	if(!ApplyHitsOnAllPlanesCut(_evt.shr_planehits_U, _evt.shr_planehits_V, _evt.shr_planehits_Y)) return false; 	// require hits on all planes
	
	if(!ApplyShowerEnergyCut(_evt.shr_energy_cali)) return false;					// shower energy

	if(!ApplyShowerScoreCut(_evt.shr_score)) return false;							// shower score

	if(!ApplyShowerSubclustersCut(_evt.shrsubclusters)) return false;				// number of shower subclusters

	if(!ApplyShowerHitRatioCut(_evt.hits_ratio)) return false;						// shower hits ratio

	return true;
}

// ------------------------------------------------------------------------------

bool Selection::ApplyGoodTrackSelection(EventContainer &_evt) {

	bool primaryTrackPasses = false;
	bool secondaryTrackPasses = false;
	bool tertiaryTrackPasses = false;

	// primary track
	if (!_evt.hasSpuriousLeadingTrack &&
		_evt.trk_pfpgeneration == 2 &&
		ApplyHitsOnAllPlanesCut(_evt.trk_planehits_U, _evt.trk_planehits_V, _evt.trk_planehits_Y) &&
		ApplyTrackLengthCut(_evt.trk_len) &&
		ApplyTrackShowerOpeningAngleCut(_evt.trk_shr_opening_angle) && 
		ApplyTrackScoreCut(_evt.trk_score) &&
		ApplyTrackContainmentCut(_evt.trk_sce_end_x, _evt.trk_sce_end_y, _evt.trk_sce_end_z) &&
	    ApplyTrackVertexDistanceCut(_evt.trk_distance)
		) primaryTrackPasses = true;  	// track length, score, distance, containment

	// secondary track	
	if (_evt.trk2_pfpgeneration == 2 &&
	    ApplyHitsOnAllPlanesCut(_evt.trk2_planehits_U, _evt.trk2_planehits_V, _evt.trk2_planehits_Y) &&
		ApplyTrackLengthCut(_evt.trk2_len) &&
		ApplyTrackShowerOpeningAngleCut(_evt.trk2_shr_opening_angle) && 
		ApplyTrackScoreCut(_evt.trk2_score) &&
		ApplyTrackContainmentCut(_evt.trk2_sce_end_x, _evt.trk2_sce_end_y, _evt.trk2_sce_end_z) &&
		ApplyTrackVertexDistanceCut(_evt.trk2_distance)
		) secondaryTrackPasses = true;		// track length, score, distance, containment

	// tertiary track
	if (_evt.trk3_pfpgeneration == 2 &&
	    ApplyHitsOnAllPlanesCut(_evt.trk3_planehits_U, _evt.trk3_planehits_V, _evt.trk3_planehits_Y) &&
		ApplyTrackLengthCut(_evt.trk3_len) &&
		ApplyTrackShowerOpeningAngleCut(_evt.trk3_shr_opening_angle) && 
		ApplyTrackScoreCut(_evt.trk3_score) &&
		ApplyTrackContainmentCut(_evt.trk3_sce_end_x, _evt.trk3_sce_end_y, _evt.trk3_sce_end_z) &&
		ApplyTrackVertexDistanceCut(_evt.trk3_distance)
		) tertiaryTrackPasses = true;		// track length, score, distance, containment
	
	// update selection status variables in event for use later
	_evt.primaryTrackValid = primaryTrackPasses;
	_evt.secondaryTrackValid = secondaryTrackPasses;
	_evt.tertiaryTrackValid = tertiaryTrackPasses;

	if (primaryTrackPasses || secondaryTrackPasses || tertiaryTrackPasses) return true;
	else return false;
}

// ------------------------------------------------------------------------------
bool Selection::ApplyReconstructionFailureChecks(const EventContainer &_evt) { 

	//if (!ApplyPrimaryTrackShowerOpeningAngleCut(_evt.n_tracks, _evt.trk_shr_opening_angle)) return false;		// primary track-shower opening angle

	if (!ApplyShrTrackFitCut(_evt.trkfit)) return false;				    // shower track fit hits fraction
    
    if (!ApplyShrTrackLengthCut(_evt.shr_trk_len)) return false;			// shower track fit length

    if (!ApplyShowerCylFractionCut(_evt.CylFrac2h_1cm)) return false;		// leading shower 1cm cylinder energy fraction	

	return true;
} 

// ------------------------------------------------------------------------------

bool Selection::ApplyCosmicRejection(const EventContainer &_evt, Utility::RunPeriodEnums runPeriod) {

	if(!ApplyTopologicalScoreCut(_evt.topological_score)) return false;				// topological score

	if(!ApplyCosmicImpactParameterCut(_evt.CosmicIPAll3D)) return false;			// cosmic impact parameter

	//if (runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
	//	if(!ApplyCRTVetoCut(_evt.crtveto)) return false;							// crt veto
	//}

	return true;
}

// ------------------------------------------------------------------------------

bool Selection::ApplyNeutralPionRejection(const EventContainer &_evt){

	if(!ApplyNumberShowersCut(_evt.n_showers_contained)) return false; 	// number of showers

	if(!ApplyNeutralPionRejectionCut(_evt.shr_trkfit_gap10_dedx_max, _evt.shr_distance)) return false;		// neutral pion rejection: shower dE/dx and vertex separation
	
	if(!ApplyMoliereAverageCut(_evt.shrmoliereavg)) return false;					// shower moliere average
	
	if(!ApplyLeadingShowerEnergyFractionCut(_evt.shr_energyFraction)) return false;	// leading shower energy fraction

	if(!ApplySecondShowerClusterCut(_evt.secondshower_Y_nhit, _evt.secondshower_Y_vtxdist, _evt.secondshower_Y_anglediff)) return false; // second shower cluster tagger

	return true;
}

bool Selection::ApplyLooseNeutralPionRejection(const EventContainer &_evt){

	//if(!ApplyLooseLeadingShowerEnergyFractionCut(_evt.shr_energyFraction)) return false;	// leading shower energy fraction (loose)

	if(!ApplyLooseShowerDistanceCut(_evt.shr_distance)) return false; 			// shower distance from vertex (loose)

	if(!ApplyLooseMoliereAverageCut(_evt.shrmoliereavg)) return false;			// shower moliere average (loose)
	
	return true;
}

// ------------------------------------------------------------------------------

bool Selection::ApplyLooseProtonRejection(EventContainer &_evt) {

	
	bool primaryTrackPasses = false;
	bool secondaryTrackPasses = false;
	bool tertiaryTrackPasses = false;

	// primary track
	if (_evt.primaryTrackValid) {
		if (ApplyLLRPIDScoreCut(_evt.trk_llr_pid_score)) primaryTrackPasses = true;   
	}

	// secondary track
	if (_evt.secondaryTrackValid) {
		if (ApplyLLRPIDScoreCut(_evt.trk2_llr_pid_score)) secondaryTrackPasses = true; 	
	}

	// tertiary track
	if (_evt.tertiaryTrackValid) {
		if (ApplyLLRPIDScoreCut(_evt.trk3_llr_pid_score)) tertiaryTrackPasses = true;
	}

	// update selection status variables in event for use later
	_evt.primaryTrackPionlikeLoose = primaryTrackPasses;
	_evt.secondaryTrackPionlikeLoose = secondaryTrackPasses;
	_evt.tertiaryTrackPionlikeLoose = tertiaryTrackPasses;

	// require at least one track passes loose cuts
	if (primaryTrackPasses || secondaryTrackPasses || tertiaryTrackPasses) return true;
	else return false;
	
}

bool Selection::ApplyProtonRejection(EventContainer &_evt) {

	
	bool primaryTrackPasses = false;
	bool secondaryTrackPasses = false;
	bool tertiaryTrackPasses = false;

	// primary track
	if (_evt.primaryTrackValid) {
		if (ApplyLLRPIDScoreCut(_evt.trk_llr_pid_score) &&
			ApplyTrackTrunkdEdxCut(_evt.trk_dEdx_trunk_max) &&
			ApplyTrackBraggPeakScoreCut(_evt.trk_bragg_mip_max, _evt.trk_bragg_pion_max)
			) primaryTrackPasses = true;   
	}

	// secondary track
	if (_evt.secondaryTrackValid) {
		if (ApplyLLRPIDScoreCut(_evt.trk2_llr_pid_score) &&
			ApplyTrackTrunkdEdxCut(_evt.trk2_dEdx_trunk_max) &&
			ApplyTrackBraggPeakScoreCut(_evt.trk2_bragg_mip_max, _evt.trk2_bragg_pion_max)
			) secondaryTrackPasses = true; 	
	}

	// tertiary track
	if (_evt.tertiaryTrackValid) {
		if (ApplyLLRPIDScoreCut(_evt.trk3_llr_pid_score) &&
			ApplyTrackTrunkdEdxCut(_evt.trk3_dEdx_trunk_max) &&
			ApplyTrackBraggPeakScoreCut(_evt.trk3_bragg_mip_max, _evt.trk3_bragg_pion_max)
			) tertiaryTrackPasses = true;
	}

	// update selection status variables in event for use later
	_evt.primaryTrackPionlike = primaryTrackPasses;
	_evt.secondaryTrackPionlike = secondaryTrackPasses;
	_evt.tertiaryTrackPionlike = tertiaryTrackPasses;

	// count number of tracks that pass
	int nPass = 0;
	if (primaryTrackPasses) nPass++;
	if (secondaryTrackPasses) nPass++;
	if (tertiaryTrackPasses) nPass++;	

	// require 1 and only 1 charged pion candidate
	//if (primaryTrackPasses || secondaryTrackPasses || tertiaryTrackPasses) return true;
	if (nPass == 1) return true;
	else return false;
	
}

// ------------------------------------------------------------------------------
// Selection cuts
// Software Trigger [MC only]
bool Selection::ApplySWTriggerCut(int swtrig){
	return swtrig;
}

// Common Optical Filter PE [MC only, not used]
bool Selection::ApplyCommonOpticalFilterPECut(float opfilter_pe_beam){
	if (opfilter_pe_beam > 0) return true;
    else return false;
}

// Common Optical Filter Michel Veto [MC only, not used]
bool Selection::ApplyCommonOpticalFilterMichelCut(float opfilter_pe_veto){
	if (opfilter_pe_veto < 20) return true;
    else return false;
}

// Slice ID
bool Selection::ApplySliceIDCut(int nslice){
	if (nslice == 1) return true;
	else return false;
}

// CC Pi+: at least two tracks, and zero showers
bool Selection::ApplySignalCanidateCut(int n_showers, int n_tracks){
	if (n_tracks >= 1 && n_showers == 0) return true;
	else return false; 
}

// Vertex within fiducial volume
bool Selection::ApplyVertexFVCut(float reco_nu_vtx_sce_x, float reco_nu_vtx_sce_y, float reco_nu_vtx_sce_z){
	bool isInFV = _utility.inFV(reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z);
	if (isInFV) return true;
	else return false;
}

// Contained fraction of hits
bool Selection::ApplyContainedFractionCut(float contained_fraction){
	if(contained_fraction >= 0.7) return true;
	else return false;
}

// Fraction of hits associated with tracks and showers
bool Selection::ApplyAssociatedHitsFractionCut(float associated_hits_fraction) {
	if(associated_hits_fraction >= 0.5) return true;
	else return false;
}

// PFP hits on all planes
bool Selection::ApplyHitsOnAllPlanesCut(int pfp_planehits_U, int pfp_planehits_V, int pfp_planehits_Y) {
	if (pfp_planehits_U > 0 && pfp_planehits_V > 0 && pfp_planehits_Y > 0) return true;
	else return false;
}

// Shower energy
bool Selection::ApplyShowerEnergyCut(float shr_energy) {
	if (shr_energy >= 0.03) return true;
	else return false;
}	

// Shower score
bool Selection::ApplyShowerScoreCut(float shr_score){
	if (shr_score <= 0.2 && shr_score >= -0.1) return true; // >= -0.1  to check is defined, rather than default value
	else return false;
}

// Shower hit ratio
bool Selection::ApplyShowerHitRatioCut(float hits_ratio){
	if (hits_ratio >= 0.5) return true; // CCincl -- to 30%? 
	else return false;
}

// Pandora topological score
bool Selection::ApplyTopologicalScoreCut(float topological_score){
	if(topological_score > 0.15) return true; // 0.4
	else return false;
}

// Cosmic impact parameter
bool Selection::ApplyCosmicImpactParameterCut(float CosmicIPAll3D) {
	if (CosmicIPAll3D > 10) return true;
	else return false;
}

// CRT Veto
bool Selection::ApplyCRTVetoCut(float crtveto) {
	if (!crtveto) return true;
	else return false;
}

// Number of showers
bool Selection::ApplyNumberShowersCut(int n_showers_contained) {
	if(n_showers_contained <= 2) return true;
	else return false;
}

// Leading shower Moliere average
bool Selection::ApplyMoliereAverageCut(float shrmoliereavg){
	if(shrmoliereavg < 8) return true;
	else return false;
}

// Leading shower Moliere average (loose)
bool Selection::ApplyLooseMoliereAverageCut(float shrmoliereavg){
	if(shrmoliereavg < 15) return true;
	else return false;
}

// Leading shower subclusters 
bool Selection::ApplyShowerSubclustersCut(unsigned int shrsubclusters) {
	if(shrsubclusters >= 5) return true;
	else return false;
}

// Leading shower energy fraction
bool Selection::ApplyLeadingShowerEnergyFractionCut(float shr_energyFraction){
	if(shr_energyFraction >= 0.8) return true;
	else return false;
}

// Loose leading shower energy fraction
bool Selection::ApplyLooseLeadingShowerEnergyFractionCut(float shr_energyFraction){
	if(shr_energyFraction >= 0.7) return true;
	else return false;
}

// Leading shower fraction of energy in 1cm cylinder from shower center [see PeLEE]
bool Selection::ApplyShowerCylFractionCut(float CylFrac2h_1cm) {
	if(CylFrac2h_1cm < 0.95)  return true;	// default value -1
	else return false;	
}

// Neutral pion rejection: shower dE/dx and vertex distance
// apply 2D cut
// note: altered to match Krishan's
bool Selection::ApplyNeutralPionRejectionCut(float dEdxMax, float shr_distance) {
	
    if (dEdxMax <= 0) {
    	return false;
    }
    else if (dEdxMax > 0 && dEdxMax < 1.75){
        if (shr_distance >= 3.0 ) return false;
        else return true;
    }
    else if (dEdxMax >= 1.75 && dEdxMax < 2.5){ 
        if (shr_distance >= 12.0 ) return false;
        else return true;
    }
    else if (dEdxMax >= 2.5 && dEdxMax < 3.5){
        if (shr_distance >= 3.0 ) return false;
        else return true;
    }
    else if (dEdxMax >= 3.5 && dEdxMax < 4.7){
        return false;
    }
    else if (dEdxMax >= 4.7){
        if (shr_distance >= 3.0 ) return false;
        else return true;
    }
    else{
        std::cout << "Error: [Selection.h] Uncaught dEdx values." << std::endl;
        return false;
    }
}

// Loose shower distance from vertex cut
bool Selection::ApplyLooseShowerDistanceCut(float shr_distance) {
	if (shr_distance < 10) return true;
	else return false;
}

// Neutral pion rejection: second shower tagger
// Only using Y plane - increased noise etc, causing issues with other planes
bool Selection::ApplySecondShowerClusterCut(int secondshower_Y_nhit, float secondshower_Y_vtxdist, float secondshower_Y_anglediff) {
	if (secondshower_Y_nhit > 25 && secondshower_Y_vtxdist < 100 && secondshower_Y_anglediff > 10 && secondshower_Y_anglediff != 9999) return false;
	return true;
}

// Neutral pion rejection BDT
bool Selection::ApplyElectronPhotonBDTCutFHC(float bdtscore_electronPhoton) {
	if (bdtscore_electronPhoton > 0.525) return true;
	
	else return false;
}
bool Selection::ApplyElectronPhotonBDTCutRHC(float bdtscore_electronPhoton) {
	if (bdtscore_electronPhoton > 0.575) return true;
	
	else return false;
}

// Track Length
bool Selection::ApplyTrackLengthCut(float trk_len) {
	if (trk_len >= 5 && trk_len <= 200) return true; // 5
	else return false;
}

// Track Energy
bool Selection::ApplyTrackEnergyCut(float trk_energy){
	if (trk_energy >= 0.04) return true; // 40 MeV
	else return false;
}

// Track Vertex Distance
bool Selection::ApplyTrackVertexDistanceCut(float trk_distance) {
	if (trk_distance <= 4) return true;
	else return false;
}

// Track Score
bool Selection::ApplyTrackScoreCut(float trk_score) {
	if(trk_score >= 0.5 && trk_score <= 1.1) return true;
	else return false;
}

// Track Containment
bool Selection::ApplyTrackContainmentCut(float trk_sce_end_x, float trk_sce_end_y, float trk_sce_end_z) {
	bool isExiting = _utility.isExiting(trk_sce_end_x, trk_sce_end_y, trk_sce_end_z);
	return isExiting;
}

// Track Shower Opening Angle
bool Selection::ApplyTrackShowerOpeningAngleCut(float tksh_angle) {
	//if (tksh_angle >= 10 && tksh_angle <= 170) return true;
	if (tksh_angle <= 170) return true;
	else return false;
}

// Primary Track Shower Opening Angle (remove back-to-back tracks, reconstruction failure)
bool Selection::ApplyPrimaryTrackShowerOpeningAngleCut(float ntracks, float tksh_angle) {
	if (ntracks > 0 && tksh_angle <= 170) return true;
	else return false;
}

// LLR PID Score
bool Selection::ApplyLLRPIDScoreCut(float trk_llr_pid_score) {	
	if(trk_llr_pid_score > 0 && trk_llr_pid_score < 1.1) return true;
	else return false;
}

// Proton rejection: BDT
bool Selection::ApplyProtonRejectionBDTCutFHC(float bdtscore_pionProton) {
	if (bdtscore_pionProton > 0.350) return true; // Old Flux 0.350
	//if (bdtscore_pionProton > 0.250) return true; // Old Flux 0.350
	else return false;
}
bool Selection::ApplyProtonRejectionBDTCutRHC(float bdtscore_pionProton) {
	if (bdtscore_pionProton > 0.300) return true; // Old Flux 0.325
	//if (bdtscore_pionProton > 0.250) return true;
	else return false;
}

// Track Trunk dE/dx
bool Selection::ApplyTrackTrunkdEdxCut(float trk_dEdx_trunk) {
	if(trk_dEdx_trunk <= 3.8) return true;
	else return false;
}

// Track Bragg Peak Score
bool Selection::ApplyTrackBraggPeakScoreCut(float trk_bragg_mip, float trk_bragg_pion) {
	if (trk_bragg_mip < 0.1 && trk_bragg_pion < 0.5) return false;
	else return true;
}

// Shower Track Fit Score
bool Selection::ApplyShrTrackFitCut(float trkfit) {
	if (trkfit < 0.5) return true;	// default value = 1 when unfilled (float)
	else return false;
}

// Shower Track Fit Length
bool Selection::ApplyShrTrackLengthCut(float shr_trk_len) {
	if (shr_trk_len <= 250) return true;
	else return false;
}

int Selection::CountProtons(EventContainer &_evt, Utility::RunPeriodEnums runPeriod) {

	int numberProtons = 0;

	// consider all tracks in the event
	for (int i = 0; i < _evt.trk_llr_pid_score_v->size(); i++) {

		// check entry does not correspond to the pion
		if ( i == _evt.trk_id-1 && _evt.primaryTrackPionlike ) continue;
		if ( i == _evt.trk2_id-1 && _evt.secondaryTrackPionlike ) continue;
		if ( i == _evt.trk3_id-1 && _evt.tertiaryTrackPionlike ) continue;

		// look at tracks only
		if (_evt.trk_score_v->at(i) <= 0.5) continue;

		// cut on vertex distance
		if (_evt.trk_distance_v->at(i) > 10) continue;

		// count proton candiates
		numberProtons++;
		
	}
	
	// fudge for plotting purposes, to put all N proton events in single bin
	if (numberProtons > 2) numberProtons = 2;

	return numberProtons;
}

void Selection::setSelectedPionInformation(EventContainer &_evt) {

    // loose
    if (_evt.primaryTrackPionlikeLoose) {
    	_evt.trk_bragg_mip_pion_loose = _evt.trk_bragg_mip_max;
     	_evt.trk_daughters_pion_loose = _evt.trk_daughters;
    	_evt.trk_dEdx_trunk_pion_loose = _evt.trk_dEdx_trunk_max;
    	_evt.trk_bragg_pion_pion_loose = _evt.trk_bragg_pion_max;
    	_evt.trk_llr_pid_score_pion_loose = _evt.trk_llr_pid_score;
    	_evt.trk_score_pion_loose = _evt.trk_score;
    	_evt.trk_end_spacepoints_pion_loose = _evt.trk_end_spacepoints;
    }
    else if (_evt.secondaryTrackPionlikeLoose) {
    	_evt.trk_bragg_mip_pion_loose = _evt.trk2_bragg_mip_max;
     	_evt.trk_daughters_pion_loose = _evt.trk2_daughters;
    	_evt.trk_dEdx_trunk_pion_loose = _evt.trk2_dEdx_trunk_max;
    	_evt.trk_bragg_pion_pion_loose = _evt.trk2_bragg_pion_max;
    	_evt.trk_llr_pid_score_pion_loose = _evt.trk2_llr_pid_score;
    	_evt.trk_score_pion_loose = _evt.trk2_score;
    	_evt.trk_end_spacepoints_pion_loose = _evt.trk2_end_spacepoints;
    }
    else {
    	_evt.trk_bragg_mip_pion_loose = _evt.trk3_bragg_mip_max;
     	_evt.trk_daughters_pion_loose = _evt.trk3_daughters;
    	_evt.trk_dEdx_trunk_pion_loose = _evt.trk3_dEdx_trunk_max;
    	_evt.trk_bragg_pion_pion_loose = _evt.trk3_bragg_pion_max;
    	_evt.trk_llr_pid_score_pion_loose = _evt.trk3_llr_pid_score;
    	_evt.trk_score_pion_loose = _evt.trk3_score;
    	_evt.trk_end_spacepoints_pion_loose = _evt.trk3_end_spacepoints;
    }

    // full
    if (_evt.primaryTrackPionlike) {
    	_evt.trk_bragg_mip_pion = _evt.trk_bragg_mip_max;
     	_evt.trk_daughters_pion = _evt.trk_daughters;
    	_evt.trk_dEdx_trunk_pion = _evt.trk_dEdx_trunk_max;
    	_evt.trk_bragg_pion_pion = _evt.trk_bragg_pion_max;
    	_evt.trk_llr_pid_score_pion = _evt.trk_llr_pid_score;
    	_evt.trk_score_pion = _evt.trk_score;
    	_evt.trk_end_spacepoints_pion = _evt.trk_end_spacepoints;
    	_evt.reco_momentum_pion = _evt.trk_momentum_pion;
    }
    else if (_evt.secondaryTrackPionlike) {
    	_evt.trk_bragg_mip_pion = _evt.trk2_bragg_mip_max;
     	_evt.trk_daughters_pion = _evt.trk2_daughters;
    	_evt.trk_dEdx_trunk_pion = _evt.trk2_dEdx_trunk_max;
    	_evt.trk_bragg_pion_pion = _evt.trk2_bragg_pion_max;
    	_evt.trk_llr_pid_score_pion = _evt.trk2_llr_pid_score;
    	_evt.trk_score_pion = _evt.trk2_score;
    	_evt.trk_end_spacepoints_pion = _evt.trk2_end_spacepoints;
    	_evt.reco_momentum_pion = _evt.trk2_momentum_pion;
	}
	else {
		_evt.trk_bragg_mip_pion = _evt.trk3_bragg_mip_max;
     	_evt.trk_daughters_pion = _evt.trk3_daughters;
    	_evt.trk_dEdx_trunk_pion = _evt.trk3_dEdx_trunk_max;
    	_evt.trk_bragg_pion_pion = _evt.trk3_bragg_pion_max;
    	_evt.trk_llr_pid_score_pion = _evt.trk3_llr_pid_score;
    	_evt.trk_score_pion = _evt.trk3_score;
    	_evt.trk_end_spacepoints_pion = _evt.trk3_end_spacepoints;
    	_evt.reco_momentum_pion = _evt.trk3_momentum_pion;
	}
}

