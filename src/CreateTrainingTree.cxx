#include "../include/CreateTrainingTree.h"

#include <iostream>

CreateTrainingTree::CreateTrainingTree(Utility::BDTEnums type) {
	
	std::cout << "Initialising Training Tree Creation Class" << std::endl;

	if (type == Utility::kElectronPhoton) {
		
		std::cout << "Creating electron photon separation BDT training tree." << std::endl;
		BDTType = Utility::kElectronPhoton;

		// create file
		outputFile = new TFile("TrainingTree_ElectronPhotonSeparation.root", "RECREATE", "Output File");
		if (outputFile->IsOpen()) std::cout << "Output file created successfully." << std::endl;
		else {
			std::cout << "Training Tree: output file could not be opened, check file name is valid." << std::endl << std::endl;
			exit(1);
		}

		// create tree
		trainingTree = new TTree("trainingTree", "trainingTree");

		// create branches
		trainingTree->Branch("isSignal", &isSignal);
		trainingTree->Branch("n_showers_contained", &n_showers_contained);
		trainingTree->Branch("shr_distance", &shr_distance);
		trainingTree->Branch("shr_trkfit_gap10_dedx_max", &shr_trkfit_gap10_dedx_max);
		trainingTree->Branch("shr_trkfit_2cm_dedx_max", &shr_trkfit_2cm_dedx_max);
		trainingTree->Branch("shrmoliereavg", &shrmoliereavg);
		trainingTree->Branch("shr_energyFraction", &shr_energyFraction);
		trainingTree->Branch("shrsubclusters", &shrsubclusters);
		trainingTree->Branch("shrPCA1CMed_5cm", &shrPCA1CMed_5cm);
		trainingTree->Branch("CylFrac2h_1cm", &CylFrac2h_1cm);
		trainingTree->Branch("DeltaRMS2h", &DeltaRMS2h);
		trainingTree->Branch("shrMCSMom", &shrMCSMom);

		trainingTree->Branch("shr_energy_second_cali", &shr_energy_second_cali);
		trainingTree->Branch("shr2_distance", &shr2_distance);
		trainingTree->Branch("shr12_p1_dstart", &shr12_p1_dstart);
		trainingTree->Branch("shr2_trackEndProximity", &shr2_trackEndProximity);
		trainingTree->Branch("shr2_pfpgeneration", &shr2_pfpgeneration);		

		trainingTree->Branch("secondshower_Y_nhit", &secondshower_Y_nhit);
		trainingTree->Branch("secondshower_Y_vtxdist", &secondshower_Y_vtxdist);
		trainingTree->Branch("secondshower_Y_dot", &secondshower_Y_dot);
		trainingTree->Branch("secondshower_Y_anglediff", &secondshower_Y_anglediff);
		trainingTree->Branch("secondshower_U_nhit", &secondshower_U_nhit);
		trainingTree->Branch("secondshower_U_vtxdist", &secondshower_U_vtxdist);
		trainingTree->Branch("secondshower_U_dot", &secondshower_U_dot);
		trainingTree->Branch("secondshower_U_anglediff", &secondshower_U_anglediff);
		trainingTree->Branch("secondshower_V_nhit", &secondshower_V_nhit);
		trainingTree->Branch("secondshower_V_vtxdist", &secondshower_V_vtxdist);
		trainingTree->Branch("secondshower_V_dot", &secondshower_V_dot);
		trainingTree->Branch("secondshower_V_anglediff", &secondshower_V_anglediff);

	}
	else if (type == Utility::kMuonPion) {
		
		std::cout << "Creating muon pion separation BDT training tree." << std::endl;
		BDTType = Utility::kMuonPion;

		// create file
		outputFile = new TFile("TrainingTree_MuonPionSeparation.root", "RECREATE", "Output File");
		if (outputFile->IsOpen()) std::cout << "Output file created successfully." << std::endl;
		else {
			std::cout << "Training Tree: output file could not be opened, check file name is valid." << std::endl << std::endl;
			exit(1);
		}

		// create tree
		trainingTree = new TTree("trainingTree", "trainingTree");

		// create branches
		trainingTree->Branch("isSignal", &isSignal);
		trainingTree->Branch("wc_track_length", &wc_track_length);
		trainingTree->Branch("wc_reco_pdg", &wc_reco_pdg);
		trainingTree->Branch("wc_reco_larpid_pidScore_mu", &wc_reco_larpid_pidScore_mu);
		trainingTree->Branch("wc_reco_larpid_pidScore_pi", &wc_reco_larpid_pidScore_pi);
		trainingTree->Branch("wc_reco_larpid_pidScore_pr", &wc_reco_larpid_pidScore_pr);
		trainingTree->Branch("wc_n_pion_candidate_daughters", &wc_n_pion_candidate_daughters);
		trainingTree->Branch("wc_pion_candidate_daughters_total_energy", &wc_pion_candidate_daughters_total_energy);
		trainingTree->Branch("wc_n_blips_25cm", &wc_n_blips_25cm);
		trainingTree->Branch("wc_n_blips_50cm", &wc_n_blips_50cm);
		trainingTree->Branch("wc_n_blips_100cm", &wc_n_blips_100cm);
	}
	else if (type == Utility::kPionProtonAlternate) {
		
		std::cout << "Creating alternate pion proton separation BDT training tree." << std::endl;
		BDTType = Utility::kPionProtonAlternate;

		// create file
		outputFile = new TFile("TrainingTree_PionProtonSeparation_Alternate.root", "RECREATE", "Output File");
		if (outputFile->IsOpen()) std::cout << "Output file created successfully." << std::endl;
		else {
			std::cout << "Training Tree: output file could not be opened, check file name is valid." << std::endl << std::endl;
			exit(1);
		}

		// create tree
		trainingTree = new TTree("trainingTree", "trainingTree");

		// create branches
		trainingTree->Branch("isSignal", &isSignal);
		trainingTree->Branch("trk_llr_pid_score", &trk_llr_pid_score);
		trainingTree->Branch("trk2_llr_pid_score", &trk2_llr_pid_score);
		trainingTree->Branch("trk3_llr_pid_score", &trk3_llr_pid_score);
		trainingTree->Branch("trk_bragg_mip_max", &trk_bragg_mip_max);
		trainingTree->Branch("trk2_bragg_mip_max", &trk2_bragg_mip_max);
		trainingTree->Branch("trk3_bragg_mip_max", &trk3_bragg_mip_max);
		trainingTree->Branch("trk_bragg_pion_max", &trk_bragg_pion_max);
		trainingTree->Branch("trk2_bragg_pion_max", &trk2_bragg_pion_max);
		trainingTree->Branch("trk3_bragg_pion_max", &trk3_bragg_pion_max);
		trainingTree->Branch("trk_dEdx_trunk_max", &trk_dEdx_trunk_max);
		trainingTree->Branch("trk2_dEdx_trunk_max", &trk2_dEdx_trunk_max);
		trainingTree->Branch("trk3_dEdx_trunk_max", &trk3_dEdx_trunk_max);
		trainingTree->Branch("trk_daughters", &trk_daughters);
		trainingTree->Branch("trk2_daughters", &trk2_daughters);
		trainingTree->Branch("trk3_daughters", &trk3_daughters);
		trainingTree->Branch("trk_end_spacepoints", &trk_end_spacepoints);
		trainingTree->Branch("trk2_end_spacepoints", &trk2_end_spacepoints);
		trainingTree->Branch("trk3_end_spacepoints", &trk3_end_spacepoints);

	}
	else {
		std::cout << "Training Tree: invalid BDT type provided." << std::endl << std::endl;
		exit(1);
	}
}

CreateTrainingTree::~CreateTrainingTree(){
	delete trainingTree;
	delete outputFile;
}

// ------------------------------------------

void CreateTrainingTree::addEvent(const EventContainer &_evt, Utility::ClassificationEnums classification) {

	if (BDTType == Utility::kElectronPhoton) {

		// populate variables, sanitising values
		n_showers_contained = _evt.n_showers_contained;

		shr_distance = _evt.shr_distance;

		if (_evt.shr_trkfit_gap10_dedx_max >= 0) shr_trkfit_gap10_dedx_max = _evt.shr_trkfit_gap10_dedx_max;
		else shr_trkfit_gap10_dedx_max = 9999;
		if (_evt.shr_trkfit_2cm_dedx_max >= 0) shr_trkfit_2cm_dedx_max = _evt.shr_trkfit_2cm_dedx_max;
		else shr_trkfit_gap10_dedx_max = 9999;

		if (_evt.shrmoliereavg >= 0) shrmoliereavg = _evt.shrmoliereavg;
		else shrmoliereavg = 9999;
		shr_energyFraction = _evt.shr_energyFraction;
		shrsubclusters = _evt.shrsubclusters;

		if (_evt.shrPCA1CMed_5cm > 0) shrPCA1CMed_5cm = _evt.shrPCA1CMed_5cm;
		else shrPCA1CMed_5cm = 9999;
		if (_evt.CylFrac2h_1cm > 0) CylFrac2h_1cm = _evt.CylFrac2h_1cm;
		else CylFrac2h_1cm = 9999;
		if (_evt.DeltaRMS2h > 0) DeltaRMS2h = _evt.DeltaRMS2h;
		else DeltaRMS2h = 9999;
		if (_evt.shrMCSMom > 0) shrMCSMom = _evt.shrMCSMom;
		else shrMCSMom = 9999;

		if (_evt.shr_energy_second_cali > 0) {
			shr_energy_second_cali = _evt.shr_energy_second_cali;
			shr2_distance = _evt.shr2_distance;
			shr12_p1_dstart = _evt.shr12_p1_dstart;
			shr2_trackEndProximity = _evt.shr2_trackEndProximity;
			shr2_pfpgeneration = _evt.shr2_pfpgeneration;
		}
		else {
			shr_energy_second_cali = 9999;
			shr2_distance = 9999;
			shr12_p1_dstart = 9999;
			shr2_trackEndProximity = 9999;
			shr2_pfpgeneration = 9999;
		}
		if (shr2_pfpgeneration == 0) shr2_pfpgeneration = 9999;		

		if (_evt.secondshower_Y_nhit >= 0) secondshower_Y_nhit = _evt.secondshower_Y_nhit;
		else secondshower_Y_nhit = 9999;
		if (_evt.secondshower_Y_vtxdist >= 0 && _evt.secondshower_Y_vtxdist < 9000) secondshower_Y_vtxdist = _evt.secondshower_Y_vtxdist;
		else secondshower_Y_vtxdist = 9999;
		if (_evt.secondshower_Y_dot >= 0 && _evt.secondshower_Y_dot < 9000) secondshower_Y_dot = _evt.secondshower_Y_dot;
		else secondshower_Y_dot = 9999;
		if (_evt.secondshower_Y_anglediff >= 0 && _evt.secondshower_Y_anglediff < 9000) secondshower_Y_anglediff = _evt.secondshower_Y_anglediff;
		else secondshower_Y_anglediff = 9999;

		if (_evt.secondshower_U_nhit >= 0) secondshower_U_nhit = _evt.secondshower_U_nhit;
		else secondshower_U_nhit = 9999;
		if (_evt.secondshower_U_vtxdist >= 0 && _evt.secondshower_U_vtxdist < 9000) secondshower_U_vtxdist = _evt.secondshower_U_vtxdist;
		else secondshower_U_vtxdist = 9999;
		if (_evt.secondshower_U_dot >= 0 && _evt.secondshower_U_dot < 9000) secondshower_U_dot = _evt.secondshower_U_dot;
		else secondshower_U_dot = 9999;
		if (_evt.secondshower_U_anglediff >= 0 && _evt.secondshower_U_anglediff < 9000) secondshower_U_anglediff = _evt.secondshower_U_anglediff;
		else secondshower_U_anglediff = 9999;

		if (_evt.secondshower_V_nhit >= 0) secondshower_V_nhit = _evt.secondshower_V_nhit;
		else secondshower_V_nhit = 9999;
		if (_evt.secondshower_V_vtxdist >= 0 && _evt.secondshower_V_vtxdist < 9000) secondshower_V_vtxdist = _evt.secondshower_V_vtxdist;
		else secondshower_V_vtxdist = 9999;
		if (_evt.secondshower_V_dot >= 0 && _evt.secondshower_V_dot < 9000) secondshower_V_dot = _evt.secondshower_V_dot;
		else secondshower_V_dot = 9999;
		if (_evt.secondshower_V_anglediff >= 0 && _evt.secondshower_V_anglediff < 9000) secondshower_V_anglediff = _evt.secondshower_V_anglediff;
		else secondshower_V_anglediff = 9999;

		// determine whether signal or background, and save
		if (classification == Utility::kNC1pi) {
				isSignal = 1;
				trainingTree->Fill();
				return;
		}
		else if (classification == Utility::kCCNue) {
			isSignal = 0;
			trainingTree->Fill();
			return;
		}
		else {
			return;
		}
	}
	else if (BDTType == Utility::kPionProton) {

		// primary track
		if (_evt.primaryTrackPionlikeLoose) {

			// populate variables
			trk_score = _evt.trk_score;
			trk_llr_pid_score = _evt.trk_llr_pid_score;
			trk_bragg_mip_max = _evt.trk_bragg_mip_max;
			trk_bragg_pion_max = _evt.trk_bragg_pion_max;
			trk_dEdx_trunk_max = _evt.trk_dEdx_trunk_max;
			trk_daughters = _evt.trk_daughters;
			trk_end_spacepoints = _evt.trk_end_spacepoints;

			// determine whether signal or background, and save
			if (_evt.trk_bkt_pdg == 211 || _evt.trk_bkt_pdg == -211) {
				isSignal = 1;
				trainingTree->Fill();
			}
			else if (_evt.trk_bkt_pdg == 2212) {
				isSignal = 0;
				trainingTree->Fill();
			}

			// otherwise don't fill
		}
		// secondary track
		if (_evt.secondaryTrackPionlikeLoose) {

			// populate variables
			trk_score = _evt.trk2_score; 
			trk_llr_pid_score = _evt.trk2_llr_pid_score;
			trk_bragg_mip_max = _evt.trk2_bragg_mip_max;
			trk_bragg_pion_max = _evt.trk2_bragg_pion_max;
			trk_dEdx_trunk_max = _evt.trk2_dEdx_trunk_max;
			trk_daughters = _evt.trk2_daughters;
			trk_end_spacepoints = _evt.trk2_end_spacepoints;

			// determine whether signal or background, and save
			if (_evt.trk2_bkt_pdg == 211 || _evt.trk2_bkt_pdg == -211) {
				isSignal = 1;
				trainingTree->Fill();
			}
			else if (_evt.trk2_bkt_pdg == 2212) {
				isSignal = 0;
				trainingTree->Fill();
			}

			// otherwise don't fill
		}
		// tertiary track
		if (_evt.tertiaryTrackPionlikeLoose) {

			// populate variables
			trk_score = _evt.trk3_score; 
			trk_llr_pid_score = _evt.trk3_llr_pid_score;
			trk_bragg_mip_max = _evt.trk3_bragg_mip_max;
			trk_bragg_pion_max = _evt.trk3_bragg_pion_max;
			trk_dEdx_trunk_max = _evt.trk3_dEdx_trunk_max;
			trk_daughters = _evt.trk3_daughters;
			trk_end_spacepoints = _evt.trk3_end_spacepoints;

			// determine whether signal or background, and save
			if (_evt.trk3_bkt_pdg == 211 || _evt.trk3_bkt_pdg == -211) {
				isSignal = 1;
				trainingTree->Fill();
			}
			else if (_evt.trk3_bkt_pdg == 2212) {
				isSignal = 0;
				trainingTree->Fill();
			}

			// otherwise don't fill
		}
	}
	else if (BDTType == Utility::kPionProtonAlternate) {

		// primary track
		if (_evt.primaryTrackPionlikeLoose) {
			// populate variables 
			trk_llr_pid_score = _evt.trk_llr_pid_score;
			trk_bragg_mip_max = _evt.trk_bragg_mip_max;
			trk_bragg_pion_max = _evt.trk_bragg_pion_max;
			trk_dEdx_trunk_max = _evt.trk_dEdx_trunk_max;
			trk_daughters = _evt.trk_daughters;
			trk_end_spacepoints = _evt.trk_end_spacepoints;
		}
		else {
			// fill with missing value
			trk_llr_pid_score = 9999;
			trk_bragg_mip_max = 9999;
			trk_bragg_pion_max = 9999;
			trk_dEdx_trunk_max = 9999;
			trk_daughters = 9999;
			trk_end_spacepoints = 9999;
		}
		// secondary track
		if (_evt.secondaryTrackPionlikeLoose) {

			// populate variables 
			trk2_llr_pid_score = _evt.trk2_llr_pid_score;
			trk2_bragg_mip_max = _evt.trk2_bragg_mip_max;
			trk2_bragg_pion_max = _evt.trk2_bragg_pion_max;
			trk2_dEdx_trunk_max = _evt.trk2_dEdx_trunk_max;
			trk2_daughters = _evt.trk2_daughters;
			trk2_end_spacepoints = _evt.trk2_end_spacepoints;
		}
		else {
			// fill with missing value
			trk2_llr_pid_score = 9999;
			trk2_bragg_mip_max = 9999;
			trk2_bragg_pion_max = 9999;
			trk2_dEdx_trunk_max = 9999;
			trk2_daughters = 9999;
			trk2_end_spacepoints = 9999;
		}
		// tertiary track
		if (_evt.tertiaryTrackPionlikeLoose) {

			// populate variables 
			trk3_llr_pid_score = _evt.trk3_llr_pid_score;
			trk3_bragg_mip_max = _evt.trk3_bragg_mip_max;
			trk3_bragg_pion_max = _evt.trk3_bragg_pion_max;
			trk3_dEdx_trunk_max = _evt.trk3_dEdx_trunk_max;
			trk3_daughters = _evt.trk3_daughters;
			trk3_end_spacepoints = _evt.trk3_end_spacepoints;
		}
		else {
			// fill with missing value
			trk3_llr_pid_score = 9999;
			trk3_bragg_mip_max = 9999;
			trk3_bragg_pion_max = 9999;
			trk3_dEdx_trunk_max = 9999;
			trk3_daughters = 9999;
			trk3_end_spacepoints = 9999;
		}

		// determine whether signal or background, and save
		if (classification == Utility::kNC1pi) {
			isSignal = 1;
			trainingTree->Fill();
			return;
		}
		else if (classification == Utility::kCCNue) {
			isSignal = 0;
			trainingTree->Fill();
			return;
		}
		else {
			return;
		}

		return;
	}
	else if (BDTType == Utility::kMuonPion) {

		// loop through tracks
		for (unsigned int wc_idx = 0; wc_idx < _evt.wc_reco_Ntrack; wc_idx++) {

			// only keep pion candidate
			if (wc_idx != _evt.wc_pion_candidate_index) continue;

			// keep only truth matched muons, pions and protons
			if (_evt.wc_particle_classification_v[wc_idx] != Utility::kMuon && 
				_evt.wc_particle_classification_v[wc_idx] != Utility::kPionDecay //&& 
				//_evt.wc_particle_classification_v[wc_idx] != Utility::kProton
			) continue;
			
			// populate variables
			wc_track_length = std::sqrt(
			std::pow(_evt.wc_reco_startXYZT[wc_idx][0] - _evt.wc_reco_endXYZT[wc_idx][0], 2) +
			std::pow(_evt.wc_reco_startXYZT[wc_idx][1] - _evt.wc_reco_endXYZT[wc_idx][1], 2) +
			std::pow(_evt.wc_reco_startXYZT[wc_idx][2] - _evt.wc_reco_endXYZT[wc_idx][2], 2) );

			wc_reco_pdg = static_cast<float>(_evt.wc_reco_pdg[wc_idx]);

			wc_reco_larpid_pidScore_mu = _evt.wc_reco_larpid_pidScore_mu[wc_idx];
			wc_reco_larpid_pidScore_pi = _evt.wc_reco_larpid_pidScore_pi[wc_idx];
			wc_reco_larpid_pidScore_pr = _evt.wc_reco_larpid_pidScore_pr[wc_idx];

			wc_n_pion_candidate_daughters = static_cast<float>(_evt.wc_n_pion_candidate_daughters);
			wc_pion_candidate_daughters_total_energy = _evt.wc_pion_candidate_daughters_total_energy;

			wc_n_blips_25cm = static_cast<float>(_evt.wc_n_blips_25cm);
			wc_n_blips_50cm = static_cast<float>(_evt.wc_n_blips_50cm);
			wc_n_blips_100cm = static_cast<float>(_evt.wc_n_blips_100cm);

			// determine whether signal or background, and save
			if (_evt.wc_particle_classification_v[wc_idx] == Utility::kPionDecay) {
				isSignal = 1;
				trainingTree->Fill();
			}
			else if (_evt.wc_particle_classification_v[wc_idx] == Utility::kMuon) { // || _evt.wc_particle_classification_v[wc_idx] == Utility::kProton) {
				isSignal = 0;
				trainingTree->Fill();
			}			
		}
	}
}

// --------------------------------------------

void CreateTrainingTree::writeOutputFile() { 
	outputFile->Write(); 
}

