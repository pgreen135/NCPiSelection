#include "../include/SelectionDriver.h"

#include <iostream>
#include <iomanip>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TLine.h>

#include "../include/EventContainer.h"
#include "../include/CreateTrainingTree.h"
#include "../include/Selection.h"
#include "../include/StackedHistTool.h"
#include "../include/BDTTool.h"


// Constructor
SelectionDriver::SelectionDriver(): _utility(false, false, false, false) {
	
	std::cout << "Initialising Selection Driver Class" << std::endl;

}

// ------------------------------------------------------------------------------

// Draw stacked histogram


// ------------------------------------------------------------------------------

// Run BDT selection on all samples
void SelectionDriver::runBDTSelectionFull() {

	// lists of file names, weights and types to run over
	
	// SURPRISE test samples
	std::vector<std::string> filename_list = {filename_mc_run4b, filename_dirt_run4b, filename_beamoff_run4b,
											  filename_mc_run4c, filename_dirt_run4c,
											  filename_mc_run4d, filename_dirt_run4d,
											  filename_mc_run5,  filename_dirt_run5
											  };
	std::vector<double> pot_weight_list = {pot_weight_mc_run4b, pot_weight_dirt_run4b, pot_weight_beamoff_run4b,
										   pot_weight_mc_run4c, pot_weight_dirt_run4c,
										   pot_weight_mc_run4d, pot_weight_dirt_run4d,
										   pot_weight_mc_run5,  pot_weight_dirt_run5
										   };
	std::vector<Utility::FileTypeEnums> file_types_list = {Utility::kMC, Utility::kDirt, Utility::kEXT,
												 		   Utility::kMC, Utility::kDirt,
												           Utility::kMC, Utility::kDirt,
												           Utility::kMC, Utility::kDirt
													    };
	std::vector<Utility::RunPeriodEnums> run_periods_list = {Utility::kRun4ab, Utility::kRun4ab, Utility::kRun4ab,
															 Utility::kRun4ab, Utility::kRun4ab,
															 Utility::kRun4ab, Utility::kRun4ab,
															 Utility::kRun4ab,    Utility::kRun4ab
															};

	// construct classes
	Selection _selection(_utility);
	CreateTrainingTree _trainingTree(Utility::kMuonPion);
	BDTTool _BDTTool(false, false, true, false, _utility);

	// construct stacks
	StackedHistTool _histStack_contained_fraction("", "", 3, -1.05, 2.10, _utility);
	StackedHistTool _histStack_single_bin_interaction("", "", 3, -1.05, 2.10, _utility, "Interaction");
	StackedHistTool _histStack_associated_hits_fraction("", "", 20, 0, 1.001, _utility);
	StackedHistTool _histStack_nu_vtx_x("", "", 100, -50, 50, _utility);
	StackedHistTool _histStack_nu_vtx_y("", "", 100, -50, 50, _utility);
	StackedHistTool _histStack_nu_vtx_z("", "", 100, -50, 50, _utility);
	StackedHistTool _histStack_trk_length("", "", 60, 0, 300, _utility);
	StackedHistTool _histStack_trk2_length("", "", 60, 0, 300, _utility);
	StackedHistTool _histStack_trk_score("", "", 20, 0, 1.001, _utility);
	StackedHistTool _histStack_trk_llrpid("", "", 20, -1, 1.001, _utility);
	StackedHistTool _histStack_trk_vertex_distance("", "", 15, 0, 15, _utility);
	StackedHistTool _histStack_topological_score("", "", 20, 0, 1.001, _utility);
	StackedHistTool _histStack_cosmic_IP("", "", 20, 0, 100, _utility);
	StackedHistTool _histStack_number_tracks("", "", 6, 0, 6, _utility);
	StackedHistTool _histStack_number_showers("", "", 6, 0, 6, _utility);
	StackedHistTool _histStack_n_non_proton_tracks("", "", 6, 0, 6, _utility);
	StackedHistTool _histStack_n_pion_candidate_tracks("", "", 6, 0, 6, _utility);

	StackedHistTool _histStack_wc_numu_cc_flag("", "", 3, -1.05, 1.05, _utility);
	StackedHistTool _histStack_vertex_discrepancy("", "", 50, 0, 50, _utility);
	StackedHistTool _histStack_wc_match_isFC("", "", 2, -0.05, 1.05, _utility);
	StackedHistTool _histStack_wc_reco_larpid_pidScore_el("", "", 40, -15, 5, _utility);
	StackedHistTool _histStack_wc_reco_larpid_pidScore_ph("", "", 40, -15, 5, _utility);
	StackedHistTool _histStack_wc_reco_larpid_pidScore_mu("", "", 40, -15, 5, _utility);
	StackedHistTool _histStack_wc_reco_larpid_pidScore_pi("", "", 40, -15, 5, _utility);
	StackedHistTool _histStack_wc_reco_larpid_pidScore_pr("", "", 40, -15, 5, _utility);
	StackedHistTool _histStack_wc_reco_larpid_procScore_prim("", "", 40, -15, 5, _utility);
	StackedHistTool _histStack_wc_reco_larpid_procScore_chgd("", "", 40, -15, 5, _utility);
	StackedHistTool _histStack_wc_reco_larpid_procScore_ntrl("", "", 40, -15, 5, _utility);
	StackedHistTool _histStack_wc_bdt_score("", "", 20, 0, 1, _utility);

	// per particle stacks
	StackedHistTool _histStack_trk_length_particle("", "", 30, 0, 300, _utility, "Particle");
	StackedHistTool _histStack_trk_score_particle("", "", 20, 0, 1.001, _utility, "Particle");
	
	StackedHistTool _histStack_trk_llrpid_particle("", "", 20, -1, 1.001, _utility, "Particle");
	StackedHistTool _histStack_trk_bragg_mu_particle("", "", 20, 0, 1.001, _utility, "Particle");
	StackedHistTool _histStack_trk_bragg_p_particle("", "", 20, 0, 1.001, _utility, "Particle");
	StackedHistTool _histStack_trk_bragg_pion_particle("", "", 20, 0, 1.001, _utility, "Particle");
	StackedHistTool _histStack_trk_bragg_mip_particle("", "", 20, 0, 1.001, _utility, "Particle");
	StackedHistTool _histStack_trk_dEdx_trunk_particle("", "", 20, 0, 10, _utility, "Particle");

	StackedHistTool _histStack_trk_daughters_particle("", "", 6, 0, 6, _utility, "Particle");
	StackedHistTool _histStack_trk_end_spacepoints_particle("", "", 50, 0, 500, _utility, "Particle");
	StackedHistTool _histStack_trk_avg_deflection_stdev_particle("", "", 50, 0, 5e-3, _utility, "Particle");
	
	// NuGraph
	StackedHistTool _histStack_pfng2mipfrac_particle("", "", 20, 0, 1.001, _utility, "Particle");
	StackedHistTool _histStack_pfng2hipfrac_particle("", "", 20, 0, 1.001, _utility, "Particle");
	StackedHistTool _histStack_pfng2shrfrac_particle("", "", 20, 0, 1.001, _utility, "Particle");
	StackedHistTool _histStack_pfng2mclfrac_particle("", "", 20, 0, 1.001, _utility, "Particle");
	StackedHistTool _histStack_pfng2dfsfrac_particle("", "", 20, 0, 1.001, _utility, "Particle");

	// LArPID
	StackedHistTool _histStack_wc_reco_larpid_pidScore_el_particle("", "", 40, -15, 5, _utility, "Particle");
	StackedHistTool _histStack_wc_reco_larpid_pidScore_ph_particle("", "", 40, -15, 5, _utility, "Particle");
	StackedHistTool _histStack_wc_reco_larpid_pidScore_mu_particle("", "", 40, -15, 5, _utility, "Particle");
	StackedHistTool _histStack_wc_reco_larpid_pidScore_pi_particle("", "", 40, -15, 5, _utility, "Particle");
	StackedHistTool _histStack_wc_reco_larpid_pidScore_pr_particle("", "", 40, -15, 5, _utility, "Particle");
	StackedHistTool _histStack_wc_reco_larpid_pdg_particle("", "", 2500, 0, 2500, _utility, "Particle");

	// Other WC
	StackedHistTool _histStack_wc_n_pion_candidate_daughters_particle("", "", 6, 0, 6, _utility, "Particle");
	StackedHistTool _histStack_wc_pion_candidate_daughters_total_energy_particle("", "", 25, 0, 0.5, _utility, "Particle");

	// Blips
	StackedHistTool _histStack_wc_n_blips_3cm("", "", 25, 0, 25, _utility, "Particle");
	StackedHistTool _histStack_wc_n_blips_5cm("", "", 25, 0, 25, _utility, "Particle");
	StackedHistTool _histStack_wc_n_blips_10cm("", "", 15, 0, 30, _utility, "Particle");

	// lantern
	StackedHistTool _histStack_lantern_vtxScore("", "", 40, 0, 1.001, _utility);
	StackedHistTool _histStack_lantern_vtxContainment("", "", 4, -1, 3, _utility);
	StackedHistTool _histStack_lantern_vtxFracHitsOnCosmic("", "", 20, 0, 1.001, _utility);
	StackedHistTool _histStack_lantern_larpid_pidScore_mu("", "", 20, 0, 1.001, _utility);
	StackedHistTool _histStack_lantern_larpid_pidScore_pi("", "", 20, 0, 1.001, _utility);
	StackedHistTool _histStack_lantern_larpid_pidScore_pr("", "", 20, 0, 1.001, _utility);
	StackedHistTool _histStack_lantern_larpid_pidScore_mu_particle("", "", 20, 0, 1.001, _utility, "Particle");
	StackedHistTool _histStack_lantern_larpid_pidScore_pi_particle("", "", 40, -15, 1, _utility, "Particle");
	StackedHistTool _histStack_lantern_larpid_pidScore_pr_particle("", "", 20, 0, 1.001, _utility, "Particle");

	StackedHistTool _histStack_lantern_larpid_mupi_llr_particle("", "", 20, -1, 1.001, _utility, "Particle");
	StackedHistTool _histStack_lantern_larpid_prpi_llr_particle("", "", 20, -1, 1.001, _utility, "Particle");
	StackedHistTool _histStack_lantern_larpid_prmu_llr_particle("", "", 20, -1, 1.001, _utility, "Particle");
	
	// Counters
	int n_data_pass = 0;
	double n_signal_all = 0;
	double n_signal_pass = 0;
	double n_total_pass = 0;

	double n_signal_decay = 0;
	double n_signal_absorption0p = 0;
	double n_signal_absorptionNp = 0;
	double n_signal_chargeexchange = 0;
	
	double n_signal_pass_decay = 0;
	double n_signal_pass_absorption0p = 0;
	double n_signal_pass_absorptionNp = 0;
	double n_signal_pass_chargeexchange = 0;

	double n_track_muon = 0;
	double n_track_decay = 0;
	double n_track_absorption0p = 0;
	double n_track_absorptionNp = 0;
	double n_track_chargeexchange = 0;
	
	double n_track_muon_cut = 0;
	double n_track_decay_cut = 0;
	double n_track_absorption0p_cut = 0;
	double n_track_absorptionNp_cut = 0;
	double n_track_chargeexchange_cut = 0;
	
	double muon_counter = 0;
	double pion_counter = 0;

	double muon_counter_cut = 0;
	double pion_counter_cut = 0;

	// driver loop
	for (unsigned int idx = 0; idx < filename_list.size(); idx++) {
		//if (idx == 0) continue;

		std::cout << "Processing file " << filename_list[idx] << ", file type: " << file_types_list[idx] << ", run period: " << run_periods_list[idx]  << std::endl;

		// read file
		TFile *f = new TFile(filename_list[idx].c_str()); 
	 	TTree *pelee_tree = (TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
		TTree *wc_BDTvars_tree = (TTree*)f->Get("wcpselection/T_BDTvars");
		TTree *wc_eval_tree = (TTree*)f->Get("wcpselection/T_eval");
		TTree *wc_PFeval_tree = (TTree*)f->Get("wcpselection/T_PFeval");
		TTree *lantern_tree = (TTree*)f->Get("lantern/EventTree");

	 	// initialise event container
	  	EventContainer _event(pelee_tree, wc_eval_tree, wc_BDTvars_tree, wc_PFeval_tree, lantern_tree, _utility);

	  	// loop through events
	  	int n_entries = pelee_tree->GetEntries();
	  	std::cout << "Initial number events: " << n_entries << std::endl;

	  	for (int e = 0; e < n_entries; e++) {
		//for (int e = 180000; e < 190000; e++) {

			// skip events with problems in test sample
			// need to load PeLEE tree first to get RSE
			pelee_tree->GetEntry(e);

			// why are these an issue?
			if (file_types_list[idx] == Utility::kMC) {
				//std::cout << "Processing event " << e << std::endl;
				// split sample, test
				//if (e == 115063) continue;
				// split sample, train
				//if (e == 275274) continue;
				
				// Run 4b
				if (_event.run == 20872 && _event.sub == 14 && _event.evt == 749) {
					std::cout << "Skipping problematic MC event: Run " << _event.run << ", Subrun " << _event.sub << ", Event " << _event.evt << std::endl;
					continue;
				}
				if (_event.run == 20607 && _event.sub == 355 && _event.evt == 17791) {
					std::cout << "Skipping problematic MC event: Run " << _event.run << ", Subrun " << _event.sub << ", Event " << _event.evt << std::endl;
					continue;
				}

				//if (e == 180530) continue;
				//if (e == 564197) continue;
			}
			
			// get event for each tree
			wc_BDTvars_tree->GetEntry(e);
			//std::cout << "Loaded WC BDT Tree" << std::endl;	
			wc_eval_tree->GetEntry(e);
			//std::cout << "Loaded WC Eval Tree" << std::endl;
			wc_PFeval_tree->GetEntry(e); 
			//std::cout << "Loaded WC PFEval Tree" << std::endl;
			lantern_tree->GetEntry(e);
			
		    if ( (e != 0) && (n_entries >= 10) &&  (e % (n_entries/10) == 0) ) {
		      std::cout << Form("%i0%% Completed...\n", e / (n_entries/10));
		    }
			
			//std::cout << "Processing event " << e << ": Run " << _event.run << ", Subrun " << _event.sub << ", Event " << _event.evt << std::endl;

			// RSE check
			if (_event.run != _event.wc_run || _event.sub != _event.wc_sub || _event.evt != _event.wc_evt) {
				std::cout << "RSE mismatch for event: " << e << ". Pandora Run: " << _event.run << ", Subrun: " << _event.sub << ", Event: " << _event.evt 
				          << ". WC Run: " << _event.wc_run << ", Subrun: " << _event.wc_sub << ", Event: " << _event.wc_evt << std::endl;
			}

		    // bool passSelection = _selection.ApplyCutBasedSelection(_event, file_types_list[idx], run_periods_list[idx]);
			// WC LArPID selection
			//bool passSelection = _selection.ApplyWCSelection(_event, _BDTTool, file_types_list[idx], run_periods_list[idx]);
			
			// Lantern selection
			bool passSelection = _selection.ApplyLanternSelection(_event, file_types_list[idx], run_periods_list[idx]);
			
		    // evaluate correct POT weight accounting for special cases
		    double pot_weight = pot_weight_list[idx];

			// increment all signal counter
			if (_event.classification == Utility::kNC1pi) {
				n_signal_all += pot_weight * _event.weight_cv;

				if (_event.truePionHasMuon) {
					n_signal_decay += pot_weight * _event.weight_cv;
				}
				else if (_event.truePionHasPi0) {
					n_signal_chargeexchange += pot_weight * _event.weight_cv;
				}
				else if (_event.truePionNProtons == 0) {
					n_signal_absorption0p += pot_weight * _event.weight_cv;
				}
				else if (_event.truePionNProtons > 0) {
					n_signal_absorptionNp += pot_weight * _event.weight_cv;
				}
				
				/*
				// truth pion process
				for (unsigned int idx_track = 0; idx_track < _event.lantern_nTracks; idx_track++) {

					if (_event.lantern_particle_classification_v[idx_track] == Utility::kPionDecay) {
						n_signal_decay += pot_weight * _event.weight_cv;
						break;
					}
					else if (_event.lantern_particle_classification_v[idx_track] == Utility::kPionAbsorption0p) {
						n_signal_absorption0p += pot_weight * _event.weight_cv;
						break;
					}
					else if (_event.lantern_particle_classification_v[idx_track] == Utility::kPionAbsorptionNp) {
						n_signal_absorptionNp += pot_weight * _event.weight_cv;
						break;
					}
					else if (_event.lantern_particle_classification_v[idx_track] == Utility::kPionChargeExchange) {
						n_signal_chargeexchange += pot_weight * _event.weight_cv;
						break;
					}
				}
				*/
			}
		    
		    if (!passSelection) continue;
		   
		    if (file_types_list[idx] == Utility::kData) {
		    	n_data_pass++;
		    	std::cout << "Passing data event -- " << "Run: " << _event.run << ", Subrun: " << _event.sub << ", Event: " << _event.evt << ", Electron energy: " << _event.shr_energy_cali << ", Visible energy: " << _event.NeutrinoEnergy2/1000 << " GeV" << std::endl;
		    	continue;
		    }

			// output event of interest
			/*
			if (_event.classification == Utility::kCCNumuNpi) {
				if (_event.npion != 1) continue;
				if (_event.nproton != 0) continue;
				if (_event.pion_e < 0.3 || _event.pion_e > 0.4) continue;
				if (_event.muon_e < 0.3 || _event.muon_e > 0.4) continue;
				std::cout << "Run: " << _event.run << ", Subrun: " << _event.sub << ", Event: " << _event.evt;
				std::cout << ", Pandora truth: muon_e = " << _event.muon_e << ", npion" << _event.npion << ", pion_e" << _event.pion_e << ", nproton" << _event.nproton << std::endl;
			}
			*/

			// increment passing signal counter
			if (_event.classification == Utility::kNC1pi) {
				n_signal_pass += pot_weight * _event.weight_cv;

				if (_event.truePionHasMuon) {
					n_signal_pass_decay += pot_weight * _event.weight_cv;
				}
				else if (_event.truePionHasPi0) {
					n_signal_pass_chargeexchange += pot_weight * _event.weight_cv;
					//std::cout << "Run: " << _event.run << ", Subrun: " << _event.sub << ", Event: " << _event.evt << " passing charged exchange" << ", NDaughterProtons: " << _event.truePionNProtons << std::endl;
					//std::cout << "Pandora truth: npion" << _event.npion << ", pion_p" << _event.pion_p << ", nproton" << _event.nproton << std::endl;
				}
				else if (_event.truePionNProtons == 0) {
					n_signal_pass_absorption0p += pot_weight * _event.weight_cv;
				}
				else if (_event.truePionNProtons > 0) {
					n_signal_pass_absorptionNp += pot_weight * _event.weight_cv;
					
				}

				if (_event.truePionNProtons > 2){
					//std::cout << "Run: " << _event.run << ", Subrun: " << _event.sub << ", Event: " << _event.evt << " passing absorptionNp" << ", NProtons: " << _event.truePionNProtons << std::endl;
				}
				/*
				// truth pion process
				for (unsigned int idx_track = 0; idx_track < _event.lantern_nTracks; idx_track++) {
					if (_event.lantern_particle_classification_v[idx_track] == Utility::kPionDecay) {
						n_signal_pass_decay += pot_weight * _event.weight_cv;
						break;
					}
					else if (_event.lantern_particle_classification_v[idx_track] == Utility::kPionAbsorption0p) {
						n_signal_pass_absorption0p += pot_weight * _event.weight_cv;
						break;
					}
					else if (_event.lantern_particle_classification_v[idx_track] == Utility::kPionAbsorptionNp) {
						n_signal_pass_absorptionNp += pot_weight * _event.weight_cv;
						break;
					}
					else if (_event.lantern_particle_classification_v[idx_track] == Utility::kPionChargeExchange) {
						n_signal_pass_chargeexchange += pot_weight * _event.weight_cv;
						break;
					}
				}
				*/
			}

			// increment total passing counter
			n_total_pass += pot_weight * _event.weight_cv;

		    // fill histogram(s) - pelee variables
		    _histStack_contained_fraction.Fill(_event.classification, _event.contained_fraction, pot_weight * _event.weight_cv);
			_histStack_associated_hits_fraction.Fill(_event.classification, _event.associated_hits_fraction, pot_weight * _event.weight_cv);
			_histStack_nu_vtx_x.Fill(_event.classification, _event.true_nu_vtx_sce_x - _event.lantern_trueVtxX, pot_weight * _event.weight_cv);
			//_histStack_nu_vtx_x.Fill(_event.classification, _event.lantern_vtxX - _event.lantern_trueVtxX, pot_weight * _event.weight_cv);
			_histStack_nu_vtx_y.Fill(_event.classification, _event.lantern_vtxY - _event.lantern_trueVtxY, pot_weight * _event.weight_cv);
			_histStack_nu_vtx_z.Fill(_event.classification, _event.lantern_vtxZ - _event.lantern_trueVtxZ, pot_weight * _event.weight_cv);	    
			_histStack_topological_score.Fill(_event.classification, _event.topological_score, pot_weight * _event.weight_cv);
			_histStack_cosmic_IP.Fill(_event.classification, _event.CosmicIPAll3D, pot_weight * _event.weight_cv);		
			_histStack_number_tracks.Fill(_event.classification, _event.n_primary_tracks, pot_weight * _event.weight_cv);
			_histStack_number_showers.Fill(_event.classification, _event.n_primary_showers, pot_weight * _event.weight_cv);
			_histStack_trk_llrpid.Fill(_event.classification, _event.trk2_llr_pid_score, pot_weight * _event.weight_cv);
			_histStack_trk_score.Fill(_event.classification, _event.trk_score, pot_weight * _event.weight_cv);

			if (_event.classification == Utility::kNC1pi) {
				_histStack_single_bin_interaction.Fill(_event.interactionMode, _event.contained_fraction, pot_weight * _event.weight_cv);
			}

			_histStack_n_non_proton_tracks.Fill(_event.classification, _event.n_non_proton_tracks, pot_weight * _event.weight_cv);
			_histStack_n_pion_candidate_tracks.Fill(_event.classification, _event.n_pion_candidate_tracks, pot_weight * _event.weight_cv);
			
			// fill histogram(s) - WC variables
			_histStack_wc_numu_cc_flag.Fill(_event.classification, _event.wc_numu_cc_flag, pot_weight * _event.weight_cv);
			_histStack_vertex_discrepancy.Fill(_event.classification, _event.vertex_discrepancy, pot_weight * _event.weight_cv);
			_histStack_wc_match_isFC.Fill(_event.classification, _event.wc_match_isFC, pot_weight * _event.weight_cv);

			_histStack_trk_length.Fill(_event.classification, _event.wc_primaryTrackLength, pot_weight * _event.weight_cv);
			_histStack_trk2_length.Fill(_event.classification, _event.trk2_len, pot_weight * _event.weight_cv);
			_histStack_trk_vertex_distance.Fill(_event.classification, _event.lantern_primaryTrackVertexSeparation, pot_weight * _event.weight_cv);

			_histStack_wc_reco_larpid_pidScore_el.Fill(_event.classification, _event.wc_reco_larpid_pidScore_el[_event.wc_pion_candidate_index], pot_weight * _event.weight_cv);
			_histStack_wc_reco_larpid_pidScore_ph.Fill(_event.classification, _event.wc_reco_larpid_pidScore_ph[_event.wc_pion_candidate_index], pot_weight * _event.weight_cv);
			_histStack_wc_reco_larpid_pidScore_mu.Fill(_event.classification, _event.wc_reco_larpid_pidScore_mu[_event.wc_pion_candidate_index], pot_weight * _event.weight_cv);
			_histStack_wc_reco_larpid_pidScore_pi.Fill(_event.classification, _event.wc_reco_larpid_pidScore_pi[_event.wc_pion_candidate_index], pot_weight * _event.weight_cv);
			_histStack_wc_reco_larpid_pidScore_pr.Fill(_event.classification, _event.wc_reco_larpid_pidScore_pr[_event.wc_pion_candidate_index], pot_weight * _event.weight_cv);
			_histStack_wc_reco_larpid_procScore_prim.Fill(_event.classification, _event.wc_reco_larpid_procScore_prim[_event.wc_pion_candidate_index], pot_weight * _event.weight_cv);
			_histStack_wc_reco_larpid_procScore_chgd.Fill(_event.classification, _event.wc_reco_larpid_procScore_chgd[_event.wc_pion_candidate_index], pot_weight * _event.weight_cv);
			_histStack_wc_reco_larpid_procScore_ntrl.Fill(_event.classification, _event.wc_reco_larpid_procScore_ntrl[_event.wc_pion_candidate_index], pot_weight * _event.weight_cv);

			// lantern
			_histStack_lantern_vtxScore.Fill(_event.classification, _event.lantern_vtxScore, pot_weight * _event.weight_cv);
			_histStack_lantern_vtxContainment.Fill(_event.classification, _event.lantern_vtxContainment, pot_weight * _event.weight_cv);
			_histStack_lantern_vtxFracHitsOnCosmic.Fill(_event.classification, _event.lantern_vtxFracHitsOnCosmic, pot_weight * _event.weight_cv);

			// BDT
			_histStack_wc_bdt_score.Fill(_event.classification, _event.wc_bdt_score, pot_weight * _event.weight_cv);

			// fill histogram(s) - per particle
			// loop through primary tracks
			for (unsigned int i_trk = 0; i_trk < _event.trk_sce_start_x_v->size(); i_trk++) {			
		
				// check if track is primary
				if (_event.pfp_generation_v->at(i_trk) != 2) continue;

				// hits on all planes
				if (_event.pfnplanehits_U->at(i_trk) == 0 || _event.pfnplanehits_V->at(i_trk) == 0 || _event.pfnplanehits_Y->at(i_trk) == 0) continue;
				
				// track score
				if (_event.trk_score_v->at(i_trk) < 0.5) continue;

				// vertex distance
				if (_event.trk_distance_v->at(i_trk) > 4) continue; // 4cm

				// track length
				if (_event.trk_len_v->at(i_trk) < 5) continue; // 10cm
				if (_event.trk_len_v->at(i_trk) > 150) continue; // 200cm

				// fill
				//_histStack_trk_length_particle.Fill(_event.pfp_classification_v[i_trk], _event.trk_len_v->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_trk_score_particle.Fill(_event.pfp_classification_v[i_trk], _event.trk_score_v->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_trk_llrpid_particle.Fill(_event.pfp_classification_v[i_trk], _event.trk_llr_pid_score_v->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_trk_bragg_mu_particle.Fill(_event.pfp_classification_v[i_trk], _event.trk_bragg_mu_v->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_trk_bragg_p_particle.Fill(_event.pfp_classification_v[i_trk], _event.trk_bragg_p_v->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_trk_bragg_pion_particle.Fill(_event.pfp_classification_v[i_trk], _event.trk_bragg_pion_v->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_trk_bragg_mip_particle.Fill(_event.pfp_classification_v[i_trk], _event.trk_bragg_mip_v->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_trk_dEdx_trunk_particle.Fill(_event.pfp_classification_v[i_trk], _event.trk_trunk_dEdx_y_v->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_trk_daughters_particle.Fill(_event.pfp_classification_v[i_trk], _event.pfp_trk_daughters_v->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_trk_end_spacepoints_particle.Fill(_event.pfp_classification_v[i_trk], _event.trk_end_spacepoints_v->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_trk_avg_deflection_stdev_particle.Fill(_event.pfp_classification_v[i_trk], _event.trk_avg_deflection_stdev_v->at(i_trk), pot_weight * _event.weight_cv);
				// NuGraph
				_histStack_pfng2mipfrac_particle.Fill(_event.pfp_classification_v[i_trk], _event.pfng2mipfrac->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_pfng2hipfrac_particle.Fill(_event.pfp_classification_v[i_trk], _event.pfng2hipfrac->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_pfng2shrfrac_particle.Fill(_event.pfp_classification_v[i_trk], _event.pfng2shrfrac->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_pfng2mclfrac_particle.Fill(_event.pfp_classification_v[i_trk], _event.pfng2mclfrac->at(i_trk), pot_weight * _event.weight_cv);
				_histStack_pfng2dfsfrac_particle.Fill(_event.pfp_classification_v[i_trk], _event.pfng2dfsfrac->at(i_trk), pot_weight * _event.weight_cv);
			}

			// fill histogram(s) - WC per particle
			for (unsigned int wc_idx = 0; wc_idx < _event.wc_reco_Ntrack; wc_idx++) {

				// only fill pion candidate
				if (wc_idx != _event.wc_pion_candidate_index) continue;

				// fill only longest track
				//if (wc_idx != _event.wc_primaryTrackIndex) continue;

				// check primary
				//if (_event.wc_reco_mother[wc_idx] != 0) continue; // require WC primary track

				// track length
				// calculate track length
				float wc_track_length = std::sqrt(
				std::pow(_event.wc_reco_startXYZT[wc_idx][0] - _event.wc_reco_endXYZT[wc_idx][0], 2) +
				std::pow(_event.wc_reco_startXYZT[wc_idx][1] - _event.wc_reco_endXYZT[wc_idx][1], 2) +
				std::pow(_event.wc_reco_startXYZT[wc_idx][2] - _event.wc_reco_endXYZT[wc_idx][2], 2) );

				// track length
			    //if (wc_track_length < 5) continue; // 10cm
				//if (wc_track_length > 150) continue; // 150cm

				// track separation from vertex
				float wc_TrackVertexSeparation = std::sqrt(
				std::pow(_event.wc_reco_startXYZT[wc_idx][0] - _event.wc_reco_nuvtxX, 2) +
				std::pow(_event.wc_reco_startXYZT[wc_idx][1] - _event.wc_reco_nuvtxY, 2) +
				std::pow(_event.wc_reco_startXYZT[wc_idx][2] - _event.wc_reco_nuvtxZ, 2) );

				//if (wc_TrackVertexSeparation > 4) continue; // 4cm

				// fill histograms
				//_histStack_trk_length_particle.Fill(_event.wc_particle_classification_v[wc_idx], wc_track_length, pot_weight * _event.weight_cv);

				// LArPID unclassified tracks
				//if (_event.wc_reco_larpid_classified[wc_idx] == 0) continue; // only fill for classified tracks
				//if (_event.wc_reco_larpid_classified[wc_idx] != 0) continue; // only fill for classified tracks

				// skip WC pseudo-particles
				//if (_event.wc_reco_pdg[wc_idx] == 22 || _event.wc_reco_pdg[wc_idx] == 2112) continue; // skip pseudo-particles

				// fill histograms
				// vertex distance
				_histStack_trk_length_particle.Fill(_event.wc_particle_classification_v[wc_idx], wc_track_length, pot_weight * _event.weight_cv);
				// Other WC
				_histStack_wc_n_pion_candidate_daughters_particle.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_n_pion_candidate_daughters, pot_weight * _event.weight_cv);
				if (_event.wc_n_pion_candidate_daughters > 0) {
					_histStack_wc_pion_candidate_daughters_total_energy_particle.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_pion_candidate_daughters_total_energy, pot_weight * _event.weight_cv);
				}
				// Blips
				_histStack_wc_n_blips_3cm.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_n_blips_25cm, pot_weight * _event.weight_cv);
				_histStack_wc_n_blips_5cm.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_n_blips_50cm, pot_weight * _event.weight_cv);
				_histStack_wc_n_blips_10cm.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_n_blips_100cm, pot_weight * _event.weight_cv);
					
				// LARPID
				_histStack_wc_reco_larpid_pidScore_el_particle.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_reco_larpid_pidScore_el[wc_idx], pot_weight * _event.weight_cv);
				_histStack_wc_reco_larpid_pidScore_ph_particle.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_reco_larpid_pidScore_ph[wc_idx], pot_weight * _event.weight_cv);
				_histStack_wc_reco_larpid_pidScore_mu_particle.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_reco_larpid_pidScore_mu[wc_idx], pot_weight * _event.weight_cv);
				_histStack_wc_reco_larpid_pidScore_pi_particle.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_reco_larpid_pidScore_pi[wc_idx], pot_weight * _event.weight_cv);
				_histStack_wc_reco_larpid_pidScore_pr_particle.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_reco_larpid_pidScore_pr[wc_idx], pot_weight * _event.weight_cv);
				_histStack_wc_reco_larpid_pdg_particle.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_reco_pdg[wc_idx], pot_weight * _event.weight_cv);
			
			}

			/*
			// output details about CC Npi events
			if (_event.classification == Utility::kCCNumuNpi) {
				std::cout << "Passing CC Numu Npi event --- " << std::endl;
				std::cout << "Candidate track index: " << _event.wc_pion_candidate_index << std::endl;
				for (unsigned int wc_idx = 0; wc_idx < _event.wc_reco_Ntrack; wc_idx++) {

					// skip WC pseudo-particles
					if (_event.wc_reco_pdg[wc_idx] == 22 || _event.wc_reco_pdg[wc_idx] == 2112) continue; // skip pseudo-particles
					// check primary
					if (_event.wc_reco_mother[wc_idx] != 0) continue; // require WC primary track

					float wc_track_length = std::sqrt(
					std::pow(_event.wc_reco_startXYZT[wc_idx][0] - _event.wc_reco_endXYZT[wc_idx][0], 2) +
					std::pow(_event.wc_reco_startXYZT[wc_idx][1] - _event.wc_reco_endXYZT[wc_idx][1], 2) +
					std::pow(_event.wc_reco_startXYZT[wc_idx][2] - _event.wc_reco_endXYZT[wc_idx][2], 2) );

					std::cout << "Track ID: " << wc_idx << ", truth matched classification: " << _event.wc_particle_classification_v[wc_idx]
							  << ", WC PDG: " << _event.wc_reco_pdg[wc_idx] << ", LArPID PDG: " << _event.wc_reco_larpid_pdg[wc_idx]
							  << ", Length: " << wc_track_length
							  << std::endl;
				}
			}
			*/

			// add event to BDT training tree,
			//if (file_types_list[idx] == Utility::kMC) {
			//	_trainingTree.addEvent(_event, _event.classification);
			//}

			// lantern candidate particle
			// lantern candidate topology
			// loop over tracks
			for (unsigned int idx_track = 0; idx_track < _event.lantern_nTracks; idx_track++) {

				// only fill pion candidate
				if (idx_track != _event.lantern_pion_candidate_index) continue;
				// skip pion candidate, look at others
				//if (idx_track == _event.lantern_pion_candidate_index) continue;
				
				// primary only
				if (_event.lantern_trackIsSecondary[idx_track] != 0) continue;

				// skip short tracks
				float track_length = std::sqrt(
					std::pow(_event.lantern_trackEndPosX[idx_track] - _event.lantern_trackStartPosX[idx_track], 2) +
					std::pow(_event.lantern_trackEndPosY[idx_track] - _event.lantern_trackStartPosY[idx_track], 2) +
					std::pow(_event.lantern_trackEndPosZ[idx_track] - _event.lantern_trackStartPosZ[idx_track], 2) );
				
				//if (track_length < 1) continue;

				//if (_event.classification != Utility::kCCNumuNpi) continue;
				//if (_event.classification != Utility::kNC1pi) continue;
				//if (_event.lantern_trackTruePID[idx_track] != 211 && _event.lantern_trackTruePID[idx_track] != -211) continue;
				// skip protons
				//if (std::abs(_event.lantern_trackTruePID[idx_track]) == 2212) continue;

				if (_event.lantern_trackClassified[idx_track] == 0) continue; // only fill for classified tracks

				/*
				std::cout << "CC Numu Npi Event RSE: " << _event.run << "-" << _event.sub << "-" << _event.evt 
						  << " -- Lantern Track Index: " << idx_track 
						  << ", True PID: " << _event.lantern_trackTruePID[idx_track]
						  << ", Length: " << track_length
						  << ", Mu Score: " << std::exp(_event.lantern_trackMuScore[idx_track])
						  << ", Pi Score: " << std::exp(_event.lantern_trackPiScore[idx_track])
						  << ", Pr Score: " << std::exp(_event.lantern_trackPrScore[idx_track])
						  << std::endl;
				*/

				double mu_pi_llr_norm = std::tanh(0.5 * (_event.lantern_trackPiScore[idx_track] - _event.lantern_trackMuScore[idx_track]));
				double pr_pi_llr_norm = std::tanh(0.5 * (_event.lantern_trackPiScore[idx_track] - _event.lantern_trackPrScore[idx_track]));
				double pr_mu_llr_norm = std::tanh(0.5 * (_event.lantern_trackMuScore[idx_track] - _event.lantern_trackPrScore[idx_track]));	


				Utility::ParticleEnums particle_class = _event.lantern_particle_classification_v[idx_track];
				// group pion final state processes, testing
				/*
				if (particle_class == Utility::kPionAbsorption0p ||
					particle_class == Utility::kPionAbsorptionNp ||
					particle_class == Utility::kPionChargeExchange) {
					
						particle_class = Utility:: kPionDecay;
				}
				*/
				
				if (particle_class == Utility::kPionDecay 
					|| particle_class == Utility::kPionAbsorption0p ||
					particle_class == Utility::kPionAbsorptionNp ||
					particle_class == Utility::kPionChargeExchange
				) {
					pion_counter += pot_weight * _event.weight_cv;
					if (pion_counter > 1250) continue;
					if (mu_pi_llr_norm > 0.0) pion_counter_cut += pot_weight * _event.weight_cv;
					
					// track counters by-type
					if (particle_class == Utility::kPionDecay) n_track_decay += pot_weight * _event.weight_cv;
					else if (particle_class == Utility::kPionAbsorption0p) n_track_absorption0p += pot_weight * _event.weight_cv;
					else if (particle_class == Utility::kPionAbsorptionNp) n_track_absorptionNp += pot_weight * _event.weight_cv;
					else if (particle_class == Utility::kPionChargeExchange) n_track_chargeexchange += pot_weight * _event.weight_cv;
					// track counters by-type passing cut
					if (mu_pi_llr_norm > 0.4) {
						if (particle_class == Utility::kPionDecay) n_track_decay_cut += pot_weight * _event.weight_cv;
						else if (particle_class == Utility::kPionAbsorption0p) n_track_absorption0p_cut += pot_weight * _event.weight_cv;
						else if (particle_class == Utility::kPionAbsorptionNp) n_track_absorptionNp_cut += pot_weight * _event.weight_cv;
						else if (particle_class == Utility::kPionChargeExchange) n_track_chargeexchange_cut += pot_weight * _event.weight_cv;
					}


				}
				else if (particle_class == Utility::kMuon) {
					muon_counter += pot_weight * _event.weight_cv;
					if (muon_counter > 1250) continue;
					if (mu_pi_llr_norm > 0.0) muon_counter_cut += pot_weight * _event.weight_cv;

					// track counters by-type
					n_track_muon += pot_weight * _event.weight_cv;
					// track counters by-type passing cut
					if (mu_pi_llr_norm > 0.4) n_track_muon_cut += pot_weight * _event.weight_cv;	
				}
				else {
					continue; // only plot pions and muons (testing)
				}
				

				_histStack_lantern_larpid_pidScore_mu_particle.Fill(particle_class, std::exp(_event.lantern_trackMuScore[idx_track]), pot_weight * _event.weight_cv);
			    _histStack_lantern_larpid_pidScore_pi_particle.Fill(particle_class, (_event.lantern_trackPiScore[idx_track]), pot_weight * _event.weight_cv);
			    _histStack_lantern_larpid_pidScore_pr_particle.Fill(particle_class, std::exp(_event.lantern_trackPrScore[idx_track]), pot_weight * _event.weight_cv);

				//std::cout << "Original: " << _event.lantern_trackPiScore[idx_track] << ", Log: " << std::log(_event.lantern_trackPiScore[idx_track]) << ", Exp: " << std::exp(_event.lantern_trackPiScore[idx_track]) << std::endl;

				_histStack_lantern_larpid_pidScore_mu.Fill(_event.classification, std::exp(_event.lantern_trackMuScore[idx_track]), pot_weight * _event.weight_cv);
				_histStack_lantern_larpid_pidScore_pi.Fill(_event.classification, std::exp(_event.lantern_trackPiScore[idx_track]), pot_weight * _event.weight_cv);
				_histStack_lantern_larpid_pidScore_pr.Fill(_event.classification, std::exp(_event.lantern_trackPrScore[idx_track]), pot_weight * _event.weight_cv);
				
				// LLR
				//double mu_pi_llr = std::exp(_event.lantern_trackMuScore[idx_track] - _event.lantern_trackPiScore[idx_track]);
				//double pr_pi_llr = std::exp(_event.lantern_trackPrScore[idx_track] - _event.lantern_trackPiScore[idx_track]);
				_histStack_lantern_larpid_mupi_llr_particle.Fill(particle_class, mu_pi_llr_norm, pot_weight * _event.weight_cv);
				_histStack_lantern_larpid_prpi_llr_particle.Fill(particle_class, pr_pi_llr_norm, pot_weight * _event.weight_cv);
				_histStack_lantern_larpid_prmu_llr_particle.Fill(particle_class, pr_mu_llr_norm, pot_weight * _event.weight_cv);
			}
			// truth pion process for signal events
			if (_event.classification == Utility::kNC1pi) {

				// loop of truth particles, WC tree
				// loop over truth particles
				/*
				for (unsigned int wc_truth_idx = 0; wc_truth_idx < _event.wc_truth_Ntrack; wc_truth_idx++) {
					
					// primary particles
					if (_event.wc_truth_mother[wc_truth_idx] != 0) continue; 

					// find primary pion
					if (std::abs(_event.wc_truth_pdg[wc_truth_idx]) != 211) continue;

					// check momentum
					double momentum = std::sqrt( 	std::pow(_event.wc_truth_startMomentum[wc_truth_idx][0],2) + 
										std::pow(_event.wc_truth_startMomentum[wc_truth_idx][1],2) + 
										std::pow(_event.wc_truth_startMomentum[wc_truth_idx][2],2) );
					if (momentum < 0.1) continue;

					// have found primary pion
					std::cout << "Truth Selected Pion PDG " << _event.wc_truth_pdg[wc_truth_idx] << " -- Start Process: " << _event.wc_truth_process->at(wc_truth_idx) 
							  << ", End Process: " << _event.wc_truth_endprocess->at(wc_truth_idx) 
							  << std::endl;

					int pionID = _event.wc_truth_id[wc_truth_idx];
					
					PrintDaughters(_event, pionID);
					

					break;

				
				}
				*/

				/*
				// loop of truth particles, Pandora tree
				for (int i = 0; i < _event.mc_pdg_v->size(); i++) {
					
					// momentum
					double mc_momentum = std::sqrt( pow(_event.mc_px_v->at(i),2) + 
													pow(_event.mc_py_v->at(i),2) + 
													pow(_event.mc_pz_v->at(i),2) );

					
						if ((_event.mc_pdg_v->at(i) == 211 || _event.mc_pdg_v->at(i) == -211) && mc_momentum > 0.1) {
						std::cout << "Truth Selected Pion -- Start Process: " << _event.mc_process_v->at(i) 
								  << ", End Process: " << _event.mc_end_process_v->at(i) 
								  << std::endl;	
					}
				
				}
				*/
			
			}
		}

		


		// clean up
		f->Close();
		delete f;
	}

	// write training tree
	std::cout << "Writing BDT training tree." << std::endl;
	_trainingTree.writeOutputFile();

	std::cout << "Passing Data Events: " << n_data_pass << std::endl;

	std::cout << "Pion Candidate Count (Lantern): " << pion_counter << std::endl;
	std::cout << "Muon Candidate Count (Lantern): " << muon_counter << std::endl;

	std::cout << "Pion Candidate Count after LLR cut (Lantern): " << pion_counter_cut << std::endl;
	std::cout << "Muon Candidate Count after LLR cut (Lantern): " << muon_counter_cut << std::endl;

	std::cout << "Pion Track Counts by Type (Lantern): " << std::endl;
	std::cout << "  Decay: " << n_track_decay << " , After LLR cut: " << n_track_decay_cut << std::endl;
	std::cout << "  Absorption 0p: " << n_track_absorption0p << " , After LLR cut: " << n_track_absorption0p_cut << std::endl;
	std::cout << "  Absorption Np: " << n_track_absorptionNp << " , After LLR cut: " << n_track_absorptionNp_cut << std::endl;
	std::cout << "  Charge Exchange: " << n_track_chargeexchange << " , After LLR cut: " << n_track_chargeexchange_cut << std::endl;
	std::cout << "Muon Track Counts (Lantern): " << n_track_muon << " , After LLR cut: " << n_track_muon_cut << std::endl;

	// print event integrals
	_histStack_contained_fraction.PrintEventIntegrals();
	double signal_integral = _histStack_contained_fraction.GetSignalIntegral();
	std::cout << "Signal integral: " << signal_integral << std::endl;

	// print signal interaction modes
	std::cout << "Signal interaction modes:" << std::endl;
	_histStack_single_bin_interaction.PrintEventIntegrals();

	std::cout << "Signal All: " << n_signal_all << std::endl;
	std::cout << "Signal Pass: " << n_signal_pass << std::endl;
	std::cout << "Total Pass: " << n_total_pass << std::endl;

	std::cout << "Signal Pion Decay: " << n_signal_decay << ", Efficiency = " << (n_signal_pass_decay / n_signal_decay) * 100 << "%" << std::endl;
	std::cout << "Signal Pion Absorption 0p: " << n_signal_absorption0p << ", Efficiency = " << (n_signal_pass_absorption0p / n_signal_absorption0p) * 100 << "%" << std::endl;
	std::cout << "Signal Pion Absorption Np: " << n_signal_absorptionNp << ", Efficiency = " << (n_signal_pass_absorptionNp / n_signal_absorptionNp) * 100 << "%" << std::endl;
	std::cout << "Signal Pion Charge Exchange: " << n_signal_chargeexchange << ", Efficiency = " << (n_signal_pass_chargeexchange / n_signal_chargeexchange) * 100 << "%" << std::endl;

	std::cout << "Efficiency: " << (n_signal_pass / n_signal_all) * 100 << std::endl;
	std::cout << "Purity: " << (n_signal_pass / n_total_pass) * 100 << std::endl;
	
	// draw histograms
	TCanvas *c1 = new TCanvas("c1", "c1", 1080, 1080);
  	_histStack_contained_fraction.DrawStack(c1, Utility::kSingleBin);
  	c1->Print("plots/plot_contained_fraction.root");

	TCanvas *c1a = new TCanvas("c1a", "c1a", 1080, 1080);
  	_histStack_single_bin_interaction.DrawStack(c1a, Utility::kSingleBin);
  	c1a->Print("plots/single_bin_interaction.root");

  	TCanvas *c2 = new TCanvas("c2", "c2", 1080, 1080);
  	c2->cd();
  	_histStack_associated_hits_fraction.DrawStack(c2, Utility::kAssociatedHitsFraction);
  	c2->Print("plots/plot_associated_hits_fraction.root");

	TCanvas *c3 = new TCanvas("c3", "c3", 1080, 1080);
  	c3->cd();
  	_histStack_nu_vtx_x.DrawStack(c3, Utility::kNuVtxX);
  	c3->Print("plots/plot_nu_vtx_x.root");

	TCanvas *c3a = new TCanvas("c3a", "c3a", 1080, 1080);
  	c3a->cd();
  	_histStack_nu_vtx_y.DrawStack(c3a, Utility::kNuVtxY);
  	c3a->Print("plots/plot_nu_vtx_y.root");

	TCanvas *c3b = new TCanvas("c3b", "c3b", 1080, 1080);
  	c3b->cd();
  	_histStack_nu_vtx_z.DrawStack(c3b, Utility::kNuVtxZ);
  	c3b->Print("plots/plot_nu_vtx_z.root");

  	TCanvas *c5 = new TCanvas("c5", "c5", 1080, 1080);
  	c5->cd();
  	_histStack_topological_score.DrawStack(c5, Utility::kTopologicalScore);
  	c5->Print("plots/plot_topological_score.root");

  	TCanvas *c6 = new TCanvas("c6", "c6", 1080, 1080);
  	c6->cd();
  	_histStack_cosmic_IP.DrawStack(c6, Utility::kCosmicImpactParameter);
  	c6->Print("plots/plot_cosmic_IP.root");

  	TCanvas *c6b = new TCanvas("c6b", "c6b", 1080, 1080);
  	c6b->cd();
  	_histStack_number_tracks.DrawStack(c6b, Utility::kNTrack);
  	c6b->Print("plots/plot_number_tracks.root");

	TCanvas *c6c = new TCanvas("c6c", "c6c", 1080, 1080);
  	c6c->cd();
  	_histStack_number_showers.DrawStack(c6c, Utility::kNShower);
  	c6c->Print("plots/plot_number_showers.root");  	

  	TCanvas *c8 = new TCanvas("c8", "c8", 1080, 1080);
  	c8->cd();
  	_histStack_trk_length.DrawStack(c8, Utility::kTrackLength);
  	c8->Print("plots/plot_trk_length.root");

	TCanvas *c8aab = new TCanvas("c8aab", "c8aab", 1080, 1080);
  	c8aab->cd();
  	_histStack_trk2_length.DrawStack(c8aab, Utility::kTrackLength);
  	c8aab->Print("plots/plot_trk2_length.root");

  	TCanvas *c8aa = new TCanvas("c8aa", "c8aa", 1080, 1080);
  	c8aa->cd();
  	_histStack_trk_score.DrawStack(c8aa, Utility::kTrackScore);
  	c8aa->Print("plots/plot_trk_score.root");  	

  	TCanvas *c8b = new TCanvas("c8b", "c8b", 1080, 1080);
  	c8b->cd();
  	_histStack_trk_llrpid.DrawStack(c8b, Utility::kLLRPID);
  	c8b->Print("plots/plot_llrpid.root");
  	
  	TCanvas *c9 = new TCanvas("c9", "c9", 1080, 1080);
  	c9->cd();
  	_histStack_trk_vertex_distance.DrawStack(c9, Utility::kTrackDistance);
  	c9->Print("plots/plot_trk_vertex_distance.root");

	TCanvas *c10 = new TCanvas("c10", "c10", 1080, 1080);
  	c10->cd();
  	_histStack_n_non_proton_tracks.DrawStack(c10, Utility::kNTrack);
  	c10->Print("plots/plot_n_non_proton_tracks.root");

	TCanvas *c11 = new TCanvas("c11", "c11", 1080, 1080);
  	c11->cd();
  	_histStack_n_pion_candidate_tracks.DrawStack(c11, Utility::kNTrack);
  	c11->Print("plots/plot_n_pion_candidate_tracks.root");
	
	TCanvas *wc1 = new TCanvas("wc1", "wc1", 1080, 1080);
	_histStack_wc_numu_cc_flag.DrawStack(wc1, Utility::kWCNumuCCFlag);
  	wc1->Print("plots/plot_wc_numu_cc_flag.root");

	TCanvas *wc2 = new TCanvas("wc2", "wc2", 1080, 1080);
	_histStack_vertex_discrepancy.DrawStack(wc2, Utility::kTrackDistance);
  	wc2->Print("plots/plot_vertex_discrepancy.root");
	
	TCanvas *wc3 = new TCanvas("wc3", "wc3", 1080, 1080);
	_histStack_wc_match_isFC.DrawStack(wc3, Utility::kisFC);
  	wc3->Print("plots/plot_wc_match_isFC.root");

	TCanvas *wc4 = new TCanvas("wc4", "wc4", 1080, 1080);
	_histStack_wc_reco_larpid_pidScore_el.DrawStack(wc4, Utility::kLArPID_el);
  	wc4->Print("plots/plot_wc_reco_larpid_pidScore_el.root");

	TCanvas *wc5 = new TCanvas("wc5", "wc5", 1080, 1080);
	_histStack_wc_reco_larpid_pidScore_ph.DrawStack(wc5, Utility::kLArPID_ph);
  	wc5->Print("plots/plot_wc_reco_larpid_pidScore_ph.root");

	TCanvas *wc6 = new TCanvas("wc6", "wc6", 1080, 1080);
	_histStack_wc_reco_larpid_pidScore_mu.DrawStack(wc6, Utility::kLArPID_mu);
  	wc6->Print("plots/plot_wc_reco_larpid_pidScore_mu.root");

	TCanvas *wc7 = new TCanvas("wc7", "wc7", 1080, 1080);
	_histStack_wc_reco_larpid_pidScore_pi.DrawStack(wc7, Utility::kLArPID_pi);
  	wc7->Print("plots/plot_wc_reco_larpid_pidScore_pi.root");

	TCanvas *wc8 = new TCanvas("wc8", "wc8", 1080, 1080);
	_histStack_wc_reco_larpid_pidScore_pr.DrawStack(wc8, Utility::kLArPID_pr);
  	wc8->Print("plots/plot_wc_reco_larpid_pidScore_pr.root");

	TCanvas *wc9 = new TCanvas("wc9", "wc9", 1080, 1080);
	_histStack_wc_reco_larpid_procScore_prim.DrawStack(wc9, Utility::kLArPID_prim);
  	wc9->Print("plots/plot_wc_reco_larpid_procScore_prim.root");

	TCanvas *wc10 = new TCanvas("wc10", "wc10", 1080, 1080);
	_histStack_wc_reco_larpid_procScore_chgd.DrawStack(wc10, Utility::kLArPID_chgd);
  	wc10->Print("plots/plot_wc_reco_larpid_procScore_chgd.root");

	TCanvas *wc11 = new TCanvas("wc11", "wc11", 1080, 1080);
	_histStack_wc_reco_larpid_procScore_ntrl.DrawStack(wc11, Utility::kLArPID_ntrl);
  	wc11->Print("plots/plot_wc_reco_larpid_procScore_ntrl.root");

	TCanvas *wc12 = new TCanvas("wc12", "wc12", 1080, 1080);
	_histStack_wc_bdt_score.DrawStack(wc12, Utility::kPionProtonBDT);
  	wc12->Print("plots/plot_wc_bdt_score.root");

	// lantern
	TCanvas *l1 = new TCanvas("l1", "l1", 1080, 1080);
	_histStack_lantern_vtxScore.DrawStack(l1, Utility::kLanternVtxScore);
  	l1->Print("plots/plot_lantern_vtxScore.root");

	TCanvas *l2 = new TCanvas("l2", "l2", 1080, 1080);
	_histStack_lantern_vtxContainment.DrawStack(l2, Utility::kLanternVtxContainment);
  	l2->Print("plots/plot_lantern_vtxContainment.root");

	TCanvas *l3 = new TCanvas("l3", "l3", 1080, 1080);
	_histStack_lantern_vtxFracHitsOnCosmic.DrawStack(l3, Utility::kLanternVtxFracHitsOnCosmic);
  	l3->Print("plots/plot_lantern_vtxFracHitsOnCosmic.root");

	// lantern per particle
	TCanvas *l4 = new TCanvas("l4", "l4", 1080, 1080);
	_histStack_lantern_larpid_pidScore_mu_particle.DrawStack(l4, Utility::kLArPID_mu);
  	l4->Print("plots/plot_lantern_larpid_pidScore_mu_particle.root");

	TCanvas *l5 = new TCanvas("l5", "l5", 1080, 1080);
	_histStack_lantern_larpid_pidScore_pi_particle.DrawStack(l5, Utility::kLArPID_pi);
  	l5->Print("plots/plot_lantern_larpid_pidScore_pi_particle.root");

	TCanvas *l6 = new TCanvas("l6", "l6", 1080, 1080);
	_histStack_lantern_larpid_pidScore_pr_particle.DrawStack(l6, Utility::kLArPID_pr);
  	l6->Print("plots/plot_lantern_larpid_pidScore_pr_particle.root");

	
	TCanvas *l6a = new TCanvas("l6a", "l6a", 1080, 1080);
	_histStack_lantern_larpid_mupi_llr_particle.DrawStack(l6a, Utility::kLArPID_pi);
  	l6a->Print("plots/plot_lantern_larpid_mupi_llr_particle.root");
	TCanvas *l6b = new TCanvas("l6b", "l6b", 1080, 1080);
	_histStack_lantern_larpid_prpi_llr_particle.DrawStack(l6b, Utility::kLArPID_pi);
  	l6b->Print("plots/plot_lantern_larpid_prpi_llr_particle.root");
	TCanvas *l6c = new TCanvas("l6c", "l6c", 1080, 1080);
	_histStack_lantern_larpid_prmu_llr_particle.DrawStack(l6c, Utility::kLArPID_pr);
  	l6c->Print("plots/plot_lantern_larpid_prmu_llr_particle.root");

	// lantern candiate particle topology
	TCanvas *l7 = new TCanvas("l7", "l7", 1080, 1080);
	_histStack_lantern_larpid_pidScore_mu.DrawStack(l7, Utility::kLArPID_mu);
  	l7->Print("plots/plot_lantern_larpid_pidScore_mu.root");

	TCanvas *l8 = new TCanvas("l8", "l8", 1080, 1080);
	_histStack_lantern_larpid_pidScore_pi.DrawStack(l8, Utility::kLArPID_pi);
  	l8->Print("plots/plot_lantern_larpid_pidScore_pi.root");

	TCanvas *l9 = new TCanvas("l9", "l9", 1080, 1080);
	_histStack_lantern_larpid_pidScore_pr.DrawStack(l9, Utility::kLArPID_pr);
  	l9->Print("plots/plot_lantern_larpid_pidScore_pr.root");
	

	// per particle
	TCanvas *p1 = new TCanvas("p1", "p1", 1080, 1080);
  	p1->cd();
  	_histStack_trk_length_particle.DrawStack(p1, Utility::kTrackLength);
  	p1->Print("plots/plot_trk_length_particle.root");
	
	TCanvas *p2 = new TCanvas("p2", "p2", 1080, 1080);
  	p2->cd();
  	_histStack_trk_score_particle.DrawStack(p2, Utility::kTrackScore);
  	p2->Print("plots/plot_trk_score_particle.root");
	TCanvas *p3 = new TCanvas("p3", "p3", 1080, 1080);
  	p3->cd();
  	_histStack_trk_llrpid_particle.DrawStack(p3, Utility::kLLRPID);
  	p3->Print("plots/plot_llrpid_particle.root");
	TCanvas *p4 = new TCanvas("p4", "p4", 1080, 1080);
  	p4->cd();
  	_histStack_trk_bragg_mu_particle.DrawStack(p4, Utility::kTrackBraggMu);
  	p4->Print("plots/plot_trk_bragg_mu_particle.root");
	TCanvas *p5 = new TCanvas("p5", "p5", 1080, 1080);
  	p5->cd();
  	_histStack_trk_bragg_p_particle.DrawStack(p5, Utility::kTrackBraggP);
  	p5->Print("plots/plot_trk_bragg_p_particle.root");
	TCanvas *p6 = new TCanvas("p6", "p6", 1080, 1080);
  	p6->cd();
  	_histStack_trk_bragg_pion_particle.DrawStack(p6, Utility::kTrackBraggPion);
  	p6->Print("plots/plot_trk_bragg_pion_particle.root");
	TCanvas *p7 = new TCanvas("p7", "p7", 1080, 1080);
  	p7->cd();
  	_histStack_trk_bragg_mip_particle.DrawStack(p7, Utility::kTrackBraggMIP);
  	p7->Print("plots/plot_trk_bragg_mip_particle.root");
	TCanvas *p8 = new TCanvas("p8", "p8", 1080, 1080);
  	p8->cd();
  	_histStack_trk_dEdx_trunk_particle.DrawStack(p8, Utility::kTrackdEdx);
  	p8->Print("plots/plot_trk_dEdx_trunk_particle.root");
	TCanvas *p9 = new TCanvas("p9", "p9", 1080, 1080);
  	p9->cd();
  	_histStack_trk_daughters_particle.DrawStack(p9, Utility::kNTrackDaughters);
  	p9->Print("plots/plot_trk_daughters_particle.root");
	TCanvas *p10 = new TCanvas("p10", "p10", 1080, 1080);
  	p10->cd();
  	_histStack_trk_end_spacepoints_particle.DrawStack(p10, Utility::kTrackEndSpacepoints);
  	p10->Print("plots/plot_trk_end_spacepoints_particle.root");
	TCanvas *p11 = new TCanvas("p11", "p11", 1080, 1080);
  	p11->cd();
  	_histStack_trk_avg_deflection_stdev_particle.DrawStack(p11, Utility::kTrackWiggliness);
  	p11->Print("plots/plot_trk_avg_deflection_stdev_particle.root");
	TCanvas *p12 = new TCanvas("p12", "p12", 1080, 1080);
  	p12->cd();
  	_histStack_pfng2mipfrac_particle.DrawStack(p12, Utility::kNuGraphMipFrac);
  	p12->Print("plots/plot_pfng2mipfrac_particle.root");
	TCanvas *p13 = new TCanvas("p13", "p13", 1080, 1080);
  	p13->cd();
  	_histStack_pfng2hipfrac_particle.DrawStack(p13, Utility::kNuGraphHipFrac);
  	p13->Print("plots/plot_pfng2hipfrac_particle.root");
	TCanvas *p14 = new TCanvas("p14", "p14", 1080, 1080);
  	p14->cd();
  	_histStack_pfng2shrfrac_particle.DrawStack(p14, Utility::kNuGraphShrFrac);
  	p14->Print("plots/plot_pfng2shrfrac_particle.root");
	TCanvas *p15 = new TCanvas("p15", "p15", 1080, 1080);
  	p15->cd();
  	_histStack_pfng2mclfrac_particle.DrawStack(p15, Utility::kNuGraphMclFrac);
  	p15->Print("plots/plot_pfng2mclfrac_particle.root");
	TCanvas *p16 = new TCanvas("p16", "p16", 1080, 1080);
  	p16->cd();
  	_histStack_pfng2dfsfrac_particle.DrawStack(p16, Utility::kNuGraphDfsFrac);
  	p16->Print("plots/plot_pfng2dfsfrac_particle.root");

	// WC LArPID per particle
	TCanvas *wc_p1 = new TCanvas("wc_p1", "wc_p1", 1080, 1080);
  	wc_p1->cd();
  	_histStack_wc_reco_larpid_pidScore_el_particle.DrawStack(wc_p1, Utility::kLArPID_el);
  	wc_p1->Print("plots/plot_wc_reco_larpid_pidScore_el_particle.root");
	TCanvas *wc_p2 = new TCanvas("wc_p2", "wc_p2", 1080, 1080);
  	wc_p2->cd();
  	_histStack_wc_reco_larpid_pidScore_ph_particle.DrawStack(wc_p2, Utility::kLArPID_ph);
  	wc_p2->Print("plots/plot_wc_reco_larpid_pidScore_ph_particle.root");
	TCanvas *wc_p3 = new TCanvas("wc_p3", "wc_p3", 1080, 1080);
  	wc_p3->cd();
  	_histStack_wc_reco_larpid_pidScore_mu_particle.DrawStack(wc_p3, Utility::kLArPID_mu);
  	wc_p3->Print("plots/plot_wc_reco_larpid_pidScore_mu_particle.root");
	TCanvas *wc_p4 = new TCanvas("wc_p4", "wc_p4", 1080, 1080);
  	wc_p4->cd();
  	_histStack_wc_reco_larpid_pidScore_pi_particle.DrawStack(wc_p4, Utility::kLArPID_pi);
  	wc_p4->Print("plots/plot_wc_reco_larpid_pidScore_pi_particle.root");
	TCanvas *wc_p5 = new TCanvas("wc_p5", "wc_p5", 1080, 1080);
  	wc_p5->cd();
  	_histStack_wc_reco_larpid_pidScore_pr_particle.DrawStack(wc_p5, Utility::kLArPID_pr);
  	wc_p5->Print("plots/plot_wc_reco_larpid_pidScore_pr_particle.root");
	TCanvas *wc_p6 = new TCanvas("wc_p6", "wc_p6", 1080, 1080);
  	wc_p6->cd();
	_histStack_wc_reco_larpid_pdg_particle.DrawStack(wc_p6, Utility::kLArPID_pdg);
  	wc_p6->Print("plots/plot_wc_reco_larpid_pdg_particle.root");
	TCanvas *wc_p7 = new TCanvas("wc_p7", "wc_p7", 1080, 1080);
  	wc_p7->cd();
	_histStack_wc_n_pion_candidate_daughters_particle.DrawStack(wc_p7, Utility::kNTrackDaughters);
  	wc_p7->Print("plots/plot_wc_n_pion_candidate_daughters_particle.root");
	TCanvas *wc_p8 = new TCanvas("wc_p8", "wc_p8", 1080, 1080);
  	wc_p8->cd();
	_histStack_wc_pion_candidate_daughters_total_energy_particle.DrawStack(wc_p8, Utility::kDaughterEnergy);
  	wc_p8->Print("plots/plot_wc_pion_candidate_daughters_total_energy_particle.root");
	TCanvas *wc_p9 = new TCanvas("wc_p9", "wc_p9", 1080, 1080);
  	wc_p9->cd();
	_histStack_wc_n_blips_3cm.DrawStack(wc_p9, Utility::kNBlips25cm);
  	wc_p9->Print("plots/plot_wc_n_blips_25cm.root");
	TCanvas *wc_p10 = new TCanvas("wc_p10", "wc_p10", 1080, 1080);
  	wc_p10	->cd();
	_histStack_wc_n_blips_5cm.DrawStack(wc_p10, Utility::kNBlips50cm);
  	wc_p10->Print("plots/plot_wc_n_blips_50cm.root");
	TCanvas *wc_p11 = new TCanvas("wc_p11", "wc_p11", 1080, 1080);
  	wc_p11->cd();
	_histStack_wc_n_blips_10cm.DrawStack(wc_p11, Utility::kNBlips10cm);
  	wc_p11->Print("plots/plot_wc_n_blips_100cm.root");	

}