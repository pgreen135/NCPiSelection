#include "../include/SelectionDriver.h"

#include <iostream>
#include <iomanip>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TLine.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include "../include/EventContainer.h"
#include "../include/CreateTrainingTree.h"
#include "../include/Selection.h"
#include "../include/StackedHistTool.h"

#include "../include/PhysdEdx.h"


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
	
	// Beam On Only
	//std::vector<std::string> filename_list= {filename_beamon_run4};
	//std::vector<double> pot_weight_list = {1.0};
	//std::vector<Utility::FileTypeEnums> file_types_list = {Utility::kData};
	//std::vector<Utility::RunPeriodEnums> run_periods_list = {Utility::kRun4ab};

	// NuWro Only
	/*
	std::vector<std::string> filename_list= {filename_nuwro_run5};
	std::vector<double> pot_weight_list = {1.0};
	std::vector<Utility::FileTypeEnums> file_types_list = {Utility::kFakeData};
    std::vector<Utility::RunPeriodEnums> run_periods_list = {Utility::kRun4ab};
	*/

	// SURPRISE MC samples
	std::vector<std::string> filename_list =  {//filename_mc_run1, //filename_dirt_run1, filename_beamoff_run1,
											  //filename_mc_run2, //filename_dirt_run2, filename_beamoff_run2,
											  //filename_mc_run3, //filename_dirt_run3, filename_beamoff_run3,
											  //filename_mc_run4, //filename_dirt_run4, filename_beamoff_run4,
											  filename_mc_run5, filename_dirt_run5, filename_beamoff_run5
											  };
	std::vector<double> pot_weight_list = { //pot_weight_mc_run1,  //pot_weight_dirt_run1, pot_weight_beamoff_run1,
											//pot_weight_mc_run2,  //pot_weight_dirt_run2, pot_weight_beamoff_run2,
											//pot_weight_mc_run3,  pot_weight_dirt_run3, pot_weight_beamoff_run3,
		                                  	//pot_weight_mc_run4, pot_weight_dirt_run4, pot_weight_beamoff_run4,
										    pot_weight_mc_run5, pot_weight_dirt_run5, pot_weight_beamoff_run5
										   };
	std::vector<Utility::FileTypeEnums> file_types_list = {//Utility::kMC, //Utility::kDirt, Utility::kEXT,
		 												   //Utility::kMC, //Utility::kDirt, Utility::kEXT,
												 		   //Utility::kMC, //Utility::kDirt, Utility::kEXT,
												           //Utility::kMC, //Utility::kDirt, Utility::kEXT,
												           Utility::kMC, Utility::kDirt, Utility::kEXT
													    };
	std::vector<Utility::RunPeriodEnums> run_periods_list = {Utility::kRun4ab, //Utility::kRun4ab, Utility::kRun4ab,
															 Utility::kRun4ab, //Utility::kRun4ab, Utility::kRun4ab,
															 Utility::kRun4ab, //Utility::kRun4ab, Utility::kRun4ab,
															 Utility::kRun4ab, //Utility::kRun4ab, Utility::kRun4ab,
															 Utility::kRun4ab, //Utility::kRun4ab, Utility::kRun4ab,
															};
	
	// construct classes
	Selection _selection(_utility);
	CreateTrainingTree _trainingTree(Utility::kMuonPion);

	// construct stacks
	StackedHistTool _histStack_contained_fraction("", "", 3, -1.05, 2.10, _utility);
	StackedHistTool _histStack_single_bin_interaction("", "", 3, -1.05, 2.10, _utility, "Interaction");
	StackedHistTool _histStack_associated_hits_fraction("", "", 20, 0, 1.001, _utility);
	StackedHistTool _histStack_nu_vtx_x("", "", 100, -50, 50, _utility);
	StackedHistTool _histStack_nu_vtx_y("", "", 100, -50, 50, _utility);
	StackedHistTool _histStack_nu_vtx_z("", "", 100, -50, 50, _utility);
	StackedHistTool _histStack_trk_length("", "", 60, 0, 300, _utility);

	StackedHistTool _histStack_wc_numu_cc_flag("", "", 3, -1.05, 1.05, _utility);
	StackedHistTool _histStack_vertex_discrepancy("", "", 50, 0, 50, _utility);
	StackedHistTool _histStack_wc_match_isFC("", "", 2, -0.05, 1.05, _utility);
	
	// Other WC
	StackedHistTool _histStack_wc_n_pion_candidate_daughters_particle("", "", 6, 0, 6, _utility, "Particle");
	StackedHistTool _histStack_wc_pion_candidate_daughters_total_energy_particle("", "", 25, 0, 0.5, _utility, "Particle");

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

	// WC Lantern Candidate Matching
	StackedHistTool _histStack_wc_lantern_candidate_start_distance("", "", 50, 0, 100, _utility);
	StackedHistTool _histStack_wc_lantern_candidate_end_distance("", "", 50, 0, 100, _utility);
	StackedHistTool _histStack_wc_lantern_candidate_angle_difference("", "", 40, 0, 180, _utility);

	// Observables
	StackedHistTool _histStack_sel_pion_theta("", "", 20, -1.0, 1.0, _utility);
	StackedHistTool _histStack_sel_pion_phi("", "", 18, -180, 180, _utility);
	StackedHistTool _histStack_sel_pion_momentum("", "", 40, 0.0, 2.0, _utility);
	StackedHistTool _histStack_sel_pion_momentum_interaction("", "", 20, 0.0, 1.0, _utility,  "Interaction");

	StackedHistTool _histStack_sel_proton_theta("", "", 20, -1.0, 1.0, _utility);													
	StackedHistTool _histStack_sel_proton_momentum("", "", 40, 0, 2, _utility);

	StackedHistTool _histStack_sel_pion_proton_opening_angle("", "", 20, 0, 180, _utility);

	StackedHistTool _histStack_sel_sideband_larpidScore("", "", 20, -1, 1.001, _utility, "Particle");

	// resolution plots
	TH1F *pionMomentum1D = new TH1F("", "", 50, -1.0, 1.0);														
	TH2F *pionMomentum2D = new TH2F("", "", 20, 0, 2.0, 20, 0, 2.0);

	TH1F *pionTheta1D = new TH1F("", "", 50, -1.0, 1.0);	
	TH2F *pionTheta2D = new TH2F("", "", 20, -1.0, 1.0, 20, -1.0, 1.0);

	// pion-proton opening angle resolution plots
	TH1F *pionProtonOpeningAngle1D = new TH1F("", "", 50, -1.0, 1.0);	
	TH2F *pionProtonOpeningAngle2D = new TH2F("", "", 20, 0, 180, 20, 0, 180);

	// dE/dx vs residual range histogram
	TH2F *dEdx_vs_ResRange = new TH2F("", "", 100, 0, 200, 100, 0, 10);
	PhysdEdx _physdedx(13); // == muon

	// create output root file for signal momentum and angle
	TFile outputFile("binningTree.root", "RECREATE");
	TTree binningTree("Events", "Example event tree");
	double sel_pion_momentum = 0;
	double sel_pion_theta = 0;
	double sel_pion_phi = 0;
	binningTree.Branch("sel_pion_momentum", &sel_pion_momentum, "sel_pion_momentum/D");
	binningTree.Branch("sel_pion_theta", &sel_pion_theta, "sel_pion_theta/D");
	binningTree.Branch("sel_pion_phi", &sel_pion_phi, "sel_pion_phi/D");
	double true_pion_momentum = 0;
	double true_pion_theta = 0;
	double true_pion_phi = 0;
	binningTree.Branch("true_pion_momentum", &true_pion_momentum, "true_pion_momentum/D");
	binningTree.Branch("true_pion_theta", &true_pion_theta, "true_pion_theta/D");
	binningTree.Branch("true_pion_phi", &true_pion_phi, "true_pion_phi/D");
	double binning_weight = 0;
	binningTree.Branch("weight", &binning_weight, "weight/D");													

	// Counters
	int n_data_pass = 0;
	double n_signal_all = 0;
	double n_signal_pass = 0;
	double n_total_pass = 0;

	double n_signal_0p = 0;
	double n_signal_Np = 0;
	double n_signal_0p_pass = 0;
	double n_signal_Np_pass = 0;

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
		TTree *wc_spacepoints_tree = (TTree*)f->Get("wcpselection/T_spacepoints");
		TTree *lantern_tree = (TTree*)f->Get("lantern/EventTree");

	 	// initialise event container
	  	EventContainer _event(pelee_tree, wc_eval_tree, wc_BDTvars_tree, wc_PFeval_tree, wc_spacepoints_tree, lantern_tree, _utility);

	  	// loop through events
	  	int n_entries = pelee_tree->GetEntries();
	  	std::cout << "Initial number events: " << n_entries << std::endl;

	  	for (int e = 0; e < n_entries; e++) {
		//for (int e = 0; e < 10000; e++) {

			// skip events with problems in test sample
			// need to load PeLEE tree first to get RSE
			pelee_tree->GetEntry(e);

			// Hypfit Testing -- Pick single track
			//if (!(_event.run == 20342 && _event.sub == 83 && _event.evt == 4177)) continue;
			//if (!(_event.run == 19925 && _event.sub == 535 && _event.evt == 26768)) continue;
			//std::cout << "Hypfit Debugging Mode" << std::endl;

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
			wc_spacepoints_tree->GetEntry(e);
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
			//bool passSelection = _selection.ApplyWCSelection(_event, file_types_list[idx], run_periods_list[idx]);

			// Lantern selection
			bool passSelection = _selection.ApplyLanternSelection(_event, file_types_list[idx], run_periods_list[idx]);
			
		    // evaluate correct POT weight accounting for special cases
		    double pot_weight = pot_weight_list[idx];
			

			// increment all signal counter
			if (_event.classification == Utility::kNC1pi) {
				n_signal_all += pot_weight * _event.weight_cv;

				if (_event.mc_NProtons_ == 0) {
					n_signal_0p += pot_weight * _event.weight_cv;
				}
				else if (_event.mc_NProtons_ > 0) {
					n_signal_Np += pot_weight * _event.weight_cv;
				}
				
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
				else {
					std::cout << "Unclassified Pion!" << std::endl;
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
			//if (!_event.sel_CC0piFarSideband_) continue;
			//if (!_event.sel_CC0piNearSideband_) continue;
			//if (!_event.sel_CC1piSideband_) continue;
			//if (!_event.sel_NCNpSideband_) continue;
		   
			
		    if (file_types_list[idx] == Utility::kData) {
		    	n_data_pass++;
		    	std::cout << "Passing data event -- " << "Run: " << _event.run << ", Subrun: " << _event.sub << ", Event: " << _event.evt << ", Electron energy: " << _event.shr_energy_cali << ", Visible energy: " << _event.NeutrinoEnergy2/1000 << " GeV" << std::endl;
		    	continue;
		    }
			

			// increment passing signal counter
			if (_event.classification == Utility::kNC1pi) {
				n_signal_pass += pot_weight * _event.weight_cv;

				if (_event.mc_NProtons_ == 0) {
					n_signal_0p_pass += pot_weight * _event.weight_cv;
				}
				else if (_event.mc_NProtons_ > 0) {
					n_signal_Np_pass += pot_weight * _event.weight_cv;
				}

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

				// Print out RSE:
				//std::cout << "Signal Event Run: " << _event.run << ", Subrun: " << _event.sub << ", Event: " << _event.evt;

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

			// fill observables for passing events
			if (file_types_list[idx] != Utility::kData) {
			
			_histStack_sel_pion_theta.Fill(_event.classification, _event.sel_pion_theta_, pot_weight * _event.weight_cv);
			_histStack_sel_proton_theta.Fill(_event.classification, _event.sel_proton_theta_, pot_weight * _event.weight_cv);

			_histStack_sel_pion_phi.Fill(_event.classification, (_event.sel_pion_phi_)*180/3.1416, pot_weight * _event.weight_cv);
			
			// [testing] truth momentum
			/*
			if (_event.classification == Utility::kNC1pi) {
			//if (_event.classification == Utility::kCCNumu0pi) {

				// find true muon, matching signal thresholds (pion with momentum > 100 MeV)
				int truth_muon_index = -1;
				double truth_muon_momentum = -1;
				for (int i = 0; i < _event.mc_pdg_v->size(); i++) {
					
					// momentum
					double mc_momentum = std::sqrt( pow(_event.mc_px_v->at(i),2) + 
													pow(_event.mc_py_v->at(i),2) + 
													pow(_event.mc_pz_v->at(i),2) );

					// pion
					//if ((_event.mc_pdg_v->at(i) == 13 || _event.mc_pdg_v->at(i) == -13) && mc_momentum > 0.1) {
					if ((_event.mc_pdg_v->at(i) == 211 || _event.mc_pdg_v->at(i) == -211) && mc_momentum > 0.1) {
						truth_muon_index = i;
						truth_muon_momentum = mc_momentum;
						break;
					}						
				}

				_histStack_sel_pion_momentum.Fill(_event.classification, truth_muon_momentum, pot_weight * _event.weight_cv);
			}
			*/
			
			
			_histStack_sel_pion_momentum.Fill(_event.classification, _event.sel_pion_momentum_, pot_weight * _event.weight_cv);
			//_histStack_sel_pion_momentum.Fill(_event.classification, _event.mc_pion_momentum_, pot_weight * _event.weight_cv);
			_histStack_sel_pion_momentum_interaction.Fill(_event.interactionMode, _event.sel_pion_momentum_, pot_weight * _event.weight_cv);

			// proton momentum
			_histStack_sel_proton_momentum.Fill(_event.classification, _event.sel_proton_momentum_, pot_weight * _event.weight_cv);

			// pion-proton opening angle
			_histStack_sel_pion_proton_opening_angle.Fill(_event.classification, _event.sel_pion_proton_opening_angle_, pot_weight * _event.weight_cv);

			_histStack_sel_sideband_larpidScore.Fill(_event.sel_CC1piSideband_track_classification_, _event.sel_CC1piSideband_track_larpidScore_, pot_weight * _event.weight_cv);
			

		    // fill histogram(s) - pelee variables
		    _histStack_contained_fraction.Fill(_event.classification, _event.contained_fraction, pot_weight * _event.weight_cv);
			_histStack_associated_hits_fraction.Fill(_event.classification, _event.associated_hits_fraction, pot_weight * _event.weight_cv);
			//_histStack_nu_vtx_x.Fill(_event.classification, _event.true_nu_vtx_sce_x - _event.lantern_trueVtxX, pot_weight * _event.weight_cv);
			_histStack_nu_vtx_x.Fill(_event.classification, _event.wc_reco_nuvtxX , pot_weight * _event.weight_cv);
			_histStack_nu_vtx_y.Fill(_event.classification, _event.wc_reco_nuvtxY, pot_weight * _event.weight_cv);
			_histStack_nu_vtx_z.Fill(_event.classification, _event.wc_reco_nuvtxZ, pot_weight * _event.weight_cv);

			if (_event.classification == Utility::kNC1pi) {
				_histStack_single_bin_interaction.Fill(_event.interactionMode, _event.contained_fraction, pot_weight * _event.weight_cv);
			}

			// fill histogram(s) - WC variables
			_histStack_wc_numu_cc_flag.Fill(_event.classification, _event.wc_numu_cc_flag, pot_weight * _event.weight_cv);
			_histStack_vertex_discrepancy.Fill(_event.classification, _event.vertex_discrepancy, pot_weight * _event.weight_cv);
			_histStack_wc_match_isFC.Fill(_event.classification, _event.wc_match_isFC, pot_weight * _event.weight_cv);

			// lantern
			_histStack_lantern_vtxScore.Fill(_event.classification, _event.lantern_vtxScore, pot_weight * _event.weight_cv);
			_histStack_lantern_vtxContainment.Fill(_event.classification, _event.lantern_vtxContainment, pot_weight * _event.weight_cv);
			_histStack_lantern_vtxFracHitsOnCosmic.Fill(_event.classification, _event.lantern_vtxFracHitsOnCosmic, pot_weight * _event.weight_cv);

			// fill histogram(s) - WC per particle
			for (unsigned int wc_idx = 0; wc_idx < _event.wc_reco_Ntrack; wc_idx++) {

				// only fill pion candidate
				if (wc_idx != _event.wc_pion_candidate_index) continue;

				// only fill non-pion candidate tracks
				//if (wc_idx == _event.wc_pion_candidate_index) continue;

				// fill only longest track
				//if (wc_idx != _event.wc_primaryTrackIndex) continue;

				// check primary
				if (_event.wc_reco_mother[wc_idx] != 0) continue; // require WC primary track

				// skip WC pseudo-particles
				if (_event.wc_reco_pdg[wc_idx] == 22 || _event.wc_reco_pdg[wc_idx] == 2112) continue; // skip pseudo-particles

				// fill histograms
				
				// Other WC
				_histStack_wc_n_pion_candidate_daughters_particle.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_n_pion_candidate_daughters, pot_weight * _event.weight_cv);
				if (_event.wc_n_pion_candidate_daughters > 0) {
					_histStack_wc_pion_candidate_daughters_total_energy_particle.Fill(_event.wc_particle_classification_v[wc_idx], _event.wc_pion_candidate_daughters_total_energy, pot_weight * _event.weight_cv);
				}
			}

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
				
				//if (track_length < 2) continue;

				//if (_event.classification != Utility::kCCNumuNpi) continue;
				//if (_event.classification != Utility::kNC1pi) continue;
				//if (_event.lantern_trackTruePID[idx_track] != 211 && _event.lantern_trackTruePID[idx_track] != -211) continue;
				// skip protons
				//if (std::abs(_event.lantern_trackTruePID[idx_track]) == 2212) continue;

				if (_event.lantern_trackClassified[idx_track] == 0) continue; // only fill for classified tracks

				double mu_pi_llr_norm = std::tanh(0.5 * (_event.lantern_trackPiScore[idx_track] - _event.lantern_trackMuScore[idx_track]));
				double pr_pi_llr_norm = std::tanh(0.5 * (_event.lantern_trackPiScore[idx_track] - _event.lantern_trackPrScore[idx_track]));
				double pr_mu_llr_norm = std::tanh(0.5 * (_event.lantern_trackMuScore[idx_track] - _event.lantern_trackPrScore[idx_track]));	


				Utility::ParticleEnums particle_class = _event.lantern_particle_classification_v[idx_track];
				// group pion final state processes, testing
	
				if (particle_class == Utility::kPionAbsorption0p ||
					particle_class == Utility::kPionAbsorptionNp ||
					particle_class == Utility::kPionChargeExchange) {
					
						particle_class = Utility:: kPionDecay;
				}
				
				/*
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
				*/

				_histStack_lantern_larpid_pidScore_mu_particle.Fill(particle_class, std::exp(_event.lantern_trackMuScore[idx_track]), pot_weight * _event.weight_cv);
			    _histStack_lantern_larpid_pidScore_pi_particle.Fill(particle_class, std::exp(_event.lantern_trackPiScore[idx_track]), pot_weight * _event.weight_cv);
			    _histStack_lantern_larpid_pidScore_pr_particle.Fill(particle_class, std::exp(_event.lantern_trackPrScore[idx_track]), pot_weight * _event.weight_cv);

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

				double truth_pion_momentum = -1;
				for (int i = 0; i < _event.mc_pdg_v->size(); i++) {
					
					// momentum
					double mc_momentum = std::sqrt( pow(_event.mc_px_v->at(i),2) + 
													pow(_event.mc_py_v->at(i),2) + 
													pow(_event.mc_pz_v->at(i),2) );

					// pion
					//if ((_event.mc_pdg_v->at(i) == 13 || _event.mc_pdg_v->at(i) == -13)) std::cout << "Truth Muon Momentum: " << mc_momentum << std::endl;
					if ((_event.mc_pdg_v->at(i) == 211 || _event.mc_pdg_v->at(i) == -211) && mc_momentum > 0.1) {
						//std::cout << "Truth Pion Momentum: " << mc_momentum << std::endl;
						truth_pion_momentum = mc_momentum;
					}						
				}

				//std::cout << "Truth Pion Momentum: " << truth_pion_momentum << ", Reco Pion Momentum (Hypfit): " << _event.pionMomentumHypfit << ", Reco Pion Momentum (Range): " << _event.pionMomentumRange << std::endl;

				// check hypfit, populate histograms if reasonable (for testing)
				if (_utility.isNumber(_event.pionMomentumHypfit) && _event.pionMomentumHypfit > 0 && _event.pionMomentumHypfit < 2) {
					pionMomentum1D->Fill((_event.pionMomentumHypfit - truth_pion_momentum) / truth_pion_momentum);
					pionMomentum2D->Fill(truth_pion_momentum, _event.pionMomentumHypfit);
				}
				else {
					pionMomentum1D->Fill((_event.pionMomentumRange - truth_pion_momentum) / truth_pion_momentum);
					pionMomentum2D->Fill(truth_pion_momentum, _event.pionMomentumRange);
				}
				

				// fill resolution plots
				pionTheta1D->Fill((_event.sel_pion_theta_ - _event.mc_pion_theta_) / _event.mc_pion_theta_);
				pionTheta2D->Fill(_event.mc_pion_theta_, _event.sel_pion_theta_, pot_weight * _event.weight_cv);

				pionProtonOpeningAngle1D->Fill((_event.sel_pion_proton_opening_angle_ - _event.mc_pion_proton_opening_angle_) / _event.mc_pion_proton_opening_angle_);
				pionProtonOpeningAngle2D->Fill(_event.mc_pion_proton_opening_angle_, _event.sel_pion_proton_opening_angle_, pot_weight * _event.weight_cv);

				// fill binning tree
				sel_pion_theta = _event.sel_pion_theta_;
				sel_pion_phi = _event.sel_pion_phi_;
				sel_pion_momentum = _event.sel_pion_momentum_;
				true_pion_theta = _event.mc_pion_theta_;
				true_pion_phi = _event.mc_pion_phi_;
				true_pion_momentum = _event.mc_pion_momentum_;
				binning_weight = pot_weight * _event.weight_cv;
				

				//binningTree.Fill();

			}

			// populate WC/Lantern candidate consistency histograms
			// only when both populated 
			if (_event.wc_pion_candidate_index != -1 && _event.lantern_pion_candidate_index != -1) {

				// calculate start distance between WC and Lantern candidates
				float wc_candidate_start[3] = {_event.wc_reco_startXYZT[_event.wc_pion_candidate_index][0], _event.wc_reco_startXYZT[_event.wc_pion_candidate_index][1], _event.wc_reco_startXYZT[_event.wc_pion_candidate_index][2]};
				float lantern_candidate_start[3] = {_event.lantern_trackStartPosX[_event.lantern_pion_candidate_index], _event.lantern_trackStartPosY[_event.lantern_pion_candidate_index], _event.lantern_trackStartPosZ[_event.lantern_pion_candidate_index]};
				float start_distance = std::sqrt(
					std::pow(wc_candidate_start[0] - lantern_candidate_start[0], 2) +
					std::pow(wc_candidate_start[1] - lantern_candidate_start[1], 2) +
					std::pow(wc_candidate_start[2] - lantern_candidate_start[2], 2) );
				
				// calculate end distance between WC and Lantern candidates
				float wc_candidate_end[3] = {_event.wc_reco_endXYZT[_event.wc_pion_candidate_index][0], _event.wc_reco_endXYZT[_event.wc_pion_candidate_index][1], _event.wc_reco_endXYZT[_event.wc_pion_candidate_index][2]};
				float lantern_candidate_end[3] = {_event.lantern_trackEndPosX[_event.lantern_pion_candidate_index], _event.lantern_trackEndPosY[_event.lantern_pion_candidate_index], _event.lantern_trackEndPosZ[_event.lantern_pion_candidate_index]};
				float end_distance = std::sqrt(
					std::pow(wc_candidate_end[0] - lantern_candidate_end[0], 2) +
					std::pow(wc_candidate_end[1] - lantern_candidate_end[1], 2) +
					std::pow(wc_candidate_end[2] - lantern_candidate_end[2], 2) );
				
			    // calculate angle between WC and Lantern candidates, in degrees
				// lantern
				TVector3 lantern_start_direction(_event.lantern_trackStartDirX[_event.lantern_pion_candidate_index], _event.lantern_trackStartDirY[_event.lantern_pion_candidate_index], _event.lantern_trackStartDirZ[_event.lantern_pion_candidate_index]);
				lantern_start_direction = lantern_start_direction.Unit();
				// WC
				TVector3 wc_start_direction(_event.wc_reco_startMomentum[_event.wc_pion_candidate_index][0], _event.wc_reco_startMomentum[_event.wc_pion_candidate_index][1], _event.wc_reco_startMomentum[_event.wc_pion_candidate_index][2]);
				wc_start_direction = wc_start_direction.Unit();
				// angle (radians)
				float angle_rad = lantern_start_direction.Angle(wc_start_direction);
				// convert to degrees
				float angle_deg = angle_rad * 180. / M_PI;

				// populate histograms
				_histStack_wc_lantern_candidate_start_distance.Fill(_event.classification, start_distance, pot_weight * _event.weight_cv);
				_histStack_wc_lantern_candidate_end_distance.Fill(_event.classification, end_distance, pot_weight * _event.weight_cv);
				_histStack_wc_lantern_candidate_angle_difference.Fill(_event.classification, angle_deg, pot_weight * _event.weight_cv);
			}

			} // End of kData Restriction

			// Energy Scale Tests
			/*
			// CC 0pi sample
			//if (_event.classification == Utility::kCCNumu0pi) {
			
				// get candidate track ID
				int recoID = _event.wc_reco_id[_event.wc_pion_candidate_index];

				// get spacepoints for candidate track
				std::vector<float> spacepoints_x; spacepoints_x.reserve(1000);
				std::vector<float> spacepoints_y; spacepoints_y.reserve(1000);
				std::vector<float> spacepoints_z; spacepoints_z.reserve(1000);
				std::vector<float> spacepoints_q; spacepoints_q.reserve(1000);

				for (int i = 0; i < _event.wc_Trecchargeblob_spacepoints_real_cluster_id->size(); i++) {
					if (_event.wc_Trecchargeblob_spacepoints_real_cluster_id->at(i) == recoID) {
						spacepoints_x.push_back(_event.wc_Trecchargeblob_spacepoints_x->at(i));
						spacepoints_y.push_back(_event.wc_Trecchargeblob_spacepoints_y->at(i));
						spacepoints_z.push_back(_event.wc_Trecchargeblob_spacepoints_z->at(i));
						spacepoints_q.push_back(_event.wc_Trecchargeblob_spacepoints_q->at(i));
					}
				}

				// check number of spacepoints for candidate, if zero return zero momentum (error case)
				int n_spacepoints = spacepoints_x.size();
				if (n_spacepoints == 0) {
					std::cout << "Error: no spacepoints found for candidate track!" << std::endl;
					continue;
				}

				// check whether first space point is start or end of track, by comparing to track start and end positions
				double start_dist = std::sqrt(std::pow(spacepoints_x.at(0) - _event.wc_reco_startXYZT[_event.wc_pion_candidate_index][0], 2) + std::pow(spacepoints_y.at(0) - _event.wc_reco_startXYZT[_event.wc_pion_candidate_index][1], 2) + std::pow(spacepoints_z.at(0) - _event.wc_reco_startXYZT[_event.wc_pion_candidate_index][2], 2));
				double end_dist = std::sqrt(std::pow(spacepoints_x.at(0) - _event.wc_reco_endXYZT[_event.wc_pion_candidate_index][0], 2) + std::pow(spacepoints_y.at(0) - _event.wc_reco_endXYZT[_event.wc_pion_candidate_index][1], 2) + std::pow(spacepoints_z.at(0) - _event.wc_reco_endXYZT[_event.wc_pion_candidate_index][2], 2));
				// reverse if filled backwards
				if (end_dist < start_dist) {
					//std::cout << "Reversing spacepoint order for candidate track..." << std::endl;
					std::reverse(spacepoints_x.begin(), spacepoints_x.end());
					std::reverse(spacepoints_y.begin(), spacepoints_y.end());
					std::reverse(spacepoints_z.begin(), spacepoints_z.end());
					std::reverse(spacepoints_q.begin(), spacepoints_q.end());
				}
		
				// calculate dx, dq, and dE for each spacepoint
				std::vector<double> dx; dx.reserve(n_spacepoints-1);
				std::vector<double> rr; rr.reserve(n_spacepoints-1);
				std::vector<double> dqdx; dqdx.reserve(n_spacepoints-1);
				std::vector<double> dEdx; dEdx.reserve(n_spacepoints-1);

				// WC values (???), unclear where these come from
				//double alpha = 1.0; 
				//double beta = 0.255;
				// Standard values
				double alpha = 0.93;
				double beta = 0.212;

				for (int i = 1; i < n_spacepoints; i++) {

					double sp_dx = std::sqrt(std::pow(spacepoints_x.at(i) - spacepoints_x.at(i-1), 2) + std::pow(spacepoints_y.at(i) - spacepoints_y.at(i-1), 2) + std::pow(spacepoints_z.at(i) - spacepoints_z.at(i-1), 2));
					double sp_q = spacepoints_q.at(i);
					double sp_dQ = (sp_q + 1000) * 10; // this is correct
					double sp_dQdx = sp_dQ / sp_dx;
					if (sp_dQdx/43e3 > 1000) sp_dQdx = 0;

					// convert to dE/dx
					// proper way
					double sp_dEdx = (std::exp(sp_dQdx * 23.6e-6*beta/1.38/0.273) - alpha)/(beta/1.38/0.273);
					if (sp_dEdx < 0) sp_dEdx = 0;
					if (sp_dEdx > 50) sp_dEdx = 50;
					
					// is there an energy scaling needed? Yes, ~1.2 factor shifts dE/dx 
					// dE/dx seems low -- to be investigated.
					sp_dEdx *= 1.2;					
					
					dx.push_back(sp_dx);
					dEdx.push_back(sp_dEdx);

					// calcuate rr, looping over remaining spacepoints to get total range
					// alternative - just use track end position to calculate rr
					double sp_rr = std::sqrt(std::pow(spacepoints_x.at(i) - _event.wc_reco_endXYZT[_event.wc_pion_candidate_index][0], 2) + std::pow(spacepoints_y.at(i) - _event.wc_reco_endXYZT[_event.wc_pion_candidate_index][1], 2) + std::pow(spacepoints_z.at(i) - _event.wc_reco_endXYZT[_event.wc_pion_candidate_index][2], 2));
					rr.push_back(sp_rr);

					// fill histogram
					dEdx_vs_ResRange->Fill(sp_rr, sp_dEdx);
				}
				
			//}
			*/
		}

		// clean up
		f->Close();
		delete f;
	}

	// write binning tree
	std::cout << "Writing binning tree." << std::endl;
	outputFile.cd();
	binningTree.Write();
	outputFile.Close();

	std::cout << "Passing Data Events: " << n_data_pass << std::endl;

	/*
	// draw dE/dx vs residual range plot
	TCanvas *res4 = new TCanvas("res4", "res4", 1080, 1080);
	res4->cd();
	dEdx_vs_ResRange->GetXaxis()->SetTitle("Residual Range [cm]");
	dEdx_vs_ResRange->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
	dEdx_vs_ResRange->Draw("COLZ");

	// create theoretical prediction
	std::vector<double> rr_vals;
	std::vector<double> mpv_dedx_vals;

	double pitch = 0.6; // cm

	for (int i = 0; i < 330; i++) {
		double rr = i*pitch;
		double KE = _physdedx.KEFromRangeSpline(rr);
		double mpv_dedx = _physdedx.MPVdEdx(KE, pitch);

		rr_vals.push_back(rr);
		mpv_dedx_vals.push_back(mpv_dedx);

	}

	// draw as TGraph
	TGraph *graph = new TGraph(mpv_dedx_vals.size(), &rr_vals[0], &mpv_dedx_vals[0]);
	graph->SetLineColor(kBlack);
	graph->SetLineStyle(9);
	graph->SetLineWidth(3);
	graph->Draw("L same");

	res4->Print("plots/plot_dEdx_vs_ResRange.root");
	*/

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

	std::cout << "Signal 0p: " << n_signal_0p << std::endl;
	std::cout << "Signal Np: " << n_signal_Np << std::endl;
	std::cout << "Signal 0p pass: " << n_signal_0p_pass << std::endl;
	std::cout << "Signal Np pass: " << n_signal_Np_pass << std::endl;

	std::cout << "Signal Pion Decay: " << n_signal_decay << ", Efficiency = " << (n_signal_pass_decay / n_signal_decay) * 100 << "%" << std::endl;
	std::cout << "Signal Pion Absorption 0p: " << n_signal_absorption0p << ", Efficiency = " << (n_signal_pass_absorption0p / n_signal_absorption0p) * 100 << "%" << std::endl;
	std::cout << "Signal Pion Absorption Np: " << n_signal_absorptionNp << ", Efficiency = " << (n_signal_pass_absorptionNp / n_signal_absorptionNp) * 100 << "%" << std::endl;
	std::cout << "Signal Pion Charge Exchange: " << n_signal_chargeexchange << ", Efficiency = " << (n_signal_pass_chargeexchange / n_signal_chargeexchange) * 100 << "%" << std::endl;

	std::cout << "Efficiency: " << (n_signal_pass / n_signal_all) * 100 << std::endl;
	std::cout << "Purity: " << (n_signal_pass / n_total_pass) * 100 << std::endl;

	std::cout << "Efficiency 0p: " << (n_signal_0p_pass / n_signal_0p) * 100 << std::endl;
	std::cout << "Efficiency Np: " << (n_signal_Np_pass / n_signal_Np) * 100 << std::endl;
	
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

	
	TCanvas *wc1 = new TCanvas("wc1", "wc1", 1080, 1080);
	_histStack_wc_numu_cc_flag.DrawStack(wc1, Utility::kWCNumuCCFlag);
  	wc1->Print("plots/plot_wc_numu_cc_flag.root");

	TCanvas *wc2 = new TCanvas("wc2", "wc2", 1080, 1080);
	_histStack_vertex_discrepancy.DrawStack(wc2, Utility::kTrackDistance);
  	wc2->Print("plots/plot_vertex_discrepancy.root");
	
	TCanvas *wc3 = new TCanvas("wc3", "wc3", 1080, 1080);
	_histStack_wc_match_isFC.DrawStack(wc3, Utility::kisFC);
  	wc3->Print("plots/plot_wc_match_isFC.root");

	// observables
	TCanvas *o1 = new TCanvas("o1", "o1", 1080, 1080);
  	o1->cd();
  	_histStack_sel_pion_theta.DrawStack(o1, Utility::kTheta);
  	o1->Print("plots/plot_pion_theta.root");

	TCanvas *o2 = new TCanvas("o2", "o2", 1080, 1080);
  	o2->cd();
  	_histStack_sel_pion_phi.DrawStack(o2, Utility::kPhi);
  	o2->Print("plots/plot_pion_phi.root");

	TCanvas *o3 = new TCanvas("o3", "o3", 1080, 1080);
  	o3->cd();
  	_histStack_sel_pion_momentum.DrawStack(o3, Utility::kTrackMomentum);
  	o3->Print("plots/plot_pion_momentum.root"); 

	TCanvas *o3b = new TCanvas("o3b", "o3b", 1080, 1080);
  	o3b->cd();
  	_histStack_sel_pion_momentum_interaction.DrawStack(o3b, Utility::kTrackMomentum);
  	o3b->Print("plots/plot_pion_momentum_interaction.root"); 

	// pion-proton opening angle
	TCanvas *o3c = new TCanvas("o3c", "o3c", 1080, 1080);
  	o3c->cd();
  	_histStack_sel_pion_proton_opening_angle.DrawStack(o3c, Utility::kOpeningAngle);
  	o3c->Print("plots/plot_pion_proton_opening_angle.root");

	TCanvas *o4 = new TCanvas("o4", "o4", 1080, 1080);
  	o4->cd();
  	_histStack_sel_sideband_larpidScore.DrawStack(o4, Utility::kLArPID_pi);
  	o4->Print("plots/plot_sideband_larpidScore.root"); 

	TCanvas *o5 = new TCanvas("o5", "o5", 1080, 1080);
  	o5->cd();
  	_histStack_sel_proton_momentum.DrawStack(o5, Utility::kTrackMomentum);
  	o5->Print("plots/plot_proton_momentum.root"); 

	TCanvas *o6 = new TCanvas("o6", "o6", 1080, 1080);
  	o6->cd();
  	_histStack_sel_proton_theta.DrawStack(o6, Utility::kTheta);
  	o6->Print("plots/plot_proton_theta.root");

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

	// lantern/wc candidate consistency
	TCanvas *l10 = new TCanvas("l10", "l10", 1080, 1080);
	_histStack_wc_lantern_candidate_start_distance.DrawStack(l10, Utility::kTrackDistance);
  	l10->Print("plots/plot_wc_lantern_candidate_start_distance.root");
	TCanvas *l11 = new TCanvas("l11", "l11", 1080, 1080);
	_histStack_wc_lantern_candidate_end_distance.DrawStack(l11, Utility::kTrackDistance);
  	l11->Print("plots/plot_wc_lantern_candidate_end_distance.root");
	TCanvas *l12 = new TCanvas("l12", "l12", 1080, 1080);
	_histStack_wc_lantern_candidate_angle_difference.DrawStack(l12, Utility::kOpeningAngle);
  	l12->Print("plots/plot_wc_lantern_candidate_angle_difference.root");

	// draw resolution plots
	TCanvas *res1 = new TCanvas("res1", "res1", 1080, 1080);
	res1->cd();
	pionMomentum1D->GetXaxis()->SetTitle("Pion Momentum ((Reco - True) / True)");
	pionMomentum1D->Draw("HIST");
  	res1->Print("plots/plot_pion_momentum_resolution_1D.root");
	
	TCanvas *res2 = new TCanvas("res2", "res2", 1080, 1080);
	res2->cd();
	pionMomentum2D->GetXaxis()->SetTitle("True Pion Momentum [GeV/c]");
    pionMomentum2D->GetYaxis()->SetTitle("Reconstructed Pion Momentum [GeV/c]");
	pionMomentum2D->Draw("COLZ");
  	res2->Print("plots/plot_pion_momentum_resolution_2D.root");	

	TCanvas *res3 = new TCanvas("res3", "res3", 1080, 1080);
	res3->cd();
	pionTheta2D->GetXaxis()->SetTitle("True Pion cos(#theta)");
	pionTheta2D->GetYaxis()->SetTitle("Reconstructed Pion cos(#theta)");
	pionTheta2D->Draw("COLZ");
  	res3->Print("plots/plot_pion_theta_resolution_2D.root");

	// draw resolution plots
	TCanvas *res3a = new TCanvas("res3a", "res3a", 1080, 1080);
	res3a->cd();
	pionTheta1D->GetXaxis()->SetTitle("Pion Theta ((Reco - True) / True)");
	pionTheta1D->Draw("HIST");
  	res3a->Print("plots/plot_pion_theta_resolution_1D.root");

	// pion-proton opening angle resolution 1D
	TCanvas *res4 = new TCanvas("res4", "res4", 1080, 1080);
	res4->cd();
	pionProtonOpeningAngle1D->GetXaxis()->SetTitle("Pion-Proton Opening Angle ((Reco - True) / True)");
	pionProtonOpeningAngle1D->Draw("HIST");
  	res4->Print("plots/plot_pion_proton_opening_angle_resolution_1D.root");

	TCanvas *res5 = new TCanvas("res5", "res5", 1080, 1080);
	res5->cd();
	pionProtonOpeningAngle2D->GetXaxis()->SetTitle("True Pion-Proton Opening Angle [deg]");
	pionProtonOpeningAngle2D->GetYaxis()->SetTitle("Reconstructed Pion-Proton Opening Angle [deg]");
	pionProtonOpeningAngle2D->Draw("COLZ");
  	res5->Print("plots/plot_pion_proton_opening_angle_resolution_2D.root");

}
// ------------------------------------------------------------------------------

// Compare Lantern vertex score for true signal events with a Lantern vertex, Genie vs NuWro (Run 5)
void SelectionDriver::compareLanternVtxScoreGenieNuWro() {

	const double data_pot_run5 = 1.425e20;

	std::vector<std::string> filename_list = {filename_mc_run5, filename_nuwro_run5};
	std::vector<Utility::FileTypeEnums> file_types_list = {Utility::kMC, Utility::kFakeData};
	std::vector<std::string> label_list = {"GENIE", "NuWro"};

	Selection _selection(_utility);

	// [sample][0 = 0p, 1 = Np]
	std::vector<std::vector<TH1F*>> hist_list;
	std::vector<std::string> proton_label_list = {"0p", "Np"};

	for (unsigned int idx = 0; idx < filename_list.size(); idx++) {

		std::cout << "Processing file " << filename_list[idx] << std::endl;

		TFile *f = new TFile(filename_list[idx].c_str());

		// sum MC POT from Pandora subrun tree
		TTree *subrun_tree = (TTree*)f->Get("nuselection/SubRun");
		float pot = 0;
		subrun_tree->SetBranchAddress("pot", &pot);
		double mc_pot = 0;
		for (int s = 0; s < subrun_tree->GetEntries(); s++) {
			subrun_tree->GetEntry(s);
			mc_pot += pot;
		}
		double pot_weight = data_pot_run5 / mc_pot;
		std::cout << label_list[idx] << " MC POT: " << mc_pot << ", POT weight: " << pot_weight << std::endl;

		TTree *pelee_tree = (TTree*)f->Get("nuselection/NeutrinoSelectionFilter");
		TTree *wc_BDTvars_tree = (TTree*)f->Get("wcpselection/T_BDTvars");
		TTree *wc_eval_tree = (TTree*)f->Get("wcpselection/T_eval");
		TTree *wc_PFeval_tree = (TTree*)f->Get("wcpselection/T_PFeval");
		TTree *wc_spacepoints_tree = (TTree*)f->Get("wcpselection/T_spacepoints");
		TTree *lantern_tree = (TTree*)f->Get("lantern/EventTree");

		EventContainer _event(pelee_tree, wc_eval_tree, wc_BDTvars_tree, wc_PFeval_tree, wc_spacepoints_tree, lantern_tree, _utility);

		std::vector<TH1F*> h;
		for (unsigned int ip = 0; ip < proton_label_list.size(); ip++) {
			h.push_back(new TH1F(Form("h_vtxScore_%s_%s", label_list[idx].c_str(), proton_label_list[ip].c_str()), "", 40, 0, 1.001));
			h[ip]->Sumw2();
		}

		int n_entries = pelee_tree->GetEntries();
		std::cout << "Initial number events: " << n_entries << std::endl;

		for (int e = 0; e < n_entries; e++) {

			pelee_tree->GetEntry(e);

			// skip known problematic MC events
			if (file_types_list[idx] == Utility::kMC) {
				if (_event.run == 20872 && _event.sub == 14 && _event.evt == 749) continue;
				if (_event.run == 20607 && _event.sub == 355 && _event.evt == 17791) continue;
			}

			wc_BDTvars_tree->GetEntry(e);
			wc_eval_tree->GetEntry(e);
			wc_PFeval_tree->GetEntry(e);
			wc_spacepoints_tree->GetEntry(e);
			lantern_tree->GetEntry(e);

			if ( (e != 0) && (n_entries >= 10) &&  (e % (n_entries/10) == 0) ) {
				std::cout << Form("%i0%% Completed...\n", e / (n_entries/10));
			}

			// populates classification and CV weight
			_selection.ApplyLanternSelection(_event, file_types_list[idx], Utility::kRun4ab);

			if (_event.classification != Utility::kNC1pi) continue;
			if (_event.lantern_foundVertex != 1) continue;

			// split by number of protons above threshold
			if (_event.mc_NProtons_ == 0) h[0]->Fill(_event.lantern_vtxScore, pot_weight * _event.weight_cv);
			else if (_event.mc_NProtons_ > 0) h[1]->Fill(_event.lantern_vtxScore, pot_weight * _event.weight_cv);
		}

		for (unsigned int ip = 0; ip < proton_label_list.size(); ip++) {
			std::cout << label_list[idx] << " signal " << proton_label_list[ip] << " events with Lantern vertex (scaled to " << data_pot_run5 << " POT): " << h[ip]->Integral() << std::endl;
			h[ip]->SetDirectory(0);
		}
		hist_list.push_back(h);
		f->Close();
	}

	// draw: one pad per proton multiplicity, GENIE vs NuWro in each
	gStyle->SetOptStat(0);
	TCanvas *c = new TCanvas("c_vtxScore", "", 1400, 600);
	c->Divide(2, 1);

	std::vector<int> colour_list = {kAzure+2, kOrange+7};

	for (unsigned int ip = 0; ip < proton_label_list.size(); ip++) {
		c->cd(ip + 1);
		gPad->SetLeftMargin(0.13);

		double ymax = 0;
		for (unsigned int idx = 0; idx < hist_list.size(); idx++) {
			hist_list[idx][ip]->SetLineColor(colour_list[idx]);
			hist_list[idx][ip]->SetLineWidth(2);
			ymax = std::max(ymax, hist_list[idx][ip]->GetMaximum());
		}

		hist_list[0][ip]->SetMaximum(1.2 * ymax);
		hist_list[0][ip]->SetMinimum(0);
		hist_list[0][ip]->SetTitle(Form("True NC1#pi %s with Lantern vertex, Run 5 (1.425 #times 10^{20} POT);Lantern vertex score;Events", proton_label_list[ip].c_str()));
		hist_list[0][ip]->Draw("hist e");
		for (unsigned int idx = 1; idx < hist_list.size(); idx++) hist_list[idx][ip]->Draw("hist e same");

		TLegend *leg = new TLegend(0.15, 0.72, 0.5, 0.88);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		for (unsigned int idx = 0; idx < hist_list.size(); idx++) {
			leg->AddEntry(hist_list[idx][ip], Form("%s (%.1f events)", label_list[idx].c_str(), hist_list[idx][ip]->Integral()), "l");
		}
		leg->Draw();
	}

	gSystem->mkdir("plots", true);
	c->Print("plots/plot_lantern_vtxScore_genie_vs_nuwro_0pNp.root");
	c->Print("plots/plot_lantern_vtxScore_genie_vs_nuwro_0pNp.png");
	c->Print("plots/plot_lantern_vtxScore_genie_vs_nuwro_0pNp.pdf");
}
