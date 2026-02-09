// Driver to post-proccess ntuples for S. Gardiner's cross-section tools

#include "../include/Utility.h"
#include "../include/EventContainer.h"
#include "../include/Selection.h"

#include <string>
#include <vector>
#include <iostream>
#include <map>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"

// Helper function to set branch addresses for the output TTree
void set_event_output_branch_addresses(TTree& out_tree, EventContainer& ev, bool create);

// Helper functions that creates a branch (or just sets a new address)
// simple variable in an output TTree
void set_output_branch_address(TTree& out_tree, const std::string& branch_name, void* address, bool create, const std::string& leaf_spec);
// pointer to an object in an output TTree
template <typename T> void set_object_output_branch_address( TTree& out_tree, const std::string& branch_name, T*& address, bool create = false);
// vector of doubles
void set_output_branch_address_vector(TTree& out_tree, const std::string& branch_name, void* address, std::vector<double>* address_ptr, bool create);

int main(int argc, char *argv[]) {

	// parse arguments
	if ( argc != 5 ) {
    	std::cout << "Usage: ProcessNTuple INPUT_PELEE_NTUPLE_FILE OUTPUT_FILE TYPE_ENUM RUN_ENUM" << std::endl;
    	return 1;
	}

	std::string input_filename( argv[1] );
	std::string output_filename( argv[2] );
	std::string type_str ( argv[3] );
	std::string runPeriod_str ( argv[4] );
	Utility::FileTypeEnums type = static_cast<Utility::FileTypeEnums>(stoi(type_str));
	Utility::RunPeriodEnums runPeriod = static_cast<Utility::RunPeriodEnums>(stoi(runPeriod_str));

	std::cout << "Processing file: " << input_filename << ", type enum: " << type << ", run period enum: " << runPeriod << std::endl;
  
	// Get the TTrees containing the event ntuples and subrun POT information
	TFile *inputFile = new TFile(input_filename.c_str()); 
	TTree *pelee_tree = (TTree*)inputFile->Get("nuselection/NeutrinoSelectionFilter");
	TTree *subrunsTree = (TTree*)inputFile->Get("nuselection/SubRun");

	TTree *wc_BDTvars_tree = (TTree*)inputFile->Get("wcpselection/T_BDTvars");
	TTree *wc_eval_tree = (TTree*)inputFile->Get("wcpselection/T_eval");
	TTree *wc_PFeval_tree = (TTree*)inputFile->Get("wcpselection/T_PFeval");
	TTree *wc_POT_tree = (TTree*)inputFile->Get("wcpselection/T_pot");
	TTree *lantern_tree = (TTree*)inputFile->Get("lantern/EventTree");

	// initialise classes
	Utility _utility(false, false, false, false);
	Selection _selection(_utility);
	EventContainer _event(pelee_tree, wc_eval_tree, wc_BDTvars_tree, wc_PFeval_tree, lantern_tree, _utility);
    
  	// Create output tree
  	TFile* outFile = new TFile( output_filename.c_str(), "recreate" );
	outFile->cd();
	TTree* outTree = new TTree( "stv_tree", "STV analysis tree" );

	// Get the total POT from the subruns TTree. Save it in the output
	// TFile as a TParameter<float>. Real data doesn't have this TTree,
	// so check that it exists first.
	double pot;
	float summed_pot = 0.;
	//bool has_pot_branch = ( subrunsTree->GetBranch("pot") != nullptr );
	bool has_pot_branch = ( wc_POT_tree->GetBranch("pot_tor875good") != nullptr );
	if ( has_pot_branch ) {
    //subrunsTree->SetBranchAddress( "pot", &pot );
	wc_POT_tree->SetBranchAddress( "pot_tor875good", &pot );
    for ( int se = 0; se < subrunsTree->GetEntries(); ++se ) {
    	//subrunsTree->GetEntry( se );
		wc_POT_tree->GetEntry( se );
    	summed_pot += pot;
   	}
	std::cout << "Total POT in file: " << summed_pot << std::endl;
	}

  	TParameter<float>* summed_pot_param = new TParameter<float>( "summed_pot", summed_pot );
	summed_pot_param->Write();
	
	bool created_output_branches = false;
	
	// Event loop
	int n_entries = pelee_tree->GetEntries();
	std::cout << "Number events: " << n_entries << std::endl;

	for (int e = 0; e < n_entries; e++) {
	//for (int e = 180500; e < 180600; e++) {

		pelee_tree->GetEntry(e);

		// skip issue MC events
		if (type == Utility::kMC) {
			// Run 4b
			if (_event.run == 20872 && _event.sub == 14 && _event.evt == 749) {
				std::cout << "Skipping problematic MC event: Run " << _event.run << ", Subrun " << _event.sub << ", Event " << _event.evt << std::endl;
				continue;
			}
			if (_event.run == 20607 && _event.sub == 355 && _event.evt == 17791) {
				std::cout << "Skipping problematic MC event: Run " << _event.run << ", Subrun " << _event.sub << ", Event " << _event.evt << std::endl;
				continue;
			}
		}
		  
		// get current entry
		wc_BDTvars_tree->GetEntry(e);
		wc_eval_tree->GetEntry(e);
		wc_PFeval_tree->GetEntry(e); 
		lantern_tree->GetEntry(e); 

		if ( (e != 0) && (n_entries >= 10) &&  (e % (n_entries/10) == 0) ) {
		std::cout << Form("%i0%% Completed...\n", e / (n_entries/10));
		}

		// apply selection
		// populates necessary variables to pass to stv tree
		bool passSelection = _selection.ApplyLanternSelection(_event, type, runPeriod);
    
		// set the output TTree branch addresses, creating the branches if needed
		// (during the first event loop iteration)
		if ( !created_output_branches ) {
			set_event_output_branch_addresses( *outTree, _event, true );
			created_output_branches = true;
		}
		else set_event_output_branch_addresses( *outTree, _event, false );

		// fill output tree
		outTree->Fill();

	} // end of event loop

	std::cout << "Finished processing events. Now writing output file: " << output_filename << std::endl;

	// write output file
	outTree->Write();
  	outFile->Close();
  	delete outFile;

	return 0;
}


// Helper function to set branch addresses for the output TTree
void set_event_output_branch_addresses(TTree& out_tree, EventContainer& ev, bool create = false) {
	
	// Signal definition flags
	set_output_branch_address(out_tree, "is_mc", &ev.is_mc_, create, "is_mc/O");
	set_output_branch_address(out_tree, "mc_is_signal", &ev.mc_is_signal_, create, "mc_is_signal/O");
	
	// MC event category
	set_output_branch_address( out_tree, "NC1pi_EventCategory", &ev.category_, create, "NC1pi_EventCategory/I");

	// CV Event weights
	set_output_branch_address( out_tree, "spline_weight", &ev.weightSpline, create, "spline_weight/F");
	set_output_branch_address( out_tree, "tuned_cv_weight", &ev.weightTune, create, "tuned_cv_weight/F");
	
	// If MC weights are available, store them in the output TTree
	if ( ev.mc_weights_map_ ) {

		// Make separate branches for the various sets of systematic variation
		// weights in the map
		for ( auto& pair : *ev.mc_weights_map_ ) {

				if (pair.first == "ppfx_all") continue; 	// skip empty PPFX flux weights

				// Prepend "weight_" to the name of the vector of weights in the map
				std::string weight_branch_name = "weight_" + pair.first;

				// Store a pointer to the vector of weights (needed to set the branch
				// address properly) in the temporary map of pointers
				ev.mc_weights_ptr_map_[ weight_branch_name ] = &pair.second;

				// Set the branch address for this vector of weights
				set_object_output_branch_address< std::vector<double> >( out_tree, weight_branch_name, ev.mc_weights_ptr_map_.at(weight_branch_name), create);
		}
	}

	// NC pi+ selection criteria
	set_output_branch_address( out_tree, "sel_passInitialSelection", &ev.sel_passInitialSelection_, create, "sel_passInitialSelection/O");
	set_output_branch_address( out_tree, "sel_passMuPiLLR", &ev.sel_passMuPiLLR_, create, "sel_passMuPiLLR/O");
	set_output_branch_address( out_tree, "sel_NC1pi", &ev.sel_NC1pi_, create, "sel_NC1pi/O");

	// Observables
	// LLR PID Scores
	set_output_branch_address( out_tree, "sel_LanternPID_llr_mu_pi", &ev.sel_LanternPID_llr_mu_pi_, create, "sel_LanternPID_llr_mu_pi/F");	// Reco
	set_output_branch_address( out_tree, "sel_LanternPID_llr_pr_pi", &ev.sel_LanternPID_llr_pr_pi_, create, "sel_LanternPID_llr_pr_pi/F");	// Reco

	// Shower/Electron energy
	//set_output_branch_address( out_tree, "sel_shr_energy_cali", &ev.shr_energy_cali, create, "sel_shr_energy_cali/F");	// Reco
	//set_output_branch_address( out_tree, "mc_shr_bkt_E", &ev.shr_bkt_E, create, "mc_shr_bkt_E/F");						// Truth - shower energy
	//set_output_branch_address( out_tree, "mc_elec_e", &ev.elec_e, create, "mc_elec_e/F");										// Truth - electron energy

	// Electron NuMI angle (beta)
	//set_output_branch_address( out_tree, "sel_reco_cos_electron_effective_angle", &ev.reco_cos_electron_effective_angle, create, "sel_reco_cos_electron_effective_angle/F");	// Reco
	//set_output_branch_address( out_tree, "mc_true_cos_electron_effective_angle", &ev.true_cos_electron_effective_angle, create, "mc_true_cos_electron_effective_angle/F");		// Truth

	// Pion Momentum
	//set_output_branch_address( out_tree, "sel_reco_momentum_pion", &ev.reco_momentum_pion, create, "sel_reco_momentum_pion/F");	// Reco - range-based
	//set_output_branch_address( out_tree, "mc_pion_p", &ev.pion_p, create, "mc_pion_p/F");											// Truth - pion energy

	// Pion NuMI angle (beta)
	//set_output_branch_address( out_tree, "sel_reco_cos_pion_effective_angle", &ev.reco_cos_pion_effective_angle, create, "sel_reco_cos_pion_effective_angle/F");	// Reco
	//set_output_branch_address( out_tree, "mc_true_cos_pion_effective_angle", &ev.true_cos_pion_effective_angle, create, "mc_true_cos_pion_effective_angle/F");		// Truth

	// Opening angle
	//set_output_branch_address( out_tree, "sel_reco_cos_electron_pion_opening_angle", &ev.reco_cos_electron_pion_opening_angle, create, "sel_reco_cos_electron_pion_opening_angle/F");	// Reco
	//set_output_branch_address( out_tree, "mc_true_cos_electron_pion_opening_angle", &ev.true_cos_electron_pion_opening_angle, create, "mc_true_cos_electron_pion_opening_angle/F");		// Truth

	// Number protons
	//set_output_branch_address( out_tree, "sel_numberProtons", &ev.numberProtons, create, "sel_numberProtons/I");	// Reco
	//set_output_branch_address( out_tree, "mc_nproton", &ev.nproton, create, "mc_nproton/I");											// Truth

	/*
	// Track BDT properties
	set_output_branch_address( out_tree, "sel_bragg_mip_pion_loose", &ev.trk_bragg_mip_pion_loose, create, "sel_bragg_mip_pion_loose/F");
	set_output_branch_address( out_tree, "sel_trk_daughters_pion_loose", &ev.trk_daughters_pion_loose, create, "sel_trk_daughters_pion_loose/I");
	set_output_branch_address( out_tree, "sel_trk_dEdx_trunk_pion_loose", &ev.trk_dEdx_trunk_pion_loose, create, "sel_trk_dEdx_trunk_pion_loose/F");
	set_output_branch_address( out_tree, "sel_trk_bragg_pion_pion_loose", &ev.trk_bragg_pion_pion_loose, create, "sel_trk_bragg_pion_pion_loose/F");
	set_output_branch_address( out_tree, "sel_trk_llr_pid_score_pion_loose", &ev.trk_llr_pid_score_pion_loose, create, "sel_trk_llr_pid_score_pion_loose/F");
	set_output_branch_address( out_tree, "sel_trk_score_pion_loose", &ev.trk_score_pion_loose, create, "sel_trk_score_pion_loose/F");
	set_output_branch_address( out_tree, "sel_trk_end_spacepoints_pion_loose", &ev.trk_end_spacepoints_pion_loose, create, "sel_trk_end_spacepoints_pion_loose/I");

	set_output_branch_address( out_tree, "sel_bragg_mip_pion", &ev.trk_bragg_mip_pion, create, "sel_bragg_mip_pion/F");
	set_output_branch_address( out_tree, "sel_trk_daughters_pion", &ev.trk_daughters_pion, create, "sel_trk_daughters_pion/I");
	set_output_branch_address( out_tree, "sel_trk_dEdx_trunk_pion", &ev.trk_dEdx_trunk_pion, create, "sel_trk_dEdx_trunk_pion/F");
	set_output_branch_address( out_tree, "sel_trk_bragg_pion_pion", &ev.trk_bragg_pion_pion, create, "sel_trk_bragg_pion_pion/F");
	set_output_branch_address( out_tree, "sel_trk_llr_pid_score_pion", &ev.trk_llr_pid_score_pion, create, "sel_trk_llr_pid_score_pion/F");
	set_output_branch_address( out_tree, "sel_trk_score_pion", &ev.trk_score_pion, create, "sel_trk_score_pion/F");
	set_output_branch_address( out_tree, "sel_trk_end_spacepoints_pion", &ev.trk_end_spacepoints_pion, create, "sel_trk_end_spacepoints_pion/I");

	// Shower BDT properties
	set_output_branch_address(out_tree,  "sel_n_showers_contained", &ev.n_showers_contained, create, "sel_n_showers_contained/I");
	set_output_branch_address( out_tree, "sel_shrmoliereavg", &ev.shrmoliereavg, create, "sel_shrmoliereavg/F");
	set_output_branch_address( out_tree, "sel_shr_distance", &ev.shr_distance, create, "sel_shr_distance/F");
	set_output_branch_address( out_tree, "sel_shr2_pfpgeneration", &ev.shr2_pfpgeneration, create, "sel_shr2_pfpgeneration/I");
	set_output_branch_address( out_tree, "sel_shr_trkfit_2cm_dedx_best", &ev.shr_trkfit_2cm_dedx_max, create, "sel_shr_trkfit_2cm_dedx_best/F");
	set_output_branch_address( out_tree, "sel_shr_trkfit_gap10_dedx_best", &ev.shr_trkfit_gap10_dedx_max, create, "sel_shr_trkfit_gap10_dedx_best/F");
	set_output_branch_address( out_tree, "sel_shr_energyFraction", &ev.shr_energyFraction, create, "sel_shr_energyFraction/F");
	set_output_branch_address( out_tree, "sel_shr2_distance", &ev.shr2_distance, create, "sel_shr2_distance/F");
  set_output_branch_address( out_tree, "sel_shrMCSMom", &ev.shrMCSMom, create, "sel_shrMCSMom/F");
  */

  // Interaction type
  set_output_branch_address( out_tree, "mc_ccnc", &ev.ccnc, create, "mc_ccnc/I");
  set_output_branch_address( out_tree, "mc_interaction", &ev.interaction, create, "mc_interaction/I");
}

// Helper function that creates a branch (or just sets a new address) for a
// simple variable in an output TTree
void set_output_branch_address(TTree& out_tree, const std::string& branch_name, void* address, bool create = false, const std::string& leaf_spec = "" ) {
  if ( create ) {
    if ( leaf_spec != "" ) {
      out_tree.Branch( branch_name.c_str(), address, leaf_spec.c_str() );
    }
    else {
      out_tree.Branch( branch_name.c_str(), address );
    }
  }
  else {
    out_tree.SetBranchAddress( branch_name.c_str(), address );
  }
}
// pointer to an object in an output TTree
template <typename T> void set_object_output_branch_address( TTree& out_tree, const std::string& branch_name, T*& address, bool create) {
  if ( create ) out_tree.Branch( branch_name.c_str(), &address );
  else out_tree.SetBranchAddress( branch_name.c_str(), &address );
}
// vector of doubles for beamline weights (to do, make template)
void set_output_branch_address_vector(TTree& out_tree, const std::string& branch_name, void* address, std::vector<double>* address_ptr, bool create) {
	 if ( create ) out_tree.Branch( branch_name.c_str(), "std::vector<double>", address );
   else {
   	out_tree.SetBranchAddress( branch_name.c_str(), address_ptr );
   }
} 

