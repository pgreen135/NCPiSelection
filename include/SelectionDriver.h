#ifndef SELECTIONDRIVER_H
#define SELECTIONDRIVER_H

// class to run selections on datasets

#include <string>

#include "../include/Utility.h"
#include "../include/EventContainer.h"


class SelectionDriver {

public:	

	// ----------------------------------
	// Constructor
	SelectionDriver();

	// Destructor
	~SelectionDriver(){};

	// Run BDT selection
    void runBDTSelectionFull();

	// ----------------------------------


protected:

	// classes
	Utility _utility;
  
	// --- NTuple Files ---
	// Beam On
	// Using newly processed ntuple files
	std::string filename_beamon_run4 = ""; 
	
	// Beam Off
	std::string filename_beamoff_run4b = "/Volumes/T9/BNBSurpriseNTuples/beam_off/checkout_MCC9.10_Run4b_v10_04_07_20_BNB_beam_off_metapatch_retuple_retuple_hist.root";
	std::string filename_beamoff_run4c = "/Volumes/T9/BNBSurpriseNTuples/beam_off/checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4c.root";
	std::string filename_beamoff_run4d = "/Volumes/T9/BNBSurpriseNTuples/beam_off/checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4d.root";
	std::string filename_beamoff_run5 = "/Volumes/T9/BNBSurpriseNTuples/beam_off/checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_5.root";
	
	// Nu Overlay
	std::string filename_mc_run4b = "/Volumes/T9/BNBSurpriseNTuples/nu_overlay/checkout_MCC9.10_Run4b_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist.root";
	std::string filename_mc_run4c = "/Volumes/T9/BNBSurpriseNTuples/nu_overlay/checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4c.root";
	std::string filename_mc_run4d = "/Volumes/T9/BNBSurpriseNTuples/nu_overlay/checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4d.root";
	std::string filename_mc_run5 =  "/Volumes/T9/BNBSurpriseNTuples/nu_overlay/checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_5.root"; 
	
	// Dirt Overlay
	std::string filename_dirt_run4b = "/Volumes/T9/BNBSurpriseNTuples/dirt/checkout_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root";
	std::string filename_dirt_run4c = "/Volumes/T9/BNBSurpriseNTuples/dirt/checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4c.root";
	std::string filename_dirt_run4d = "/Volumes/T9/BNBSurpriseNTuples/dirt/checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4d.root";
	std::string filename_dirt_run5 =  "/Volumes/T9/BNBSurpriseNTuples/dirt/checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_5.root"; 

	// --- POT Weights ---
	// Beam On
    double pot_weight_beamon_run4 = 1; // 1.3e21 POT, ~294236984 Triggers
   
	// Beam Off
	//double pot_weight_beamoff_run4b = (double)294236984.0 / (double)89010180.0; // 89010180.0 Triggers
	double pot_weight_beamoff_run4b = (double)32305463.0 / (double)89010180.0; // 89010180.0 Triggers
	double pot_weight_beamoff_run4c = (double)20273291.0 / (double)53659787.0; // 53659787.0 Triggers
	double pot_weight_beamoff_run4d = (double)11192660.0 / (double)76563108.0; // 76563108.0 Triggers
	double pot_weight_beamoff_run5 = (double)35265730.0 / (double)111457148; // 89010180.0 Triggers
	//double pot_weight_beamoff_run4b = (double)(32305463.0 + 20273291.0 + 11192660.0 + 35265730.0) / (double)88445969.0;

	// Nu Overlay
	//double pot_weight_mc_run4b = (double)1.3e21 / (double)7.87801e+20;
	double pot_weight_mc_run4b = (double)1.36e20 / (double)7.87801e+20;
	double pot_weight_mc_run4c = (double)0.895e20 / (double)4.70313e+20; 
	double pot_weight_mc_run4d = (double)0.493e20 / (double)8.94895e+20; 
	double pot_weight_mc_run5 =  (double)1.48e20 / (double)1.00406e+21;

	// Dirt Overlay
	//double pot_weight_dirt_run4b = (double)1.3e21 / (double)3.05858e+20;
	double pot_weight_dirt_run4b = (double)1.36e20 / (double)3.05858e+20;
	double pot_weight_dirt_run4c = (double)0.895e20 / (double)1.79887e+20;
	double pot_weight_dirt_run4d = (double)0.493e20 / (double)3.47696e+20;
	double pot_weight_dirt_run5 =  (double)1.48e20 / (double)3.52038e+20;

};

#endif