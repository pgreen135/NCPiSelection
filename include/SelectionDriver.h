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

	// Compare Lantern vtxScore for signal, Genie vs NuWro (Run 5)
	void compareLanternVtxScoreGenieNuWro();

	// ----------------------------------


protected:

	// classes
	Utility _utility;
  
	// --- NTuple Files ---
	// NuWro
	std::string filename_nuwro_run1 = "/Volumes/T9/BNBSurpriseNTuples/nuwro/checkout_MCC9.10_Run123_v10_04_07_23_BNB_nuwro_overlay_surprise_reco2_hist_1_filtered_slim.root";
	std::string filename_nuwro_run2 = "/Volumes/T9/BNBSurpriseNTuples/nuwro/checkout_MCC9.10_Run123_v10_04_07_23_BNB_nuwro_overlay_surprise_reco2_hist_2_filtered_slim.root";
	std::string filename_nuwro_run3 = "/Volumes/T9/BNBSurpriseNTuples/nuwro/checkout_MCC9.10_Run123_v10_04_07_23_BNB_nuwro_overlay_surprise_reco2_hist_3_filtered_slim.root";
	std::string filename_nuwro_run4 = "/Volumes/T9/BNBSurpriseNTuples/nuwro/checkout_MCC9.10_Run45_v10_04_07_23_BNB_nuwro_overlay_surprise_reco2_hist_4_slim.root";
	std::string filename_nuwro_run5 = "/Volumes/T9/BNBSurpriseNTuples/nuwro/checkout_MCC9.10_Run45_v10_04_07_23_BNB_nuwro_overlay_surprise_reco2_hist_5_slim.root";

	// Beam On
	// Using newly processed ntuple files
	std::string filename_beamon_run4 = "/Volumes/T9/BNBSurpriseNTuples/beam_on/checkout_MCC9.10_Run4_v10_04_07_23_BNB_beam_on_retuple_retuple_hist_4_slim.root"; 

	// Beam Off
	std::string filename_beamoff_run1 = "/Volumes/T9/BNBSurpriseNTuples/beam_off/checkout_MCC9.10_Run123_v10_04_07_20_BNB_beam_off_data_surprise_reco2_hist_1_slim.root";
	std::string filename_beamoff_run2 = "/Volumes/T9/BNBSurpriseNTuples/beam_off/checkout_MCC9.10_Run123_v10_04_07_20_BNB_beam_off_data_surprise_reco2_hist_2_slim.root";
	std::string filename_beamoff_run3 = "/Volumes/T9/BNBSurpriseNTuples/beam_off/checkout_MCC9.10_Run123_v10_04_07_20_BNB_beam_off_data_surprise_reco2_hist_3_slim.root";
	std::string filename_beamoff_run4 = "/Volumes/T9/BNBSurpriseNTuples/beam_off/checkout_MCC9.10_Run4_v10_04_07_23_BNB_beam_off_retuple_retuple_hist_4_slim.root";
	std::string filename_beamoff_run5 = "/Volumes/T9/BNBSurpriseNTuples/beam_off/checkout_MCC9.10_Run4acd5_v10_04_07_23_BNB_beam_off_retuple_retuple_hist_5_slim.root";
	
	// Nu Overlay 
	std::string filename_mc_run1 = "/Volumes/T9/BNBSurpriseNTuples/nu_overlay/checkout_MCC9.10_Run123_v10_04_07_20_BNB_nu_overlay_surprise_reco2_hist_1_FSIrwgt.root";
	std::string filename_mc_run2 = "/Volumes/T9/BNBSurpriseNTuples/nu_overlay/checkout_MCC9.10_Run123_v10_04_07_20_BNB_nu_overlay_surprise_reco2_hist_2_FSIrwgt.root";
	std::string filename_mc_run3 = "/Volumes/T9/BNBSurpriseNTuples/nu_overlay/checkout_MCC9.10_Run123_v10_04_07_20_BNB_nu_overlay_surprise_reco2_hist_3_FSIrwgt.root";
	std::string filename_mc_run4 = "/Volumes/T9/BNBSurpriseNTuples/nu_overlay/checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4_FSIrwgt.root";
	std::string filename_mc_run5 = "/Volumes/T9/BNBSurpriseNTuples/nu_overlay/checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_5_FSIrwgt.root"; 
	
	// Dirt Overlay
	std::string filename_dirt_run1 = "/Volumes/T9/BNBSurpriseNTuples/dirt/checkout_MCC9.10_Run123_v10_04_07_23_BNB_dirt_overlay_surprise_reco2_hist_1_FSIrwgt.root";
	std::string filename_dirt_run2 = "/Volumes/T9/BNBSurpriseNTuples/dirt/checkout_MCC9.10_Run123_v10_04_07_23_BNB_dirt_overlay_surprise_reco2_hist_2_FSIrwgt.root";
	std::string filename_dirt_run3 = "/Volumes/T9/BNBSurpriseNTuples/dirt/checkout_MCC9.10_Run123_v10_04_07_23_BNB_dirt_overlay_surprise_reco2_hist_3_FSIrwgt.root";
	std::string filename_dirt_run4 = "/Volumes/T9/BNBSurpriseNTuples/dirt/checkout_MCC9.10_Run4b_v10_04_07_20_BNB_dirt_overlay_retuple_retuple_hist_FSIrwgt.root";
	std::string filename_dirt_run5 = "/Volumes/T9/BNBSurpriseNTuples/dirt/checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_5_FSIrwgt.root"; 

	// --- POT Weights ---
	// Beam On
	// 1.3e21 POT, ~294236984 Triggers
	// Run 1: 1.545e+20, 34344460.0
	// Run 2: 2.589e+20, 61296355.0
	// Run 3: 2.545e+20, 61112238.0
	// Run 4:  1.313e+20 + 8.961e+19 + 4.906e+19; 31160144.0 + 20316612.0 + 11150094.0 [4b, 4c, 4d]
	// Run 5:  1.425e+20;  34060578.0
    double pot_weight_beamon_run4 = 1; 

	// Nuwro (don't scale)
	double pot_weight_nuwro_run12345 = 1; 
   
	// Beam Off [ Up to date ]
	// Run 1: 49551151.0
	// Run 2: 152116269.0
	// Run 3: 143426455.0
	// Run 4: 84966908.0 + 52102363.0 + 73918592.0 [4b, 4c, 4d]
	// Run 5: 105586073.0
	double pot_weight_beamoff_run1 = (double)34344460.0 / (double)49551151.0;
	double pot_weight_beamoff_run2 = (double)61296355.0 / (double)152116269.0;
	double pot_weight_beamoff_run3 = (double)61112238.0 / (double)143426455.0;
	double pot_weight_beamoff_run4 = (double)62626850.0 / (double)210987863.0;
	double pot_weight_beamoff_run5 = (double)34060578.0 / (double)105586073.0;

	// Nu Overlay [Up to date]
	double pot_weight_mc_run1 = (double)1.545e20 / (double)7.38236e20;
	double pot_weight_mc_run2 = (double)2.589e20 / (double)9.06254e20;
	double pot_weight_mc_run3 = (double)2.545e20 / (double)8.6286e20;
	double pot_weight_mc_run4 = (double)2.700e20 / (double)2.21635e21;
	double pot_weight_mc_run5 =  (double)1.425e20 / (double)9.46919e20;

	// Dirt Overlay [ Up to date ]
	double pot_weight_dirt_run1 = (double)1.545e20 / (double)3.7306e20;
	double pot_weight_dirt_run2 = (double)2.589e20 / (double)7.43595e20;
	double pot_weight_dirt_run3 = (double)2.545e20 / (double)3.62877e20;
	double pot_weight_dirt_run4 = (double)2.700e20 / (double)8.95758e20; // Run 4b only until retupled samples available
	double pot_weight_dirt_run5 =  (double)1.425e20 / (double)3.26707e20;

};

#endif