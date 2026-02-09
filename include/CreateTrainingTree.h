#ifndef CREATETRAININGTREE_H
#define CREATETRAININGTREE_H

// class for creating root trees used for BDT training
// selects required sanitized variables from event and writes to flat root file

#include "../include/Utility.h"
#include "../include/EventContainer.h"

#include <TFile.h>
#include <TTree.h>

class CreateTrainingTree {

public:	

	// ----------------------------------
	// Constructor
	CreateTrainingTree(Utility::BDTEnums type);

	// Destructor
	~CreateTrainingTree();

	// Functions to add an event
	void addEvent(const EventContainer &_evt, Utility::ClassificationEnums classification);

	// Function to write output file
	void writeOutputFile();

	// ----------------------------------

private:

	// type
	Utility::BDTEnums BDTType;

	// output file
	TFile *outputFile;

	// trees
	TTree *trainingTree;

	// branches
	bool isSignal;

	// --- electron-photon separation tree ---	
	float n_showers_contained;
	
	float shr_distance;
	float shr_trkfit_gap10_dedx_max;
	float shr_trkfit_2cm_dedx_max;

	float shrmoliereavg;
	float shr_energyFraction;
	float shrsubclusters;

	float shrPCA1CMed_5cm;
    float CylFrac2h_1cm;
    float DeltaRMS2h;
    float shrMCSMom;

	float shr_energy_second_cali;
	float shr2_distance;
	float shr12_p1_dstart;
	float shr2_trackEndProximity;
	float shr2_pfpgeneration;

	float secondshower_Y_nhit;
	float secondshower_Y_vtxdist;
	float secondshower_Y_dot;
	float secondshower_Y_anglediff;
	
	float secondshower_U_nhit;
	float secondshower_U_vtxdist;
	float secondshower_U_dot;
	float secondshower_U_anglediff;

	float secondshower_V_nhit;
	float secondshower_V_vtxdist;
	float secondshower_V_dot;
	float secondshower_V_anglediff;

	// --- pion-proton separation tree --
	float trk_score;
    float trk_llr_pid_score;
    float trk_bragg_mip_max;
    float trk_bragg_pion_max;
    float trk_dEdx_trunk_max;
    float trk_daughters;
    float trk_end_spacepoints;

	// --- additional variables for alternate pion-proton separation tree -- 
	float trk2_llr_pid_score;
    float trk2_bragg_mip_max;
    float trk2_bragg_pion_max;
    float trk2_dEdx_trunk_max;
    float trk2_daughters;
    float trk2_end_spacepoints;
    float trk3_llr_pid_score;
    float trk3_bragg_mip_max;
    float trk3_bragg_pion_max;
    float trk3_dEdx_trunk_max;
    float trk3_daughters;
    float trk3_end_spacepoints;

	// WC muon pion separation tree
	float wc_track_length;
	float wc_reco_pdg;
	float wc_reco_larpid_pidScore_mu;
	float wc_reco_larpid_pidScore_pi;
	float wc_reco_larpid_pidScore_pr;
	float wc_n_pion_candidate_daughters;
	float wc_pion_candidate_daughters_total_energy;
	float wc_n_blips_25cm;
	float wc_n_blips_50cm;
	float wc_n_blips_100cm;

}; 

#endif