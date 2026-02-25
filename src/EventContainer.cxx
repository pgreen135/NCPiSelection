#include "../include/EventContainer.h"

#include <iostream>
#include <algorithm>
#include <numeric>

#include "TVector3.h"
#include "TRotation.h"

// Constructor
EventContainer::EventContainer(TTree *pelee_tree, TTree *wc_eval_tree, TTree *wc_BDTvars_tree, TTree *wc_PFeval_tree, TTree *wc_spacepoints_tree, TTree *lantern_tree, const Utility &utility): _utility{ utility } {
	
	std::cout << "Initialising Event Container Class" << std::endl;

	// set wc_BDTvars_tree addresses
	wc_BDTvars_tree->SetBranchAddress("run", &wc_run);
	wc_BDTvars_tree->SetBranchAddress("subrun", &wc_sub);
	wc_BDTvars_tree->SetBranchAddress("event", &wc_evt);

	wc_BDTvars_tree->SetBranchAddress("numu_cc_flag", &wc_numu_cc_flag);

	// set wc_Eval_tree addresses
	wc_eval_tree->SetBranchAddress("match_isFC", &wc_match_isFC);

	// set wc_PFEval_tree addresses
	wc_PFeval_tree->SetBranchAddress("reco_nuvtxX", &wc_reco_nuvtxX);
	wc_PFeval_tree->SetBranchAddress("reco_nuvtxY", &wc_reco_nuvtxY);
	wc_PFeval_tree->SetBranchAddress("reco_nuvtxZ", &wc_reco_nuvtxZ);

	wc_PFeval_tree->SetBranchAddress("reco_Ntrack", &wc_reco_Ntrack);
	wc_PFeval_tree->SetBranchAddress("reco_id", wc_reco_id);
	wc_PFeval_tree->SetBranchAddress("reco_pdg", wc_reco_pdg);
	wc_PFeval_tree->SetBranchAddress("reco_mother", wc_reco_mother);
	wc_PFeval_tree->SetBranchAddress("reco_startXYZT", wc_reco_startXYZT);
	wc_PFeval_tree->SetBranchAddress("reco_endXYZT", wc_reco_endXYZT);
	wc_PFeval_tree->SetBranchAddress("reco_startMomentum", wc_reco_startMomentum);
	wc_PFeval_tree->SetBranchAddress("reco_larpid_classified", wc_reco_larpid_classified);
	wc_PFeval_tree->SetBranchAddress("reco_larpid_pdg", wc_reco_larpid_pdg);
	wc_PFeval_tree->SetBranchAddress("reco_larpid_pidScore_el", wc_reco_larpid_pidScore_el);
	wc_PFeval_tree->SetBranchAddress("reco_larpid_pidScore_ph", wc_reco_larpid_pidScore_ph);
	wc_PFeval_tree->SetBranchAddress("reco_larpid_pidScore_mu", wc_reco_larpid_pidScore_mu);
	wc_PFeval_tree->SetBranchAddress("reco_larpid_pidScore_pi", wc_reco_larpid_pidScore_pi);
	wc_PFeval_tree->SetBranchAddress("reco_larpid_pidScore_pr", wc_reco_larpid_pidScore_pr);
	wc_PFeval_tree->SetBranchAddress("reco_larpid_procScore_prim", wc_reco_larpid_procScore_prim);
	wc_PFeval_tree->SetBranchAddress("reco_larpid_procScore_chgd", wc_reco_larpid_procScore_chgd);
	wc_PFeval_tree->SetBranchAddress("reco_larpid_procScore_ntrl", wc_reco_larpid_procScore_ntrl);

	wc_PFeval_tree->SetBranchAddress("reco_truthMatch_pdg", wc_reco_truthMatch_pdg);
	wc_PFeval_tree->SetBranchAddress("reco_truthMatch_id", wc_reco_truthMatch_id);

	wc_PFeval_tree->SetBranchAddress("truth_Ntrack", &wc_truth_Ntrack);
	wc_PFeval_tree->SetBranchAddress("truth_id", wc_truth_id);
	wc_PFeval_tree->SetBranchAddress("truth_pdg", wc_truth_pdg);
	wc_PFeval_tree->SetBranchAddress("truth_mother", wc_truth_mother);
	wc_PFeval_tree->SetBranchAddress("truth_startMomentum", wc_truth_startMomentum);

	wc_PFeval_tree->SetBranchAddress("truth_process", &wc_truth_process);
	wc_PFeval_tree->SetBranchAddress("truth_endprocess", &wc_truth_endprocess);

	// spacepoints tree addresses
	wc_spacepoints_tree->SetBranchAddress("Trecchargeblob_spacepoints_x", &wc_Trecchargeblob_spacepoints_x);
	wc_spacepoints_tree->SetBranchAddress("Trecchargeblob_spacepoints_y", &wc_Trecchargeblob_spacepoints_y);
	wc_spacepoints_tree->SetBranchAddress("Trecchargeblob_spacepoints_z", &wc_Trecchargeblob_spacepoints_z);
	wc_spacepoints_tree->SetBranchAddress("Trecchargeblob_spacepoints_real_cluster_id", &wc_Trecchargeblob_spacepoints_real_cluster_id);
	wc_spacepoints_tree->SetBranchAddress("Trecchargeblob_spacepoints_q", &wc_Trecchargeblob_spacepoints_q);

	// set lantern_tree addresses
	lantern_tree->SetBranchAddress("run", &lantern_run);
	lantern_tree->SetBranchAddress("subrun", &lantern_sub);
	lantern_tree->SetBranchAddress("event", &lantern_evt);

	lantern_tree->SetBranchAddress("trueNuPDG", &lantern_trueNuPDG);
	lantern_tree->SetBranchAddress("trueNuCCNC", &lantern_trueNuCCNC);
	lantern_tree->SetBranchAddress("trueVtxX", &lantern_trueVtxX);
	lantern_tree->SetBranchAddress("trueVtxY", &lantern_trueVtxY);
	lantern_tree->SetBranchAddress("trueVtxZ", &lantern_trueVtxZ);

	lantern_tree->SetBranchAddress("foundVertex", &lantern_foundVertex);
	lantern_tree->SetBranchAddress("vtxX", &lantern_vtxX);
	lantern_tree->SetBranchAddress("vtxY", &lantern_vtxY);
	lantern_tree->SetBranchAddress("vtxZ", &lantern_vtxZ);
	lantern_tree->SetBranchAddress("vtxScore", &lantern_vtxScore);
	lantern_tree->SetBranchAddress("vtxContainment", &lantern_vtxContainment);

	lantern_tree->SetBranchAddress("vtxFracHitsOnCosmic", &lantern_vtxFracHitsOnCosmic);

	lantern_tree->SetBranchAddress("nTracks", &lantern_nTracks);
	lantern_tree->SetBranchAddress("trackIsSecondary", lantern_trackIsSecondary);
	lantern_tree->SetBranchAddress("trackStartPosX", lantern_trackStartPosX);
	lantern_tree->SetBranchAddress("trackStartPosY", lantern_trackStartPosY);
	lantern_tree->SetBranchAddress("trackStartPosZ", lantern_trackStartPosZ);
	lantern_tree->SetBranchAddress("trackEndPosX", lantern_trackEndPosX);
	lantern_tree->SetBranchAddress("trackEndPosY", lantern_trackEndPosY);
	lantern_tree->SetBranchAddress("trackEndPosZ", lantern_trackEndPosZ);

	lantern_tree->SetBranchAddress("trackClassified", lantern_trackClassified);
	lantern_tree->SetBranchAddress("trackPID", lantern_trackPID);
	lantern_tree->SetBranchAddress("trackMuScore", lantern_trackMuScore);
	lantern_tree->SetBranchAddress("trackPiScore", lantern_trackPiScore);
	lantern_tree->SetBranchAddress("trackPrScore", lantern_trackPrScore);

	lantern_tree->SetBranchAddress("trackTruePID", lantern_trackTruePID);

	lantern_tree->SetBranchAddress("nShowers", &lantern_nShowers);
	lantern_tree->SetBranchAddress("showerIsSecondary", lantern_showerIsSecondary);

	// set pelee_tree addresses
	pelee_tree->SetBranchAddress("run", &run);
	pelee_tree->SetBranchAddress("sub", &sub);
	pelee_tree->SetBranchAddress("evt", &evt);

	pelee_tree->SetBranchAddress("nu_pdg", &nu_pdg);
	pelee_tree->SetBranchAddress("ccnc", &ccnc);
	pelee_tree->SetBranchAddress("interaction", &interaction);
	pelee_tree->SetBranchAddress("nu_e", &nu_e);
	
	pelee_tree->SetBranchAddress("nmuon", &nmuon);
	pelee_tree->SetBranchAddress("muon_e", &muon_e);
	pelee_tree->SetBranchAddress("nelec", &nelec);
	pelee_tree->SetBranchAddress("elec_e", &elec_e);
	pelee_tree->SetBranchAddress("elec_px", &elec_px);
    pelee_tree->SetBranchAddress("elec_py", &elec_py);
    pelee_tree->SetBranchAddress("elec_pz", &elec_pz);
	pelee_tree->SetBranchAddress("npion", &npion);
	pelee_tree->SetBranchAddress("pion_e", &pion_e);
	pelee_tree->SetBranchAddress("npi0", &npi0);
	pelee_tree->SetBranchAddress("pi0_e", &pi0_e);
	pelee_tree->SetBranchAddress("nproton", &nproton);
	pelee_tree->SetBranchAddress("proton_e", &proton_e);
	pelee_tree->SetBranchAddress("neta", &neta);
    pelee_tree->SetBranchAddress("eta_e", &eta_e);

	pelee_tree->SetBranchAddress("true_nu_vtx_x", &true_nu_vtx_x);
	pelee_tree->SetBranchAddress("true_nu_vtx_y", &true_nu_vtx_y);
	pelee_tree->SetBranchAddress("true_nu_vtx_z", &true_nu_vtx_z);
	pelee_tree->SetBranchAddress("true_nu_vtx_sce_x", &true_nu_vtx_sce_x);
	pelee_tree->SetBranchAddress("true_nu_vtx_sce_y", &true_nu_vtx_sce_y);
	pelee_tree->SetBranchAddress("true_nu_vtx_sce_z", &true_nu_vtx_sce_z);
	pelee_tree->SetBranchAddress("true_nu_px", &true_nu_px);
    pelee_tree->SetBranchAddress("true_nu_py", &true_nu_py);
    pelee_tree->SetBranchAddress("true_nu_pz", &true_nu_pz);

	pelee_tree->SetBranchAddress("swtrig", &swtrig);    

    pelee_tree->SetBranchAddress("_opfilter_pe_beam", &opfilter_pe_beam);
    pelee_tree->SetBranchAddress("_opfilter_pe_veto", &opfilter_pe_veto);

    pelee_tree->SetBranchAddress("weightSplineTimesTune", &weightSplineTimesTune);
    pelee_tree->SetBranchAddress("weightSpline", &weightSpline);
    pelee_tree->SetBranchAddress("weightTune", &weightTune);
    pelee_tree->SetBranchAddress("ppfx_cv", &ppfx_cv);
    pelee_tree->SetBranchAddress("weights", &mc_weights_map_);

	pelee_tree->SetBranchAddress("nslice", &nslice);
	pelee_tree->SetBranchAddress("n_tracks", &n_tracks);
    pelee_tree->SetBranchAddress("n_showers", &n_showers);
    pelee_tree->SetBranchAddress("n_showers_contained", &n_showers_contained);
    pelee_tree->SetBranchAddress("n_tracks_contained", &n_tracks_contained);
    pelee_tree->SetBranchAddress("slpdg", &slpdg);
	pelee_tree->SetBranchAddress("nu_purity_from_pfp", &nu_purity_from_pfp);
	pelee_tree->SetBranchAddress("slice_orig_pass_id", &slice_orig_pass_id);

	pelee_tree->SetBranchAddress("shr_hits_tot", &shr_hits_tot);
	pelee_tree->SetBranchAddress("trk_hits_tot", &trk_hits_tot);
	pelee_tree->SetBranchAddress("trk_hits_max", &trk_hits_max);
	pelee_tree->SetBranchAddress("trk_hits_2nd", &trk_hits_2nd);

	pelee_tree->SetBranchAddress("total_hits_y", &total_hits_y);
	pelee_tree->SetBranchAddress("shr_hits_y_tot", &shr_hits_y_tot);
	pelee_tree->SetBranchAddress("trk_hits_y_tot", &trk_hits_y_tot);
	pelee_tree->SetBranchAddress("extra_energy_y", &extra_energy_y);

	pelee_tree->SetBranchAddress("crtveto", &crtveto);

	pelee_tree->SetBranchAddress("shr_hits_u_tot", &shr_hits_u_tot);
	pelee_tree->SetBranchAddress("trk_hits_u_tot", &trk_hits_u_tot);

	pelee_tree->SetBranchAddress("shr_hits_v_tot", &shr_hits_v_tot);
	pelee_tree->SetBranchAddress("trk_hits_v_tot", &trk_hits_v_tot);

	pelee_tree->SetBranchAddress("reco_nu_vtx_sce_x", &reco_nu_vtx_sce_x);
    pelee_tree->SetBranchAddress("reco_nu_vtx_sce_y", &reco_nu_vtx_sce_y);
    pelee_tree->SetBranchAddress("reco_nu_vtx_sce_z", &reco_nu_vtx_sce_z);

    pelee_tree->SetBranchAddress("contained_fraction", &contained_fraction);

    pelee_tree->SetBranchAddress("NeutrinoEnergy2", &NeutrinoEnergy2);

    pelee_tree->SetBranchAddress("topological_score", &topological_score);
    pelee_tree->SetBranchAddress("CosmicIPAll3D", &CosmicIPAll3D);
    pelee_tree->SetBranchAddress("CosmicDirAll3D", &CosmicDirAll3D);

    pelee_tree->SetBranchAddress("shr_id", &shr_id);
    pelee_tree->SetBranchAddress("shr2_id", &shr2_id);
    pelee_tree->SetBranchAddress("shr3_id", &shr3_id);
    pelee_tree->SetBranchAddress("shr_distance", &shr_distance);
    pelee_tree->SetBranchAddress("shr_score", &shr_score);
    pelee_tree->SetBranchAddress("shr_energy_cali", &shr_energy_cali);
    pelee_tree->SetBranchAddress("shr_energy_second_cali", &shr_energy_second_cali);
    pelee_tree->SetBranchAddress("shr_energy_third_cali", &shr_energy_third_cali);
    pelee_tree->SetBranchAddress("shr_energy_tot_cali", &shr_energy_tot_cali);    
    pelee_tree->SetBranchAddress("hits_ratio", &hits_ratio);
    pelee_tree->SetBranchAddress("shrmoliereavg", &shrmoliereavg);
    pelee_tree->SetBranchAddress("shr1shr2moliereavg", &shr1shr2moliereavg);
    pelee_tree->SetBranchAddress("shr1trk2moliereavg", &shr1trk2moliereavg);

    pelee_tree->SetBranchAddress("shr_llrpid_dedx", &shr_llrpid_dedx);

    pelee_tree->SetBranchAddress("trkfit", &trkfit);
 
    pelee_tree->SetBranchAddress("shr_tkfit_dedx_Y", &shr_tkfit_dedx_Y);
    pelee_tree->SetBranchAddress("shr_tkfit_dedx_V", &shr_tkfit_dedx_V);
    pelee_tree->SetBranchAddress("shr_tkfit_dedx_U", &shr_tkfit_dedx_U);

    pelee_tree->SetBranchAddress("shr_tkfit_nhits_Y", &shr_tkfit_nhits_Y);
    pelee_tree->SetBranchAddress("shr_tkfit_nhits_V", &shr_tkfit_nhits_V);
    pelee_tree->SetBranchAddress("shr_tkfit_nhits_U", &shr_tkfit_nhits_U);

    pelee_tree->SetBranchAddress("shr_tkfit_gap05_dedx_Y", &shr_tkfit_gap05_dedx_Y);
    pelee_tree->SetBranchAddress("shr_tkfit_gap05_dedx_V", &shr_tkfit_gap05_dedx_V);
    pelee_tree->SetBranchAddress("shr_tkfit_gap05_dedx_U", &shr_tkfit_gap05_dedx_U);

    pelee_tree->SetBranchAddress("shr_tkfit_gap05_nhits_Y", &shr_tkfit_gap05_nhits_Y);
    pelee_tree->SetBranchAddress("shr_tkfit_gap05_nhits_V", &shr_tkfit_gap05_nhits_V);
    pelee_tree->SetBranchAddress("shr_tkfit_gap05_nhits_U", &shr_tkfit_gap05_nhits_U);

    pelee_tree->SetBranchAddress("shr_tkfit_gap10_dedx_Y", &shr_tkfit_gap10_dedx_Y);
    pelee_tree->SetBranchAddress("shr_tkfit_gap10_dedx_V", &shr_tkfit_gap10_dedx_V);
    pelee_tree->SetBranchAddress("shr_tkfit_gap10_dedx_U", &shr_tkfit_gap10_dedx_U);

    pelee_tree->SetBranchAddress("shr_tkfit_gap10_nhits_Y", &shr_tkfit_gap10_nhits_Y);
    pelee_tree->SetBranchAddress("shr_tkfit_gap10_nhits_V", &shr_tkfit_gap10_nhits_V);
    pelee_tree->SetBranchAddress("shr_tkfit_gap10_nhits_U", &shr_tkfit_gap10_nhits_U);

    pelee_tree->SetBranchAddress("shr_tkfit_2cm_dedx_Y", &shr_tkfit_2cm_dedx_Y);
    pelee_tree->SetBranchAddress("shr_tkfit_2cm_dedx_V", &shr_tkfit_2cm_dedx_V);
    pelee_tree->SetBranchAddress("shr_tkfit_2cm_dedx_U", &shr_tkfit_2cm_dedx_U);

    pelee_tree->SetBranchAddress("shr_tkfit_2cm_nhits_Y", &shr_tkfit_2cm_nhits_Y);
    pelee_tree->SetBranchAddress("shr_tkfit_2cm_nhits_V", &shr_tkfit_2cm_nhits_V);
    pelee_tree->SetBranchAddress("shr_tkfit_2cm_nhits_U", &shr_tkfit_2cm_nhits_U);

    pelee_tree->SetBranchAddress("shrclusdir0", &shrclusdir0);
    pelee_tree->SetBranchAddress("shrclusdir1", &shrclusdir1);
    pelee_tree->SetBranchAddress("shrclusdir2", &shrclusdir2);

    pelee_tree->SetBranchAddress("shrsubclusters0", &shrsubclusters0);
    pelee_tree->SetBranchAddress("shrsubclusters1", &shrsubclusters1);
    pelee_tree->SetBranchAddress("shrsubclusters2", &shrsubclusters2);

    pelee_tree->SetBranchAddress("shrPCA1CMed_5cm", &shrPCA1CMed_5cm);
    pelee_tree->SetBranchAddress("CylFrac1h_1cm", &CylFrac1h_1cm);
    pelee_tree->SetBranchAddress("CylFrac2h_1cm", &CylFrac2h_1cm);
    pelee_tree->SetBranchAddress("DeltaRMS2h", &DeltaRMS2h);
    pelee_tree->SetBranchAddress("shrMCSMom", &shrMCSMom);

    pelee_tree->SetBranchAddress("shr_bkt_E", &shr_bkt_E);

    pelee_tree->SetBranchAddress("secondshower_U_nhit", &secondshower_U_nhit);
    pelee_tree->SetBranchAddress("secondshower_V_nhit", &secondshower_V_nhit);
    pelee_tree->SetBranchAddress("secondshower_Y_nhit", &secondshower_Y_nhit);
    pelee_tree->SetBranchAddress("secondshower_U_vtxdist", &secondshower_U_vtxdist);
    pelee_tree->SetBranchAddress("secondshower_V_vtxdist", &secondshower_V_vtxdist);
    pelee_tree->SetBranchAddress("secondshower_Y_vtxdist", &secondshower_Y_vtxdist);
    pelee_tree->SetBranchAddress("secondshower_U_dot", &secondshower_U_dot);
    pelee_tree->SetBranchAddress("secondshower_V_dot", &secondshower_V_dot);
    pelee_tree->SetBranchAddress("secondshower_Y_dot", &secondshower_Y_dot);	
    pelee_tree->SetBranchAddress("secondshower_U_dir", &secondshower_U_dir);
    pelee_tree->SetBranchAddress("secondshower_V_dir", &secondshower_V_dir);
    pelee_tree->SetBranchAddress("secondshower_Y_dir", &secondshower_Y_dir);
    
    pelee_tree->SetBranchAddress("trk_id", &trk_id);
    pelee_tree->SetBranchAddress("trk2_id", &trk2_id);
    pelee_tree->SetBranchAddress("trk3_id", &trk3_id);
    pelee_tree->SetBranchAddress("trk_len", &trk_len);
    pelee_tree->SetBranchAddress("trk_distance", &trk_distance);
    pelee_tree->SetBranchAddress("trk_score", &trk_score);
    pelee_tree->SetBranchAddress("trk_theta", &trk_theta);
    pelee_tree->SetBranchAddress("trk_phi", &trk_phi);
    pelee_tree->SetBranchAddress("trk_bkt_pdg", &trk_bkt_pdg);

    pelee_tree->SetBranchAddress("trk_bragg_p", &trk_bragg_p);
    pelee_tree->SetBranchAddress("trk_bragg_mu", &trk_bragg_mu);
    pelee_tree->SetBranchAddress("trk_bragg_pion", &trk_bragg_pion);
    pelee_tree->SetBranchAddress("trk_bragg_mip", &trk_bragg_mip);

    pelee_tree->SetBranchAddress("tksh_distance", &tksh_distance);
    pelee_tree->SetBranchAddress("tksh_angle", &tksh_angle);

    pelee_tree->SetBranchAddress("pfnplanehits_U",          &pfnplanehits_U);
    pelee_tree->SetBranchAddress("pfnplanehits_V",          &pfnplanehits_V);
    pelee_tree->SetBranchAddress("pfnplanehits_Y",          &pfnplanehits_Y);

    pelee_tree->SetBranchAddress("trk_score_v",          &trk_score_v);
    pelee_tree->SetBranchAddress("trk_start_x_v",        &trk_start_x_v);
    pelee_tree->SetBranchAddress("trk_start_y_v",        &trk_start_y_v);
    pelee_tree->SetBranchAddress("trk_start_z_v",        &trk_start_z_v);
    pelee_tree->SetBranchAddress("trk_end_x_v",          &trk_end_x_v);
    pelee_tree->SetBranchAddress("trk_end_y_v",          &trk_end_y_v);
    pelee_tree->SetBranchAddress("trk_end_z_v",          &trk_end_z_v);
    pelee_tree->SetBranchAddress("trk_sce_start_x_v",        &trk_sce_start_x_v);
    pelee_tree->SetBranchAddress("trk_sce_start_y_v",        &trk_sce_start_y_v);
    pelee_tree->SetBranchAddress("trk_sce_start_z_v",        &trk_sce_start_z_v);
    pelee_tree->SetBranchAddress("trk_sce_end_x_v",          &trk_sce_end_x_v);
    pelee_tree->SetBranchAddress("trk_sce_end_y_v",          &trk_sce_end_y_v);
    pelee_tree->SetBranchAddress("trk_sce_end_z_v",          &trk_sce_end_z_v);
    pelee_tree->SetBranchAddress("trk_dir_x_v",          &trk_dir_x_v);
    pelee_tree->SetBranchAddress("trk_dir_y_v",          &trk_dir_y_v);
    pelee_tree->SetBranchAddress("trk_dir_z_v",          &trk_dir_z_v);

    pelee_tree->SetBranchAddress("trk_len_v",                &trk_len_v);
    pelee_tree->SetBranchAddress("trk_distance_v",           &trk_distance_v);
	pelee_tree->SetBranchAddress("trk_calo_energy_y_v",      &trk_calo_energy_y_v);
	pelee_tree->SetBranchAddress("trk_energy_proton_v",      &trk_energy_proton_v);
	pelee_tree->SetBranchAddress("trk_energy_muon_v",      &trk_energy_muon_v);
	
    pelee_tree->SetBranchAddress("trk_bragg_p_v", &trk_bragg_p_v);
    pelee_tree->SetBranchAddress("trk_bragg_mu_v", &trk_bragg_mu_v);
    pelee_tree->SetBranchAddress("trk_bragg_pion_v", &trk_bragg_pion_v);
    pelee_tree->SetBranchAddress("trk_bragg_pion_u_v", &trk_bragg_pion_u_v);
    pelee_tree->SetBranchAddress("trk_bragg_pion_v_v", &trk_bragg_pion_v_v);
    pelee_tree->SetBranchAddress("trk_bragg_mip_v", &trk_bragg_mip_v);
    pelee_tree->SetBranchAddress("trk_bragg_mip_u_v", &trk_bragg_mip_u_v);
    pelee_tree->SetBranchAddress("trk_bragg_mip_v_v", &trk_bragg_mip_v_v);
    pelee_tree->SetBranchAddress("trk_llr_pid_score_v", &trk_llr_pid_score_v);

    pelee_tree->SetBranchAddress("all_trk_hits", &all_trk_hits);
    
    pelee_tree->SetBranchAddress("all_shr_energies",          &all_shr_energies);
    pelee_tree->SetBranchAddress("shr_energy_u_v",          &shr_energy_u_v);
    pelee_tree->SetBranchAddress("shr_energy_v_v",          &shr_energy_v_v);
    pelee_tree->SetBranchAddress("shr_energy_y_v",          &shr_energy_y_v);
    pelee_tree->SetBranchAddress("shr_px_v",          &shr_px_v);
    pelee_tree->SetBranchAddress("shr_py_v",          &shr_py_v);
    pelee_tree->SetBranchAddress("shr_pz_v",          &shr_pz_v);
    pelee_tree->SetBranchAddress("shr_dist_v",        &shr_dist_v);
    pelee_tree->SetBranchAddress("shr_moliere_avg_v", &shr_moliere_avg_v);
    
    pelee_tree->SetBranchAddress("shr_start_x_v",     &shr_start_x_v);
    pelee_tree->SetBranchAddress("shr_start_y_v",     &shr_start_y_v);
    pelee_tree->SetBranchAddress("shr_start_z_v",     &shr_start_z_v);

    pelee_tree->SetBranchAddress("shr_tkfit_dedx_u_v",     &shr_tkfit_dedx_u_v);
    pelee_tree->SetBranchAddress("shr_tkfit_dedx_v_v",     &shr_tkfit_dedx_v_v);
    pelee_tree->SetBranchAddress("shr_tkfit_dedx_y_v",     &shr_tkfit_dedx_y_v);
    pelee_tree->SetBranchAddress("shr_tkfit_dedx_nhits_u_v",    &shr_tkfit_dedx_nhits_u_v);
    pelee_tree->SetBranchAddress("shr_tkfit_dedx_nhits_v_v",    &shr_tkfit_dedx_nhits_v_v);
    pelee_tree->SetBranchAddress("shr_tkfit_dedx_nhits_y_v",    &shr_tkfit_dedx_nhits_y_v);

    pelee_tree->SetBranchAddress("pfp_trk_daughters_v", &pfp_trk_daughters_v);
 	pelee_tree->SetBranchAddress("pfp_shr_daughters_v", &pfp_shr_daughters_v); 

 	pelee_tree->SetBranchAddress("trk_end_spacepoints_v", &trk_end_spacepoints_v);
 	 
 	pelee_tree->SetBranchAddress("pfpplanesubclusters_U", &pfpplanesubclusters_U_v);
 	pelee_tree->SetBranchAddress("pfpplanesubclusters_V", &pfpplanesubclusters_V_v);
 	pelee_tree->SetBranchAddress("pfpplanesubclusters_Y", &pfpplanesubclusters_Y_v);
    
    pelee_tree->SetBranchAddress("pfp_generation_v", &pfp_generation_v);
 	  
 	pelee_tree->SetBranchAddress("backtracked_pdg", &backtracked_pdg_v);
	pelee_tree->SetBranchAddress("backtracked_px", &backtracked_px_v);
	pelee_tree->SetBranchAddress("backtracked_py", &backtracked_py_v);
	pelee_tree->SetBranchAddress("backtracked_pz", &backtracked_pz_v); 	
    
    pelee_tree->SetBranchAddress("pi0_mass_U",         &pi0_mass_U);

    pelee_tree->SetBranchAddress("mc_pdg",       &mc_pdg_v);
    pelee_tree->SetBranchAddress("mc_E",         &mc_E_v);
    pelee_tree->SetBranchAddress("mc_px",           &mc_px_v);
    pelee_tree->SetBranchAddress("mc_py",           &mc_py_v);
    pelee_tree->SetBranchAddress("mc_pz",           &mc_pz_v); 

	pelee_tree->SetBranchAddress("mc_process",       &mc_process_v);
	pelee_tree->SetBranchAddress("mc_end_process",   &mc_end_process_v);

    // track trunk dE/dx
    pelee_tree->SetBranchAddress("trk_nhits_u_v",		&trk_nhits_u_v);
    pelee_tree->SetBranchAddress("trk_nhits_v_v",		&trk_nhits_v_v);
    pelee_tree->SetBranchAddress("trk_nhits_y_v",		&trk_nhits_y_v);
    pelee_tree->SetBranchAddress("trk_trunk_dEdx_u_v",	&trk_trunk_dEdx_u_v);
    pelee_tree->SetBranchAddress("trk_trunk_dEdx_v_v",	&trk_trunk_dEdx_v_v);
    pelee_tree->SetBranchAddress("trk_trunk_dEdx_y_v",	&trk_trunk_dEdx_y_v);

	pelee_tree->SetBranchAddress("trk_avg_deflection_stdev_v", &trk_avg_deflection_stdev_v);

	// NuGraph
	pelee_tree->SetBranchAddress("pfng2mipfrac", &pfng2mipfrac);
	pelee_tree->SetBranchAddress("pfng2hipfrac", &pfng2hipfrac);
	pelee_tree->SetBranchAddress("pfng2shrfrac", &pfng2shrfrac);
	pelee_tree->SetBranchAddress("pfng2mclfrac", &pfng2mclfrac);
	pelee_tree->SetBranchAddress("pfng2dfsfrac", &pfng2dfsfrac);

	// Blips
	pelee_tree->SetBranchAddress("nblips_saved", &nblips_saved);
	pelee_tree->SetBranchAddress("blip_energy", &blip_energy_v);
	pelee_tree->SetBranchAddress("blip_x", &blip_x_v);
	pelee_tree->SetBranchAddress("blip_y", &blip_y_v);
	pelee_tree->SetBranchAddress("blip_z", &blip_z_v);
	pelee_tree->SetBranchAddress("blip_proxtrkdist", &blip_proxtrkdist_v);
	pelee_tree->SetBranchAddress("blip_touchtrk", &blip_touchtrk_v);

	// Load FSI weight histograms
	TFile *FSI_NC_piminus = new TFile("../fsiweights/fsi_weights_nc_piminus.root");
    TFile *FSI_NC_piplus = new TFile("../fsiweights/fsi_weights_nc_piplus.root");
	TFile *FSI_CC_piminus = new TFile("../fsiweights/fsi_weights_cc_piminus.root");
    TFile *FSI_CC_piplus = new TFile("../fsiweights/fsi_weights_cc_piplus.root");
   
    if (FSI_NC_piminus->IsZombie() || FSI_NC_piplus->IsZombie() || FSI_CC_piminus->IsZombie() || FSI_CC_piplus->IsZombie()) {
      std::cout << "Error opening ratio files\n";
      exit(111);
    }

    h_fsi_nc_piminus  = (TH2D*)FSI_NC_piminus->Get("h_ratio");
    h_fsi_nc_piplus = (TH2D*)FSI_NC_piplus->Get("h_ratio");
    h_fsi_cc_piminus  = (TH2D*)FSI_CC_piminus->Get("h_ratio");
    h_fsi_cc_piplus = (TH2D*)FSI_CC_piplus->Get("h_ratio");

    if (!h_fsi_nc_piminus || !h_fsi_nc_piplus || !h_fsi_cc_piminus || !h_fsi_cc_piplus) {
      std::cout << "Error retrieving ratio histograms\n";
      exit(222);
    }

    h_fsi_nc_piminus->SetDirectory(0);
    h_fsi_nc_piplus->SetDirectory(0);
    h_fsi_cc_piminus->SetDirectory(0);
    h_fsi_cc_piplus->SetDirectory(0);

    FSI_NC_piminus->Close();
    FSI_NC_piplus->Close();
    FSI_CC_piminus->Close();
    FSI_CC_piplus->Close();

	// Hypfit initialization
	hypfit = new Hypfit();
}

// Destructor
EventContainer::~EventContainer(){
	
}

// Function to classify the event
void EventContainer::EventClassifier(Utility::FileTypeEnums type){

	// --- MC classification --
	if (type == Utility::kMC || type == Utility::kDetVar || type == Utility::kIntrinsic || type == Utility::kCCNCPiZero || type == Utility::kFakeData) {

		// check whether vertex is within fiducial volume
		bool isInFV = _utility.inFV(true_nu_vtx_x, true_nu_vtx_y, true_nu_vtx_z);

		// Outside of FV
		if (!isInFV) {
			// classify as outFV
			classification = Utility::kOutFV;
			category_ = static_cast<int>(classification);
			return;
		}

		// Inside of FV
		else {
			// Charged Current
			if (ccnc == 0) {
				// numu
				if (std::abs(nu_pdg) == 14) {
					// count muons and pions above threshold
					int nmuon_threshold = 0;
					int npion_threshold = 0;
					// count other mesons present
					int nkaon = 0;
					int nrho = 0;
					int nomega = 0;
					int nphi = 0;
						
					for (int i = 0; i < mc_pdg_v->size(); i++) {
						
						// momentum
						double mc_momentum = std::sqrt( pow(mc_px_v->at(i),2) + 
													 pow(mc_py_v->at(i),2) + 
													 pow(mc_pz_v->at(i),2) );

						if (mc_pdg_v->at(i) == 13 && mc_momentum > 0.1) nmuon_threshold++; // 100 MeV/c threshold
						if ((mc_pdg_v->at(i) == 211 || mc_pdg_v->at(i) == -211) && mc_momentum > 0.1) npion_threshold++; // 100 MeV/c threshold						
						if ((mc_pdg_v->at(i) == 321 || mc_pdg_v->at(i) == -321 || mc_pdg_v->at(i) == 311)) nkaon++;
						if ((mc_pdg_v->at(i) == 213 || mc_pdg_v->at(i) == -213 || mc_pdg_v->at(i) == 113)) nrho++;
						if (mc_pdg_v->at(i) == 223) nomega++;
						if (mc_pdg_v->at(i) == 333) nphi++;   	
					}

					if (nmuon_threshold == 1 && npion_threshold == 0 && npi0 == 0 && nkaon == 0 && nrho == 0 && nomega == 0 && nphi == 0 && neta == 0) {
						// classify as CC numu 0pi
						classification = Utility::kCCNumu0pi;
						category_ = static_cast<int>(classification);
						return;
					}
					else if (nmuon_threshold == 1 && npion_threshold == 1 && npi0 == 0 && nkaon == 0 && nrho == 0 && nomega == 0 && nphi == 0 && neta == 0) {
						// classify as CC numu 1pi
						classification = Utility::kCCNumu1pi;
						category_ = static_cast<int>(classification);
						return;
					}
					else if (nmuon_threshold == 1 && npion_threshold >= 1 && npi0 == 0 && nkaon == 0 && nrho == 0 && nomega == 0 && nphi == 0 && neta == 0) {
						// classify as CC numu Npi
						classification = Utility::kCCNumuNpi;
						category_ = static_cast<int>(classification);
						return;
					}
					else {
						// classify as CC numu other
						classification = Utility::kCCNumuOther;
						category_ = static_cast<int>(classification);
						return;
					}
				}				
				// nue
				if (std::abs(nu_pdg) == 12) {
					classification = Utility::kCCNue;
					category_ = static_cast<int>(classification);
					return;						
				}
			}
			// Neutral Current
			else {
				// NC pi+ (signal)
				// count pions above threshold
				int npion_threshold = 0;
				// count other mesons present
				int nkaon = 0;
				int nrho = 0;
				int nomega = 0;
				int nphi = 0;
						
				for (int i = 0; i < mc_pdg_v->size(); i++) {
					
					// momentum
					double mc_momentum = std::sqrt( pow(mc_px_v->at(i),2) + 
													pow(mc_py_v->at(i),2) + 
													pow(mc_pz_v->at(i),2) );

					if ((mc_pdg_v->at(i) == 211 || mc_pdg_v->at(i) == -211) && mc_momentum > 0.1) npion_threshold++; // 100 MeV/c threshold						
					if ((mc_pdg_v->at(i) == 321 || mc_pdg_v->at(i) == -321 || mc_pdg_v->at(i) == 311)) nkaon++;
					if ((mc_pdg_v->at(i) == 213 || mc_pdg_v->at(i) == -213 || mc_pdg_v->at(i) == 113)) nrho++;
					if (mc_pdg_v->at(i) == 223) nomega++;
					if (mc_pdg_v->at(i) == 333) nphi++;   	
				}

				if (npion_threshold == 1 && npi0 == 0 && nkaon == 0 && nrho == 0 && nomega == 0 && nphi == 0 && neta == 0) {
					// classify as NC 1pi
					classification = Utility::kNC1pi;
					category_ = static_cast<int>(classification);
					mc_is_signal_ = true;

					// determine pion reinteraction type 
					// loop of truth particles, WC tree
					std::vector<int> daughtersPDG;
					for (unsigned int wc_truth_idx = 0; wc_truth_idx < wc_truth_Ntrack; wc_truth_idx++) {
						
						// primary particles
						if (wc_truth_mother[wc_truth_idx] != 0) continue; 

						// find primary pion
						if (std::abs(wc_truth_pdg[wc_truth_idx]) != 211) continue;

						// check momentum
						double momentum = std::sqrt( 	std::pow(wc_truth_startMomentum[wc_truth_idx][0],2) + 
											std::pow(wc_truth_startMomentum[wc_truth_idx][1],2) + 
											std::pow(wc_truth_startMomentum[wc_truth_idx][2],2) );
						
						if (momentum <= 0.1) continue;

						int pionID = wc_truth_id[wc_truth_idx];

						// get daughters of pion
						PrintDaughters(pionID, daughtersPDG);

						break;
					}

					// classify based on daughters
					bool hasMuon = false;
					bool hasPi0 = false;
					int nProtons = 0;
					for (int i = 0; i < daughtersPDG.size(); i++) {
						if (std::abs(daughtersPDG[i]) == 13) hasMuon = true;
						if (daughtersPDG[i] == 111) hasPi0 = true;
						if (daughtersPDG[i] == 2212) nProtons++;
						
					}

					truePionHasMuon = hasMuon;
					truePionHasPi0 = hasPi0;
					truePionNProtons = nProtons;

					return;		
				}			
				else {
					// classify as NC other
					classification = Utility::kNCOther;
					category_ = static_cast<int>(classification);
					return;
				}				
			}
		}
	}
	// Beam Off Classification
	else if (type == Utility::kEXT) {
		classification = Utility::kBeamOff;
		category_ = static_cast<int>(classification);
		return;
	}
	// Dirt
	else if (type == Utility::kDirt) {
		classification = Utility::kOutOfCryo;
		category_ = static_cast<int>(classification);
		return;
	}
	// Data
	else if (type == Utility::kData) {
		classification = Utility::kBeamOn;
		category_ = static_cast<int>(classification);
		return;
	}

	// Catch any events without classification
	// Shouldn't be any of these! Will cause issues.
	std::cout << "Warning: Unknown Event" << std::endl;
	std::cout << "Type: " << type << std::endl;
	exit(10);
	return;
}

// Function to get event classification
Utility::ClassificationEnums EventContainer::getEventClassification(Utility::FileTypeEnums type) {
	EventClassifier(type);
	return classification;
}

// Function to set interaction type enum
void EventContainer::InteractionClassifier(Utility::FileTypeEnums type){
	// --- MC classification --
	if (type == Utility::kMC || type == Utility::kDetVar || type == Utility::kIntrinsic || type == Utility::kCCNCPiZero || type == Utility::kFakeData) {
		
		if (interaction == 0) {
			interactionMode = Utility::kQE;
		}
		else if (interaction == 1) {
			interactionMode = Utility::kRES;
		}
		else if (interaction == 2) {
			interactionMode = Utility::kDIS;
		}
		else if (interaction == 3) {
			interactionMode = Utility::kCOH;
		}
		else if (interaction == 10) {
			interactionMode = Utility::kMEC;
		}
		else {
			interactionMode = Utility::kOther;
		}
	}
	else {
		// undefined for data
		interactionMode = Utility::kOther;
	}
}

// Function to set particle classification enum
void EventContainer::ParticleClassifier(Utility::FileTypeEnums type){
	
	if (type == Utility::kMC || type == Utility::kIntrinsic || type == Utility::kDirt) {

		// loop through all reconstructed particles
		for (unsigned int i_trk = 0; i_trk < trk_sce_start_x_v->size(); i_trk++) {	
			// get the PDG code of the particle
			int pdg = backtracked_pdg_v->at(i_trk);

			// momentum
			double momentum = std::sqrt( pow(backtracked_px_v->at(i_trk),2) + 
										 pow(backtracked_px_v->at(i_trk),2) + 
										 pow(backtracked_px_v->at(i_trk),2) );


			// classify the particle
			if (std::abs(pdg) == 13) {
				pfp_classification_v[i_trk] = Utility::kMuon;
			}
			else if (pdg == 211 || pdg == -211) {
				// check threshold
				if (momentum > 0.1) { // 100 MeV/c threshold
					pfp_classification_v[i_trk] = Utility::kPionDecay;
				}
				else {
					pfp_classification_v[i_trk] = Utility::kPionBelowThreshold;
				}
			}
			else if (pdg == 2212) {
				pfp_classification_v[i_trk] = Utility::kProton;
			}
			else {
				pfp_classification_v[i_trk] = Utility::kCosmic;
			}
		}
	}
	// otherwise leave as default "kOther"
	else {
	}
}

// Function to set particle classification enum for WC 
void EventContainer::ParticleClassifierWC(Utility::FileTypeEnums type){
	
	if (type == Utility::kMC || type == Utility::kIntrinsic || type == Utility::kDirt) {

		// loop through all reconstructed particles
		for (unsigned int wc_idx = 0; wc_idx < wc_reco_Ntrack; wc_idx++) {
		
			// get the PDG code of the particle
			int pdg = wc_reco_truthMatch_pdg[wc_idx];

			// get the backtracker truth match index
			int backtracked_id = wc_reco_truthMatch_id[wc_idx];

			// find correct particle to calculate truth momentum
			double momentum = 0.1;
			for (unsigned int wc_truth_idx = 0; wc_truth_idx < wc_truth_Ntrack; wc_truth_idx++) {
				
				if (wc_truth_id[wc_truth_idx] != backtracked_id) continue;
				
				// set momentum
				momentum = std::sqrt( 	pow(wc_truth_startMomentum[wc_truth_idx][0],2) + 
										pow(wc_truth_startMomentum[wc_truth_idx][1],2) + 
										pow(wc_truth_startMomentum[wc_truth_idx][2],2) );
			}

			// classify the particle
			if (std::abs(pdg) == 13) {
				wc_particle_classification_v[wc_idx] = Utility::kMuon;
			}
			else if (pdg == 211 || pdg == -211) {
				// check threshold
				if (momentum > 0.1) { // 100 MeV/c threshold
					wc_particle_classification_v[wc_idx] = Utility::kPionDecay;

					// check creation process
					//std::cout << "Pion process = " << wc_truth_process->at(wc_idx) << std::endl;

				}
				else {
					wc_particle_classification_v[wc_idx] = Utility::kPionBelowThreshold;
				}
			}
			else if (pdg == 2212) {
				wc_particle_classification_v[wc_idx] = Utility::kProton;
			}
			else {
				wc_particle_classification_v[wc_idx] = Utility::kCosmic;
			}
		}
	}
	// otherwise leave as default "kOther"
	else {
	}
}

// Function to set particle classification enum for Lantern 
void EventContainer::ParticleClassifierLantern(Utility::FileTypeEnums type) {
	
	if (type == Utility::kMC || type == Utility::kIntrinsic || type == Utility::kDirt) {

		// loop through all lantern tracks 
		for (unsigned int track_idx = 0; track_idx < lantern_nTracks; track_idx++) {
		
			// get the truth matched PDG code of the particle
			int pdg = lantern_trackTruePID[track_idx];

			// classify the particle
			if (std::abs(pdg) == 13) {
				lantern_particle_classification_v[track_idx] = Utility::kMuon;
			}
			else if (pdg == 211 || pdg == -211) {	
		
				// found pion, classify interaction modes based on best-guess truth particle
				// works for single-pion events, will confuse multi-pion events
				// loop of truth particles, WC tree
				std::vector<int> daughtersPDG;
				for (unsigned int wc_truth_idx = 0; wc_truth_idx < wc_truth_Ntrack; wc_truth_idx++) {
					
					// primary particles
					if (wc_truth_mother[wc_truth_idx] != 0) continue; 

					// find primary pion
					if (std::abs(wc_truth_pdg[wc_truth_idx]) != 211) continue;

					// check momentum
					double momentum = std::sqrt( 	std::pow(wc_truth_startMomentum[wc_truth_idx][0],2) + 
										std::pow(wc_truth_startMomentum[wc_truth_idx][1],2) + 
										std::pow(wc_truth_startMomentum[wc_truth_idx][2],2) );
					
					if (momentum <= 0.1) continue;

					// have found primary pion
					//std::cout << "Truth Selected Pion PDG " << wc_truth_pdg[wc_truth_idx] << " -- Start Process: " << wc_truth_process->at(wc_truth_idx) 
					//		<< ", End Process: " << wc_truth_endprocess->at(wc_truth_idx) 
					//		<< std::endl;

					int pionID = wc_truth_id[wc_truth_idx];

					// get daughters of pion
					PrintDaughters(pionID, daughtersPDG);

					break;
				}

				// print out full list of daughters
				/*
				std::cout << "Full Daughters List: ";
				for (unsigned int i = 0; i < daughtersPDG.size(); i++) {
					std::cout << daughtersPDG[i] << ", ";
				}
				std::cout << std::endl;
				*/

				// classify based on daughters
				bool hasMuon = false;
				bool hasPi0 = false;
				int nProtons = 0;
				for (int i = 0; i < daughtersPDG.size(); i++) {
					if (std::abs(daughtersPDG[i]) == 13) hasMuon = true;
					if (daughtersPDG[i] == 111) hasPi0 = true;
					if (daughtersPDG[i] == 2212) nProtons++;
					
				}
				if (hasMuon) {
					//std::cout << "Classification: Pion Decay" << std::endl;
					lantern_particle_classification_v[track_idx] = Utility::kPionDecay;
				}
				else if (hasPi0) {
					//std::cout << "Classification: Pion Charge Exchange" << std::endl;
					lantern_particle_classification_v[track_idx] = Utility::kPionChargeExchange;
				}
				else {
					//std::cout << "Classification: Pion Absorption" << std::endl;
					if (nProtons == 0) lantern_particle_classification_v[track_idx] = Utility::kPionAbsorption0p;
					else lantern_particle_classification_v[track_idx] = Utility::kPionAbsorptionNp;
				} 
				
			}
			else if (pdg == 2212) {
				lantern_particle_classification_v[track_idx] = Utility::kProton;
			}
			else {
				lantern_particle_classification_v[track_idx] = Utility::kCosmic;
			}
		}
	}
	// otherwise leave as default "kOther"
	else {
	}
}

// Function to calculate CV event weight
void EventContainer::calculateCVEventWeight(Utility::FileTypeEnums type, Utility::RunPeriodEnums runPeriod) {

	float weight = 1.0;

	// check weights are sensible
	weightSpline = checkWeight(weightSpline);
	weightTune = checkWeight(weightTune);

	// overlay MC events
	if (type == Utility::kMC || type == Utility::kDetVar || type == Utility::kIntrinsic || type == Utility::kCCNCPiZero || type == Utility::kFakeData) {			
		
		// CV weight
		weight = weightSpline * weightTune;

		// set normalisation weight
		normalisation_weight = 1;

		// apply FSI reweighting for pions in signal events
		double fsi_weight = GetFSIWeight();
		weight *= fsi_weight;

	}
	
	// dirt events
	if (type == Utility::kDirt) {

		weight = weightSpline * weightTune;

		normalisation_weight = 1;
		
	}

	// beam off events
	if (type == Utility::kEXT) {
		weight = 1;

		// set normalisation weight
	    normalisation_weight = 1;
	}

	// beam on
	if (type == Utility::kData) {
		weight = 1;

		// set normalisation weight
	    normalisation_weight = 1;
	}	

	// check weight is sensible
	weight = checkWeight(weight);

	// set weight in event
	weight_cv = weight;

}

// Function to check event weight is sensible
float EventContainer::checkWeight(float weight) {


	// infinite or nan weight
	if (!_utility.isNumber(weight)) {
		//std::cout << "Warning: infinite/nan event weight" << std::endl;
		weight = 1.0;
	}

	// overly large weight
	else if (weight > 15.0) {
		//std::cout << "Warning: overly large event weight, " << weight << std::endl;
		weight = 1.0;
	}

	// negative weight
	else if (weight < 0.0) {
		//std::cout << "Warning: negative event weight" << std::endl;
		weight = 0.0;
	}

	// approximately zero weight
	else if (weight > -1.0e-4 && weight < 1e-4) {
		weight = 0.0;
	}

	return weight;
}

// Function to populate derived variables
void EventContainer::populateDerivedVariables(Utility::FileTypeEnums type, Utility::RunPeriodEnums runPeriod){

	// Lantern and WC vertex comparison
	vertex_discrepancy = std::sqrt(
		std::pow(lantern_vtxX - wc_reco_nuvtxX, 2) +
		std::pow(lantern_vtxY - wc_reco_nuvtxY, 2) +
		std::pow(lantern_vtxZ - wc_reco_nuvtxZ, 2)
	);

	// reset Lantern selection variables 
	sel_NC1pi_ = false;
	sel_passInitialSelection_ = false;
	sel_passMuPiLLR_ = false;
	sel_LanternPID_llr_mu_pi_ = -2;
	sel_LanternPID_llr_pr_pi_ = -2;

	// fraction of hits associated with tracks or showers
	associated_hits_fraction = ((float)trk_hits_y_tot + (float)shr_hits_y_tot) / (float)total_hits_y;
	
	// number of showers above threshold
	n_showers_alt = all_shr_energies->size();
	n_showers_above_threshold = 0;

	//std::cout << "Number of all_shr_energies: " << all_shr_energies->size() << std::endl;

	for (unsigned int i = 0; i < all_shr_energies->size(); i++) {
		if (all_shr_energies->at(i) > 0.07) n_showers_above_threshold++; // 70 MeV, michels
	}

	// Shower dE/dx
	GetdEdxMax(false);	// first 2cm only 
	GetdEdxMax(true); 	// with 10mm gap from start
 
	// Primary track information
	// if information primary track present -- vectors are all the same size by construction
	if (trk_id-1 < trk_sce_start_x_v->size() && trk_id != 0) {

		trk_llr_pid_score = trk_llr_pid_score_v->at(trk_id-1); 	// LLR PID Score
		trk_daughters = pfp_trk_daughters_v->at(trk_id-1); 		// Track daughters
		trk_start_x = trk_start_x_v->at(trk_id-1);				// Track start, without SCE to allow comparison with shower start
		trk_start_y = trk_start_y_v->at(trk_id-1);				// Track start, without SCE to allow comparison with shower start
		trk_start_z = trk_start_z_v->at(trk_id-1);				// Track start, without SCE to allow comparison with shower start
		trk_end_x = trk_end_x_v->at(trk_id-1); 					// Track end, without SCE to allow comparison with shower start
        trk_end_y = trk_end_y_v->at(trk_id-1); 					// Track end, without SCE to allow comparison with shower start
        trk_end_z = trk_end_z_v->at(trk_id-1); 					// Track end, without SCE to allow comparison with shower start
		trk_sce_end_x = trk_sce_end_x_v->at(trk_id-1); 			// Track end
        trk_sce_end_y = trk_sce_end_y_v->at(trk_id-1); 			// Track end
        trk_sce_end_z = trk_sce_end_z_v->at(trk_id-1); 			// Track end
        trk_dir_x = trk_dir_x_v->at(trk_id-1);
        trk_dir_y = trk_dir_y_v->at(trk_id-1);  	
        trk_dir_z = trk_dir_z_v->at(trk_id-1);
        trk_calo_energy = trk_calo_energy_y_v->at(trk_id-1);
        trk_energy_proton = trk_energy_proton_v->at(trk_id-1);
        trk_energy_muon = trk_energy_muon_v->at(trk_id-1);
        trk_end_spacepoints = trk_end_spacepoints_v->at(trk_id-1);
        trk_planehits_U = pfnplanehits_U->at(trk_id-1);
        trk_planehits_V = pfnplanehits_U->at(trk_id-1);
        trk_planehits_Y = pfnplanehits_U->at(trk_id-1);
        trk_pfpgeneration = pfp_generation_v->at(trk_id-1);
	}
	else {
		trk_llr_pid_score = 9999;
		trk_daughters = 9999;
		trk_start_x = 9999;
		trk_start_y = 9999;
		trk_start_z = 9999;
		trk_end_x = 9999;
		trk_end_y = 9999;
		trk_end_z = 9999;
		trk_sce_end_x = 9999;
		trk_sce_end_y = 9999;
		trk_sce_end_z = 9999;
		trk_dir_x = 9999;
		trk_dir_y = 9999;
		trk_dir_z = 9999;
		trk_calo_energy = 9999;
		trk_energy_proton = 0;
        trk_energy_muon = 0;
        trk_end_spacepoints = 9999;
        trk_planehits_U = 0;
        trk_planehits_V = 0;
        trk_planehits_Y = 0;
        trk_pfpgeneration = 0;
	}

	// Secondary track information
	if (trk2_id-1 < trk_sce_start_x_v->size() && trk2_id != 0) {
		
		trk2_bragg_p = trk_bragg_p_v->at(trk2_id-1); 				// Bragg Proton
		trk2_bragg_mu = trk_bragg_mu_v->at(trk2_id-1); 				// Bragg Muon
		trk2_bragg_mip = trk_bragg_mip_v->at(trk2_id-1); 			// Bragg MIP
		trk2_llr_pid_score = trk_llr_pid_score_v->at(trk2_id-1); 	// LLR PID Score
		trk2_len = trk_len_v->at(trk2_id-1); 						// Track length
		trk2_score = trk_score_v->at(trk2_id-1); 					// Track score
		trk2_distance = trk_distance_v->at(trk2_id-1);				// Track distance
		trk2_daughters = pfp_trk_daughters_v->at(trk2_id-1); 		// Track daughters
		trk2_start_x = trk_start_x_v->at(trk2_id-1);				// Track start, without SCE to allow comparison with shower start
		trk2_start_y = trk_start_y_v->at(trk2_id-1);				// Track start, without SCE to allow comparison with shower start
		trk2_start_z = trk_start_z_v->at(trk2_id-1);				// Track start, without SCE to allow comparison with shower start
		trk2_end_x = trk_end_x_v->at(trk2_id-1); 					// Track end, without SCE to allow comparison with shower start
        trk2_end_y = trk_end_y_v->at(trk2_id-1); 					// Track end, without SCE to allow comparison with shower start
        trk2_end_z = trk_end_z_v->at(trk2_id-1); 					// Track end, without SCE to allow comparison with shower start
		trk2_sce_end_x = trk_sce_end_x_v->at(trk2_id-1); 			// Track end
        trk2_sce_end_y = trk_sce_end_y_v->at(trk2_id-1); 			// Track end
        trk2_sce_end_z = trk_sce_end_z_v->at(trk2_id-1); 			// Track end
        trk2_dir_x = trk_dir_x_v->at(trk2_id-1);
        trk2_dir_y = trk_dir_y_v->at(trk2_id-1);  	
        trk2_dir_z = trk_dir_z_v->at(trk2_id-1);
        trk2_tkfit_dedx_u = shr_tkfit_dedx_u_v->at(trk2_id-1);
        trk2_tkfit_dedx_v = shr_tkfit_dedx_v_v->at(trk2_id-1); 
        trk2_tkfit_dedx_y = shr_tkfit_dedx_y_v->at(trk2_id-1);
        trk2_tkfit_nhits_u = shr_tkfit_dedx_nhits_u_v->at(trk2_id-1);
        trk2_tkfit_nhits_v = shr_tkfit_dedx_nhits_v_v->at(trk2_id-1);
        trk2_tkfit_nhits_y = shr_tkfit_dedx_nhits_y_v->at(trk2_id-1);
        trk2subclusters0 = pfpplanesubclusters_U_v->at(trk2_id-1);
        trk2subclusters1 = pfpplanesubclusters_V_v->at(trk2_id-1);
        trk2subclusters2 = pfpplanesubclusters_Y_v->at(trk2_id-1);
		trk2_bragg_pion = trk_bragg_pion_v->at(trk2_id-1); 			// Bragg Pion
		trk2_calo_energy = trk_calo_energy_y_v->at(trk2_id-1);
		trk2_energy_proton = trk_energy_proton_v->at(trk2_id-1);
        trk2_energy_muon = trk_energy_muon_v->at(trk2_id-1);
        trk2_end_spacepoints = trk_end_spacepoints_v->at(trk2_id-1);
		if (type != Utility::kEXT && type != Utility::kData) trk2_bkt_pdg = backtracked_pdg_v->at(trk2_id-1);			// Backtracked PDG
		else trk2_bkt_pdg = 9999;
		trk2_planehits_U = pfnplanehits_U->at(trk2_id-1);
        trk2_planehits_V = pfnplanehits_U->at(trk2_id-1);
        trk2_planehits_Y = pfnplanehits_U->at(trk2_id-1);
        trk2_pfpgeneration = pfp_generation_v->at(trk2_id-1);
	}
	else {
		trk2_len = 0;
		trk2_score = 9999;
		trk2_distance =9999;
		trk2_bragg_p = 9999;
		trk2_bragg_mu = 9999;
		trk2_bragg_pion = 9999;
		trk2_bragg_mip = 9999;
		trk2_llr_pid_score = 9999;
		trk2_daughters = 9999;
		trk2_bkt_pdg = 9999;
		trk2_start_x = 9999;
		trk2_start_y = 9999;
		trk2_start_z = 9999;
		trk2_end_x = 9999;
		trk2_end_y = 9999;
		trk2_end_z = 9999;
		trk2_sce_end_x = 9999;
		trk2_sce_end_y = 9999;
		trk2_sce_end_z = 9999;
		trk2_dir_x = 9999;
		trk2_dir_y = 9999;
		trk2_dir_z = 9999;
		trk2_tkfit_dedx_u = 9999;
		trk2_tkfit_dedx_v = 9999;
		trk2_tkfit_dedx_y = 9999;
		trk2_tkfit_nhits_u = 0;
		trk2_tkfit_nhits_v = 0;
		trk2_tkfit_nhits_y = 0;
		trk2subclusters0 = 9999;
		trk2subclusters1 = 9999;
		trk2subclusters2 = 9999;
		trk2_calo_energy = 9999;
		trk2_energy_proton = 0;
        trk2_energy_muon = 0;
        trk2_end_spacepoints = 9999;
        trk2_planehits_U = 0;
        trk2_planehits_V = 0;
        trk2_planehits_Y = 0;
        trk2_pfpgeneration = 0;
	}

	// tertiary track information
	if (trk3_id-1 < trk_sce_start_x_v->size() && trk3_id != 0) {

		trk3_bragg_p = trk_bragg_p_v->at(trk3_id-1); 				// Bragg Proton
		trk3_bragg_mu = trk_bragg_mu_v->at(trk3_id-1); 				// Bragg Muon
		trk3_bragg_mip = trk_bragg_mip_v->at(trk3_id-1); 			// Bragg MIP
		trk3_llr_pid_score = trk_llr_pid_score_v->at(trk3_id-1); 	// LLR PID Score
		trk3_len = trk_len_v->at(trk3_id-1); 						// Track length
		trk3_score = trk_score_v->at(trk3_id-1); 					// Track score
		trk3_distance = trk_distance_v->at(trk3_id-1);				// Track distance
		trk3_daughters = pfp_trk_daughters_v->at(trk3_id-1); 		// Track daughters
		trk3_end_x = trk_end_x_v->at(trk3_id-1); 					// Track end, without SCE to allow comparison with shower start
        trk3_end_y = trk_end_y_v->at(trk3_id-1); 					// Track end, without SCE to allow comparison with shower start
        trk3_end_z = trk_end_z_v->at(trk3_id-1); 					// Track end, without SCE to allow comparison with shower start
		trk3_sce_end_x = trk_sce_end_x_v->at(trk3_id-1); 			// Track end
        trk3_sce_end_y = trk_sce_end_y_v->at(trk3_id-1); 			// Track end
        trk3_sce_end_z = trk_sce_end_z_v->at(trk3_id-1); 			// Track end
        trk3_dir_x = trk_dir_x_v->at(trk3_id-1);
        trk3_dir_y = trk_dir_y_v->at(trk3_id-1);  	
        trk3_dir_z = trk_dir_z_v->at(trk3_id-1);
		trk3_bragg_pion = trk_bragg_pion_v->at(trk3_id-1);
		trk3_calo_energy = trk_calo_energy_y_v->at(trk3_id-1);
		trk3_energy_proton = trk_energy_proton_v->at(trk3_id-1);
        trk3_energy_muon = trk_energy_muon_v->at(trk3_id-1); 									
        trk3_end_spacepoints = trk_end_spacepoints_v->at(trk3_id-1);
		if (type != Utility::kEXT && type != Utility::kData) trk3_bkt_pdg = backtracked_pdg_v->at(trk3_id-1);			// Backtracked PDG
		else trk3_bkt_pdg = 9999;
		trk3_planehits_U = pfnplanehits_U->at(trk3_id-1);
        trk3_planehits_V = pfnplanehits_U->at(trk3_id-1);
        trk3_planehits_Y = pfnplanehits_U->at(trk3_id-1);
        trk3_pfpgeneration = pfp_generation_v->at(trk3_id-1);
	}
	else {
		trk3_len = 0;
		trk3_score = 9999;
		trk3_distance =9999;
		trk3_bragg_p = 9999;
		trk3_bragg_mu = 9999;
		trk3_bragg_pion = 9999;
		trk3_bragg_mip = 9999;
		trk3_llr_pid_score = 9999;
		trk3_daughters = 9999;
		trk3_bkt_pdg = 9999;
		trk3_end_x = 9999;
		trk3_end_y = 9999;
		trk3_end_z = 9999;
		trk3_sce_end_x = 9999;
		trk3_sce_end_y = 9999;
		trk3_sce_end_z = 9999;
		trk3_calo_energy = 0;
		trk3_energy_proton = 0;
        trk3_energy_muon = 0;
        trk3_end_spacepoints = 9999;
        trk3_planehits_U = 0;
        trk3_planehits_V = 0;
        trk3_planehits_Y = 0;
        trk3_pfpgeneration = 0;
        trk3_dir_x = 9999;
		trk3_dir_y = 9999;
		trk3_dir_z = 9999;
	}


	// track PID variables, best plane 
	// primary track
	if (trk_id-1 < trk_sce_start_x_v->size() && trk_id != 0) {
		trk_dEdx_trunk_max = GetTrackTrunkdEdxBestPlane(trk_id);
		trk_bragg_pion_max = GetTrackBraggPionBestPlane(trk_id);
		trk_bragg_mip_max = GetTrackBraggMIPBestPlane(trk_id);
	}
	else {
		trk_dEdx_trunk_max = 9999;
		trk_bragg_pion_max = 9999;
		trk_bragg_mip_max = 9999;
	}
	// secondary track
	if (trk2_id-1 < trk_sce_start_x_v->size() && trk2_id != 0) {
		trk2_dEdx_trunk_max = GetTrackTrunkdEdxBestPlane(trk2_id);
		trk2_bragg_pion_max = GetTrackBraggPionBestPlane(trk2_id);
		trk2_bragg_mip_max = GetTrackBraggMIPBestPlane(trk2_id);
	}
	else {
		trk2_dEdx_trunk_max = 9999;
		trk2_bragg_pion_max = 9999;
		trk2_bragg_mip_max = 9999;
	}
	// tertiary track
	if (trk3_id-1 < trk_sce_start_x_v->size() && trk3_id != 0) {
		trk3_dEdx_trunk_max = GetTrackTrunkdEdxBestPlane(trk3_id);
		trk3_bragg_pion_max = GetTrackBraggPionBestPlane(trk3_id);
		trk3_bragg_mip_max = GetTrackBraggMIPBestPlane(trk3_id);
	}
	else {
		trk3_dEdx_trunk_max = 9999;
		trk3_bragg_pion_max = 9999;
		trk3_bragg_mip_max = 9999;
	}

	// Shower energy correction
	shr_energy_cali = shr_energy_cali / 0.83;
	shr_energy_tot_cali = shr_energy_tot_cali / 0.83;

	// Primary shower energy fraction
    shr_energyFraction = shr_energy_cali / shr_energy_tot_cali;

    // primary shower subclusters
    shrsubclusters = shrsubclusters0 + shrsubclusters1 + shrsubclusters2;

    // primary shower back-tracked truth, set to zero when not filled
    if (shr_bkt_E < -1 || shr_bkt_E > 100) shr_bkt_E=0; 

    // Primary shower information
    if (shr_id-1 < shr_energy_y_v->size() && shr_id != 0) {
    	shr_start_x = shr_start_x_v->at(shr_id-1);
    	shr_start_y = shr_start_y_v->at(shr_id-1);
    	shr_start_z = shr_start_z_v->at(shr_id-1);
    	shr_dir_x = shr_px_v->at(shr_id-1);
    	shr_dir_y = shr_py_v->at(shr_id-1);
    	shr_dir_z = shr_pz_v->at(shr_id-1);
    	shr_planehits_U = pfnplanehits_U->at(shr_id-1);
        shr_planehits_V = pfnplanehits_U->at(shr_id-1);
        shr_planehits_Y = pfnplanehits_U->at(shr_id-1);
        shr_pfpgeneration = pfp_generation_v->at(shr_id-1);	
    } 
    else {
    	shr_start_x = 9999;
    	shr_start_y = 9999;
    	shr_start_z = 9999;
    	shr_dir_x = 9999;
    	shr_dir_y = 9999;
    	shr_dir_z = 9999;
    	shr_planehits_U = 0;
    	shr_planehits_V = 0;
    	shr_planehits_Y = 0;
    	shr_pfpgeneration = 0;
    }

    // Secondary shower information
    if (shr2_id-1 < shr_energy_y_v->size() && shr2_id != 0) {
    	shr2_energy = shr_energy_y_v->at(shr2_id-1)/1000;	// shower energy, in GeV
    	shr2_score = trk_score_v->at(shr2_id-1); 			// note track-shower score, full vector contains all PFPs not just "track-like" ones
    	shr2_start_x = shr_start_x_v->at(shr2_id-1);
    	shr2_start_y = shr_start_y_v->at(shr2_id-1);
    	shr2_start_z = shr_start_z_v->at(shr2_id-1);
    	shr2_dir_x = shr_px_v->at(shr2_id-1);
    	shr2_dir_y = shr_py_v->at(shr2_id-1);
    	shr2_dir_z = shr_pz_v->at(shr2_id-1);
    	shr2_distance = shr_dist_v->at(shr2_id-1);
    	shr2moliereavg = shr_moliere_avg_v->at(shr2_id-1); 	
    	shr2subclusters = pfpplanesubclusters_U_v->at(shr2_id-1) + pfpplanesubclusters_V_v->at(shr2_id-1) + pfpplanesubclusters_Y_v->at(shr2_id-1);
    	shr2pid = trk_llr_pid_score_v->at(shr2_id-1);
    	shr2_pfpgeneration = pfp_generation_v->at(shr2_id-1);	
    } 
    else {
    	shr2_energy = 9999;
    	shr2_score = 9999;
    	shr2_dir_x = 9999;
    	shr2_dir_y = 9999;
    	shr2_dir_z = 9999;
    	shr2_start_x = 9999;
    	shr2_start_y = 9999;
    	shr2_start_z = 9999;
    	shr2_distance = 9999;
    	shr2moliereavg = 9999;
    	shr2subclusters = 9999;
    	shr2pid = 9999;
    	shr2_pfpgeneration = 0;
    }

    // Tertiary shower information
    if (shr3_id-1 < shr_start_x_v->size() && shr3_id != 0) {
    	shr3_energy = shr_energy_y_v->at(shr3_id-1)/1000;	// shower energy, in GeV
    	shr3_start_x = shr_start_x_v->at(shr3_id-1);
    	shr3_start_y = shr_start_y_v->at(shr3_id-1);
    	shr3_start_z = shr_start_z_v->at(shr3_id-1);
    	shr3subclusters = pfpplanesubclusters_U_v->at(shr3_id-1) + pfpplanesubclusters_V_v->at(shr3_id-1) + pfpplanesubclusters_Y_v->at(shr3_id-1);
    	shr3_distance = shr_dist_v->at(shr3_id-1);	
    }
    else {
    	shr3_energy = 9999;
    	shr3_start_x = 9999;
    	shr3_start_y = 9999;
    	shr3_start_z = 9999;
    	shr3subclusters = 9999;
    	shr3_distance = 9999;
    }
    

    // convert tksh angle
    tksh_angle = std::acos(tksh_angle) * 180 / 3.14159;

    // Second shower opening angle
    if (n_showers > 1 && shr_dir_x != 9999 && shr2_start_x != 9999) {

    	TVector3 showerDir(shr_dir_x, shr_dir_y, shr_dir_z);
    	TVector3 shower12StartDist(shr2_start_x - shr_start_x, shr2_start_y - shr_start_y, shr2_start_z - shr_start_z);
    	shr12_p1_dstart = showerDir.Angle(shower12StartDist) * 180 / 3.14159;
	}
	else shr12_p1_dstart = 9999;

	// Tertiary shower opening angle
	if (n_showers > 1 && shr_dir_x != 9999 && shr3_start_x != 9999) {

    	TVector3 showerDir(shr_dir_x, shr_dir_y, shr_dir_z);
    	TVector3 shower13StartDist(shr3_start_x - shr_start_x, shr3_start_y - shr_start_y, shr3_start_z - shr_start_z);
    	shr13_p1_dstart = showerDir.Angle(shower13StartDist) * 180 / 3.14159;
	}
	else shr13_p1_dstart = 9999;	

	// second shower distance from primary track
	if (n_showers_contained > 1 && n_tracks_contained > 0 && trk_start_x != 9999 & shr2_start_x != 9999) {
		tk1sh2_distance = std::sqrt( std::pow(shr2_start_x - trk_start_x, 2) + std::pow(shr2_start_y - trk_start_y, 2) + std::pow(shr2_start_z - trk_start_z, 2) );
	}
	else {
		tk1sh2_distance = 9999;
	}

	// tertiary shower distance from primary track
	if (n_showers_contained > 1 && n_tracks_contained > 0 && trk_start_x != 9999 & shr3_start_x != 9999) {
		tk1sh3_distance = std::sqrt( std::pow(shr3_start_x - trk_start_x, 2) + std::pow(shr3_start_y - trk_start_y, 2) + std::pow(shr3_start_z - trk_start_z, 2) );
	}
	else {
		tk1sh3_distance = 9999;
	}

	// distance between primary track and secondary track
	if (n_tracks_contained > 1 && trk_start_x != 9999 && trk2_start_x != 9999) {
		tk1tk2_distance = std::sqrt( std::pow(trk2_start_x - trk_start_x, 2) + std::pow(trk2_start_y - trk_start_y, 2) + std::pow(trk2_start_z - trk_start_z, 2) );
	}
	else {
		tk1tk2_distance = 9999;
	}

	// primary shower and second track opening angle
	if (n_tracks_contained > 1 && n_showers_contained > 0 && shr_dir_x != 9999 && trk2_dir_x != 9999) {

		TVector3 showerDir(shr_dir_x, shr_dir_y, shr_dir_z);
		TVector3 trk2Dir(trk2_dir_x, trk2_dir_y, trk2_dir_z);
		tk2sh1_angle = trk2Dir.Angle(showerDir) * 180 / 3.14159;
	}
	
	// primary track and second track opening angle
	if (n_tracks_contained > 1 && trk_dir_x != 9999 && trk2_dir_x != 9999) {

		TVector3 trkDir(trk_dir_x, trk_dir_x, trk_dir_x);
		TVector3 trk2Dir(trk2_dir_x, trk2_dir_y, trk2_dir_z);
		tk1tk2_angle = trkDir.Angle(trk2Dir) * 180 / 3.14159;
	}
	

    // Second shower cluster opening angle
    GetSecondShowerClusterAngleDiff();	

    // Shower track fit
    if (shr_id-1 < trk_sce_start_x_v->size() && shr_id != 0) {
    	shr_trk_sce_start_x = trk_sce_start_x_v->at(shr_id-1);
    	shr_trk_sce_start_y = trk_sce_start_y_v->at(shr_id-1);
    	shr_trk_sce_start_z = trk_sce_start_z_v->at(shr_id-1);
    	shr_trk_sce_end_x = trk_sce_end_x_v->at(shr_id-1);
    	shr_trk_sce_end_y = trk_sce_end_y_v->at(shr_id-1);
    	shr_trk_sce_end_z = trk_sce_end_z_v->at(shr_id-1);

    	shr_trk_len = std::sqrt( std::pow(shr_trk_sce_end_x - shr_trk_sce_start_x, 2) + std::pow(shr_trk_sce_end_y - shr_trk_sce_start_y, 2) + std::pow(shr_trk_sce_end_z - shr_trk_sce_start_z, 2) );
    }
    else {
    	shr_trk_sce_start_x = 9999;
    	shr_trk_sce_start_y = 9999;
    	shr_trk_sce_start_z = 9999;
    	shr_trk_sce_end_x = 9999;
    	shr_trk_sce_end_y = 9999;
    	shr_trk_sce_end_z = 9999;
    	shr_trk_len = 9999;
    }

    // Pion-hypothesis track energies (A. Smith)
    if (trk_len > 0) trk_momentum_pion = CalculatePionMomentumRange(trk_len);
    else trk_momentum_pion = 0;
    if (trk2_len > 0) trk2_momentum_pion = CalculatePionMomentumRange(trk2_len);
    else trk2_momentum_pion = 0;
    if (trk3_len > 0) trk3_momentum_pion = CalculatePionMomentumRange(trk3_len);
    else trk3_momentum_pion = 0;

    // second shower track proximity
    shr2_trackEndProximity = GetShowerTrackEndProximity(shr2_id);
    shr3_trackEndProximity = GetShowerTrackEndProximity(shr3_id);

    // truth pion momentum
    pion_px = 9999;
    pion_py = 9999;
    pion_pz = 9999;

	//std::cout << "Number of truth particles: " << mc_pdg_v->size() << std::endl;

    // loop through truth particles and find pion
    for (int i = 0; i < mc_pdg_v->size(); i++) {
    	if (mc_pdg_v->at(i) == 211 || mc_pdg_v->at(i) == -211) {
    		// check energy threshold, need to get correct pion if ones below threshold
    		if (mc_E_v->at(i) > 0.04 + (139.57039/1000)) {
    			pion_px = mc_px_v->at(i);
    			pion_py = mc_py_v->at(i); 
    			pion_pz = mc_pz_v->at(i); 

                pion_p = std::sqrt ( std::pow(pion_px, 2) + std::pow(pion_py, 2) + std::pow(pion_pz, 2));

    			break; 
    		}
    	}
    }

	//std::cout << "Here 1a" << std::endl;

    // track - shower opening angles
    trk_shr_opening_angle = 9999;
    trk2_shr_opening_angle = 9999;
    trk3_shr_opening_angle = 9999;

    if (n_showers_contained >= 1) {
	    TVector3 shr_dir(shr_dir_x, shr_dir_y, shr_dir_z); // Shower direction
	    shr_dir.Unit();

	    if (trk_id != 0) {
		    TVector3 trk_dir(trk_dir_x, trk_dir_y, trk_dir_z);	
		    trk_dir.Unit();
		    trk_shr_opening_angle = trk_dir.Angle(shr_dir) * 180 / 3.14159;
		}
		if (trk2_id != 0) {
			TVector3 trk2_dir(trk2_dir_x, trk2_dir_y, trk2_dir_z);	
		    trk2_dir.Unit();
		    trk2_shr_opening_angle = trk2_dir.Angle(shr_dir) * 180 / 3.14159;
		}
		if (trk3_id != 0) {
			TVector3 trk3_dir(trk3_dir_x, trk3_dir_y, trk3_dir_z);	
		    trk3_dir.Unit();
		    trk3_shr_opening_angle = trk3_dir.Angle(shr_dir) * 180 / 3.14159;
		}
	}

    // re-define truth nproton with alternate threshold
    // loop through truth particles and find protons
    nprotonalternate = 0; 
    for (int i = 0; i < mc_pdg_v->size(); i++) {
    	if (mc_pdg_v->at(i) == 2212) {
    		// check energy threshold
    		if (mc_E_v->at(i) > 0.04 + (938.27208816/1000)) { // 50 MeV
    			nprotonalternate++;
    		}
    	}
    }
	//std::cout << "Here 1b" << std::endl;


    // initialise reconstruction failure condition variables
    hasSplitPrimaryShower = false;
    hasSpuriousLeadingTrack = false;
    hasSplitTrackShower = false;
    hasSecondShowerProton = false;

	// Initialise selection condition variables
	primaryTrackValid = false;
	secondaryTrackValid = false;
	tertiaryTrackValid = false;

	primaryTrackPionlike = false;
	secondaryTrackPionlike = false;

	tertiaryTrackPionlike = false;

	primaryTrackPionlikeLoose = false;
	secondaryTrackPionlikeLoose = false;
	tertiaryTrackPionlikeLoose = false;


	n_non_proton_tracks = 0;
    n_pion_candidate_tracks = 0;

	wc_pion_candidate_index = -1;
	n_wc_reco_larpid_unclassified = 0;

	// number protons
	numberProtons = 0;
	if (nproton > 2) nproton = 2; // fudge to make N proton bin
	if (nprotonalternate > 2) nprotonalternate = 2; // fudge to make N proton bin

	// subtract electron mass from electron energy
	elec_e = elec_e - (0.51099895/1000);

	// initialise BDT score variables
	BDTScoreElectronPhoton = -1;
    primaryTrackBDTScorePionProton = -1;
    secondaryTrackBDTScorePionProton = -1;
    tertiaryTrackBDTScorePionProton = -1;
    BDTScorePionProtonAlternate = -1;

	wc_bdt_score = -1;

    // initialise STV pelee_tree variables
    if (type == Utility::kMC || type == Utility::kDetVar || type == Utility::kIntrinsic || type == Utility::kDirt || type == Utility::kFakeData) is_mc_ = true;
    else is_mc_ = false;
    mc_is_signal_ = false;
    mc_is_inclusive_signal_ = false;
    category_ = 9999;

    // set horn current
    // FHC
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) hornCurrent_ = 0;
	// RHC
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) hornCurrent_ = 1;
	// Unknown
	else {
		std::cout << "[EventContainer::populateDerivedVariables] Error: unknown run period." << std::endl;
		exit(1);
	}

	//std::cout << "Here 1c" << std::endl;

    // clear weights
    beamlineVarWeightsPresent = false;
	
    Horn_2kA.clear();
    Horn1_x_3mm.clear();
    Horn1_y_3mm.clear();
    Beam_spot_1_1mm.clear();
    Beam_spot_1_5mm.clear();
    Horn2_x_3mm.clear();
    Horn2_y_3mm.clear();
    Horns_0mm_water.clear();
    Horns_2mm_water.clear();
    Beam_shift_x_1mm.clear();
    Beam_shift_y_1mm.clear();
    Target_z_7mm.clear();
	

    fluggWeightsPresent = false;
    FluggCV.clear();

    reweightRatioPresent = false;	
	reweightRatio = 1;

	n_primary_tracks = 0;
    n_primary_track_exiting = 0;
    n_primary_showers = 0;

	wc_n_pion_candidate_daughters = 0;
	wc_pion_candidate_daughters_total_energy = 0.0;
	wc_n_blips_25cm = 0;
	wc_n_blips_50cm = 0;
	wc_n_blips_100cm = 0;

	//std::cout << "trk_sce_start_x_v->size() = " << trk_sce_start_x_v->size() << std::endl;
	//std::cout << "wc_reco_Ntrack = " << wc_reco_Ntrack << std::endl;

	pfp_classification_v = std::vector<Utility::ParticleEnums>(trk_sce_start_x_v->size(), Utility::kCosmic);
	//std::cout << "Here 1d" << std::endl;
	wc_particle_classification_v = std::vector<Utility::ParticleEnums>(wc_reco_Ntrack, Utility::kCosmic);
	lantern_particle_classification_v = std::vector<Utility::ParticleEnums>(lantern_nTracks, Utility::kCosmic);
	//std::cout << "Here 1e" << std::endl;

	// WC variables
	wc_n_primary_particles = 0;	
	wc_primaryTrackIndex = -1;
	wc_primaryTrackLength = -1;
	wc_primaryTrackVertexSeparation = -1;

	// Lantern variables
	lantern_primaryTrackIndex = -1;
    lantern_primaryTrackLength = 0.0; 
    lantern_primaryTrackVertexSeparation = 0.0;
  	lantern_pion_candidate_index = -1;

	truePionHasMuon = false;
	truePionHasPi0 = false;
	truePionNProtons = 0;

	//std::cout << "Here 1f" << std::endl;

}

// get shower dE/dx on plane with the most hits
void EventContainer::GetdEdxMax(bool includeGap = false) {

	float temp_shr_trkfit_dedx_max = -1;

    // We want to also use the dedx when it is defined properly. Sometimes, the plane can have hits but an undefined dedx
    // use the dedx where we get the max number of hits and the dedx > 0
    int temp_shr_hits_u_tot, temp_shr_hits_v_tot, temp_shr_hits_y_tot;    
	float temp_shr_tkfit_dedx_U, temp_shr_tkfit_dedx_V, temp_shr_tkfit_dedx_Y;

	// Case where gap is included
	if (includeGap) {
		temp_shr_hits_u_tot = shr_tkfit_gap10_nhits_U;
    	temp_shr_hits_v_tot = shr_tkfit_gap10_nhits_V;
    	temp_shr_hits_y_tot = shr_tkfit_gap10_nhits_Y;
		temp_shr_tkfit_dedx_U = shr_tkfit_gap10_dedx_U;
		temp_shr_tkfit_dedx_V = shr_tkfit_gap10_dedx_V;
		temp_shr_tkfit_dedx_Y = shr_tkfit_gap10_dedx_Y;
	}
	else {
		temp_shr_hits_u_tot = shr_tkfit_2cm_nhits_U;
    	temp_shr_hits_v_tot = shr_tkfit_2cm_nhits_V;
    	temp_shr_hits_y_tot = shr_tkfit_2cm_nhits_Y;
		temp_shr_tkfit_dedx_U = shr_tkfit_2cm_dedx_U;
		temp_shr_tkfit_dedx_V = shr_tkfit_2cm_dedx_V;
		temp_shr_tkfit_dedx_Y = shr_tkfit_2cm_dedx_Y;   
	}
	// If the dedx is undefined, set the hits to zero
    if (temp_shr_tkfit_dedx_U <= 0 || !_utility.isNumber(temp_shr_tkfit_dedx_U)) temp_shr_hits_u_tot = 0;
    if (temp_shr_tkfit_dedx_V <= 0 || !_utility.isNumber(temp_shr_tkfit_dedx_V)) temp_shr_hits_v_tot = 0;
    if (temp_shr_tkfit_dedx_Y <= 0 || !_utility.isNumber(temp_shr_tkfit_dedx_Y)) temp_shr_hits_y_tot = 0;

    // Collection plane is the largest
    if (temp_shr_hits_y_tot > temp_shr_hits_u_tot && temp_shr_hits_y_tot > temp_shr_hits_v_tot ){
        temp_shr_trkfit_dedx_max = temp_shr_tkfit_dedx_Y;
    }
    // V Plane is the largest
    else if (temp_shr_hits_v_tot > temp_shr_hits_u_tot && temp_shr_hits_v_tot > temp_shr_hits_y_tot) {
        temp_shr_trkfit_dedx_max = temp_shr_tkfit_dedx_V;        
    }
    // U Plane is the largest
    else if (temp_shr_hits_u_tot > temp_shr_hits_v_tot && temp_shr_hits_u_tot > temp_shr_hits_y_tot){
        temp_shr_trkfit_dedx_max = temp_shr_tkfit_dedx_U;
    }
    // One plane was equal, so need to prioritise planes in preference of y, v, u
    else {

        // If y == any other plane, then y wins
        if (temp_shr_hits_y_tot == temp_shr_hits_u_tot || temp_shr_hits_y_tot == temp_shr_hits_v_tot ){
            temp_shr_trkfit_dedx_max = temp_shr_tkfit_dedx_Y;           
        }
        // U == V, ALL Y cases have been used up, so default to v
        else if (temp_shr_hits_u_tot == temp_shr_hits_v_tot ){
            temp_shr_trkfit_dedx_max = temp_shr_tkfit_dedx_V;            
        }
        else {
            temp_shr_trkfit_dedx_max = temp_shr_tkfit_dedx_U;
        }
    }

    if (temp_shr_trkfit_dedx_max == -1) {
        std::cout << shr_tkfit_dedx_U << " " << shr_tkfit_dedx_V << " " << shr_tkfit_dedx_Y<< std::endl;
        std::cout << "Error [Selection.h]: Edge case of dEdx comparisons." << std::endl;
    }

    // set variable
    if (includeGap) shr_trkfit_gap10_dedx_max = temp_shr_trkfit_dedx_max;
    else shr_trkfit_2cm_dedx_max = temp_shr_trkfit_dedx_max;

}

// get second shower opening angle, if present
void EventContainer::GetSecondShowerClusterAngleDiff() {
	
	// U plane, if defined
	if (secondshower_U_nhit > 0 && shrclusdir0 > 0) {

		// check each combination to determine correct direction - to account for cases where showers are on different sides of the wire direction
		float anglediff_case_1 = std::abs(secondshower_U_dir - shrclusdir0);
		float anglediff_case_2 = std::abs((360 - secondshower_U_dir) + shrclusdir0);
		float anglediff_case_3 = std::abs(secondshower_U_dir + (360 - shrclusdir0));

		// set as minimum of the three cases
		if (anglediff_case_1 <= anglediff_case_2 && anglediff_case_1 <= anglediff_case_3) secondshower_U_anglediff = anglediff_case_1;
		else if (anglediff_case_2 <= anglediff_case_1 && anglediff_case_2 <= anglediff_case_3) secondshower_U_anglediff = anglediff_case_2;
		else if (anglediff_case_3 <= anglediff_case_1 && anglediff_case_3 <= anglediff_case_2) secondshower_U_anglediff = anglediff_case_3; 
	}
	else secondshower_U_anglediff = 9999;

	// V plane, if defined
	if (secondshower_V_nhit > 0 && shrclusdir1 > 0) {

		// check each combination to determine correct direction - to account for cases where showers are on different sides of the wire direction
		float anglediff_case_1 = std::abs(secondshower_V_dir - shrclusdir1);
		float anglediff_case_2 = std::abs((360 - secondshower_V_dir) + shrclusdir1);
		float anglediff_case_3 = std::abs(secondshower_V_dir + (360 - shrclusdir1));

		// set as minimum of the three cases
		if (anglediff_case_1 <= anglediff_case_2 && anglediff_case_1 <= anglediff_case_3) secondshower_V_anglediff = anglediff_case_1;
		else if (anglediff_case_2 <= anglediff_case_1 && anglediff_case_2 <= anglediff_case_3) secondshower_V_anglediff = anglediff_case_2;
		else if (anglediff_case_3 <= anglediff_case_1 && anglediff_case_3 <= anglediff_case_2) secondshower_V_anglediff = anglediff_case_3; 
	}
	else secondshower_V_anglediff = 9999;

	// Y plane, if defined
	if (secondshower_Y_nhit > 0 && shrclusdir2 > 0) {

		// check each combination to determine correct direction - to account for cases where showers are on different sides of the wire direction
		float anglediff_case_1 = std::abs(secondshower_Y_dir - shrclusdir2);
		float anglediff_case_2 = std::abs((360 - secondshower_Y_dir) + shrclusdir2);
		float anglediff_case_3 = std::abs(secondshower_Y_dir + (360 - shrclusdir2));

		// set as minimum of the three cases
		//secondshower_Y_anglediff = anglediff_case_1;
		if (anglediff_case_1 <= anglediff_case_2 && anglediff_case_1 <= anglediff_case_3) secondshower_Y_anglediff = anglediff_case_1;
		else if (anglediff_case_2 <= anglediff_case_1 && anglediff_case_2 <= anglediff_case_3) secondshower_Y_anglediff = anglediff_case_2;
		else if (anglediff_case_3 <= anglediff_case_1 && anglediff_case_3 <= anglediff_case_2) secondshower_Y_anglediff = anglediff_case_3; 
	}
	else secondshower_Y_anglediff = 9999;
}

void EventContainer::isAlongWire(unsigned int trackID, bool &isAlongWire_u, bool &isAlongWire_v, bool &isAlongWire_y) {
	float tolerance = 25; 	// degrees

	// calculate theta_yz
	float delta_y = trk_sce_end_y_v->at(trackID-1) - trk_sce_start_y_v->at(trackID-1);
	float delta_z = trk_sce_end_z_v->at(trackID-1) - trk_sce_start_z_v->at(trackID-1);
	float theta_yz = std::atan2(delta_z, delta_y) * 180 / 3.1415;

	// check whether track along wire directions
	if ((theta_yz > 60 - tolerance && theta_yz < 60 + tolerance) || (theta_yz > -120 - tolerance && theta_yz < -120 + tolerance)) isAlongWire_u = true; // +60
	if ((theta_yz > 120 - tolerance && theta_yz < 120 + tolerance) || (theta_yz > -60 - tolerance && theta_yz < -60 + tolerance)) isAlongWire_v = true; // - 60	
	if (std::abs(theta_yz) < tolerance || std::abs(theta_yz) > 180 - tolerance) isAlongWire_y = true; // vertical
}

// get track trunk dE/dx - use collection plane when available, unless track along wire direction
float EventContainer::GetTrackTrunkdEdxBestPlane(unsigned int trackID) {

	// check whether track is along wire
	bool isAlongWire_u = false; bool isAlongWire_v = false; bool isAlongWire_y = false;
	isAlongWire(trackID, isAlongWire_u, isAlongWire_v, isAlongWire_y);

	// get number of hits on each plane
	int Nhits_Y = trk_nhits_y_v->at(trackID-1);
	int Nhits_U = trk_nhits_u_v->at(trackID-1);
	int Nhits_V = trk_nhits_v_v->at(trackID-1);
	 
	// prefer collection plane where available, otherwise take average of u,v planes if available
	// track trunk dE/dx
	if (trk_trunk_dEdx_y_v->at(trackID-1) != 9999 && _utility.isNumber(trk_trunk_dEdx_y_v->at(trackID-1)) && trk_trunk_dEdx_y_v->at(trackID-1) > 0.0001 && !isAlongWire_y && Nhits_Y > 0 && trk_trunk_dEdx_y_v->at(trackID-1) < 100) {
		return trk_trunk_dEdx_y_v->at(trackID-1);
	}
	else if (trk_trunk_dEdx_u_v->at(trackID-1) != 9999 && _utility.isNumber(trk_trunk_dEdx_u_v->at(trackID-1)) && trk_trunk_dEdx_u_v->at(trackID-1) > 0.0001 && !isAlongWire_u && trk_trunk_dEdx_v_v->at(trackID-1) != 9999 && _utility.isNumber(trk_trunk_dEdx_v_v->at(trackID-1)) && trk_trunk_dEdx_v_v->at(trackID-1) > 0.0001 && !isAlongWire_v && Nhits_U > 0 && Nhits_V > 0 && trk_trunk_dEdx_u_v->at(trackID-1) < 100 && trk_trunk_dEdx_v_v->at(trackID-1) < 100) {
		return (Nhits_U * trk_trunk_dEdx_u_v->at(trackID-1) + Nhits_V * trk_trunk_dEdx_v_v->at(trackID-1)) / (Nhits_U + Nhits_V);
	}
	else if (trk_trunk_dEdx_v_v->at(trackID-1) != 9999 && _utility.isNumber(trk_trunk_dEdx_v_v->at(trackID-1)) && trk_trunk_dEdx_v_v->at(trackID-1) > 0.0001 && !isAlongWire_v && Nhits_V > 0 && trk_trunk_dEdx_v_v->at(trackID-1) < 100) {
		return trk_trunk_dEdx_v_v->at(trackID-1);
	} 
	else if (trk_trunk_dEdx_u_v->at(trackID-1) != 9999 && _utility.isNumber(trk_trunk_dEdx_u_v->at(trackID-1)) && trk_trunk_dEdx_u_v->at(trackID-1) > 0.0001 && !isAlongWire_u && Nhits_U > 0 && trk_trunk_dEdx_u_v->at(trackID-1) < 100) {
		return trk_trunk_dEdx_u_v->at(trackID-1);
	}
	else return 9999;
}

// bragg peak pion
float EventContainer::GetTrackBraggPionBestPlane(unsigned int trackID) {

	// check whether track is along wire
	bool isAlongWire_u = false; bool isAlongWire_v = false; bool isAlongWire_y = false;
	isAlongWire(trackID, isAlongWire_u, isAlongWire_v, isAlongWire_y);

	// get number of hits on each plane
	int Nhits_Y = trk_nhits_y_v->at(trackID-1);
	int Nhits_U = trk_nhits_u_v->at(trackID-1);
	int Nhits_V = trk_nhits_v_v->at(trackID-1);

	// track bragg peak pion
	if (trk_bragg_pion_v->at(trackID-1) != 9999 && !isAlongWire_y && Nhits_Y > 0 && trk_bragg_pion_v->at(trackID-1) >= 0.0001) {
		return trk_bragg_pion_v->at(trackID-1);
	}
	else if (trk_bragg_pion_u_v->at(trackID-1) != 9999 && !isAlongWire_u && trk_bragg_pion_v_v->at(trackID-1) != 9999 && !isAlongWire_v && Nhits_U > 0 && Nhits_V > 0 && trk_bragg_pion_u_v->at(trackID-1) >= 0.0001 && trk_bragg_pion_v_v->at(trackID-1) >= 0.0001) {
		return (Nhits_U * trk_bragg_pion_u_v->at(trackID-1) + Nhits_V * trk_bragg_pion_v_v->at(trackID-1)) / (Nhits_U + Nhits_V);
	}
	else if (trk_bragg_pion_v_v->at(trackID-1) != 9999 && !isAlongWire_v && Nhits_V > 0 && trk_bragg_pion_v_v->at(trackID-1) >= 0.0001) {
		return trk_bragg_pion_v_v->at(trackID-1);
	} 
	else if (trk_bragg_pion_u_v->at(trackID-1) != 9999 && !isAlongWire_u && Nhits_U > 0 && trk_bragg_pion_u_v->at(trackID-1) >= 0.0001) {
		return trk_bragg_pion_u_v->at(trackID-1);
	}
	else return 9999;
}

// bragg peak mip
float EventContainer::GetTrackBraggMIPBestPlane(unsigned int trackID) {

	// check whether track is along wire
	bool isAlongWire_u = false; bool isAlongWire_v = false; bool isAlongWire_y = false;
	isAlongWire(trackID, isAlongWire_u, isAlongWire_v, isAlongWire_y);

	// get number of hits on each plane
	int Nhits_Y = trk_nhits_y_v->at(trackID-1);
	int Nhits_U = trk_nhits_u_v->at(trackID-1);
	int Nhits_V = trk_nhits_v_v->at(trackID-1);

	// track bragg peak pion
	if (trk_bragg_mip_v->at(trackID-1) != 9999 && !isAlongWire_y && Nhits_Y > 0 && trk_bragg_mip_v->at(trackID-1) >= 0.0001) {
		return trk_bragg_mip_v->at(trackID-1);
	}
	else if (trk_bragg_mip_u_v->at(trackID-1) != 9999 && !isAlongWire_u && trk_bragg_mip_v_v->at(trackID-1) != 9999 && !isAlongWire_v && Nhits_U > 0 && Nhits_V > 0 && trk_bragg_mip_u_v->at(trackID-1) >= 0.0001 && trk_bragg_mip_v_v->at(trackID-1) >= 0.0001) {
		return (Nhits_U * trk_bragg_mip_u_v->at(trackID-1) + Nhits_V * trk_bragg_mip_v_v->at(trackID-1)) / (Nhits_U + Nhits_V);
	}
	else if (trk_bragg_mip_v_v->at(trackID-1) != 9999 && !isAlongWire_v && Nhits_V > 0 && trk_bragg_mip_v_v->at(trackID-1) >= 0.0001) {
		return trk_bragg_mip_v_v->at(trackID-1);
	} 
	else if (trk_bragg_mip_u_v->at(trackID-1) != 9999 && !isAlongWire_u && Nhits_U > 0 && trk_bragg_mip_u_v->at(trackID-1) >= 0.0001) {
		return trk_bragg_mip_u_v->at(trackID-1);
	}
	else return 9999;
}

// pion energy range-based
// from A. Smith analysis -- docid=33809
float EventContainer::CalculatePionMomentumRange(float R) {

	float p_range = 0.25798 + (0.0024088 * R) - (0.18828 * std::pow(R, - 0.11687));	// GeV
	float e_range = std::sqrt(std::pow(p_range, 2) + std::pow(0.13957, 2)) - 0.13957; // GeV

	if (p_range > 0) return p_range;
	else return 0;
}

float EventContainer::CalculatePionMomentumHypfit(int candidateID) {

	//std::vector<double> test_rr = {2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30.};
	//std::vector<double> test_dedx = {6., 5.8, 5.6, 5.4, 5.2, 5.0, 4.8, 4.6, 4.4, 4.2, 4.0, 3.8, 3.6, 3.4, 3.2};

	// get candidate track ID
	int recoID = wc_reco_id[wc_pion_candidate_index];

	//std::cout << "Candidate track ID: " << recoID << std::endl;

	// get spacepoints for candidate track
	std::vector<float> spacepoints_x; spacepoints_x.reserve(1000);
	std::vector<float> spacepoints_y; spacepoints_y.reserve(1000);
	std::vector<float> spacepoints_z; spacepoints_z.reserve(1000);
	std::vector<float> spacepoints_q; spacepoints_q.reserve(1000);

	//std::cout << "Total number of spacepoints: " << wc_Trecchargeblob_spacepoints_real_cluster_id->size() << std::endl;

	for (int i = 0; i < wc_Trecchargeblob_spacepoints_real_cluster_id->size(); i++) {
		//std::cout << "Spacepoint " << i << ": cluster ID = " << wc_Trecchargeblob_spacepoints_real_cluster_id->at(i) << std::endl;
		if (wc_Trecchargeblob_spacepoints_real_cluster_id->at(i) == recoID) {
			spacepoints_x.push_back(wc_Trecchargeblob_spacepoints_x->at(i));
			spacepoints_y.push_back(wc_Trecchargeblob_spacepoints_y->at(i));
			spacepoints_z.push_back(wc_Trecchargeblob_spacepoints_z->at(i));
			spacepoints_q.push_back(wc_Trecchargeblob_spacepoints_q->at(i));
		}
	}
	
	//std::cout << "Number of spacepoints for candidate track: " << spacepoints_x.size() << std::endl;



	//double pion_likelihood_KE = hypfit->Likelihood(test_dedx, test_rr, 211);
	//double pion_gaus_KE = hypfit->Gaussian(test_dedx, test_rr, 211);
	//std::cout << "pion_likelihood_KE: " << pion_likelihood_KE << std::endl;
	//std::cout << "pion_gaus_KE: " << pion_gaus_KE << std::endl;

	return 0;
}

// shower track end proximity
float EventContainer::GetShowerTrackEndProximity(unsigned int shrID) {

	float shr_trackEndProximity = 9999;

	if (shrID < shr_start_x_v->size() && shrID != 0) {

		// check against each track end to find nearest
		// primary
		if (trk_id < trk_end_x_v->size() && trk_id != 0) { 
			float distance = std::sqrt(	std::pow(shr_start_x_v->at(shrID-1) - trk_end_x, 2) + std::pow(shr_start_y_v->at(shrID-1) - trk_end_y, 2) + std::pow(shr_start_z_v->at(shrID-1) - trk_end_z, 2));
			if (distance < shr_trackEndProximity) shr_trackEndProximity = distance; 
		}
		// secondary
		if (trk2_id < trk_end_x_v->size() && trk2_id != 0) { 
			float distance = std::sqrt(	std::pow(shr_start_x_v->at(shrID-1) - trk2_end_x, 2) + std::pow(shr_start_y_v->at(shrID-1) - trk2_end_y, 2) + std::pow(shr_start_z_v->at(shrID-1) - trk2_end_z, 2));
			if (distance < shr_trackEndProximity) shr_trackEndProximity = distance; 
		}
		// tertiary
		if (trk3_id < trk_end_x_v->size() && trk3_id != 0) { 
			float distance = std::sqrt(	std::pow(shr_start_x_v->at(shrID-1) - trk3_end_x, 2) + std::pow(shr_start_y_v->at(shrID-1) - trk3_end_y, 2) + std::pow(shr_start_z_v->at(shrID-1) - trk3_end_z, 2));
			if (distance < shr_trackEndProximity) shr_trackEndProximity = distance; 
		}
	}

	return shr_trackEndProximity;
}

// -----------------------------------------------------------------------------
// NuMI angle calculate, from Krishan
double EventContainer::GetNuMIAngle(double px, double py, double pz, std::string direction){

    // Variables
    TRotation RotDet2Beam;             // Rotations
    TVector3  detxyz, BeamCoords;      // Translations
    std::vector<double> rotmatrix;     // Inputs

    // input detector coordinates to translate
    detxyz = {px, py, pz};     

    // From beam to detector rotation matrix
    rotmatrix = {
        0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021,
        4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359,
        -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291 };

    // Return the TRotation
    TVector3 newX, newY, newZ;
    newX = TVector3(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
    newY = TVector3(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
    newZ = TVector3(rotmatrix[6], rotmatrix[7], rotmatrix[8]);

    RotDet2Beam.RotateAxes(newX, newY, newZ); // Return the TRotation now det to beam
    // RotDet2Beam.Invert(); // Invert back to the beam to det

    // Rotate to beam coords
    BeamCoords = RotDet2Beam * detxyz;

    TVector3 beamdir = {0 , 0 , 1};;
    
    // Get the angle wrt to the beam
    if (direction == "beam") beamdir = {0 , 0 , 1};
    
    // Get the angle wrt to the target to detector direction
    else if (direction == "target") {
        beamdir = {5502, 7259, 67270};
        beamdir = beamdir.Unit(); // Get the direction
    }
    else {
        std::cout << "Warning unknown angle type specified, you should check this" << std::endl;
    }
    
    double angle = BeamCoords.Angle(beamdir) * 180 / 3.1415926;

    // Create vectors to get the angle in the yz and xz planes
    TVector3 BeamCoords_yz = { 0, BeamCoords.Y(), BeamCoords.Z() }; // Angle upwards
    TVector3 BeamCoords_xz = { BeamCoords.X(), 0, BeamCoords.Z() }; // Angle across

    return angle;
}

// -----------------------------------------------------------------------------
// NuMI angular variables, from Krishan
void EventContainer::setRecoNuMIAngularVariables(bool isSideband = false){
	
	// Effective angle - beta (see K. Mistry Thesis)
    // Electron
	// Calculate the angle between the shower direction and the vector from the target to the nu vtx
    TVector3 shr_dir(shr_dir_x, shr_dir_y, shr_dir_z); // Shower direction
    shr_dir.Unit();
    
    TVector3 v_targ_uboone(-31387.58422, -3316.402543, -60100.2414);
    TVector3 v_nu_vtx(reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z);
    TVector3 v_targ_to_vtx = (-1*v_targ_uboone + v_nu_vtx).Unit(); // -1 because the vector points from uboone to tgt, we need the other way around

    // Set the values
    reco_electron_effective_angle = shr_dir.Angle(v_targ_to_vtx) * 180 / 3.14159;
    reco_cos_electron_effective_angle = std::cos(shr_dir.Angle(v_targ_to_vtx));

    // Pion
    if (isSideband) {
    	// for sidebands use primary track
    	TVector3 piontrk_dir;	
	    piontrk_dir.SetXYZ(trk_dir_x, trk_dir_y, trk_dir_z);
	    piontrk_dir.Unit();

	    // Set the values
	    reco_pion_effective_angle = piontrk_dir.Angle(v_targ_to_vtx) * 180 / 3.14159;
	    reco_cos_pion_effective_angle = std::cos(piontrk_dir.Angle(v_targ_to_vtx));	

	    // opening angle
	    reco_electron_pion_opening_angle = piontrk_dir.Angle(shr_dir) * 180 / 3.14159;
	    reco_cos_electron_pion_opening_angle = std::cos(piontrk_dir.Angle(shr_dir));

        // for the sidebands also need to set momentum of pion to primary track, otherwise undefined
        reco_momentum_pion = trk_momentum_pion;	
    }
    else {
        if (primaryTrackPionlike || secondaryTrackPionlike || tertiaryTrackPionlike) {
	    	TVector3 piontrk_dir;
	    	if (primaryTrackPionlike) {
	    		piontrk_dir.SetXYZ(trk_dir_x, trk_dir_y, trk_dir_z);
	    		
	    	}
	    	else if (secondaryTrackPionlike) {
	    		piontrk_dir.SetXYZ(trk2_dir_x, trk2_dir_y, trk2_dir_z);
	    	}
	    	else {
	    		piontrk_dir.SetXYZ(trk3_dir_x, trk3_dir_y, trk3_dir_z);
	    	}
	    	piontrk_dir.Unit();

	    	// Set the values
		    reco_pion_effective_angle = piontrk_dir.Angle(v_targ_to_vtx) * 180 / 3.14159;
		    reco_cos_pion_effective_angle = std::cos(piontrk_dir.Angle(v_targ_to_vtx));

		    // opening angle
	    	reco_electron_pion_opening_angle = piontrk_dir.Angle(shr_dir) * 180 / 3.14159;
	    	reco_cos_electron_pion_opening_angle = std::cos(piontrk_dir.Angle(shr_dir));
	    }
	}
}

void EventContainer::setTrueNuMIAngularVariables(){
	
	// neutrino
    TVector3 nu_dir(true_nu_px, true_nu_py, true_nu_pz); 
    nu_dir.Unit();

    // electron
    TVector3 elec_dir(elec_px, elec_py, elec_pz); 
    elec_dir.Unit();
    true_electron_effective_angle = elec_dir.Angle(nu_dir) * 180 / 3.14159; 
    true_cos_electron_effective_angle = std::cos(elec_dir.Angle(nu_dir));

    // pion
   	TVector3 pion_dir(pion_px, pion_py, pion_pz); 
    pion_dir.Unit();

	true_pion_effective_angle = pion_dir.Angle(nu_dir) * 180 / 3.14159;
	true_cos_pion_effective_angle = std::cos(pion_dir.Angle(nu_dir));

    // opening angle
    true_electron_pion_opening_angle = pion_dir.Angle(elec_dir) * 180 / 3.14159; 
    true_cos_electron_pion_opening_angle = std::cos(pion_dir.Angle(elec_dir));
}

// Print Daughters 
void EventContainer::PrintDaughters(int pionID, std::vector<int>& daughters) {

	// find daughters of primary pion
	for (unsigned int daughter_idx = 0; daughter_idx < wc_truth_Ntrack; daughter_idx++) {
		
		if (wc_truth_mother[daughter_idx] != pionID) continue; 

		// skip neutrons
		if (std::abs(wc_truth_pdg[daughter_idx]) == 2112) continue;
		
		// skip Ar fragments
		if (std::abs(wc_truth_pdg[daughter_idx]) > 10000) continue;

		//std::cout << "    Daughter PDG: " << wc_truth_pdg[daughter_idx] 
		//			<< ", Start Process: " << wc_truth_process->at(daughter_idx) 
		//			<< ", End Process: " << wc_truth_endprocess->at(daughter_idx) 
		//			<< std::endl;

		// add threshold
		// find correct particle to calculate truth momentum		
		double momentum = std::sqrt( pow(wc_truth_startMomentum[daughter_idx][0],2) + 
									pow(wc_truth_startMomentum[daughter_idx][1],2) + 
									pow(wc_truth_startMomentum[daughter_idx][2],2) );
		
		// proton threshold, 300 MeV/c
		if (std::abs(wc_truth_pdg[daughter_idx]) == 2212 && momentum < 0.3) continue;						

		daughters.push_back(wc_truth_pdg[daughter_idx]);
		
		// is the daughter a pion?
		if (std::abs(wc_truth_pdg[daughter_idx]) == 211) {
			// recursively print granddaughters
			//std::cout << "    Granddaughters: " << std::endl;
			int daughterID = wc_truth_id[daughter_idx];
			PrintDaughters(daughterID, daughters);
		}
	}
}

double EventContainer::GetFSIWeight() {

	double fsi_weight = 1.0;

	// only consider NC 1pi and CC 1pi events, for now
	if (classification == Utility::kNC1pi || classification == Utility::kCCNumu1pi) {

		// find pions and calculate KE / angle + flag whether pi+ or pi-
		for (int i = 0; i < mc_pdg_v->size(); i++) {

			// skip non-pions
			if (!(mc_pdg_v->at(i) == 211 || mc_pdg_v->at(i) == -211)) continue;

			// momentum
			double pz = mc_pz_v->at(i);
			double mc_momentum = std::sqrt( pow(mc_px_v->at(i),2) + 
											pow(mc_py_v->at(i),2) + 
											pow(mc_pz_v->at(i),2) );
			
			// skip below threshold pions
			if (mc_momentum <= 0.1) continue;

			// calculate KE and angle
			double mass = 0.13957; // GeV
			double mc_energy = std::sqrt( pow(mc_momentum, 2) + mass*mass );
			double mc_ke = (mc_energy - mass) * 1000; // subtract pion mass, convert to MeV
			double mc_cos_theta = pz / mc_momentum;

			// get correct histogram 
			TH2D* h;
			if (classification == Utility::kNC1pi) {
				if (mc_pdg_v->at(i) == 211) h = h_fsi_nc_piplus;
				else if (mc_pdg_v->at(i) == -211) h = h_fsi_nc_piminus;
			}
			else if (classification == Utility::kCCNumu1pi) {
				if (mc_pdg_v->at(i) == 211) h = h_fsi_cc_piplus;
				else if (mc_pdg_v->at(i) == -211) h = h_fsi_cc_piminus;
			}
			else {
				std::cout << "Warning: FSI weight calculation only implemented for NC 1pi and CC 1pi events, check this!" << std::endl;
				return 1.0;
			}
			
			// get weight
			if (mc_ke <= 1000) {
				int binx = h->GetXaxis()->FindBin(mc_ke);
				binx = std::max(1, std::min(binx, h->GetXaxis()->GetNbins()));
				int biny = h->GetYaxis()->FindBin(mc_cos_theta);
				biny = std::max(1, std::min(biny, h->GetYaxis()->GetNbins()));

				fsi_weight = h->GetBinContent(binx, biny);
			}
			else fsi_weight = 1.0; // above 1 GeV, no FSI reweighting available, so set to 1

			if (fsi_weight == 0) {
				std::cout << "Warning: FSI weight is zero, "; 
				std::cout << "KE: " << mc_ke << " MeV, cos(theta): " << mc_cos_theta << std::endl;
				fsi_weight = 1.0; // reset to 1 to avoid zeroing out events
			}

			break; // only consider first pion in event
		}

		//std::cout << "FSI weight: " << fsi_weight << std::endl;

		fsi_weight = checkWeight(fsi_weight);
	}

	return fsi_weight;
}