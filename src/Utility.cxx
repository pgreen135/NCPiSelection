#include "../include/Utility.h"

#include <iostream>
#include <numeric>
#include <string>
#include <sstream>

#include <TFile.h>

// Constructor
Utility::Utility(bool loadBeamLineWeights, bool loadPPFXCVWeights, bool loadFluggWeights, bool loadReweightRatios) {
	
	std::cout << "Initialising Utilities Class" << std::endl;

	// load beamline variation histograms
	if (loadBeamLineWeights) {
		std::cout << "Loading NuMI beamline variation weights" << std::endl;
		loadBeamLineVariationHistograms();
	}

	// load PPFX CV histograms
	if (loadPPFXCVWeights) {
		std::cout << "Loading PPFX CV weights for fake data" << std::endl;
		loadPPFXCVHistograms();
	}

	// load Flugg CV histograms
	if (loadFluggWeights) {
		std::cout << "Loading Flugg weights for fake data" << std::endl;
		loadFluggHistograms();
	}

	// load re-weighting angle ratio histograms
	if (loadReweightRatios) {
		std::cout << "Loading re-weighting ratios" << std::endl;
		loadReweightRatioHistograms();
	}
}

// ------------------------------------------------------------------------------

// Function to check if point is within FV
bool Utility::inFV(const double x, const double y, const double z) const {
	// check the point is in the boundary 										
	// NTuples definition
	if ( x   >= 10 && x <= 246
		 && y   >= -101 && y <= 101
	  	 && z   >= 10 && z <= 986  ) {
	  	return true;
	}
	else return false;
}

// ------------------------------------------------------------------------------

// Function to check whether a track is exiting
bool Utility::isExiting(const double x, const double y, const double z) const {
	// containment definition (10 cm from borders)
	if ( x   >= 10 && x <= 246
		 && y   >= -106 && y <= 106
	  	 && z   >= 10 && z <= 1026  ) {
	  	return true;
	}
	else return false;
}

// ------------------------------------------------------------------------------

// Function to check whether number is valid
bool Utility::isNumber(float input) const {

	// infinite 
	if (std::isinf(input)) {
		return false;
	}

	// nan 
	if (std::isnan(input)) {
		return false;
	}

	return true;
}

// ------------------------------------------------------------------------------

void Utility::scaleByBinWidth(TH1F* hist) {
    int num_bins = hist->GetNbinsX();
    for ( int b = 1; b <= num_bins; b++ ) {
      double width = hist->GetBinWidth(b);
      double value = hist->GetBinContent(b);
      hist->SetBinContent(b, value/width);
    }
}

// ------------------------------------------------------------------------------

// Functions to get beam line variation weights
void Utility::loadBeamLineVariationHistograms() {

	// load beamline variations TH2D
	TFile *beamlineVariationsFile = NULL;
	beamlineVariationsFile = new TFile("/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/NuMIBeamLineWeights/NuMI_Geometry_Weights_Histograms.root");
	
	// loop over beamline variations
	std::vector<std::string> horncur {"FHC", "RHC"};

	for (int current = 0; current < horncur.size(); current++) {
		for (int i = 1; i <= 20; i++) {
			// construct names
			std::stringstream name_nue_ss; name_nue_ss << "EnergyTheta2D/ratio_run" << i <<  "_" << horncur[current] << "_nue_CV_AV_TPC_2D";
			std::stringstream name_nuebar_ss; name_nuebar_ss << "EnergyTheta2D/ratio_run" << i <<  "_" << horncur[current] << "_nuebar_CV_AV_TPC_2D";
			std::stringstream name_numu_ss; name_numu_ss << "EnergyTheta2D/ratio_run" << i <<  "_" << horncur[current] << "_numu_CV_AV_TPC_2D";
			std::stringstream name_numubar_ss; name_numubar_ss << "EnergyTheta2D/ratio_run" << i <<  "_" << horncur[current] << "_numubar_CV_AV_TPC_2D";
			std::string name_nue = name_nue_ss.str();
			std::string name_nuebar = name_nuebar_ss.str();
			std::string name_numu = name_numu_ss.str();
			std::string name_numubar = name_numubar_ss.str();		

			// load histograms
			if (horncur[current] == "FHC") {
				h_nue_fhc.push_back(*(TH2F*)beamlineVariationsFile->Get(name_nue.c_str()));
				h_nuebar_fhc.push_back(*(TH2F*)beamlineVariationsFile->Get(name_nuebar.c_str()));
				h_numu_fhc.push_back(*(TH2F*)beamlineVariationsFile->Get(name_numu.c_str()));
				h_numubar_fhc.push_back(*(TH2F*)beamlineVariationsFile->Get(name_numubar.c_str()));
			}
			else if (horncur[current] == "RHC") {
				h_nue_rhc.push_back(*(TH2F*)beamlineVariationsFile->Get(name_nue.c_str()));
				h_nuebar_rhc.push_back(*(TH2F*)beamlineVariationsFile->Get(name_nuebar.c_str()));
				h_numu_rhc.push_back(*(TH2F*)beamlineVariationsFile->Get(name_numu.c_str()));
				h_numubar_rhc.push_back(*(TH2F*)beamlineVariationsFile->Get(name_numubar.c_str()));	
			}
			else {
				std::cout << "Error [Utility::loadBeamLineVariationHistograms]: invalid horn current." << std::endl;
				exit(1);
			}
		}
	}

	beamlineVariationsFile->Close();
}

std::vector<double> Utility::getWeightsNue(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	std::vector<double> weights; weights.reserve(20);

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		for (int i = 0; i < h_nue_fhc.size(); i++) {
			int binx = h_nue_fhc[i].GetXaxis()->FindBin(energy);
			int biny = h_nue_fhc[i].GetYaxis()->FindBin(angle);
			float weight = h_nue_fhc[i].GetBinContent(binx,biny);
			weights.push_back(weight);
		}
	}
	// rhc
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) {
		for (int i = 0; i < h_nue_rhc.size(); i++) {
			int binx = h_nue_rhc[i].GetXaxis()->FindBin(energy);
			int biny = h_nue_rhc[i].GetYaxis()->FindBin(angle);
			float weight = h_nue_rhc[i].GetBinContent(binx,biny);
			weights.push_back(weight);
		}
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period." << std::endl;
	}

	return weights;
}

std::vector<double> Utility::getWeightsNuebar(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	std::vector<double> weights; weights.reserve(20);

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		for (int i = 0; i < h_nuebar_fhc.size(); i++) {
			int binx = h_nuebar_fhc[i].GetXaxis()->FindBin(energy);
			int biny = h_nuebar_fhc[i].GetYaxis()->FindBin(angle);
			float weight = h_nuebar_fhc[i].GetBinContent(binx,biny);
			weights.push_back(weight);
		}
	}
	// rhc
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) {
		for (int i = 0; i < h_nuebar_rhc.size(); i++) {
			int binx = h_nuebar_rhc[i].GetXaxis()->FindBin(energy);
			int biny = h_nuebar_rhc[i].GetYaxis()->FindBin(angle);
			float weight = h_nuebar_rhc[i].GetBinContent(binx,biny);
			weights.push_back(weight);
		}
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period." << std::endl;
	}

	return weights;
}

std::vector<double> Utility::getWeightsNumu(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	std::vector<double> weights; weights.reserve(20);

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		for (int i = 0; i < h_numu_fhc.size(); i++) {
			int binx = h_numu_fhc[i].GetXaxis()->FindBin(energy);
			int biny = h_numu_fhc[i].GetYaxis()->FindBin(angle);
			float weight = h_numu_fhc[i].GetBinContent(binx,biny);
			weights.push_back(weight);
		}
	}
	// rhc
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) {
		for (int i = 0; i < h_numu_rhc.size(); i++) {
			int binx = h_numu_rhc[i].GetXaxis()->FindBin(energy);
			int biny = h_numu_rhc[i].GetYaxis()->FindBin(angle);
			float weight = h_numu_rhc[i].GetBinContent(binx,biny);
			weights.push_back(weight);
		}
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period." << std::endl;
	}

	return weights;
}

std::vector<double> Utility::getWeightsNumubar(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	std::vector<double> weights; weights.reserve(20);

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		for (int i = 0; i < h_numubar_fhc.size(); i++) {
			int binx = h_numubar_fhc[i].GetXaxis()->FindBin(energy);
			int biny = h_numubar_fhc[i].GetYaxis()->FindBin(angle);
			float weight = h_numubar_fhc[i].GetBinContent(binx,biny);
			weights.push_back(weight);
		}
	}
	// rhc
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) {
		for (int i = 0; i < h_numubar_rhc.size(); i++) {
			int binx = h_numubar_rhc[i].GetXaxis()->FindBin(energy);
			int biny = h_numubar_rhc[i].GetYaxis()->FindBin(angle);
			float weight = h_numubar_rhc[i].GetBinContent(binx,biny);
			weights.push_back(weight);
		}
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period." << std::endl;
	}

	return weights;
}

// ------------------------------------------------------------------------------

// Functions to get PPFX CV weights for NuWro samples
void Utility::loadPPFXCVHistograms() {

	// load ratios
	// Old Flux
	/*
	TFile *PPFXCVFileFHC = new TFile("/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/NuMIPPFXCVWeights/ppfx_maps_fhc_new.root");	// FHC
	TFile *PPFXCVFileRHC = new TFile("/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/NuMIPPFXCVWeights/ppfx_maps_rhc_new.root");	// RHC

	ppfx_nue_fhc = *(TH2F*)PPFXCVFileFHC->Get("nue_ratio");
	ppfx_nuebar_fhc = *(TH2F*)PPFXCVFileFHC->Get("nuebar_ratio");
	ppfx_numu_fhc = *(TH2F*)PPFXCVFileFHC->Get("numu_ratio");
	ppfx_numubar_fhc = *(TH2F*)PPFXCVFileFHC->Get("numubar_ratio");

	ppfx_nue_rhc = *(TH2F*)PPFXCVFileRHC->Get("nue_ratio");
	ppfx_nuebar_rhc = *(TH2F*)PPFXCVFileRHC->Get("nuebar_ratio");
	ppfx_numu_rhc = *(TH2F*)PPFXCVFileRHC->Get("numu_ratio");
	ppfx_numubar_rhc = *(TH2F*)PPFXCVFileRHC->Get("numubar_ratio");

	PPFXCVFileFHC->Close();
	PPFXCVFileRHC->Close();
	*/
	// New Flux (direct re-weighting)
	TFile *PPFXCVFile = new TFile("/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/NuMIPPFXCVWeights/ratios_g4_10_4_cvonly_fornuwrofakedata.root"); // old flux without ppfx --> new flux with ppfx, E/Theta binning

	ppfx_nue_fhc = *(TH2F*)PPFXCVFile->Get("ratio_cv_fhc_enu_angle_nue");
	ppfx_nuebar_fhc = *(TH2F*)PPFXCVFile->Get("ratio_cv_fhc_enu_angle_nuebar");
	ppfx_numu_fhc = *(TH2F*)PPFXCVFile->Get("ratio_cv_fhc_enu_angle_numu");
	ppfx_numubar_fhc = *(TH2F*)PPFXCVFile->Get("ratio_cv_fhc_enu_angle_numubar");

	ppfx_nue_rhc = *(TH2F*)PPFXCVFile->Get("ratio_cv_rhc_enu_angle_nue");
	ppfx_nuebar_rhc = *(TH2F*)PPFXCVFile->Get("ratio_cv_rhc_enu_angle_nuebar");
	ppfx_numu_rhc = *(TH2F*)PPFXCVFile->Get("ratio_cv_rhc_enu_angle_numu");
	ppfx_numubar_rhc = *(TH2F*)PPFXCVFile->Get("ratio_cv_rhc_enu_angle_numubar");

	PPFXCVFile->Close();

}

double Utility::getCVWeightNue(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	double weight;

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		int binx = ppfx_nue_fhc.GetXaxis()->FindBin(energy);
		int biny = ppfx_nue_fhc.GetYaxis()->FindBin(angle);
		weight = ppfx_nue_fhc.GetBinContent(binx,biny);
	}
	// rhc
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) {
		int binx = ppfx_nue_rhc.GetXaxis()->FindBin(energy);
		int biny = ppfx_nue_rhc.GetYaxis()->FindBin(angle);
		weight = ppfx_nue_rhc.GetBinContent(binx,biny);
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period." << std::endl;
	}

	return weight;
}

double Utility::getCVWeightNuebar(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	double weight;

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		int binx = ppfx_nuebar_fhc.GetXaxis()->FindBin(energy);
		int biny = ppfx_nuebar_fhc.GetYaxis()->FindBin(angle);
		weight = ppfx_nuebar_fhc.GetBinContent(binx,biny);
	}
	// rhc
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) {
		int binx = ppfx_nuebar_rhc.GetXaxis()->FindBin(energy);
		int biny = ppfx_nuebar_rhc.GetYaxis()->FindBin(angle);
		weight = ppfx_nuebar_rhc.GetBinContent(binx,biny);
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period." << std::endl;
	}

	return weight;
}

double Utility::getCVWeightNumu(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	double weight;

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		int binx = ppfx_numu_fhc.GetXaxis()->FindBin(energy);
		int biny = ppfx_numu_fhc.GetYaxis()->FindBin(angle);
		weight = ppfx_numu_fhc.GetBinContent(binx,biny);
	}
	// rhc
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) {
		int binx = ppfx_numu_rhc.GetXaxis()->FindBin(energy);
		int biny = ppfx_numu_rhc.GetYaxis()->FindBin(angle);
		weight = ppfx_numu_rhc.GetBinContent(binx,biny);
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period." << std::endl;
	}

	return weight;
}

double Utility::getCVWeightNumubar(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	double weight;

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		int binx = ppfx_numubar_fhc.GetXaxis()->FindBin(energy);
		int biny = ppfx_numubar_fhc.GetYaxis()->FindBin(angle);
		weight = ppfx_numubar_fhc.GetBinContent(binx,biny);
	}
	// rhc
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) {
		int binx = ppfx_numubar_rhc.GetXaxis()->FindBin(energy);
		int biny = ppfx_numubar_rhc.GetYaxis()->FindBin(angle);
		weight = ppfx_numubar_rhc.GetBinContent(binx,biny);
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period." << std::endl;
	}

	return weight;
}

// ------------------------------------------------------------------------------

// Functions to get Flugg CV weights for Flugg fake data
void Utility::loadFluggHistograms() {

	// load beamline variations TH2D
	TFile *FluggFile = NULL;
	FluggFile = new TFile("/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/FluggWeights/flugg_ratios_ppfx.root");
	
	flugg_nue_fhc = *(TH1F*)FluggFile->Get("ratio_nue_FHC");
	flugg_nuebar_fhc = *(TH1F*)FluggFile->Get("ratio_nuebar_FHC");
	flugg_numu_fhc = *(TH1F*)FluggFile->Get("ratio_numu_FHC");
	flugg_numubar_fhc = *(TH1F*)FluggFile->Get("ratio_numubar_FHC");

	flugg_nue_rhc = *(TH1F*)FluggFile->Get("ratio_nue_RHC");
	flugg_nuebar_rhc = *(TH1F*)FluggFile->Get("ratio_nuebar_RHC");
	flugg_numu_rhc = *(TH1F*)FluggFile->Get("ratio_numu_RHC");
	flugg_numubar_rhc = *(TH1F*)FluggFile->Get("ratio_numubar_RHC");

	FluggFile->Close();
}

double Utility::getFluggWeightNue(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	double weight;

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		int binx = flugg_nue_fhc.GetXaxis()->FindBin(energy);
		//int biny = flugg_nue_fhc.GetYaxis()->FindBin(angle);
		//weight = flugg_nue_fhc.GetBinContent(binx,biny);
		weight = flugg_nue_fhc.GetBinContent(binx);
	}
	// rhc
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) {
		int binx = flugg_nue_rhc.GetXaxis()->FindBin(energy);
		//int biny = flugg_nue_rhc.GetYaxis()->FindBin(angle);
		//weight = flugg_nue_rhc.GetBinContent(binx,biny);
		weight = flugg_nue_rhc.GetBinContent(binx);
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period." << std::endl;
	}

	return weight;
}

double Utility::getFluggWeightNuebar(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	double weight;

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		int binx = flugg_nuebar_fhc.GetXaxis()->FindBin(energy);
		//int biny = flugg_nuebar_fhc.GetYaxis()->FindBin(angle);
		//weight = flugg_nuebar_fhc.GetBinContent(binx,biny);
		weight = flugg_nuebar_fhc.GetBinContent(binx);
	}
	// rhc
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) {
		int binx = flugg_nuebar_rhc.GetXaxis()->FindBin(energy);
		//int biny = flugg_nuebar_rhc.GetYaxis()->FindBin(angle);
		//weight = flugg_nuebar_rhc.GetBinContent(binx,biny);
		weight = flugg_nuebar_rhc.GetBinContent(binx);
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period." << std::endl;
	}

	return weight;
}

double Utility::getFluggWeightNumu(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	double weight;

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		int binx = flugg_numu_fhc.GetXaxis()->FindBin(energy);
		//int biny = flugg_numu_fhc.GetYaxis()->FindBin(angle);
		//weight = flugg_numu_fhc.GetBinContent(binx,biny);
		weight = flugg_numu_fhc.GetBinContent(binx);
	}
	// rhc
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) {
		int binx = flugg_numu_rhc.GetXaxis()->FindBin(energy);
		//int biny = flugg_numu_rhc.GetYaxis()->FindBin(angle);
		//weight = flugg_numu_rhc.GetBinContent(binx,biny);
		weight = flugg_numu_rhc.GetBinContent(binx);
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period." << std::endl;
	}

	return weight;
}

double Utility::getFluggWeightNumubar(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	double weight;

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		int binx = flugg_numubar_fhc.GetXaxis()->FindBin(energy);
		//int biny = flugg_numubar_fhc.GetYaxis()->FindBin(angle);
		//weight = flugg_numubar_fhc.GetBinContent(binx,biny);
		weight = flugg_numubar_fhc.GetBinContent(binx);
	}
	// rhc
	else if (runPeriod == Utility::kRun1b || runPeriod == Utility::kRun2b || runPeriod == Utility::kRun3b || runPeriod == Utility::kRun4ab) {
		int binx = flugg_numubar_rhc.GetXaxis()->FindBin(energy);
		//int biny = flugg_numubar_rhc.GetYaxis()->FindBin(angle);
		//weight = flugg_numubar_rhc.GetBinContent(binx,biny);
		weight = flugg_numubar_rhc.GetBinContent(binx);
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period." << std::endl;
	}

	return weight;
}

// ------------------------------------------------------------------------------

// Functions to get CV re-weight ratios
void Utility::loadReweightRatioHistograms() {

	// load beamline variations TH2D
	TFile *ReweightFile = NULL;
	ReweightFile = new TFile("/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/AngleWeightCheck/ratios_set1fhc_sanitycheck_angle.root");
	
	reweight_ratio_nue_fhc = *(TH2D*)ReweightFile->Get("ratio_cv_fhc_enu_angle_nue");
	reweight_ratio_nuebar_fhc = *(TH2D*)ReweightFile->Get("ratio_cv_fhc_enu_angle_nuebar");
	reweight_ratio_numu_fhc = *(TH2D*)ReweightFile->Get("ratio_cv_fhc_enu_angle_numu");
	reweight_ratio_numubar_fhc = *(TH2D*)ReweightFile->Get("ratio_cv_fhc_enu_angle_numubar");

	ReweightFile->Close();
}

double Utility::getReweightRatioNue(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	double weight;

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		int binx = reweight_ratio_nue_fhc.GetXaxis()->FindBin(energy);
		int biny = reweight_ratio_nue_fhc.GetYaxis()->FindBin(angle);
		weight = reweight_ratio_nue_fhc.GetBinContent(binx,biny);
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period. Reweight ratios only available for FHC." << std::endl;
	}

	return weight;
}

double Utility::getReweightRatioNuebar(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	double weight;

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		int binx = reweight_ratio_nuebar_fhc.GetXaxis()->FindBin(energy);
		int biny = reweight_ratio_nuebar_fhc.GetYaxis()->FindBin(angle);
		weight = reweight_ratio_nuebar_fhc.GetBinContent(binx,biny);
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period. Reweight ratios only available for FHC." << std::endl;
	}

	return weight;
}

double Utility::getReweightRatioNumu(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	double weight;

	// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		int binx = reweight_ratio_numu_fhc.GetXaxis()->FindBin(energy);
		int biny = reweight_ratio_numu_fhc.GetYaxis()->FindBin(angle);
		weight = reweight_ratio_numu_fhc.GetBinContent(binx,biny);
	}
	else {
		std::cout << "Error [Utility::getWeights]: invalid run period. Reweight ratios only available for FHC." << std::endl;
	}

	return weight;
}

double Utility::getReweightRatioNumubar(float energy, float angle, Utility::RunPeriodEnums runPeriod) {
	
	double weight;

	/// fhc
	if (runPeriod == Utility::kRun1a || runPeriod == Utility::kRun2a || runPeriod == Utility::kRun4cd || runPeriod == Utility::kRun5) {
		int binx = reweight_ratio_numubar_fhc.GetXaxis()->FindBin(energy);
		int biny = reweight_ratio_numubar_fhc.GetYaxis()->FindBin(angle);
		weight = reweight_ratio_numubar_fhc.GetBinContent(binx,biny);
	}
	else {
		std::cout << "Error [Utility::getReweightRatio]: invalid run period. CV reweight ratios only available for FHC." << std::endl;
	}

	return weight;
}