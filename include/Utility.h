#ifndef UTILITY_H
#define UTILITY_H

// class containing utility functions and enum definitions
// allows suite of utilities to be included where required

#include <vector>

#include <TH1F.h>
#include <TH2F.h>

class Utility {

public:	
	
	// ----------------------------------
	// Constructor
	Utility(bool loadBeamLineWeights = false, bool loadPPFXCVWeights = false, bool loadFluggWeights = false, bool loadReweightRatios = false);

	// Destructor
	~Utility(){};

	// ----------------------------------
	// Enum definitions

	// File types
	enum FileTypeEnums{
		kIntrinsic,
		kCCNCPiZero,
		kMC,
		kEXT,
		kDirt,
		kDetVar,
		kData,
		kFakeData
	};

	// Run periods
	enum RunPeriodEnums {
		kRun1a,
		kRun1b,
		kRun2a,
		kRun2b,
		kRun3b,
		kRun4ab,
		kRun4cd,
		kRun5
	};

	// BDT types
	enum BDTEnums {
		kElectronPhoton,
		kPionProton,
		kPionProtonAlternate,
		kMuonPion
	};

	// Event classification enums
	enum ClassificationEnums {
		kNC1pi,
		kNCOther,
		kCCNumu0pi,
		kCCNumu1pi,
		kCCNumuNpi,
		kCCNumuOther,
		kCCNue,
		kOutFV,
		kOutOfCryo,
		kBeamOff,
		kBeamOn
	};

	// Interaction type enums
	enum InteractionEnums {
		kQE,
		kCOH,
		kRES,
		kDIS,
		kMEC,
		kOther
	};

	// Particle type enums
	enum ParticleEnums {		
		kMuon,
		kPion,
		kPionDecay,
		kPionAbsorption0p,
		kPionAbsorptionNp,
		kPionChargeExchange,
		kPionBelowThreshold,
		kProton,
		kCosmic
	};
	

	// PlotVariable enums
	enum PlotVariableEnums {
		kSliceID,
		kNuE,
		kNShower,
		kNTrack,
		kVertexInFV,
		kContainedFraction,
		kAssociatedHitsFraction,
		kTopologicalScore,
		kCosmicImpactParameter,
		kShowerScore,
		kShowerEnergy,
		kHitRatio,
		kMoliereAverage,
		kShowerEnergyRatio,
		kSecondShowerEnergy,
		kSecondShowerPFPGeneration,
		kShowerDistance,
		kShowerdEdxMax,
		kNeutralPionInvariantMass,
		kTrackLength,
		kTrackScore,
		kTrackDistance,
		kTrackTheta,
		kTrackPhi,
		kTrackBraggPion,
		kTrackBraggP,
		kTrackBraggMu,
		kTrackBraggMIP,
		kLLRPID,
		kTrackdEdx,
		kTrackWiggliness,
		kPi0InvMassDifference,
		kTrackShowerAngle,
		kNTrackDaughters,
		kDaughterEnergy,
		kTrackEndSpacepoints,
		kTrackEnergy,
		kTrackMomentum,
		kSecondShowerNHit,
		kShowerLLRPID,
		kPionETrue,
		kNeutrinoETrue,
		kElectronETrue,
		kElectronNeutralPionBDT,
		kPionProtonBDT,
		kFHC,
		kRHC,
		kNumberProtons,
		kNumberProtonsTrue,
		kElectronBeta,
		kPionBeta,
		kOpeningAngle,
		kNuVtxX,
		kNuVtxY,
		kNuVtxZ,
		kWCNumuCCFlag,
		kisFC,
		kLArPID_el,
		kLArPID_ph,
		kLArPID_mu,
		kLArPID_pi,
		kLArPID_pr,
		kLArPID_prim,
		kLArPID_chgd,
		kLArPID_ntrl,
		kLArPID_pdg,
		kNuGraphMipFrac,
		kNuGraphHipFrac,
		kNuGraphShrFrac,
		kNuGraphMclFrac,
		kNuGraphDfsFrac,
		kNBlips25cm,
		kNBlips50cm,
		kNBlips10cm,
		kSingleBin,
		kLanternVtxScore,
		kLanternVtxContainment,
		kLanternVtxFracHitsOnCosmic	
	};

	// ---------------------------------- 

	// Function to check if point is within FV
	bool inFV(const double x, const double y, const double z) const;

	// Function to check if track is exiting
	bool isExiting(const double x, const double y, const double z) const;

	// Function to check whether number is valid
	bool isNumber(float input) const;

	// Function to scale histogram by bin width
	void scaleByBinWidth(TH1F* hist);

	// Functions to get beamline variation weights
	void loadBeamLineVariationHistograms();
	std::vector<double> getWeightsNue(float energy, float angle, Utility::RunPeriodEnums runPeriod);
	std::vector<double> getWeightsNuebar(float energy, float angle, Utility::RunPeriodEnums runPeriod);
	std::vector<double> getWeightsNumu(float energy, float angle, Utility::RunPeriodEnums runPeriod);
	std::vector<double> getWeightsNumubar(float energy, float angle, Utility::RunPeriodEnums runPeriod);

	// Functions to get PPFX CV weights for fake data
	void loadPPFXCVHistograms();
	double getCVWeightNue(float energy, float angle, Utility::RunPeriodEnums runPeriod);
	double getCVWeightNuebar(float energy, float angle, Utility::RunPeriodEnums runPeriod);
	double getCVWeightNumu(float energy, float angle, Utility::RunPeriodEnums runPeriod);
	double getCVWeightNumubar(float energy, float angle, Utility::RunPeriodEnums runPeriod);

	// Functions to get Flugg weights for fake data
	void loadFluggHistograms();
	double getFluggWeightNue(float energy, float angle, Utility::RunPeriodEnums runPeriod);
	double getFluggWeightNuebar(float energy, float angle, Utility::RunPeriodEnums runPeriod);
	double getFluggWeightNumu(float energy, float angle, Utility::RunPeriodEnums runPeriod);
	double getFluggWeightNumubar(float energy, float angle, Utility::RunPeriodEnums runPeriod);

	// Functions to get CV re-weight ratios [testing]
	void loadReweightRatioHistograms();
	double getReweightRatioNue(float energy, float angle, Utility::RunPeriodEnums runPeriod);
	double getReweightRatioNuebar(float energy, float angle, Utility::RunPeriodEnums runPeriod);
	double getReweightRatioNumu(float energy, float angle, Utility::RunPeriodEnums runPeriod);
	double getReweightRatioNumubar(float energy, float angle, Utility::RunPeriodEnums runPeriod);


private:

	// beamline variation weight histograms
	// fhc
	std::vector<TH2F> h_nue_fhc;
	std::vector<TH2F> h_nuebar_fhc;
	std::vector<TH2F> h_numu_fhc;
	std::vector<TH2F> h_numubar_fhc;
	// rhc
	std::vector<TH2F> h_nue_rhc;
	std::vector<TH2F> h_nuebar_rhc;
	std::vector<TH2F> h_numu_rhc;
	std::vector<TH2F> h_numubar_rhc;

	// PPFX CV weight histograms
	// fhc
	TH2F ppfx_nue_fhc;
	TH2F ppfx_nuebar_fhc;
	TH2F ppfx_numu_fhc;
	TH2F ppfx_numubar_fhc;
	// rhc
	TH2F ppfx_nue_rhc;
	TH2F ppfx_nuebar_rhc;
	TH2F ppfx_numu_rhc;
	TH2F ppfx_numubar_rhc;

	// Flugg CV histograms
	// fhc
	TH1F flugg_nue_fhc;
	TH1F flugg_nuebar_fhc;
	TH1F flugg_numu_fhc;
	TH1F flugg_numubar_fhc;
	// rhc
	TH1F flugg_nue_rhc;
	TH1F flugg_nuebar_rhc;
	TH1F flugg_numu_rhc;
	TH1F flugg_numubar_rhc;

	// Re-weight ratio histograms
	// fhc
	TH2D reweight_ratio_nue_fhc;
	TH2D reweight_ratio_nuebar_fhc;
	TH2D reweight_ratio_numu_fhc;
	TH2D reweight_ratio_numubar_fhc;
}; 

#endif