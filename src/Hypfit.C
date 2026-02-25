#include "Hypfit.h"

Hypfit::Hypfit(){
  map_PhysdEdx[13] = new PhysdEdx(13); // == muon
  map_PhysdEdx[211] = new PhysdEdx(211); // == charged pion
  map_PhysdEdx[321] = new PhysdEdx(321); // == charged kaon
  map_PhysdEdx[2212] = new PhysdEdx(2212); // == proton
}

Hypfit::~Hypfit(){
}

double Hypfit::Gaussian(const std::vector<double> & dEdx, const std::vector<double> & ResRange, int PID){

  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }
  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  double dEdx_truncate_upper = dEdx_truncate_upper_pion;
  double dEdx_truncate_bellow = dEdx_truncate_bellow_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    dEdx_truncate_upper = dEdx_truncate_upper_proton;
    dEdx_truncate_bellow = dEdx_truncate_bellow_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialize
  double best_additional_res_length = -0.1;
  double best_chi2 = 99999.;

  int this_N_calo = dEdx.size();
  if(this_N_calo <= 10){
    return -8888.; // == Too small number of hits
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;

  // == Fit
  for(int i = 0; i < res_length_trial; i++){
    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_chi2 = 0.;
    for(int j = N_skip; j < this_N_hits - N_skip; j++){ // == Do not use first and last N_skip hits
      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_PhysdEdx[PID]->KEFromRangeSpline(this_res_length);
      double dEdx_theory = map_PhysdEdx[PID]->meandEdx(this_KE);
      double dEdx_measured = dEdx.at(j);
      if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue; // == Truncate
      // == Gaussian approx.
      //double dEdx_theory_err = dEdx_theory * 0.02;
      this_chi2 += pow(dEdx_measured - dEdx_theory, 2);
    }
    this_chi2 = this_chi2 / (this_N_hits + 0.); // == chi2 / n.d.f
    if(this_chi2 < best_chi2){
      best_chi2 = this_chi2;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
  }

  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length;

  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){ // == Fit failed : no mimumum
    return -7777.;
  }
  else if(best_chi2 > 99990.){ // == Fit failed : best_chi2 > 99990."
    return -6666.;
  }
  else if(best_chi2 < 1.0e-11){ // == Fit failed : best_chi2 < 1.0e-11
    return -5555.;
  }

  double best_KE = map_PhysdEdx[PID] -> KEFromRangeSpline(best_total_res_length);
  return best_KE;
}

double Hypfit::Likelihood(const std::vector<double> & dEdx, const std::vector<double> & ResRange, int PID){

  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }
  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  double dEdx_truncate_upper = dEdx_truncate_upper_pion;
  double dEdx_truncate_bellow = dEdx_truncate_bellow_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    dEdx_truncate_upper = dEdx_truncate_upper_proton;
    dEdx_truncate_bellow = dEdx_truncate_bellow_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialize
  double best_additional_res_length = -0.1;
  double best_m2lnL = 99999.;

  int this_N_calo = dEdx.size();
  if(this_N_calo <= 10){
    return -8888.; // == Too small number of hits
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;

  double default_m2lnL = -1.;
  // == Fit
  for(int i = 0; i < res_length_trial; i++){
    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_m2lnL = 0.;
    for(int j = N_skip; j < this_N_hits - N_skip; j++){ // == Do not use first and last N_skip this
      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_PhysdEdx[PID]->KEFromRangeSpline(this_res_length);
      double dEdx_measured = dEdx.at(j);
      if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue; // == Truncate
      // == Likelihood
      double this_pitch = fabs(ResRange.at(j - 1) - ResRange.at(j + 1)) / 2.0;
      double this_likelihood = map_PhysdEdx[PID] -> dEdx_PDF(this_KE, this_pitch, dEdx_measured);
      if(this_likelihood > 1e-6) this_m2lnL += (-2.0) * log(this_likelihood);
    }
    if(this_m2lnL < best_m2lnL){
      best_m2lnL = this_m2lnL;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
    if(i == 0){
      default_m2lnL = this_m2lnL;
    }
  }

  // == Result
  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length; // == [cm]
  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){ // == Fit failed : no mimumum
    return -7777.;
  }
  else if(best_m2lnL > 99990.){ // == Fit failed : best likelihood > 99990.
    return -6666.;
  }
  else if(fabs(best_m2lnL) < 1.0e-11){ // == Fit failed : |best_m2lnL| < 1.0e-11, no valid likelihood value
    return -5555.;
  }

  // == Return
  double best_KE = map_PhysdEdx[PID] -> KEFromRangeSpline(best_total_res_length);
  return best_KE;
}
