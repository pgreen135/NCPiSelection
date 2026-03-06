#include "Hypfit.h"

#include <TH1F.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TStyle.h>
#include <TGraph.h>

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
  
  // Make a plot for debugging 
  TH1F *likelihoodPlot = new TH1F("", "", 300, 0, 300);
  
  // Store actual used dEdx and ResRange for debugging
  std::vector<double> used_dEdx;
  std::vector<double> used_ResRange;

  // Store actual likelihood values and range for debugging
  std::vector<double> likelihood_values;
  std::vector<double> likelihood_ResRange;

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

      if (i == 0) {
        used_dEdx.push_back(dEdx_measured);
        used_ResRange.push_back(this_res_length);
      }

      // == Likelihood
      double this_pitch = fabs(ResRange.at(j - 1) - ResRange.at(j + 1)) / 2.0;
      double this_likelihood = map_PhysdEdx[PID] -> dEdx_PDF(this_KE, this_pitch, dEdx_measured);
      double this_likelihood_max = map_PhysdEdx[PID] -> dEdx_PDF_max(this_KE, this_pitch);

      
      if(this_likelihood < 1e-6) this_likelihood = 1e-6; // == Avoid zero likelihood
      //if(this_likelihood_max < 1e-6) this_likelihood_max = 1e-6;
      //if (i == 10) std::cout << "(i,j) = (" << i << "," << j << "), this_KE = " << this_KE << ", dEdx_measured = " << dEdx_measured << ", this_likelihood: " << this_likelihood << ", this_likelihood_max: " << this_likelihood_max << ", m2lnL: " << (-2.0) * (log(this_likelihood) - log(this_likelihood_max)) << std::endl;
      this_m2lnL += (-2.0) * (log(this_likelihood) - log(this_likelihood_max));
    }
    
    // scale to number of hits
    this_m2lnL = this_m2lnL / (this_N_hits - 2*N_skip);

    //std::cout << "Res Length = " << ResRange.at(this_N_calo - 1) + this_additional_res_length << ", likelihood = " << this_m2lnL << std::endl;
    likelihoodPlot->Fill(ResRange.at(this_N_calo - 1) + this_additional_res_length, this_m2lnL);

    likelihood_values.push_back(this_m2lnL);

    likelihood_ResRange.push_back(ResRange.at(this_N_calo - 1) + this_additional_res_length);


    if(this_m2lnL < best_m2lnL){
      best_m2lnL = this_m2lnL;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
    if(i == 0){
      default_m2lnL = this_m2lnL;
    }
  }

  // draw gra
  /*
  TCanvas *c1 = new TCanvas();
  c1->cd();
  
  likelihoodPlot->GetYaxis()->SetTitle("Output -2lnL");
  likelihoodPlot->GetXaxis()->SetTitle("Residual Range [cm]");
  likelihoodPlot->Draw("HIST");

  // draw line for default hypothesis
  TLine *line_default = new TLine(best_additional_res_length + ResRange.at(this_N_calo - 1), 0, best_additional_res_length + ResRange.at(this_N_calo - 1), likelihoodPlot->GetMaximum());
  line_default->SetLineColor(kOrange+7);
  line_default->SetLineStyle(9);
  line_default->SetLineWidth(3);
  line_default->Draw("SAME");

  // draw line for truth KE
  double truth_p = 0.429131 * 1000; // MeV/c 0.429131
  double truth_KE = map_PhysdEdx[PID] -> MomentumtoKE(truth_p);
  double truth_res_range = map_PhysdEdx[PID] -> RangeFromKESpline(truth_KE);
  TLine *line_truth = new TLine(truth_res_range, 0, truth_res_range, likelihoodPlot->GetMaximum());
  line_truth->SetLineColor(kGreen+2);
  line_truth->SetLineStyle(9);
  line_truth->SetLineWidth(3);
  line_truth->Draw("SAME");
 
  c1->SaveAs("plots/likelihoodPlot.root");

  // Draw scatter plot of used dEdx and ResRange for debugging
  TCanvas *c2 = new TCanvas();
  c2->cd();
  TGraph *graph = new TGraph(used_ResRange.size(), &used_ResRange[0], &used_dEdx[0]);
  graph->GetXaxis()->SetTitle("Residual Range [cm]");
  graph->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
  graph->SetTitle("Used dE/dx Profile");
  graph->SetMarkerStyle(4);
  graph->SetMarkerSize(2);
  graph->Draw("AP");
  c2->SaveAs("plots/dEdx_vs_ResRange.root");
  */

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
  double best_p = map_PhysdEdx[PID] -> KEtoMomentum(best_KE);
  return best_p;
}
