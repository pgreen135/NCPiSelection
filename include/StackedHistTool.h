#ifndef STACKEDHISTTOOL
#define STACKEDHISTTOOL

// class to contain and draw stacked histogram based on event classification

#include "../include/Utility.h"

#include "TH1F.h"
#include "THStack.h"
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>

#include <vector>
#include <string>
#include <iostream>

class StackedHistTool{

 public:
  StackedHistTool(std::string histname, std::string title, int nbins, double lowlimit, double highlimit, const Utility &utility, std::string mode = "Topology"); // Constructor
  void Fill(Utility::ClassificationEnums topology, double value);
  void Fill(Utility::InteractionEnums topology, double value);
  void Fill(Utility::ParticleEnums particle, double value);
  void Fill(Utility::ClassificationEnums topology, double value, double weight);
  void Fill(Utility::InteractionEnums topology, double value, double weight);
  void Fill(Utility::ParticleEnums particle, double value, double weight);
  
  std::string topologyenum2str(Utility::ClassificationEnums topology);
  std::string interactionenum2str(Utility::InteractionEnums interaction);
  std::string particleenum2str(Utility::ParticleEnums particle);
  std::string PlotVariableEnum2str(Utility::PlotVariableEnums plotvariable);
  
  void DrawStack(TCanvas *c1, Utility::PlotVariableEnums plotvariable);
 
  void PrintHistIntegrals();
  double GetTotalIntegral();
  double GetSignalIntegral();
  double GetBinValue(const int true_idx, const int reco_idx) const;
  void PrintEventIntegrals();
 
  void SaveHistStack();

 protected:
  Utility _utility;

  std::string configMode;

  int nHists;
  std::vector<Utility::ClassificationEnums> hist_order; // For keeping track of which hist goes with which topology
  std::vector<Utility::InteractionEnums> hist_order_interaction; // For keeping track of which hist goes with which topology
  std::vector<Utility::ParticleEnums> hist_order_particle; // For keeping track of which hist goes with which particle

  THStack *stack;
  TH1F *hists[25];

  std::vector<double> eventCount;  // to keep track of number of events of each classification 

  double invalid_total_x;
  double invalid_total_y;

  void InstantiateHistOrder();
  void StyleHistsStack();
  unsigned int GetHistN(Utility::ClassificationEnums topology);
  unsigned int GetHistN(Utility::InteractionEnums topology);
  unsigned int GetHistN(Utility::ParticleEnums particle);
  Utility::ClassificationEnums GetTopologyFromHistN(unsigned int hist_n);
  Utility::InteractionEnums GetInteractionFromHistN(unsigned int hist_n);
  Utility::ParticleEnums GetParticleFromHistN(unsigned int hist_n);
};

#endif