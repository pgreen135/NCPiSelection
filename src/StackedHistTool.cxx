#include "../include/StackedHistTool.h"

// -------------------------- Constructor -------------------------- //
// Intended to be very general so you can make histograms of whatever you like
StackedHistTool::StackedHistTool(std::string histname, std::string title, int nbins, double lowlimit, double highlimit, const Utility &utility, std::string mode): _utility{ utility }
{

  // Check mode configuration
  if (mode != "Topology" && mode != "Interaction" && mode != "Particle") {
    std::cout << "Error [StackedHistTool]: Invalid mode configured." << std::endl;
    exit(1);
  }
  else {
    configMode = mode;
  }

  InstantiateHistOrder();
    
  stack = new THStack(histname.c_str(),title.c_str());
  for (int i_hist=0; i_hist < nHists; i_hist++){
    std::string histname_i = std::string(histname)+std::string("_")+std::to_string(i_hist);
    hists[i_hist] = new TH1F(histname_i.c_str(),"",nbins,lowlimit,highlimit);
  }

  // Style the histograms
  StyleHistsStack();

  invalid_total_x = 0.;

  // set event counts to zero
  eventCount = std::vector<double>(nHists, 0.0);
}

// Set histogram order (moved to be in a separate function from the constructor to allow for overload of constructor for 1D or 2D histograms)
void StackedHistTool::InstantiateHistOrder(){

  if (configMode == "Topology" ) {
    // Topology categories for the histograms  
    hist_order.push_back(_utility.kNC1pi);
    hist_order.push_back(_utility.kNCOther);
    hist_order.push_back(_utility.kCCNumu0pi);
    hist_order.push_back(_utility.kCCNumu1pi);
    hist_order.push_back(_utility.kCCNumuNpi);
    hist_order.push_back(_utility.kCCNumuOther);
    hist_order.push_back(_utility.kCCNue);    
    hist_order.push_back(_utility.kOutFV);
    hist_order.push_back(_utility.kOutOfCryo);
    hist_order.push_back(_utility.kBeamOff);
    nHists = hist_order.size();
  }
  else if (configMode == "Interaction") {
    hist_order_interaction.push_back(_utility.kQE);
    hist_order_interaction.push_back(_utility.kCOH);
    hist_order_interaction.push_back(_utility.kRES);
    hist_order_interaction.push_back(_utility.kDIS);
    hist_order_interaction.push_back(_utility.kMEC);
    hist_order_interaction.push_back(_utility.kOther);

    nHists = hist_order_interaction.size();
  }
  else if (configMode == "Particle") {
    hist_order_particle.push_back(_utility.kMuon);
    hist_order_particle.push_back(_utility.kPionDecay);
    hist_order_particle.push_back(_utility.kPionAbsorption0p);
    hist_order_particle.push_back(_utility.kPionAbsorptionNp);
    hist_order_particle.push_back(_utility.kPionChargeExchange);
    hist_order_particle.push_back(_utility.kPionBelowThreshold);
    hist_order_particle.push_back(_utility.kProton);
    hist_order_particle.push_back(_utility.kCosmic);

    nHists = hist_order_particle.size();
  }
}

// -------------------------- Function to fill the correct histogram -------------------------- //
void StackedHistTool::Fill(Utility::ClassificationEnums topology, double value)
{
  unsigned int n_hist = StackedHistTool::GetHistN(topology);
  hists[n_hist]->Fill(value);
  eventCount[n_hist] += 1;
}
void StackedHistTool::Fill(Utility::InteractionEnums interaction, double value)
{
  unsigned int n_hist = StackedHistTool::GetHistN(interaction);
  hists[n_hist]->Fill(value);
  eventCount[n_hist] += 1;
}
void StackedHistTool::Fill(Utility::ParticleEnums particle, double value)
{
  unsigned int n_hist = StackedHistTool::GetHistN(particle);
  hists[n_hist]->Fill(value);
  eventCount[n_hist] += 1;
}

// -------------------------- Function to fill the correct histogram -------------------------- //
// Overloaded to allow for weights
void StackedHistTool::Fill(Utility::ClassificationEnums topology, double value, double weight)
{
  unsigned int n_hist = StackedHistTool::GetHistN(topology);
  hists[n_hist]->Fill(value, weight);
  eventCount[n_hist] += weight;
}
void StackedHistTool::Fill(Utility::InteractionEnums interaction, double value, double weight)
{
  unsigned int n_hist = StackedHistTool::GetHistN(interaction);
  hists[n_hist]->Fill(value, weight);
  eventCount[n_hist] += weight;
}
void StackedHistTool::Fill(Utility::ParticleEnums particle, double value, double weight)
{
  unsigned int n_hist = StackedHistTool::GetHistN(particle);
  hists[n_hist]->Fill(value, weight);
  eventCount[n_hist] += weight;
}

// -------------------------- Function to draw the histograms -------------------------- //
void StackedHistTool::DrawStack(TCanvas *c1, Utility::PlotVariableEnums plotvariable)
{

	// Next: add histograms to the stack and make TLegend
	// Only do this for histograms that have entries
	TLegend *leg = new TLegend(0.1,0.88,0.9,0.99);
	// leg->SetTextFont(132);
	leg->SetLineColor(kWhite);
	leg->SetTextAlign(12);
	leg->SetNColumns(6);

	TPaveText *pt = new TPaveText(0.4,0.88,0.90,0.91,"NDC NB");
	//pt->SetLineColor(kWhite);
	//pt->SetFillColor(kWhite);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
	pt->SetTextAlign(12);

	double underflow_total = 0.;
	double overflow_total = 0.;

	for (int i_hist = 0; i_hist < nHists; i_hist++) {
	  if (hists[i_hist]->GetEntries() == 0)
		continue;

	  stack->Add(hists[i_hist]);
    if (configMode == "Topology") {
	   Utility::ClassificationEnums topology_for_legend = StackedHistTool::GetTopologyFromHistN((unsigned int)i_hist);
	   leg->AddEntry(hists[i_hist], topologyenum2str(topology_for_legend).c_str(), "f");
    }
    else if (configMode == "Interaction") {
     Utility::InteractionEnums interaction_for_legend = StackedHistTool::GetInteractionFromHistN((unsigned int)i_hist);
     leg->AddEntry(hists[i_hist], interactionenum2str(interaction_for_legend).c_str(), "f");
    }
    else if (configMode == "Particle") {
      Utility::ParticleEnums particle_for_legend = StackedHistTool::GetParticleFromHistN((unsigned int)i_hist);
      leg->AddEntry(hists[i_hist], particleenum2str(particle_for_legend).c_str(), "f");
    }

	  underflow_total += hists[i_hist]->GetBinContent(0);
	  overflow_total += hists[i_hist]->GetBinContent(hists[i_hist]->GetXaxis()->GetNbins()+1);
	}

	//pt->AddText(TString::Format("Underflow (Invalid): %.2f (%.2f) \t Overflow: %.2f",underflow_total,invalid_total_x, overflow_total).Data());
	//pt->AddText(TString::Format("Overflow: %.2f",overflow_total).Data());

	c1->cd();
	c1->SetTopMargin(0.14);

	stack->Draw("hist");

	stack->GetXaxis()->SetTitle(PlotVariableEnum2str(plotvariable).c_str());
	stack->GetYaxis()->SetTitle("Entries");
	
  stack->GetXaxis()->SetTitleSize(0.045);
  stack->GetYaxis()->SetTitleSize(0.045);

  c1->Modified();
	
	leg->Draw();
	pt->Draw();
}

// -------------------------- Function to style the histograms -------------------------- //
// Private: only called by DrawStack function in this file
void StackedHistTool::StyleHistsStack()
{ 
  if (configMode == "Topology") {

    // Set fill colors for all histograms
    hists[0] ->SetFillColor(kCyan+3); // kNC1pi
    hists[1] ->SetFillColor(kMagenta+2); // kNCOther
    hists[2] ->SetFillColor(kBlue-3); // kCCNumu0pi
    hists[3] ->SetFillColor(kBlue-6); // kCCNumu1pi
    hists[4] ->SetFillColor(kBlue-5); // kCCNumuNpi
    hists[5] ->SetFillColor(kBlue-10); // kCCNumuOther
    hists[6] ->SetFillColor(kRed); // kCCNue
    hists[7]->SetFillColor(kGray+1); // Out FV
    hists[8]->SetFillColor(kRed-3); // Out of Cryo
    hists[9]->SetFillColor(kGray+3); // Beam Off

    // Set line color for all histograms
    hists[0] ->SetLineColor(kCyan+3); // kNC1pi
    hists[1] ->SetLineColor(kMagenta+2); // kNCOther
    hists[2] ->SetLineColor(kBlue-3); // kCCNumu0pi
    hists[3] ->SetLineColor(kBlue-6); // kCCNumu1pi
    hists[4] ->SetLineColor(kBlue-5); // kCCNumuNpi
    hists[5] ->SetLineColor(kBlue-10); // kCCNumuOther
    hists[6] ->SetLineColor(kRed); // kCCNue
    hists[7]->SetLineColor(kGray+1); // Out FV
    hists[8]->SetLineColor(kRed-3); // Out of Cryo
    hists[9]->SetLineColor(kGray+3); // Beam Off

    // Set fill style
    hists[8]->SetFillStyle(3004); // Out of Cryo
    hists[9]->SetFillStyle(3005); // Beam Off
  }
  else if (configMode == "Interaction") {
    hists[0] ->SetFillColor(kBlue-3); // QE
    hists[1] ->SetFillColor(kMagenta); // COH
    hists[2] ->SetFillColor(kGreen+2); // RES
    hists[3] ->SetFillColor(kRed-3); // DIS
    hists[4] ->SetFillColor(kOrange-3); // MEC
    hists[5] ->SetFillColor(kGray); // Other

    hists[0] ->SetLineColor(kBlue-3); // QE
    hists[1] ->SetLineColor(kMagenta); // COH
    hists[2] ->SetLineColor(kGreen+2); // RES
    hists[3] ->SetLineColor(kRed-3); // DIS
    hists[4] ->SetLineColor(kOrange-3); // MEC
    hists[5] ->SetLineColor(kGray); // Other
  }
  else if (configMode == "Particle") {
    hists[0] ->SetFillColor(kBlue-3); // muon
    hists[1] ->SetFillColor(kCyan+3); // pion decay
    hists[2] ->SetFillColor(kGreen+2); // pion absorb 0p
    hists[3] ->SetFillColor(kGreen-5); // pion absorb Np
    hists[4] ->SetFillColor(kViolet-3); // pion charge exchange
    hists[5] ->SetFillColor(kCyan-8); // pion below threshold
    hists[6] ->SetFillColor(kRed-3); // proton
    hists[7] ->SetFillColor(kGray+3); // other

    // Set line color for all histograms
    hists[0] ->SetLineColor(kBlue-3); // muon
    hists[1] ->SetLineColor(kCyan+3); // pion decay
    hists[2] ->SetLineColor(kGreen+2); // pion absorb
    hists[3] ->SetLineColor(kGreen-5); // pion absorb Np
    hists[4] ->SetLineColor(kViolet-3); // pion charge exchange
    hists[5] ->SetLineColor(kCyan-8); // pion below threshold
    hists[6] ->SetLineColor(kRed-3); // proton
    hists[7] ->SetLineColor(kGray+3); // other

    // Set fill style
    hists[7]->SetFillStyle(3004); // other
  }
  
}


// ---------------------- Function to get histogram number for given topology ---------------------- //
// Private: only called by functions in this class
unsigned int StackedHistTool::GetHistN(Utility::ClassificationEnums topology)
{
  unsigned int HistN;
  bool found_hist=false;

  for (int i=0; i<nHists; i++){
    if (hist_order.at(i) == topology){
      HistN = i;
      found_hist = true;
      break;
    }
  }

  if (!found_hist){
    std::cout << "[ERROR: StackedHistTool.h] Could not find histogram for topology " << topology << std::endl;
    HistN = nHists-1;
  }

  return HistN;
}
unsigned int StackedHistTool::GetHistN(Utility::InteractionEnums interaction)
{
  unsigned int HistN;
  bool found_hist=false;

  for (int i=0; i<nHists; i++){
    if (hist_order_interaction.at(i) == interaction){
      HistN = i;
      found_hist = true;
      break;
    }
  }

  if (!found_hist){
    std::cout << "[ERROR: StackedHistTool.h] Could not find histogram for interaction " << interaction << std::endl;
    HistN = nHists-1;
  }

  return HistN;
}

unsigned int StackedHistTool::GetHistN(Utility::ParticleEnums particle)
{
  unsigned int HistN;
  bool found_hist=false;

  for (int i=0; i<nHists; i++){
    if (hist_order_particle.at(i) == particle){
      HistN = i;
      found_hist = true;
      break;
    }
  }

  if (!found_hist){
    std::cout << "[ERROR: StackedHistTool.h] Could not find histogram for particle " << particle << std::endl;
    HistN = nHists-1;
  }

  return HistN;
}


// ---------------------- Function to get topology code for given histogram number ---------------------- //
// Private: only called by functions in this class
Utility::ClassificationEnums StackedHistTool::GetTopologyFromHistN(unsigned int hist_n)
{
  return hist_order.at(hist_n);
}
Utility::InteractionEnums StackedHistTool::GetInteractionFromHistN(unsigned int hist_n)
{
  return hist_order_interaction.at(hist_n);
}
Utility::ParticleEnums StackedHistTool::GetParticleFromHistN(unsigned int hist_n)
{
  return hist_order_particle.at(hist_n);
}

// ---------------------- Function to get histogram integrals ---------------------- //
void StackedHistTool::PrintHistIntegrals()
{
  // Calculate and print out relative integrals (percentage of events that are each topology)
    // Compute total integral
    double total_integral = 0;
    for (int i_hist=0; i_hist < nHists; i_hist++){
      if (hists[i_hist]->GetEntries() == 0) continue;

      double integral = hists[i_hist]->Integral();
      total_integral += integral;
    }
    std::cout << "Total number of events: " << total_integral << std::endl;
    for (int i_hist=0; i_hist < nHists; i_hist++){
      if (hists[i_hist]->GetEntries() == 0) continue;

      double integral = hists[i_hist]->Integral();

      if (configMode == "Topology") {
        Utility::ClassificationEnums topology = StackedHistTool::GetTopologyFromHistN((unsigned int)i_hist);
        std::cout << "Number of events for toplogy " << topologyenum2str(topology) << ": " << integral << std::endl;// << "\t\t Fraction of total: " << integral/total_integral << std::endl;
      }
      else if (configMode == "Interaction") {
        Utility::InteractionEnums interaction = StackedHistTool::GetInteractionFromHistN((unsigned int)i_hist);
        std::cout << "Number of events for interaction type " << interactionenum2str(interaction) << ": " << integral << std::endl;// << "\t\t Fraction of total: " << integral/total_integral << std::endl;
      }
      else if (configMode == "Particle") {
        Utility::ParticleEnums particle = StackedHistTool::GetParticleFromHistN((unsigned int)i_hist);
        std::cout << "Number of particles of type " << particleenum2str(particle) << ": " << integral << std::endl;// << "\t\t Fraction of total: " << integral/total_integral << std::endl;
      }

    }
}
void StackedHistTool::PrintEventIntegrals()
{
  // Calculate and print out relative integrals (percentage of events that are each topology)
	// Compute total integral
	double total_integral = 0;
	for (int i_hist=0; i_hist < nHists; i_hist++){
	  total_integral += eventCount[i_hist];
	}
  double nue_integral = 0;

	std::cout << "Total number of events: " << total_integral << std::endl;
	for (int i_hist=0; i_hist < nHists; i_hist++){
	  
    if (configMode == "Topology") {
      Utility::ClassificationEnums topology = StackedHistTool::GetTopologyFromHistN((unsigned int)i_hist);
  	  std::cout << "Number of events for toplogy " << topologyenum2str(topology) << ": " << eventCount[i_hist] << std::endl;// << "\t\t Fraction of total: " << integral/total_integral << std::endl;
      
      if (topology <= 4) nue_integral += eventCount[i_hist];

    }
    else if (configMode == "Interaction") {
      Utility::InteractionEnums interaction = StackedHistTool::GetInteractionFromHistN((unsigned int)i_hist);
      std::cout << "Number of events for interaction type " << interactionenum2str(interaction) << ": " << eventCount[i_hist] << std::endl;// << "\t\t Fraction of total: " << integral/total_integral << std::endl;
    }
    else if (configMode == "Particle") {
      Utility::ParticleEnums particle = StackedHistTool::GetParticleFromHistN((unsigned int)i_hist);
      std::cout << "Number entries of particles of type " << particleenum2str(particle) << ": " << eventCount[i_hist] << std::endl;// << "\t\t Fraction of total: " << integral/total_integral << std::endl;
    }
	}

  //std::cout << "Total Nue/Nuebar events: " << nue_integral << std::endl;
}

// ---------------------- Function to get histogram integrals ---------------------- //
double StackedHistTool::GetTotalIntegral()
{
  // Compute total integral
  double total_integral = 0;
  for (int i_hist=0; i_hist < nHists; i_hist++){
    if (hists[i_hist]->GetEntries() == 0) continue;

    double integral = hists[i_hist]->Integral();
    total_integral += integral;
  }
  return total_integral; 
}
double StackedHistTool::GetSignalIntegral()
{
 // Compute signal integral
  double signal_integral = 0;
  for (int i_hist=0; i_hist < 1 ; i_hist++){
    if (hists[i_hist]->GetEntries() == 0) continue;

    double integral = eventCount[i_hist];
    signal_integral += integral;
  }
  return signal_integral; 
}
double StackedHistTool::GetBinValue(const int true_idx, const int reco_idx) const
{
  return hists[true_idx-1]->GetBinContent(reco_idx);
} 

std::string StackedHistTool::topologyenum2str(Utility::ClassificationEnums topology)
{
  std::string returnString = "";

  switch(topology){

  case Utility::kNC1pi:
    returnString = "NC 1#pi^{#pm} (Signal)";
    break;
  case Utility::kNCOther:
    returnString = "NC Other";
    break;
  case Utility::kCCNumu0pi:
    returnString = "#nu_{#mu} CC 0#pi^{#pm}";
    break;
  case Utility::kCCNumu1pi:
    returnString = "#nu_{#mu} CC 1#pi^{#pm}";
    break;
  case Utility::kCCNumuNpi:
    returnString = "#nu_{#mu} CC N#pi^{#pm}";
    break;
  case Utility::kCCNumuOther:
    returnString = "#nu_{#mu} CC Other";
    break;  
  case Utility::kCCNue:
    returnString = "#nu_{e} CC";
    break;
  case Utility::kOutFV:
    returnString = "Out of FV";
    break;
  case Utility::kBeamOff:
    returnString = "Beam Off";
    break;
  case Utility::kOutOfCryo:
    returnString = "Out of Cryo";
    break;
  default:
    std::cout << "[ERROR: StackedHistTool] Could not find string conversion for " << topology << std::endl;
    returnString = "Unknown";
    break; 
  }

  return returnString;
}

std::string StackedHistTool::interactionenum2str(Utility::InteractionEnums interaction)
{
  std::string returnString = "";

  switch(interaction){

  case Utility::kQE:
    returnString = "QE";
    break;
  case Utility::kCOH:
    returnString = "COH";
    break;
  case Utility::kRES:
    returnString = "RES";
    break;
  case Utility::kDIS:
    returnString = "DIS";
    break;
  case Utility::kMEC:
    returnString = "MEC";
    break;
  case Utility::kOther:
    returnString = "Other";
    break;
  default:
    std::cout << "[ERROR: StackedHistTool] Could not find string conversion for " << interaction << std::endl;
    returnString = "Unknown";
    break;
  }

  return returnString;
}

std::string StackedHistTool::particleenum2str(Utility::ParticleEnums particle)
{
  std::string returnString = "";

  switch(particle){

  case Utility::kMuon:
    returnString = "Muon";
    break;
  case Utility::kPionDecay:
    returnString = "Pion (Decay)";
    break;
  case Utility::kPionAbsorption0p:
    returnString = "Pion (Absorption 0p)";
    break;
  case Utility::kPionAbsorptionNp:
    returnString = "Pion (Absorption Np)";
    break;
  case Utility::kPionChargeExchange:
    returnString = "Pion (Charge Exchange)";
    break;
  case Utility::kPionBelowThreshold:
    returnString = "Pion (non-signal)";
    break;
  case Utility::kProton:
    returnString = "Proton";
    break;
  case Utility::kCosmic:
    returnString = "Other/Cosmic";
    break;
  default:
    std::cout << "[ERROR: StackedHistTool] Could not find string conversion for " << particle << std::endl;
    returnString = "Unknown";
    break;
  }

  return returnString;
} 

std::string StackedHistTool::PlotVariableEnum2str(Utility::PlotVariableEnums plotvariable)
{
  std::string returnString = "";

  switch(plotvariable){

  case Utility::kSliceID:
    returnString = "SliceID";
    break;
  case Utility::kNuE:
    returnString = "True Neutrino Energy [GeV]";
    break;
  case Utility::kNShower:
  	returnString = "Shower Multiplicity";
  	break;
  case Utility::kNTrack:
  	returnString = "Track Multiplicity";
  	break;
  case Utility::kVertexInFV:
    returnString = "Vertex in Fiducial Volume";
    break;
  case Utility::kContainedFraction:
  	returnString = "Contained Fraction [Hits in FV / Hits in Slice]";
  	break;
  case Utility::kAssociatedHitsFraction:
    returnString = "Fraction of hits in slice associated with tracks or showers";
    break;
  case Utility::kTopologicalScore:
  	returnString = "Pandora Topological Score";
  	break;
  case Utility::kCosmicImpactParameter:
  	returnString = "Cosmic Impact Parameter [cm]";
  	break;
  case Utility::kShowerScore:
  	returnString = "Shower Score";
  	break;
  case Utility::kShowerEnergy:
    returnString = "Shower Energy [GeV]";
    break;
  case Utility::kSecondShowerEnergy:
    returnString = "Second Shower Energy [GeV]";
    break;
  case Utility::kSecondShowerPFPGeneration:
    returnString = "Second Shower PFP Generation";
    break;
  case Utility::kHitRatio:
  	returnString = "Shower Hit Fracton";
  	break;
  case Utility::kMoliereAverage:
    returnString = "Leading Shower Moliere Average [deg]";
    break;
  case Utility::kShowerEnergyRatio:
    returnString = "Leading Shower Energy / Total Shower Energy";
    break;
  case Utility::kShowerDistance:
    returnString = "Shower Distance Vertex [cm]";
    break;
  case Utility::kShowerdEdxMax:
    returnString = "Leading Shower dE/dx [MeV/cm]";
    break;
  case Utility::kNeutralPionInvariantMass:
    returnString = "Neutral Pion Invariant Mass Difference [MeV]";
    break;  
  case Utility::kTrackLength:
    returnString = "Track Length [cm]";
    break;
  case Utility::kTrackBraggPion:
    returnString = "Pion Bragg Peak Score";
    break;
  case Utility::kTrackBraggP:
    returnString = "Proton Bragg Peak Score";
    break;
  case Utility::kTrackBraggMu:
    returnString = "Muon Bragg Peak Score";
    break;
  case Utility::kTrackBraggMIP:
    returnString = "MIP Bragg Peak Score";
    break;
  case Utility::kLLRPID:
    returnString = "LLR PID Score";
    break; 
  case Utility::kTrackScore:
    returnString = "Pandora Track Score";
    break;
  case Utility::kTrackdEdx:
    returnString = "Track dE/dx [MeV/cm]";
    break;
  case Utility::kTrackWiggliness:
    returnString = "Track Average Deflections";
    break;  
  case Utility::kTrackDistance:
    returnString = "Track Distance Vertex [cm]";
    break; 
  case Utility::kTrackTheta:
    returnString = "Track #theta [deg]";
    break;
  case Utility::kTrackPhi:
    returnString = "Track #phi [deg]";
    break;
  case Utility::kPi0InvMassDifference:
    returnString = "#pi^{0} Invariant Mass Difference [MeV]";
    break;
  case Utility::kTrackShowerAngle:
    returnString = "Track-Shower Opening Angle [deg]";
    break;
  case Utility::kNTrackDaughters:
    returnString = "Number Daughters";
    break;
  case Utility::kDaughterEnergy:
    returnString = "Daughter Total Energy [GeV]";
    break;
  case Utility::kTrackEndSpacepoints:
    returnString = "Track End Spacepoints";
    break;
  case Utility::kTrackEnergy:
    returnString = "Track Energy [GeV]";
    break;
  case Utility::kTrackMomentum:
    returnString = "Track Momentum [GeV/c]";
    break;
  case Utility::kSecondShowerNHit:
    returnString = "Second Shower NHits";
    break;
  case Utility::kShowerLLRPID:
    returnString = "Shower LLR PID Score";
    break;
  case Utility::kPionETrue:
    returnString = "Pion True KE [MeV]";
    break;
  case Utility::kNeutrinoETrue:
    returnString = "True Neutrino Energy [GeV]";
    break;
  case Utility::kElectronETrue:
    returnString = "True Electron Energy [GeV]";
    break;
  case Utility::kElectronNeutralPionBDT:
    returnString = "Electron - Neutral Pion BDT Score";
    break;
  case Utility::kPionProtonBDT:
    returnString = "Pion - Proton BDT Score";
    break;
  case Utility::kFHC:
    returnString = "Forward Horn Current";
    break;
  case Utility::kRHC:
    returnString = "Reverse Horn Current";
    break;
  case Utility::kNumberProtons:
    returnString = "Proton Multiplicity";
    break;
  case Utility::kNumberProtonsTrue:
    returnString = "True Proton Multiplicity";
    break;
  case Utility::kElectronBeta:
    returnString = "Electron cos(#beta)";
    break;
  case Utility::kPionBeta:
    returnString = "Pion cos(#beta)";
    break;
  case Utility::kOpeningAngle:
    returnString = "cos(#beta_{e#pi})";
    break;
  case Utility::kNuVtxX:
    returnString = "Neutrino Vertex X [cm]";
    break;
  case Utility::kNuVtxY:
    returnString = "Neutrino Vertex Y [cm]";
    break;
  case Utility::kNuVtxZ:
    returnString = "Neutrino Vertex Z [cm]";
    break;
  case Utility::kWCNumuCCFlag:
    returnString = "WC Numu CC Flag";
    break;
  case Utility::kisFC:
    returnString = "WC Is FC";
    break;
  case Utility::kLArPID_el:
    returnString = "LAr PID Electron Score";
    break;
  case Utility::kLArPID_ph:
    returnString = "LAr PID Photon Score";
    break;
  case Utility::kLArPID_mu:
    returnString = "LAr PID Muon Score";
    break;
  case Utility::kLArPID_pi:
    returnString = "WC LArPID Muon-Pion Score";
    break;
  case Utility::kLArPID_pr:
    returnString = "WC LArPID Proton-Pion Score";
    break;
  case Utility::kLArPID_prim:
    returnString = "LAr PID Primary Process Score";
    break;
  case Utility::kLArPID_chgd:
    returnString = "LAr PID Charged Process Score";
    break;
  case Utility::kLArPID_ntrl:
    returnString = "LAr PID Neutral Process Score";
    break;
  case Utility::kLArPID_pdg:
    returnString = "LAr PID PDG";
    break;
  case Utility::kNuGraphMipFrac:
    returnString = "Nu Graph Mip Fraction";
    break;
  case Utility::kNuGraphHipFrac:
    returnString = "Nu Graph Hip Fraction";
    break;
  case Utility::kNuGraphShrFrac:
    returnString = "Nu Graph Shower Fraction";
    break;
  case Utility::kNuGraphMclFrac:
    returnString = "Nu Graph Michel Fraction";
    break;
  case Utility::kNuGraphDfsFrac:
    returnString = "Nu Graph Diffuse Fraction";
    break;
  case Utility::kNBlips25cm:
    returnString = "Number of Blips 25cm";
    break;
  case Utility::kNBlips50cm:
    returnString = "Number of Blips 50cm";
    break;
  case Utility::kNBlips10cm:
    returnString = "Number of Blips 100cm";
    break;
  case Utility::kSingleBin:
    returnString = "Single Bin";
    break;
  case Utility::kLanternVtxScore:
    returnString = "Lantern Vertex Score";
    break;
  case Utility::kLanternVtxContainment:
    returnString = "Lantern Vertex Containment";
    break;
  case Utility::kLanternVtxFracHitsOnCosmic:
    returnString = "Lantern Vertex Fraction Cosmic Hits";
    break;
  
  default:
    std::cout << "[ERROR: StackedHistTool] Could not find string conversion for enum " << plotvariable << std::endl;
    returnString = "Unknown";
    break;
  }

  return returnString;
}

void StackedHistTool::SaveHistStack() {
  // TO DO
}