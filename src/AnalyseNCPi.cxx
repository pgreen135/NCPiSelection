// Driver script to launch analysis code

#include "../include/Utility.h"
#include "../include/EventContainer.h"
#include "../include/Selection.h"
#include "../include/SelectionDriver.h"
#include "../include/StackedHistTool.h"
#include "../include/CreateTrainingTree.h"

#include <iostream>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TROOT.h>
#include <TSystem.h>


int main() {

	//ROOT::EnableImplicitMT();
	
	SelectionDriver _selectionDriver;
	
    _selectionDriver.runBDTSelectionFull();
    
	//_selectionDriver.optimiseBDTSelection();

	//_selectionDriver.createElectronPhotonBDTTrainingTreeFHC();
	//_selectionDriver.createElectronPhotonBDTTrainingTreeRHC();

	//_selectionDriver.createPionProtonBDTTrainingTreeFHC();
	//_selectionDriver.createPionProtonBDTTrainingTreeRHC();
	
	return 0;
}