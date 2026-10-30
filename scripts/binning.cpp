// simple binning optimisation

void binning() {

    // read input file "binningTree.root"
    TFile inputFile("binningTree.root", "READ");
    TTree* binningTree = (TTree*)inputFile.Get("Events");

    // load branches
    double sel_pion_momentum;
    double sel_pion_theta;
    double sel_pion_phi;
    double true_pion_momentum;
    double true_pion_theta;
    double true_pion_phi;
    double weight;
    binningTree->SetBranchAddress("sel_pion_momentum", &sel_pion_momentum);
    binningTree->SetBranchAddress("sel_pion_theta", &sel_pion_theta);
    binningTree->SetBranchAddress("sel_pion_phi", &sel_pion_phi);
    binningTree->SetBranchAddress("true_pion_momentum", &true_pion_momentum);
    binningTree->SetBranchAddress("true_pion_theta", &true_pion_theta);
    binningTree->SetBranchAddress("true_pion_phi", &true_pion_phi);
    binningTree->SetBranchAddress("weight", &weight);

    double total_weighted_events = 0;

    std::vector<double> momentum_reco_vals;
    std::vector<double> theta_reco_vals;
    std::vector<double> weight_vals;

    // loop over events
    int current_bin = 0;
    for (Long64_t i = 0; i < binningTree->GetEntries(); i++) {
        binningTree->GetEntry(i);

        // count total events
        total_weighted_events += weight;

        // store reco values for binning optimisation
        momentum_reco_vals.push_back(sel_pion_momentum);
        theta_reco_vals.push_back(sel_pion_theta);
        weight_vals.push_back(weight);

    }

    std::cout << "Total weighted events = " << total_weighted_events << std::endl;

    // binning loop
    //std::vector<double> variable_reco_vec = momentum_reco_vals;
	std::vector<double> variable_reco_vec = theta_reco_vals;


	//double lower_edge = 0.1; // pion momentum 100 MeV, from signal definition
	double lower_edge = 1.0; // angles (work other way with angles)
	
	std::vector<double> edges;
	edges.push_back(lower_edge);
    std::vector<double> bin_counts;

	double new_edge = 0.95; // first new edge to check
	//double new_edge = -0.95; // first new edge to check

	for (int bin = 0; bin < 7; bin++) {

		for (int b = 0; b < 200; b++ ) {

			double event_counter = 0;
			
			new_edge = new_edge + -0.005;

			for (int i = 0; i < variable_reco_vec.size(); i++) {
				//if ( variable_reco_vec[i] >= lower_edge && variable_reco_vec[i] < new_edge ) event_counter += weight_vals[i];
				if ( variable_reco_vec[i] <= lower_edge && variable_reco_vec[i] > new_edge ) event_counter += weight_vals[i];
			}

			std::cout << "Bin " << bin << ": [" << lower_edge << ", " << new_edge << "] --> Signal Events: " << event_counter << std::endl;

			/* // momentum tuning
			if ((bin == 0) && event_counter >= 99) break;
			else if (bin == 1 && event_counter >= 95) break;
			else if (bin == 2 && event_counter >= 95) break;
			else if (bin == 3 && event_counter >= 132) break;
			else if (bin == 4 && event_counter >= 180) break;
			else if (event_counter >= 180) break;
			*/
			if (event_counter >= 100) break;

			// angle tuning
			if ((bin == 0) && event_counter >= 99) break;
			else if(event_counter >= 100) break;
			//if (new_edge >= 5) break;

		}

		// store edge
		edges.push_back(new_edge);
        

		// set up for next loop
		lower_edge = new_edge;

		//if (new_edge >= 1) break;
	}

	std::cout << "Bin edges: " << std::endl;
	for (int i = 0; i < edges.size(); i++) {
		std::cout << edges[i] << ", ";
	}

	// compute and print final bin counts
	std::cout << std::endl << "Final bin counts:" << std::endl;
	for (int bin = 0; bin + 1 < (int)edges.size(); bin++) {
		double bin_count = 0;
		double bin_low = edges[bin];
		double bin_high = edges[bin + 1];
		for (int i = 0; i < (int)variable_reco_vec.size(); i++) {
			//if (variable_reco_vec[i] >= bin_low && variable_reco_vec[i] < bin_high) {
			if (variable_reco_vec[i] <= bin_low && variable_reco_vec[i] > bin_high) {
				bin_count += weight_vals[i];
			}
		}
		std::cout << "Bin " << bin << ": [" << bin_low << ", " << bin_high << ") --> Count = " << bin_count << std::endl;
	}
}
