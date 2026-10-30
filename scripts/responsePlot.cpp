// respones matrix plot

void responsePlot() {

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

    Double_t bin_edges_pionMomentum[7] = {0.1, 0.175, 0.215, 0.275, 0.41, 0.635, 3};
	TH1F *pionMomentum1D = new TH1F("", "", 6, bin_edges_pionMomentum);
    TH2F *pionMomentum2D = new TH2F("", "", 6, bin_edges_pionMomentum, 6, bin_edges_pionMomentum);

    Double_t bin_edges_pionAngle[8] = {-1.0, 0.0, 0.42, 0.64, 0.77, 0.865, 0.935, 1.0};
    TH1F *pionAngle1D = new TH1F("pionAngle1D", "Selected pion cos(#theta);Reconstructed cos(#theta);Weighted events", 7, bin_edges_pionAngle);
    TH2F *pionAngle2D = new TH2F("pionAngle2D", "Pion cos(#theta) response;True cos(#theta);Reconstructed cos(#theta)",
                                 7, bin_edges_pionAngle, 7, bin_edges_pionAngle);
	
    // loop over events
    for (Long64_t i = 0; i < binningTree->GetEntries(); i++) {
        binningTree->GetEntry(i);

        // count total events
        total_weighted_events += weight;

        // add to histogram
        pionMomentum1D->Fill(sel_pion_momentum, weight);
        pionAngle1D->Fill(sel_pion_theta, weight);

        // build response matrix
        pionMomentum2D->Fill(true_pion_momentum, sel_pion_momentum, weight);
        pionAngle2D->Fill(true_pion_theta, sel_pion_theta, weight);

    }

    std::cout << "Total weighted events = " << total_weighted_events << std::endl;

    // create an equal-visual-size version of the normalized confusion matrix (uniform bin indices)
    {
        int nX = pionMomentum2D->GetNbinsX();
        int nY = pionMomentum2D->GetNbinsY();
        TH2F *pionMomentum2D_equal = new TH2F("pionMomentum2D_equal", "Confusion matrix (% of true events in each reconstructed bin)",
                                              nX, 0.5, nX + 0.5, nY, 0.5, nY + 0.5);

        for (int ix = 1; ix <= nX; ++ix) {
            double true_bin_sum = 0.0;
            for (int iy = 1; iy <= nY; ++iy) {
                true_bin_sum += pionMomentum2D->GetBinContent(ix, iy);
            }

            for (int iy = 1; iy <= nY; ++iy) {
                double val = 0.0;
                double err = 0.0;
                if (true_bin_sum > 0.0) {
                    val = 100.0 * pionMomentum2D->GetBinContent(ix, iy) / true_bin_sum;
                    err = 100.0 * pionMomentum2D->GetBinError(ix, iy) / true_bin_sum;
                }
                pionMomentum2D_equal->SetBinContent(ix, iy, val);
                pionMomentum2D_equal->SetBinError(ix, iy, err);
            }

            char label[64];
            snprintf(label, sizeof(label), "%.3f-%.3f", bin_edges_pionMomentum[ix-1], bin_edges_pionMomentum[ix]);
            pionMomentum2D_equal->GetXaxis()->SetBinLabel(ix, label);
            pionMomentum2D_equal->GetYaxis()->SetBinLabel(ix, label);
        }
        pionMomentum2D_equal->GetZaxis()->SetTitle("% of true events");
        pionMomentum2D_equal->SetMinimum(0);
        pionMomentum2D_equal->SetMaximum(100);
        pionMomentum2D_equal->SetStats(false);
        
        TCanvas *c2 = new TCanvas("c2", "c2", 1080, 1080);
        c2->cd();
        gPad->SetFixedAspectRatio(1.0);
        c2->SetLeftMargin(0.16);
        c2->SetRightMargin(0.16);
        c2->SetTopMargin(0.08);
        c2->SetBottomMargin(0.13);

        pionMomentum2D_equal->GetXaxis()->SetTitle("True Pion Momentum [GeV/c]");
        pionMomentum2D_equal->GetYaxis()->SetTitle("Reconstructed Pion Momentum [GeV/c]");
        pionMomentum2D_equal->GetXaxis()->SetTitleOffset(1.5);

        // improve label readability
        pionMomentum2D_equal->GetXaxis()->LabelsOption("d"); // 45-degree X labels
        pionMomentum2D_equal->GetXaxis()->SetLabelSize(0.035);
        pionMomentum2D_equal->GetYaxis()->SetLabelSize(0.035);
        pionMomentum2D_equal->GetZaxis()->SetLabelSize(0.035);

        // draw color map and numeric values in each (equal) bin
        pionMomentum2D_equal->Draw("COLZ TEXT0");

        c2->SaveAs("pionMomentum2D_equal.root");
        c2->SaveAs("pionMomentum2D_equal.png");
    }

    // plot the selected pion-angle distribution
    pionAngle1D->SetStats(false);
    TCanvas *cAngle1D = new TCanvas("cAngle1D", "Pion angle distribution", 850, 650);
    cAngle1D->SetLeftMargin(0.14);
    cAngle1D->SetBottomMargin(0.14);
    pionAngle1D->SetLineWidth(2);
    pionAngle1D->Draw("HIST E");
    cAngle1D->SaveAs("pionAngle1D.root");
    cAngle1D->SaveAs("pionAngle1D.png");

    // create an equal-visual-size, true-bin-normalized pion-angle response matrix
    {
        int nX = pionAngle2D->GetNbinsX();
        int nY = pionAngle2D->GetNbinsY();
        TH2F *pionAngle2D_equal = new TH2F("pionAngle2D_equal",
                                           "Pion cos(#theta) response (% of true events in each reconstructed bin)",
                                           nX, 0.5, nX + 0.5, nY, 0.5, nY + 0.5);

        for (int ix = 1; ix <= nX; ++ix) {
            double true_bin_sum = 0.0;
            for (int iy = 1; iy <= nY; ++iy) {
                true_bin_sum += pionAngle2D->GetBinContent(ix, iy);
            }

            for (int iy = 1; iy <= nY; ++iy) {
                double val = 0.0;
                double err = 0.0;
                if (true_bin_sum > 0.0) {
                    val = 100.0 * pionAngle2D->GetBinContent(ix, iy) / true_bin_sum;
                    err = 100.0 * pionAngle2D->GetBinError(ix, iy) / true_bin_sum;
                }
                pionAngle2D_equal->SetBinContent(ix, iy, val);
                pionAngle2D_equal->SetBinError(ix, iy, err);
            }

            char label[64];
            snprintf(label, sizeof(label), "%.3f-%.3f", bin_edges_pionAngle[ix - 1], bin_edges_pionAngle[ix]);
            pionAngle2D_equal->GetXaxis()->SetBinLabel(ix, label);
            pionAngle2D_equal->GetYaxis()->SetBinLabel(ix, label);
        }

        pionAngle2D_equal->GetZaxis()->SetTitle("% of true events");
        pionAngle2D_equal->SetMinimum(0);
        pionAngle2D_equal->SetMaximum(100);
        pionAngle2D_equal->SetStats(false);

        TCanvas *cAngle2D = new TCanvas("cAngle2D", "Pion angle response matrix", 1080, 1080);
        cAngle2D->cd();
        gPad->SetFixedAspectRatio(1.0);
        cAngle2D->SetLeftMargin(0.16);
        cAngle2D->SetRightMargin(0.16);
        cAngle2D->SetTopMargin(0.08);
        cAngle2D->SetBottomMargin(0.13);

        pionAngle2D_equal->GetXaxis()->SetTitle("True Pion cos(#theta)");
        pionAngle2D_equal->GetYaxis()->SetTitle("Reconstructed Pion cos(#theta)");
        pionAngle2D_equal->GetXaxis()->SetTitleOffset(1.5);
        pionAngle2D_equal->GetYaxis()->SetTitleOffset(2.2);
        pionAngle2D_equal->GetZaxis()->SetTitleOffset(1.4);
        pionAngle2D_equal->GetXaxis()->LabelsOption("d");
        pionAngle2D_equal->GetXaxis()->SetLabelSize(0.03);
        pionAngle2D_equal->GetYaxis()->SetLabelSize(0.03);
        pionAngle2D_equal->GetZaxis()->SetLabelSize(0.03);
        gStyle->SetPaintTextFormat("4.1f");
        pionAngle2D_equal->Draw("COLZ TEXT0");

        cAngle2D->SaveAs("pionAngle2D_equal.root");
        cAngle2D->SaveAs("pionAngle2D_equal.png");
    }

}