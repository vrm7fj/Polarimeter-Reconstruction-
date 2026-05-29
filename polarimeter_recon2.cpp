void polarimeter_recon_singleloop() {

  TString filename = "/volatile/halla/sbs/vidura/GEP_REPLAYS/GEP3/LH2/KIN3/rootfiles/gep5_fullreplay_5969.root";
  TChain *c = new TChain("T");
  c->Add(filename);

  // Cuts
  TCut globalcut_thetafpp =
    "sqrt(pow((heep.dxECAL-0.01+0.025*earm.ecal.x)/0.0125,2)"
    "+pow((heep.dyECAL-(0.0019+0.00365*earm.ecal.x+0.0171*pow(earm.ecal.x,2)"
    "+0.01448*pow(earm.ecal.x,3)))/0.0157,2))<=3.5"
    "&&abs(heep.dt_ADC-0.5)<9"
    "&&abs(heep.dpp-0.0164)<0.04"
    "&&abs(sbs.tr.vz+0.1)<0.175"
    "&&(sbs.gemFT.track.nhits>4||sbs.gemFT.track.ngoodhits>2)"
    "&&sqrt(pow((sbs.gemFPP.track.y+sbs.gemFPP.track.yp*6.7-sbs.hcal.y+0.0132)/0.0425,2)"
    "+pow((sbs.gemFPP.track.x+sbs.gemFPP.track.xp*6.7-sbs.hcal.x-0.192)/0.0497,2))<=3.5"
    "&&sbs.gemFPP.track.ngoodhits>1"
    "&&earm.ecal.nblk>2"
    "&&sbs.hcal.nblk>1";

  TCut globalcut_zclose = globalcut_thetafpp + "&& sbs.gemFPP.track.sclose < 0.005";

  TCut cut_sAng = "sbs.gemFPP.track.theta*TMath::RadToDeg() <= 1.1";
  TCut cut_lAng = "sbs.gemFPP.track.theta*TMath::RadToDeg() > 1.1";

  // Histograms
  TH1D *htheta_fpp = new TH1D("htheta_fpp", "theta_fpp;theta[deg]; count", 100, 0, 10);
  TH1D *hdoca      = new TH1D("hdoca", "DOCA;DOCA[m];count", 100, 0, 0.05);

  TH1D *hzclose_all = new TH1D("hzclose", "zclose; zclose[m]; count", 50, 0, 3.5);
  TH1D *hzclose_sAng = new TH1D("hzclose_sAng", "zclose small angle; zclose[m]; count", 50, 0, 3.5);
  TH1D *hzclose_lAng = new TH1D("hzclose_lAng", "zclose large angle; zclose[m]; count", 50, 0, 3.5);

  // Branches
  const int MAX = 10000;
  vector<double> theta(MAX), sclose(MAX), zclose(MAX);
  int besttrack, ntracks;

  c->SetBranchStatus("*",0);
  c->SetBranchStatus("sbs.gemFPP.track.theta",1);
  c->SetBranchStatus("sbs.gemFPP.track.sclose",1);
  c->SetBranchStatus("sbs.gemFPP.track.zclose",1);
  c->SetBranchStatus("sbs.gemFPP.track.besttrack",1);
  c->SetBranchStatus("sbs.gemFPP.track.ntrack",1);

  // Also enable all branches used in cuts
  c->SetBranchStatus("heep.*",1);
  c->SetBranchStatus("earm.*",1);
  c->SetBranchStatus("sbs.*",1);

  c->SetBranchAddress("sbs.gemFPP.track.theta", &theta[0]);
  c->SetBranchAddress("sbs.gemFPP.track.sclose", &sclose[0]);
  c->SetBranchAddress("sbs.gemFPP.track.zclose", &zclose[0]);
  c->SetBranchAddress("sbs.gemFPP.track.besttrack", &besttrack);
  c->SetBranchAddress("sbs.gemFPP.track.ntrack", &ntracks);

  // Tree formulas
  TTreeFormula f_cut_theta("f_cut_theta", globalcut_thetafpp, c);
  TTreeFormula f_cut_zclose("f_cut_zclose", globalcut_zclose, c);
  TTreeFormula f_cut_sAng("f_cut_sAng", cut_sAng, c);
  TTreeFormula f_cut_lAng("f_cut_lAng", cut_lAng, c);

  Long64_t N = c->GetEntries();
  for (Long64_t i=0; i<N; i++) {
    c->GetEntry(i);

    // Validate track index
    if (besttrack < 0 || besttrack >= ntracks) continue;

    // Fill theta_fpp and doca
    if (f_cut_theta.EvalInstance()) {
      htheta_fpp->Fill(theta[besttrack] * TMath::RadToDeg());
      hdoca->Fill(sclose[besttrack]);
    }

    // zclose ALL
    if (f_cut_zclose.EvalInstance()) {
      hzclose_all->Fill(zclose[besttrack]);

      // small angle
      if (f_cut_sAng.EvalInstance())
        hzclose_sAng->Fill(zclose[besttrack]);

      // large angle
      if (f_cut_lAng.EvalInstance())
        hzclose_lAng->Fill(zclose[besttrack]);
    }
  }

    // -------------------------
    // Drawing Section
    // -------------------------

     TCanvas *canvas = new TCanvas("polarimeter_recon", "polarimeter_recon", 1200, 900);
     canvas->Divide(2,2);

     // -------------------------
     // PAD 1: theta_fpp
     // -------------------------
     canvas->cd(1);
     gPad->SetLogy();

     htheta_fpp->SetLineColor(kBlue+1);
     htheta_fpp->SetLineWidth(2);
     htheta_fpp->Draw();

     double ymax = htheta_fpp->GetMaximum();

     TLine *lfpp_theta1 = new TLine(1.1, 0, 1.1, ymax);
     TLine *lfpp_theta2 = new TLine(6.0, 0, 6.0, ymax);

     lfpp_theta1->SetLineColor(kRed);
     lfpp_theta1->SetLineStyle(2);
     lfpp_theta2->SetLineColor(kRed);
     lfpp_theta2->SetLineStyle(2);

     lfpp_theta1->Draw("same");
     lfpp_theta2->Draw("same");


     // -------------------------
     // PAD 2: DOCA
     // -------------------------
     canvas->cd(2);

     hdoca->SetLineColor(kBlack);
     hdoca->SetLineWidth(2);
     hdoca->Draw();


     // -------------------------
     // PAD 3: zclose stacked
     // -------------------------
     canvas->cd(3);

     THStack *stack_zclose = new THStack("stack_zclose", "Q^{2} = 11 GeV^{2}; zclose [m]; count");

     // Colors and fill styles
     hzclose_all->SetLineColor(kBlack);
     hzclose_all->SetFillColor(kGray+1);
     hzclose_all->SetFillStyle(3004);

     hzclose_sAng->SetLineColor(kBlue+1);
     hzclose_sAng->SetFillColor(kBlue-9);
     hzclose_sAng->SetFillStyle(3005);

     hzclose_lAng->SetLineColor(kRed+1);
     hzclose_lAng->SetFillColor(kRed-7);
     hzclose_lAng->SetFillStyle(3003);

     // Add in correct order (background first)
     stack_zclose->Add(hzclose_all);
     stack_zclose->Add(hzclose_sAng);
     stack_zclose->Add(hzclose_lAng);

     stack_zclose->Draw("hist");

     // Legend
     TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.88);
     leg->AddEntry(hzclose_all, "All", "f");
     leg->AddEntry(hzclose_sAng, "#theta_{FPP} < 1.1^{#circ}", "f");
     leg->AddEntry(hzclose_lAng, "#theta_{FPP} > 1.1^{#circ}", "f");
     leg->Draw();


     // -------------------------
     // PAD 4: optional ratio or overlay
     // -------------------------
     canvas->cd(4);

     hzclose_lAng->SetLineWidth(2);
     hzclose_sAng->SetLineWidth(2);

     hzclose_lAng->Draw("hist");
     hzclose_sAng->Draw("hist same");

     TLegend *leg2 = new TLegend(0.65, 0.75, 0.88, 0.88);
     leg2->AddEntry(hzclose_lAng, "Large angle", "l");
     leg2->AddEntry(hzclose_sAng, "Small angle", "l");
     leg2->Draw();

     canvas->Update();

}
