struct XSPlot {
  XSPlot(TString file, TString sample, TLegend* l, TString _title, int color) : title(_title) {
    TFile* f = TFile::Open(file);
    data = (TH1D*) f->Get(sample + TString("_data"));
    data->SetDirectory(0);
    data->SetLineWidth(2);
    data->SetMarkerStyle(20);
    data->SetMarkerSize(0.75);
    mc = (TH1D*) f->Get(sample + TString("_MC"));
    mc->SetDirectory(0);
    mc->SetLineColor(color);
    mc->SetFillColorAlpha(color, 0.5);
    mc->SetLineWidth(2);
    TH1D* hll = (TH1D*) f->Get("likelihood_hist");
    TH1D* hndof = (TH1D*) f->Get("ndof_hist");

    for (int i=1; i<hll->GetNbinsX()+1; i++) {
      if (TString(hll->GetXaxis()->GetBinLabel(i)).Contains(sample)) {
        float likelihood = hll->GetBinContent(i);
        float ndof = hndof->GetBinContent(i);
        char extitle[100];
        snprintf(extitle, 100, ", #chi^{2}/dof=%1.1f/%1.0f", likelihood, ndof);
        title = title + TString(extitle);
      }
    }
    l->AddEntry(mc, title);
  }

  TH1D* data;
  TH1D* mc;
  TString title;
};


struct GenDef {
  TString filename;
  TString title;
  int color;
};


struct PlotDef {
  TString var;
  std::string suffix;
  float ymax;
  float lx1;
  float ly1;
  float lx2;
  float ly2;
};


TH1D* setup(TH1D* h) {
  h->SetDirectory(0);
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
  h->SetLineStyle(2);
  h->Scale(1e-38);
  return h;

}

void compare_numinue() {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadLeftMargin(0.15);

  TCanvas* c = new TCanvas("c", "", 500, 500);

  std::string outdir = ".";

  std::vector<GenDef> gens = {
      { "numinue.gv3.root", "GENIE v3.00.06", kBlue+1   }
     ,{ "numinue.nuwro.root", "NuWro 19.02.1",    kRed+1    }
     //,{ "../neut/ccinc-svd.neut.root",   "NEUT 5.4.0.1",     kViolet   }
     ,{ "numinue.gibuu.root", "GiBUU 2019.08",    kOrange+1 }
  };

  std::vector<PlotDef> defs {
     { "MicroBooNE_CCInc_NuMI_XSec_1DEe_nue", "ee", 7.4e-39, 0.35, 0.5, 0.88, 0.88 },
     { "MicroBooNE_CCInc_NuMI_XSec_1DCosBeta_nue",  "cosbeta",  25e-39, 0.2, 0.60, 0.73, 0.88 },
  };

  TFile* inputRootFile = TFile::Open("xsec_result_run1_paper_v3.root");

  for (auto const& def : defs) {
    TLegend l(def.lx1, def.ly1, def.lx2, def.ly2);
    l.SetFillStyle(0);

    std::vector<XSPlot*> plots;
    bool first = true;
    for (auto const& gen : gens) {
      plots.push_back(new XSPlot(gen.filename, def.var, &l, gen.title, gen.color));
      plots.back()->mc->GetYaxis()->SetRangeUser(0, def.ymax);
      plots.back()->mc->DrawClone(first ? "hist e3 l" : "hist e3 l same");
      first = false;
    }

    for (auto const& plot : plots) {
      plot->mc->SetFillColor(0);
      plot->mc->Draw("hist l same");
    }

    //TH1D* fMCHist = (TH1D*) inputRootFile->Get(def.suffix == "ee" ? "mc_xsec_true_genie_v3_0_6_energy_new" : "mc_xsec_true_genie_v3_0_6_angle_new");
    TH1D* fMCHist = (TH1D*) inputRootFile->Get(def.suffix == "ee" ? "mc_xsec_true_genie_v3_0_6_energy" : "mc_xsec_true_genie_v3_0_6_angle");
    std::string objSuffix = (def.suffix == "ee" ? "energy" : "angle");

    // Set up the smearing matrix
    assert(inputRootFile && inputRootFile->IsOpen());
    TH2D* hsmear = (TH2D*) inputRootFile->Get(("ac_" + objSuffix).c_str());
    assert(hsmear);
    int nrows = hsmear->GetNbinsX();
    int ncols = hsmear->GetNbinsY();
    TMatrixD* fSmearingMatrix = new TMatrixD(nrows, ncols);
    for (int i=0; i<nrows; i++) {
      for (int j=0; j<ncols; j++) {
        fSmearingMatrix->operator()(i,j) = hsmear->GetBinContent(i+1, j+1);
      }
    }

    int n = fMCHist->GetNbinsX();
    TVectorD v(n);
    for (int i=0; i<n; i++) {
      v(i) = fMCHist->GetBinContent(i+1) * fMCHist->GetBinWidth(i+1);
    }
    TVectorD vs = (*fSmearingMatrix) * v;
    for (int i=0; i<n; i++) {
      fMCHist->SetBinContent(i+1, vs(i) / fMCHist->GetBinWidth(i+1));
    }

    fMCHist->Scale(1e-39);
    fMCHist->SetLineColor(kBlack);
    fMCHist->SetLineStyle(7);
    fMCHist->SetLineWidth(2);
    fMCHist->Draw("hist l same");
    l.AddEntry(fMCHist, "GENIE v3 MC");

    plots[0]->data->Draw("e1 same");
    l.AddEntry(plots[0]->data, "Data");

    l.Draw();
    c->SaveAs(Form("%s/numinue_%s.pdf", outdir.c_str(), def.suffix.c_str()));
    c->SaveAs(Form("%s/numinue_%s.C", outdir.c_str(), def.suffix.c_str()));
  }
}

