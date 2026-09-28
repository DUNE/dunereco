// n-nbar score plots for the signal detector-variation samples.
//   root -l -b -q 'plots/plot_nnbar_scores.C(<cut>, "<note>")'   (run from the directory holding data/plots/, see make_score_table.py)
// Input: data/plots/scores.csv (detvar mode score), mode_labels.txt, detvar_labels.txt
// Output: plots/*.pdf. Style follows forKevin/plot_nnbar_cut_efficiency.py.
#include <fstream>
#include <map>
#include <sstream>
#include <vector>

void plot_nnbar_scores(double cut = 0.9668, const char* cutnote = "provisional cut (technote value, to be updated)") {
  gROOT->SetBatch(true); gStyle->SetOptStat(0); gStyle->SetEndErrorSize(4);
  gStyle->SetTitleFont(42, "XYZ"); gStyle->SetLabelFont(42, "XYZ"); gStyle->SetLegendFont(42);
  const int NV = 4; const int NM = 44;
  std::vector<std::string> vname(NV), vlabel(NV); std::map<int, std::string> mlabel;
  { std::ifstream in("data/plots/detvar_labels.txt"); int i; std::string n, l;
    while (in >> i && std::getline(in, n, '\t') && std::getline(in, n, '\t') && std::getline(in, l)) { vname[i] = n; vlabel[i] = l; } }
  { std::ifstream in("data/plots/mode_labels.txt"); std::string line;
    while (std::getline(in, line)) { auto t = line.find('\t'); mlabel[std::stoi(line.substr(0, t))] = line.substr(t + 1); } }
  TTree tree("t", "scores"); tree.ReadFile("data/plots/scores.csv", "detvar/I:mode/I:score/F", ',');
  int detvar, mode; float score; tree.SetBranchAddress("detvar", &detvar); tree.SetBranchAddress("mode", &mode); tree.SetBranchAddress("score", &score);
  const int colors[NV] = {kBlue + 1, kRed + 1, kGreen + 2, kOrange + 7}; const int markers[NV] = {20, 21, 22, 23}; const double offsets[NV] = {-0.24, -0.08, 0.08, 0.24};
  std::vector<TH1D*> hv(NV); std::vector<TH1D*> hm(NM + 1, nullptr);
  double total[NV][NM + 1] = {{0}}, passed[NV][NM + 1] = {{0}}; double mtotal[NM + 1] = {0};
  for (int i = 0; i < NV; ++i) { hv[i] = new TH1D(Form("hv%d", i), ";n#bar{n} score;fraction of events / 0.01", 100, 0, 1); hv[i]->SetLineColor(colors[i]); hv[i]->SetLineWidth(2); }
  for (int m = 1; m <= NM; ++m) { hm[m] = new TH1D(Form("hm%d", m), ";n#bar{n} score;fraction of events / 0.02", 50, 0, 1); }
  for (Long64_t k = 0; k < tree.GetEntries(); ++k) {
    tree.GetEntry(k); hv[detvar]->Fill(score); hm[mode]->Fill(score); mtotal[mode] += 1;
    total[detvar][mode] += 1; if (score > cut) passed[detvar][mode] += 1;
  }
  TLatex note; note.SetNDC(true); note.SetTextFont(42); note.SetTextSize(0.027);

  // ---- 1. score distributions per detector variation
  { TCanvas c("scores_detvar", "", 1100, 800); c.SetLogy(); c.SetLeftMargin(0.10); c.SetRightMargin(0.03); c.SetTopMargin(0.07); c.SetBottomMargin(0.11);
    TLegend leg(0.13, 0.70, 0.52, 0.90); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.032);
    double ymax = 0;
    for (int i = 0; i < NV; ++i) { hv[i]->Scale(1.0 / hv[i]->Integral()); ymax = std::max(ymax, hv[i]->GetMaximum()); }
    for (int i = 0; i < NV; ++i) { hv[i]->SetMaximum(ymax * 3); hv[i]->SetMinimum(1e-5); hv[i]->GetXaxis()->SetTitleSize(0.045); hv[i]->GetYaxis()->SetTitleSize(0.045); hv[i]->GetYaxis()->SetTitleOffset(1.0);
      hv[i]->Draw(i ? "HIST SAME" : "HIST"); leg.AddEntry(hv[i], Form("%s  (%.0f events)", vlabel[i].c_str(), hv[i]->GetEntries()), "l"); }
    leg.Draw(); TLine l(cut, 1e-5, cut, ymax * 3); l.SetLineStyle(2); l.SetLineColor(kGray + 2); l.Draw();
    note.DrawLatex(0.13, 0.945, "n#bar{n} hA-BR signal, WireMod detector variations, events after precut, nnbar_best.ckpt");
    note.DrawLatex(0.13, 0.655, Form("dashed line: score = %.4g, %s", cut, cutnote)); c.RedrawAxis(); c.SaveAs("plots/nnbar_scores_by_detvar.pdf"); }

  // ---- 2. efficiency at the cut per decay mode, per detector variation
  { double maximum = 0; std::vector<TGraphAsymmErrors*> gs;
    for (int i = 0; i < NV; ++i) { auto g = new TGraphAsymmErrors(NM); g->SetName(Form("eff_%s", vname[i].c_str())); g->SetTitle(vlabel[i].c_str());
      g->SetMarkerColor(colors[i]); g->SetLineColor(colors[i]); g->SetMarkerStyle(markers[i]); g->SetMarkerSize(1.05); g->SetLineWidth(2);
      for (int m = 1; m <= NM; ++m) { double n = total[i][m], p = passed[i][m]; double e = n > 0 ? p / n : 0;
        double lo = n > 0 ? TEfficiency::ClopperPearson(n, p, 0.682689492, false) : 0, hi = n > 0 ? TEfficiency::ClopperPearson(n, p, 0.682689492, true) : 0;
        g->SetPoint(m - 1, m + offsets[i], e); g->SetPointError(m - 1, 0, 0, e - lo, hi - e); maximum = std::max(maximum, hi); }
      gs.push_back(g); }
    TCanvas c("eff_mode", "", 1700, 900); c.SetLeftMargin(0.085); c.SetRightMargin(0.025); c.SetBottomMargin(0.105); c.SetTopMargin(0.075); c.SetGridy(true);
    TH1D frame("frame", Form(";GENIE NNBarOsc decay mode;CVN efficiency  N(score > %.4g) / N(after precut)", cut), NM, 0.5, NM + 0.5);
    frame.SetMinimum(0); frame.SetMaximum(std::min(1.0, std::max(0.25, 1.22 * maximum)));
    for (int m = 1; m <= NM; ++m) frame.GetXaxis()->SetBinLabel(m, Form("%d", m));
    frame.GetXaxis()->SetLabelSize(0.027); frame.GetXaxis()->SetTitleSize(0.045); frame.GetYaxis()->SetLabelSize(0.038); frame.GetYaxis()->SetTitleSize(0.043); frame.GetYaxis()->SetTitleOffset(0.90);
    frame.Draw("AXIS"); TLegend leg(0.105, 0.765, 0.40, 0.915); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.032);
    for (auto g : gs) { g->Draw("P SAME"); leg.AddEntry(g, g->GetTitle(), "pe"); } leg.Draw();
    note.DrawLatex(0.50, 0.925, Form("68.27%% Clopper-Pearson intervals;  %s", cutnote)); c.RedrawAxis(); c.SaveAs("plots/nnbar_cvn_efficiency_by_mode.pdf");
    std::ofstream csv("plots/nnbar_cvn_efficiency_by_mode.csv"); csv << "detvar,decay_mode,decay_label,n_total,n_pass,efficiency\n";
    for (int i = 0; i < NV; ++i) for (int m = 1; m <= NM; ++m) csv << vname[i] << "," << m << ",\"" << mlabel[m] << "\"," << total[i][m] << "," << passed[i][m] << "," << (total[i][m] > 0 ? passed[i][m] / total[i][m] : 0) << "\n"; }

  // ---- 3. score distributions for the most populous decay modes (all four samples together)
  { std::vector<int> order; for (int m = 1; m <= NM; ++m) order.push_back(m);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return mtotal[a] > mtotal[b]; });
    TCanvas c("scores_mode", "", 1100, 800); c.SetLogy(); c.SetLeftMargin(0.10); c.SetRightMargin(0.03); c.SetTopMargin(0.07); c.SetBottomMargin(0.11);
    TLegend leg(0.12, 0.60, 0.60, 0.90); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.028);
    const int cols[8] = {kBlue + 1, kRed + 1, kGreen + 2, kOrange + 7, kMagenta + 1, kCyan + 2, kGray + 2, kYellow + 3}; double ymax = 0;
    for (int j = 0; j < 8; ++j) { auto h = hm[order[j]]; h->Scale(1.0 / h->Integral()); ymax = std::max(ymax, h->GetMaximum()); }
    for (int j = 0; j < 8; ++j) { auto h = hm[order[j]]; h->SetLineColor(cols[j]); h->SetLineWidth(2); h->SetMaximum(ymax * 3); h->SetMinimum(1e-4);
      h->GetXaxis()->SetTitleSize(0.045); h->GetYaxis()->SetTitleSize(0.045); h->GetYaxis()->SetTitleOffset(1.0); h->Draw(j ? "HIST SAME" : "HIST");
      std::string lab = mlabel[order[j]]; auto p = lab.find("-->"); if (p != std::string::npos) lab.replace(p, 3, "#rightarrow");
      leg.AddEntry(h, Form("%d: %s  (%.0f)", order[j], lab.c_str(), mtotal[order[j]]), "l"); }
    leg.Draw(); TLine l(cut, 1e-4, cut, ymax * 3); l.SetLineStyle(2); l.SetLineColor(kGray + 2); l.Draw();
    note.DrawLatex(0.12, 0.945, "n#bar{n} hA-BR signal, four detector variations combined, eight most populous decay modes");
    c.RedrawAxis(); c.SaveAs("plots/nnbar_scores_by_mode.pdf"); }

  // ---- 4. decay-mode key
  { TCanvas c("decay_mode_key", "", 1700, 1000); c.SetMargin(0.03, 0.03, 0.03, 0.06); TLatex t; t.SetNDC(true); t.SetTextFont(42);
    t.SetTextSize(0.030); t.DrawLatex(0.04, 0.955, "GENIE NNBarOsc decay-mode definitions (events after precut, four variations combined)");
    t.SetTextSize(0.021);
    for (int m = 1; m <= NM; ++m) { int col = m <= 22 ? 0 : 1, row = m <= 22 ? m - 1 : m - 23; std::string lab = mlabel[m]; auto p = lab.find("-->"); if (p != std::string::npos) lab.replace(p, 3, "#rightarrow");
      t.DrawLatex(0.04 + 0.49 * col, 0.91 - 0.0395 * row, Form("%2d: %s   (%.0f)", m, lab.c_str(), mtotal[m])); }
    c.SaveAs("plots/nnbar_decay_mode_key.pdf"); }
  printf("done: cut %.4g\n", cut);
}
