#include "PPSCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TVirtualFFT.h"
#include <iostream>

void ref_timing_gastof_combination() {

  TString files[2] = { "run359_gastof_2p6kV.root", "run359_timingref.root" };
  TFile* f[2];
  TTree* t[2];

  Int_t fNumMeasurements;
  Int_t fETTT;
  Double_t fLeadingEdge[1000], fTrailingEdge[1000];

  TH1D *h_ettt[2], *h_mult[2];
  TH1D *h_abstime_lead[2], *h_abstime_trail[2];
  TH1D *h_abstime_lead_zoom[2], *h_abstime_trail_zoom[2];
  TH1D *h_abstime_lead_merged[2], *h_abstime_trail_merged[2];

  for (UInt_t fn=0; fn<2; fn++) {
    f[fn] = new TFile(files[fn]);
    t[fn] = (TTree*)(f[fn]->Get("tdc"));
    t[fn]->SetBranchAddress("num_measurements", &fNumMeasurements);
    t[fn]->SetBranchAddress("ettt", &fETTT);
    t[fn]->SetBranchAddress("leading_edge", fLeadingEdge);
    t[fn]->SetBranchAddress("trailing_edge", fTrailingEdge);

    h_ettt[fn] = new TH1D("abs_time_ettt", "", 5000, 0., 1000.);
    h_mult[fn] = new TH1D("hits_multiplicity", "", 10, 0., 10.);
    h_abstime_lead[fn] = new TH1D("abs_time_lead", "", 1000, 0., 2000.);
    h_abstime_trail[fn] = new TH1D("abs_time_trail", "", 1000, 0., 2000.);
    h_abstime_lead_zoom[fn] = new TH1D("abs_time_lead_zoom", "", 1000, 500., 600.);
    h_abstime_trail_zoom[fn] = new TH1D("abs_time_trail_zoom", "", 1000, 500., 600.);
    h_abstime_lead_merged[fn] = new TH1D("abs_time_lead_merged", "", 1000, 0., 10.);
    h_abstime_trail_merged[fn] = new TH1D("abs_time_trail_merged", "", 1000, 0., 10.);

    Int_t max = -1, last_ettt = 0;
    Double_t ettt;
    Int_t num_overfl = 0;
    Double_t first_trig = 0., last_trig = 0.;
    const Double_t window_width = 10000.; // in ns
    //const Double_t mult = 1./t[fn]->GetEntries();
    const Double_t mult = 1.;
    for (Int_t i=0; i<t[fn]->GetEntries(); i++) {
      t[fn]->GetEntry(i);

      // first compute "proper" ETTT
      if (fETTT<last_ettt) num_overfl++;
      double shift = num_overfl*(1UL<<32)*25.;
      ettt = shift+fETTT*25.; // in ns
      if (i==0) first_trig = ettt;

      h_ettt[fn]->Fill(ettt/1.e9, mult);
      if (ettt-last_trig>2.e9) { // 2 seconds separation between triggers
        first_trig = ettt;
      }
      h_mult[fn]->Fill(fNumMeasurements, mult);
      for (Int_t j=0; j<fNumMeasurements; j++) {
        Double_t lead_abs = ettt-window_width+fLeadingEdge[j], trail_abs = ettt-window_width+fTrailingEdge[j]; 
        h_abstime_lead[fn]->Fill(lead_abs/1.e9, mult);
        h_abstime_trail[fn]->Fill(trail_abs/1.e9, mult);
        h_abstime_lead_zoom[fn]->Fill(lead_abs/1.e9, mult);
        h_abstime_trail_zoom[fn]->Fill(trail_abs/1.e9, mult);
        h_abstime_lead_merged[fn]->Fill((lead_abs-first_trig)/1.e9, mult);
        h_abstime_trail_merged[fn]->Fill((trail_abs-first_trig)/1.e9, mult);
      }
      max = TMath::Max(max, fETTT);
      last_ettt = fETTT;
      last_trig = ettt;
    }
  }

  PPSCanvas* c_ettt = new PPSCanvas("ettt_run359_absolute_timing", "ETTT");
  h_ettt[0]->Draw();
  h_ettt[0]->GetXaxis()->SetTitle("Extended trigger time tag (s)");
  h_ettt[0]->GetYaxis()->SetTitle("Triggers");
  c_ettt->AddLegendEntry(h_ettt[0], "GasToF");
  h_ettt[1]->Draw("same");
  h_ettt[1]->SetLineColor(kRed);
  h_ettt[1]->SetLineStyle(2);
  c_ettt->AddLegendEntry(h_ettt[1], "Ref. timing");
  c_ettt->Save("pdf", "ref_timing/");

  PPSCanvas* c_mult = new PPSCanvas("hits_multiplicities_run359", "Hits multiplicities / trigger");
  h_mult[0]->Draw();
  h_mult[0]->GetXaxis()->SetTitle("Hits multiplicities / trigger");
  h_mult[0]->GetYaxis()->SetTitle("Triggers fraction");
  h_mult[0]->GetYaxis()->SetRangeUser(0.1, h_mult[0]->GetMaximum()*1.3);
  c_mult->AddLegendEntry(h_mult[0], "GasToF");
  h_mult[1]->Draw("same");
  h_mult[1]->SetLineColor(kRed);
  h_mult[1]->SetLineStyle(2);
  c_mult->AddLegendEntry(h_mult[1], "Ref. timing");
  c_mult->Save("pdf", "ref_timing/");

  PPSCanvas* c_abstime_allbursts = new PPSCanvas("allhits_allbursts_run359_absolute_timing", "Lead. edges (all bursts combined)");
  h_abstime_lead_merged[0]->Draw();
  h_abstime_lead_merged[0]->GetXaxis()->SetTitle("Absolute time (s)");
  c_abstime_allbursts->AddLegendEntry(h_abstime_lead_merged[0], "GasToF");
  h_abstime_lead_merged[1]->Draw("same");
  h_abstime_lead_merged[1]->SetLineColor(kRed);
  h_abstime_lead_merged[1]->SetLineStyle(2);
  c_abstime_allbursts->AddLegendEntry(h_abstime_lead_merged[1], "Ref. timing");
  c_abstime_allbursts->Save("pdf", "ref_timing/");

  PPSCanvas* c_abstime = new PPSCanvas("allhits_run359_absolute_timing", "Leading edges (all indiv. bursts)");
  h_abstime_lead[0]->Draw();
  h_abstime_lead[0]->GetYaxis()->SetRangeUser(0.1, h_abstime_lead[0]->GetMaximum()*1.3);
  h_abstime_lead[0]->GetXaxis()->SetTitle("Absolute time (s)");
  c_abstime->AddLegendEntry(h_abstime_lead[0], "GasToF");
  h_abstime_lead[1]->Draw("same");
  h_abstime_lead[1]->SetLineColor(kRed);
  h_abstime_lead[1]->SetLineStyle(2);
  c_abstime->AddLegendEntry(h_abstime_lead[1], "Ref. timing");
  c_abstime->Save("pdf", "ref_timing/");

  PPSCanvas* c_abstime_zoom = new PPSCanvas("allhits_run359_absolute_timing_zoom", "Leading edges (all indiv. bursts)");
  h_abstime_lead_zoom[0]->Draw();
  h_abstime_lead_zoom[0]->GetYaxis()->SetRangeUser(0.1, h_abstime_lead_zoom[0]->GetMaximum()*1.3);
  h_abstime_lead_zoom[0]->GetXaxis()->SetTitle("Absolute time (s)");
  c_abstime_zoom->AddLegendEntry(h_abstime_lead_zoom[0], "GasToF");
  h_abstime_lead_zoom[1]->Draw("same");
  h_abstime_lead_zoom[1]->SetLineColor(kRed);
  h_abstime_lead_zoom[1]->SetLineStyle(2);
  c_abstime_zoom->AddLegendEntry(h_abstime_lead_zoom[1], "Ref. timing");
  c_abstime_zoom->Save("pdf", "ref_timing/");

  //cout << "max ettt = " << max << endl;

}
