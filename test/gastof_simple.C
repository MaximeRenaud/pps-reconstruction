#include "GastofCommons.h"
#include "GastofGrid.h"
#include "GastofCanvas.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>

void gastof_simple()
{
  const unsigned int run_id = 725;
  TFile* f = new TFile(Form("../../gastof_inner_run%d.root", run_id)); //FIXME
  TTree* t = (TTree*)f->Get("tdc");

  const unsigned int channels_to_probe = 32; //64;
  const double tot_min = 70., tot_max = -1.;
  const double lead_min = 6.8, lead_max = 6.95; // in \mus

  int num_hits;
  unsigned int ettt;
  double leading_edge[5000], tot[5000];
  int channel_id[5000];
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("num_measurements", 1); t->SetBranchAddress("num_measurements", &num_hits);
  t->SetBranchStatus("leading_edge", 1); t->SetBranchAddress("leading_edge", leading_edge);
  t->SetBranchStatus("channel_id", 1); t->SetBranchAddress("channel_id", channel_id);
  t->SetBranchStatus("tot", 1); t->SetBranchAddress("tot", tot);

  const unsigned int num_triggers = t->GetEntries();

  TH1D* h_tot_all_channels = new TH1D("h_tot_all_channels", "", 250, 0., 1000.);
  TH1D* h_lead_all_channels = new TH1D("h_lead_all_channels", "", 500, 0., 10.);

  double last_ettt = -1.;
  unsigned int num_overfl = 0;
  unsigned int mult[channels_to_probe];
  for (unsigned int i=0; i<num_triggers; i++) {
    t->GetEntry(i);
    for (unsigned int j=0; j<channels_to_probe; j++) {
//       lead[j].clear();    //lead[] not defined in gastof_simple
      mult[j] = 0;
    }
    if (i%10000==0) cerr << "Processing event " << i << " / " << num_triggers << endl;

    for (int j=0; j<num_hits; j++) {
      h_lead_all_channels->Fill(leading_edge[j]/1.e3, 1./num_triggers);
      if (leading_edge[j]<lead_min*1.e3 or leading_edge[j]>lead_max*1.e3) continue;
      h_tot_all_channels->Fill(tot[j]);
    }
  }

  GastofCanvas c_tot_all_chan("tot_all-channels", Form("Run %d, all channels", run_id));
  h_tot_all_channels->Draw();
  h_tot_all_channels->GetXaxis()->SetTitle("Time over threshold (ns)");
  h_tot_all_channels->GetYaxis()->SetTitle("Hits");

  {
    TPaveText* cuts = new TPaveText(0.4, 0.8, 0.5, 0.85, "ndc");
    cuts->SetTextAlign(13);
    cuts->SetTextFont(43);
    cuts->SetTextSize(18);
    cuts->SetFillColor(kWhite);
    cuts->SetLineColor(kWhite);
    //cuts->AddText(Form("ToT > %.1f ns", tot_min));
    cuts->AddText(Form("Lead. edge #in [%.1f-%.2f]  #mus", lead_min, lead_max));
    cuts->Draw("same");
  }

  c_tot_all_chan.Pad()->SetLogy();
  {
    TLine* line = new TLine(tot_min, 0., tot_min, h_tot_all_channels->GetMaximum()*0.8); //constructor for line: TLine( x1, y1, x2, y2)
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");
  }
//   if (tot_max>0.) {  //<- condition never met in gastof_simple
//     TLine* line = new TLine(tot_max, 0., tot_max, h_tot_all_channels->GetMaximum()*0.8);
//     line->SetLineColor(kBlack);
//     line->SetLineStyle(2);
//     line->SetLineWidth(2);
//     line->Draw("same");
  }
  c_tot_all_chan.Prettify(h_tot_all_channels);
  h_tot_all_channels->SetLineColor(kGray+1);
  c_tot_all_chan.Save("png", "time_diff_betw_chan");
  c_tot_all_chan.Save("pdf", "time_diff_betw_chan");

  {
    GastofCanvas c("leading_all-channels", Form("Run %d, all channels", run_id));
    h_lead_all_channels->Draw();

    {
      TLine* line1 = new TLine(lead_min, 0., lead_min, h_lead_all_channels->GetMaximum()*1.05);
      line1->SetLineColor(kBlack);
      line1->SetLineStyle(2);
      line1->SetLineWidth(2);
      line1->Draw("same");
      TLine* line2 = new TLine(lead_max, 0., lead_max, h_lead_all_channels->GetMaximum()*1.05);
      line2->SetLineColor(kBlack);
      line2->SetLineStyle(2);
      line2->SetLineWidth(2);
      line2->Draw("same");
    }

    c.Prettify(h_lead_all_channels);
    h_lead_all_channels->GetXaxis()->SetTitle("Leading edge ( #mus)");
    h_lead_all_channels->GetYaxis()->SetTitle("Hits / Trigger");
    h_lead_all_channels->SetLineColor(kGray+1);
    h_lead_all_channels->GetXaxis()->SetLabelSize(18);
    c.Save("png", "time_diff_betw_chan");
    c.Save("pdf", "time_diff_betw_chan");
  }

}
