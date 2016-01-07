#include "GastofCommons.h"
#include "GastofGrid.h"
#include "GastofCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TFitResult.h"

#include <iostream>

double doublegauss(double* xv, double* par)
{
  double x = xv[0];
  double f1 = par[0]*exp(-pow((x-par[1])/par[2], 2)/2.);
  double f2 = par[3]*exp(-pow((x-par[4])/par[5], 2)/2.);
  return f1+f2+par[6];
}

void gastof_onebyone_channel_study()
{
  const unsigned int run_id = 705;
  TFile* f = new TFile(Form("../datasets/gastof_full_run%d.root", run_id)); //FIXME
  TTree* t = (TTree*)f->Get("tdc");

  const unsigned int channels_to_probe = 64;
  const double tot_min = 70., tot_max = -1.;
  const double lead_min = 6.8, lead_max = 6.95; // in us

  int num_hits;
  unsigned int ettt;
  double leading_edge[5000], tot[5000];
  int channel_id[5000];
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("num_measurements", 1); t->SetBranchAddress("num_measurements", &num_hits);
  t->SetBranchStatus("ettt", 1); t->SetBranchAddress("ettt", &ettt);
  t->SetBranchStatus("leading_edge", 1); t->SetBranchAddress("leading_edge", leading_edge);
  t->SetBranchStatus("channel_id", 1); t->SetBranchAddress("channel_id", channel_id);
  t->SetBranchStatus("tot", 1); t->SetBranchAddress("tot", tot);

  //GastofCoordinatesMap::GetCoordinates(channel_to_probe).Print();

  const unsigned int num_triggers = t->GetEntries();
  vector<double> lead[channels_to_probe];
  vector<double> tots[channels_to_probe];

  const double range = .5;
  
  TH1D* h_diff[channels_to_probe], *h_tot_diff[channels_to_probe];
  TH1D* h_tot[channels_to_probe];
  TH1D* h_mult[channels_to_probe];
  TH1D* h_tot_all_channels = new TH1D("h_tot_all_channels", "", 250, 0., 1000.);
  TH1D* h_lead_all_channels = new TH1D("h_lead_all_channels", "", 500, 0., 10.);
  TH1D* h_lead_all_channels_zoom = new TH1D("h_lead_all_channels_zoom", "", 500, 6.75, 7.0);
  //TH1D* h_diff_all_channels = new TH1D("h_diff_all_channels", "", 1000, -50., 50.);
  TH1D* h_diff_all_channels = new TH1D("h_diff_all_channels", "", (int)range/0.1, -range, range);
  TH1D* h_mult_all_channels = new TH1D("h_mult_all_channels", "", 20, 0, 20);
  for (unsigned int i=0; i<channels_to_probe; i++) {
    h_diff[i] = new TH1D(Form("h_diff_%d", i), "", 1000, -50., 50.);
    h_tot_diff[i] = new TH1D(Form("h_tot_diff_%d", i), "", 1000, -250., 250.);
    h_tot[i] = new TH1D(Form("h_tot_%d", i), "", 100, 0., 1000.);
    h_mult[i] = new TH1D(Form("h_mult_%d", i), "", 20, 0, 20);
  }

  double last_ettt = -1.;
  unsigned int num_overfl = 0;
  unsigned int mult[channels_to_probe];
  for (unsigned int i=0; i<num_triggers; i++) {
    t->GetEntry(i);
    for (unsigned int j=0; j<channels_to_probe; j++) {
      lead[j].clear();
      mult[j] = 0;
    }
    if (i%10000==0) cerr << "Processing event " << i << " / " << num_triggers << endl;

    // compute "physical" ETTT
    if (ettt<last_ettt) num_overfl++;
    last_ettt = ettt;
    double shift = num_overfl*(1UL<<32)*25.;
    double ettt_phys = shift+ettt*25.; // in ns
   
    if (num_hits>1000) continue;
    for (int j=0; j<num_hits; j++) {
      h_lead_all_channels->Fill(leading_edge[j]/1.e3, 1./num_triggers);
      h_lead_all_channels_zoom->Fill(leading_edge[j]/1.e3, 1./num_triggers);
      if (leading_edge[j]<lead_min*1.e3 or leading_edge[j]>lead_max*1.e3) continue;
      h_tot[channel_id[j]]->Fill(tot[j]);
      h_tot_all_channels->Fill(tot[j]);
      mult[channel_id[j]]++;
      if (tot_min>0. and tot[j]<tot_min) continue;
      if (tot_max>0. and tot[j]>tot_max) continue;
      lead[channel_id[j]].push_back(ettt_phys+leading_edge[j]);
      tots[channel_id[j]].push_back(tot[j]);
    }
    for (unsigned int j=0; j<channels_to_probe; j++) {
      h_mult[j]->Fill(mult[j], 1./num_triggers);
      h_mult_all_channels->Fill(mult[j], 1./num_triggers/channels_to_probe);
    }
    for (int j=0; j<num_hits; j++) {
      if (leading_edge[j]<lead_min*1.e3 or leading_edge[j]>lead_max*1.e3) continue;
      if (tot_min>0. and tot[j]<tot_min) continue;
      if (tot_max>0. and tot[j]>tot_max) continue;
      //if (mult[channel_id[j]]!=1) continue; //FIXME
      for (unsigned int k=0; k<channels_to_probe; k++) {
        if (!GastofCoordinatesMap::IsNeighbour(k, channel_id[j])) continue;
        for (unsigned int l=0; l<lead[k].size(); l++) {
        //for (vector<double>::const_iterator l=lead[k].begin(); l!=lead[k].end(); l++) {
          //cout << leading_edge[j]-(*l) << endl;
          //if (fabs(tot[j]-tots[k][l])>15.) continue;
          h_diff[k]->Fill(ettt_phys+leading_edge[j]-lead[k][l]);
          h_tot_diff[k]->Fill(tot[j]-tots[k][l]);
          h_diff_all_channels->Fill(ettt_phys+leading_edge[j]-lead[k][l]);
        }
      }
    }
  }

  {
    GastofCanvas c("time_diff_all-channels", Form("Run %d, all channels", run_id));
    h_diff_all_channels->Draw();
    h_diff_all_channels->GetXaxis()->SetTitle("Neigh. hits time difference (ns)");
    h_diff_all_channels->GetYaxis()->SetTitle("Hits");
    c.SetLegendX(0.16);

    TPaveText* cuts = new TPaveText(0.15, 0.68, 0.2, 0.78, "ndc");
    cuts->SetTextAlign(13);
    cuts->SetTextFont(43);
    cuts->SetTextSize(18);
    cuts->SetFillColor(kWhite);
    cuts->SetLineColor(kWhite);
    if (tot_max<0.) cuts->AddText(Form("ToT > %.1f ns", tot_min));
    else            cuts->AddText(Form("ToT #in [%.1f-%.1f] ns", tot_min, tot_max));
    cuts->AddText(Form("Lead. edges #in [%.1f-%.2f]  #mus", lead_min, lead_max));
    cuts->Draw("same");

    TF1* fit = new TF1("g2", doublegauss, -range, range, 7);
    fit->SetParameter(0, h_diff_all_channels->GetMaximum()); 
    fit->SetParameter(1, 0.); 
    fit->SetParameter(3, h_diff_all_channels->GetMaximum()/2.); 
    fit->SetParameter(4, 0.);
    fit->SetParameter(6, 400.);
    
    if (h_diff_all_channels->Integral()>0.) {
      TFitResultPtr fr = h_diff_all_channels->Fit(fit,"sr0");
      if (fr->IsValid() and !fr->IsEmpty()) {
        fit->Draw("same");
        fit->SetLineColor(kBlack);
        fit->SetLineStyle(2);
        fit->SetRange(-10., 10.);
        c.AddLegendEntry(fit, Form("Gaussian fit, #mu = %.2f ns, #sigma = %.2f ps", fr->Parameter(1), fr->Parameter(2)*1.e3), "l");
        c.AddLegendEntry(fit, Form("(#chi^{2}/ndf = %.1f / %d)", fr->Chi2(), fr->Ndf()), "");
      }
    }
    c.Prettify(h_diff_all_channels);
    h_diff_all_channels->SetLineColor(kGray+1);
    h_diff_all_channels->GetYaxis()->SetTitleOffset(1.);
    double max = h_diff_all_channels->GetMaximum()*15;
    c.Pad()->SetLogy();
    h_diff_all_channels->SetMaximum(max);
    h_diff_all_channels->SetMinimum(10.);
    c.Save("png", "time_diff_betw_chan");
    c.Save("pdf", "time_diff_betw_chan");
    c.Save("root", "time_diff_betw_chan");
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
    TLine* line = new TLine(tot_min, 0., tot_min, h_tot_all_channels->GetMaximum()*0.8);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");
  }
  if (tot_max>0.) {
    TLine* line = new TLine(tot_max, 0., tot_max, h_tot_all_channels->GetMaximum()*0.8);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");
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

    TPad* zoom = new TPad("zoom", "", 0.14, 0.45, 0.62, 0.87);
    zoom->SetLeftMargin(0.15);
    zoom->Draw();
    zoom->SetTicks(1,1);
    zoom->cd();
    h_lead_all_channels_zoom->Draw();
    h_lead_all_channels_zoom->SetMinimum(0.001);
    h_lead_all_channels_zoom->SetLabelFont(43, "xy");
    h_lead_all_channels_zoom->SetLabelSize(18, "xy");
    h_lead_all_channels_zoom->SetLineColor(kGray+1);
    h_lead_all_channels_zoom->SetLineWidth(2);
    h_lead_all_channels_zoom->GetYaxis()->SetNdivisions(505);
    h_lead_all_channels_zoom->GetXaxis()->SetNdivisions(505, true);
    {
      TLine* line1 = new TLine(lead_min, 0., lead_min, h_lead_all_channels_zoom->GetMaximum()*1.05);
      line1->SetLineColor(kBlack);
      line1->SetLineStyle(2);
      line1->SetLineWidth(2);
      line1->Draw("same");
      TLine* line2 = new TLine(lead_max, 0., lead_max, h_lead_all_channels_zoom->GetMaximum()*1.05);
      line2->SetLineColor(kBlack);
      line2->SetLineStyle(2);
      line2->SetLineWidth(2);
      line2->Draw("same");
    }

    c.cd();
    c.Prettify(h_lead_all_channels);
    h_lead_all_channels->GetXaxis()->SetTitle("Leading edge ( #mus)");
    h_lead_all_channels->GetYaxis()->SetTitle("Hits / Trigger");
    h_lead_all_channels->SetLineColor(kGray+1);
    h_lead_all_channels->GetXaxis()->SetLabelSize(18);
    c.Save("png", "time_diff_betw_chan");
    c.Save("pdf", "time_diff_betw_chan");
  }

  {
    GastofCanvas c("mult_one-allchannels", Form("Run %d, all channels", run_id));
    h_mult_all_channels->Draw();

    TPaveText* cuts = new TPaveText(0.4, 0.68, 0.5, 0.78, "ndc");
    cuts->SetTextAlign(13);
    cuts->SetTextFont(43);
    cuts->SetTextSize(18);
    cuts->SetFillColor(kWhite);
    cuts->SetLineColor(kWhite);
    cuts->AddText(Form("Lead. edges #in [%.1f-%.2f]  #mus", lead_min, lead_max));
    cuts->Draw("same");

    //c.Pad()->SetLogy();
    c.Prettify(h_mult_all_channels);
    h_mult_all_channels->GetXaxis()->SetNdivisions(110);
    h_mult_all_channels->GetXaxis()->SetTitle("Hits multiplicity");
    h_mult_all_channels->GetYaxis()->SetTitle("#LTHits#GT / Trigger");
    c.Save("png", "time_diff_betw_chan");
  }

}
