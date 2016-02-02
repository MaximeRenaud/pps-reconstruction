#include "PPSCanvas.h"
// #src/#include "FileReader.cpp" (?) 
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
// #include "TVirtualFFT.h"
#include "GastofCanvas.h"
#include "GastofGrid.h"
#include "GastofCommons.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>


#define MAX_MEAS 5000

using namespace std;


int main(int argc, char** argv) {  //NoiseIssue()
  
  
  //Def des variables globales
  
  int fNumMeasurements; //fNumErrors;
//   int fRunId;
//   int fETTT, fTriggerNumber;
  int fChannelId[MAX_MEAS];
  double fLeadingEdge[MAX_MEAS]/*, fTrailingEdge[MAX_MEAS]*/, fToT[MAX_MEAS];
  double NoiseSteps[32];
//   double summing;
  
  const unsigned int run_id = (argc >1 ?atoi(argv[1]):725); //if nothing is given, default run_id will be 725 (atoi()=turn into int)
//   const unsigned int channels_to_probe = 32; //64;
  const double lead_min = 6.8, lead_max = 7.1; // in \mus
//   int Channel = 4;
  
  
  
  
  
  //Def des objets + Initialisation
  
  TFile* file = new TFile(Form("../../gastof_inner_run%d.root", run_id));
  if (file == 0) {
      // if we cannot open the file, print an error message and return immediatly
      printf("Error: cannot open gastof_inner_run725.root !\n");
      return 0;
  } cout <<endl<<"Yay ! It worked !"<<endl<<endl;
  
  TTree* tree = (TTree*)file->Get("tdc");
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("num_measurements", 1); tree->SetBranchAddress("num_measurements", &fNumMeasurements);
  tree->SetBranchStatus("leading_edge", 1); tree->SetBranchAddress("leading_edge", fLeadingEdge);
  tree->SetBranchStatus("channel_id", 1); tree->SetBranchAddress("channel_id", fChannelId);
  tree->SetBranchStatus("tot", 1); tree->SetBranchAddress("tot", fToT);
  
  const unsigned int num_triggers = tree->GetEntries();
  
//   TH1D* test = new TH1D("test", "Leading edge time for channel 4", 100,0,10100);
  TH1D* h_tot_all_channels = new TH1D("h_tot_all_channels", "", 250, 0., 1000.);
  TH1D* h_lead_all_channels = new TH1D("h_lead_all_channels", "", 500, 0., 10.);
  
  for (unsigned int i=0; i<num_triggers; i++) {
    tree->GetEntry(i);
    if (i%10000==0) cerr << "Processing event " << i << " / " << num_triggers << endl;
    for (int j=0; j<fNumMeasurements; j++) {
      h_lead_all_channels->Fill(fLeadingEdge[j]/1.e3, 1./num_triggers);  //Fill(x,w) <- fill \w weight; here: renormalise (1./1000. for beauty)
      if (fLeadingEdge[j]<lead_min*1.e3 or fLeadingEdge[j]>lead_max*1.e3) continue;
      h_tot_all_channels->Fill(fToT[j]);
    }
  }
  
  
  //Initialisation of Canvases
  
//   GastofCanvas c_tot_all_chan("tot_all-channels", Form("Run %d, all channels", run_id));
//   h_tot_all_channels->Draw();
//   h_tot_all_channels->GetXaxis()->SetTitle("Time over threshold (ns)");
//   h_tot_all_channels->GetYaxis()->SetTitle("Hits");
//   
//   c_tot_all_chan.Pad()->SetLogy();
//   
// //   c_tot_all_chan.Prettify(h_tot_all_channels);
// //   h_tot_all_channels->SetLineColor(kGray+1);
// //   c_tot_all_chan.Save("png", "time_diff_betw_chan");
// //   c_tot_all_chan.Save("pdf", "time_diff_betw_chan");
//   
//   
//   GastofCanvas c("leading_all-channels", Form("Run %d, all channels", run_id));
//   h_lead_all_channels->Draw();
//   
//   {
//       TLine* line1 = new TLine(lead_min, 0., lead_min, h_lead_all_channels->GetMaximum()*1.05);
//       line1->SetLineColor(kBlack);
//       line1->SetLineStyle(2);
//       line1->SetLineWidth(2);
//       line1->Draw("same");
//       TLine* line2 = new TLine(lead_max, 0., lead_max, h_lead_all_channels->GetMaximum()*1.05);
//       line2->SetLineColor(kBlack);
//       line2->SetLineStyle(2);
//       line2->SetLineWidth(2);
//       line2->Draw("same");
//     }
//   
//   c.Prettify(h_lead_all_channels);
//   h_lead_all_channels->GetXaxis()->SetTitle("Leading edge ( #mus)");
//   h_lead_all_channels->GetYaxis()->SetTitle("Hits / Trigger");
//   h_lead_all_channels->SetLineColor(kGray+1);
//   h_lead_all_channels->GetXaxis()->SetLabelSize(18);
  
  
  
  
  return 0;
  
  
  
}