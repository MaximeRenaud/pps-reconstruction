#include "FileReader.h"
#include "QuarticCanvas.h"
#include <dirent.h>

#include "TFile.h"
#include "TTree.h"

#define MAX_MEAS 50000

using namespace std;

int main(int argc, char* argv[]) {
  if (argc<3) { cerr << "Usage: " << argv[0] << " <run id> <trigger start> [trigger stop=-1]" << endl; exit(0); }
  const unsigned int num_channels = 64;
  int run_id = atoi(argv[1]);
  int trigger_start = atoi(argv[2]), trigger_stop = -1;
  if (argc>3) trigger_stop = atoi(argv[3]);
  
  string output = "output_test_sorted.root";
  if (argc>4) output = argv[4];

  int fNumMeasurements; 
  int fNumErrors = 0;
  int fRunId = run_id, fBurstId = 0;
  unsigned int fETTT[2], fTriggerNumber; // unsigned int should have all 32 bits meaningfull
  int fChannelId[MAX_MEAS];
  double fLeadingEdge[MAX_MEAS], fTrailingEdge[MAX_MEAS], fToT[MAX_MEAS];
  unsigned int fRawLeadingEdge[MAX_MEAS], fRawTrailingEdge[MAX_MEAS], fRawToT[MAX_MEAS];

  const unsigned int boards[2] = { 4, 3 }; // GasToF boards
  DIR* dir; struct dirent* ent;
  cout << "Search in directory: " << getenv("PPS_DATA_PATH") << endl;

  const unsigned int max_files = 500;
  std::vector<VME::TDCTrigger> trig[2];
  VME::TDCTrigger tr;

  for (unsigned int i=0; i<2; i++) {
    unsigned int num_trig = 0;
    bool finished = false;
    for (unsigned int sp=1; sp<max_files; sp++) { // we loop over all files
      bool file_found = false; string filename;
      // first we search for the proper file to open
      if ((dir=opendir(getenv("PPS_DATA_PATH")))==NULL) return -1;    //PPS_DATA_PATH = environment variable, stored in .bashrc of machine executing
      while ((ent=readdir(dir))!=NULL) {
        if (string(ent->d_name).find(Form("events_%d_%d_", run_id, sp))!=string::npos and
            string(ent->d_name).find(Form("_board%d", boards[i]))!=string::npos) {
            file_found = true;
            filename = ent->d_name;
          break;
        }
      }
      closedir(dir);
      if (!file_found) { cout << "Found " << sp << " files in this run" << endl; break; }

      // then we open it
      cout << "Opening file " << filename << endl;
      try {
        FileReader f(Form("%s/%s", getenv("PPS_DATA_PATH"), filename.c_str()));
        while (true) {
          if (!f.GetNextTrigger(&tr)) break;
          if (trigger_stop>0 and num_trig>=(unsigned int)trigger_stop) { finished = true; break; }
          trig[i].push_back(tr);
          num_trig++;
        }
      //      f.Clear(); // we return to beginning of the file
      } catch (Exception& e) { e.Dump(); }
      if (finished) break;
    }
  }
  cout << "Number of triggers: " << trig[0].size() << " / " << trig[1].size() << std::endl;

  /*for (unsigned int i=0; i<trig[0].size(); i++) {
    cout << i << ": " << trig[0][i].GetETTT() << " -- " << trig[1][i].GetETTT() << " %%% " << trig[1][i].GetETTT()-trig[0][i].GetETTT() << endl;
  }*/

  unsigned int tot_num_triggers = trig[0].size();
  if (trig[0].size()!=trig[1].size()) {
    cerr << "... not the same number!" << endl;
    tot_num_triggers = TMath::Min(trig[0].size(), trig[1].size());
    //return -1;
  }

  // Start to fill the tree

  TFile* f = new TFile(output.c_str(), "recreate");
  TTree* t = new TTree("tdc", "List of TDC measurements");
  t->Branch("num_measurements", &fNumMeasurements, "num_measurements/I");
  t->Branch("num_errors", &fNumErrors, "num_errors/I");
  t->Branch("run_id", &fRunId, "run_id/I");
  t->Branch("burst_id", &fBurstId, "burst_id/I");
  t->Branch("channel_id", fChannelId, "channel_id[num_measurements]/I");
  t->Branch("leading_edge", fLeadingEdge, "leading_edge[num_measurements]/D");
  t->Branch("trailing_edge", fTrailingEdge, "trailing_edge[num_measurements]/D");
  t->Branch("tot", fToT, "tot[num_measurements]/D");

  /*t->Branch("raw_leading_edge", fRawLeadingEdge,   "raw_leading_edge[num_measurements]/i");
  t->Branch("raw_trailing_edge", fRawTrailingEdge, "raw_trailing_edge[num_measurements]/i");
  t->Branch("raw_tot", fRawToT, "raw_tot[num_measurements]/i");*/

  t->Branch("ettt", fETTT, "ettt[2]/i"); // i for UInt_32
  t->Branch("trigger_number",&fTriggerNumber,"trigger_number/i");
 
  //const unsigned int ettt_diff_0 = trig[1][0].GetETTT()-trig[0][0].GetETTT();
  unsigned int last_ettt = 0;
  for (unsigned int i=0; i<tot_num_triggers; i++) { // loop over the triggers
    if (trigger_start>0 and i<(unsigned int)trigger_start) continue;
    if (trigger_stop>0 and i>=(unsigned int)(trigger_stop-trigger_start)) break;

    fNumMeasurements = 0;
    for (unsigned int j=0; j<2; j++) { // loop over the boards
      fETTT[j] = trig[j][i].GetETTT();
      if (j==0) {
        if ((fETTT[j]-last_ettt)*25.e9>2.) fBurstId++;
        last_ettt = fETTT[j];
      }
      fNumErrors += trig[j][i].GetErrors().size();
      for (unsigned int k=0; k<num_channels; k++) {
        std::vector<VME::TDCTrigger::Hit> hits = trig[j][i].GetHits(k);
        for (std::vector<VME::TDCTrigger::Hit>::const_iterator it=hits.begin(); it!=hits.end(); it++) {
          fChannelId[fNumMeasurements] = k+j*32;

          fRawLeadingEdge[fNumMeasurements] = it->first;
          fRawTrailingEdge[fNumMeasurements] = it->second;
          fRawToT[fNumMeasurements] = fRawTrailingEdge[fNumMeasurements]-fRawLeadingEdge[fNumMeasurements];

          fLeadingEdge[fNumMeasurements] = fRawLeadingEdge[fNumMeasurements]*25./1024.;
          fTrailingEdge[fNumMeasurements] = fRawTrailingEdge[fNumMeasurements]*25./1024.;
          fToT[fNumMeasurements] = fTrailingEdge[fNumMeasurements]-fLeadingEdge[fNumMeasurements];

          fNumMeasurements++;
        }
      }
    }
    fTriggerNumber = i;
    t->Fill();
  }
  t->Write();
  f->Close();
  
  return 0;
}
