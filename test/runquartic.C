

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <list>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include <math.h>
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLorentzVector.h"


double runquartic(TString,int,TString);
double getmeasurements(TString,int, int,TString);


double getmeasurements(TString fileinput, int triggernumber , TString globalvariable, int index ) {



  //Define the input and output files
  TFile finput(fileinput);
  


  //store the tdc tree in the inputtree variable
  TTree *inputtree = (TTree*)finput.Get("tdc");
  
  //Number of entries (events) in the input tree
  Int_t nentries = (Int_t)inputtree->GetEntries();

  // Read Branches you are interested in from the input tree
  double trailing_edge[nentries]; inputtree->SetBranchAddress("trailing_edge",&trailing_edge);
  double leading_edge[nentries]; inputtree->SetBranchAddress("leading_edge",&leading_edge);
  double tot[nentries]; inputtree->SetBranchAddress("tot",&tot);
  Int_t channel_id[nentries]; inputtree->SetBranchAddress("channel_id",&channel_id);
  
  
  
    
   inputtree->GetEntry(triggernumber);

 

if (globalvariable=="trailing_edge") return trailing_edge[index];
if (globalvariable=="leading_edge") return leading_edge[index];
if (globalvariable=="channel_id") return (double)channel_id[index];
if (globalvariable=="tot") return tot[index];
}




double runquartic(TString fileinput, int triggernumber , TString globalvariable ) {


  //Define the input and output files
  TFile finput(fileinput);


  //store the tdc tree in the inputtree variable
  TTree *inputtree = (TTree*)finput.Get("tdc");
  

  // Read Branches you are interested in from the input tree
  UInt_t num_measurements; inputtree->SetBranchAddress("num_measurements",&num_measurements);
  ULong64_t ettt; inputtree->SetBranchAddress("ettt",&ettt);
  

    
    
 inputtree->GetEntry(triggernumber);

  

if (globalvariable=="num_measurements") return (double)num_measurements;
if (globalvariable=="ettt") return (double) ettt;

}


