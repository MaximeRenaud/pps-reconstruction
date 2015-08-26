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
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "runquartic.C"

/*
 * Simple macro for 1-D fitting of tracks from RD51 GEM tracker cluster ntuples and extraploations to Quartic detector
 *
 * Fits 1 track per event in x-z and y-z planes. In case of multiple clusters per layer, the x- and y- clusters with the highest charge are 
 * used (following RD51 macros). Uncertainties on cluster position in x- and y- are currently set to 100micron. Offset corrections 
 * are applied to align the peak of the cluster position distributions of layers 1 and 3 with layer 2. 
 *
 * Usage: [0] .L runtracker_trackfits.C+
 *        [1] runtracker_trackfits(fileinput,quarticfileinput,event,plotvar,chicut,maxclusterperlayer)
 * 
 * Options:
 *     fileinput: tracker input ntuple filename (should contain a TCluster tree) 
 *     quarticfileinput: Quartic input ntuple filename (should contain a tdc tree)
 *     event: if >0, only that evtID will be processed. If =-1, all events in the ntuple will be processed
 *     plotvar: variable to plot
 *              1: "event display" of clusters and tracks
 *              2: chi2 of fitted tracks
 *              3: 1D x- and y- positions of tracks extrapolated to quartic#1 at 1330mm from last GEM plane
 *              4: multiplicities of x- and y- clusters in each GEM layer
 *              5: 2D x- and y- positions of clusters in each GEM layer
 *              6: 2D x- and y- positions of tracks extrapolated to quartic#1 at 1330mm from last GEM plane
 *              7: Quartic hit-map for leading-trailing edge pairs passing cuts in events with tracks passing track cuts
 *     chicut: select events passing this chi^2 cut applied to fitted tracks
 *     maxclusterperlayer: select events with at most this number (counting x-y combinations) of clusters per GEM layer
 *     minquartictot: minimum time-over-threshold for quartic
 *     maxquartictot: maximum time-over-threshold for quartic 
 *     minquarticleadedge: minimum leading-edge time for quartic
 *     maxquarticleadedge: maximum leading-edge time for quartic
 */

void runtracker_trackfits(TString fileinput, 
			  TString quarticfileinput, 
			  Int_t event=-1, 
			  Int_t plotvar=1, 
			  Float_t chicut=3.0, 
			  Int_t maxclustperlayer=1,
			  Float_t minquartictot=5.0,
			  Float_t maxquartictot=30.0,
			  Float_t minquarticleadedge=6800.0,
			  Float_t maxquarticleadedge=6900.0) {

  // Declare histograms
  TH2F *h1 = new TH2F("h1","h1",100,-50,50,100,-50,50); 
  TH2F *h2 = new TH2F("h2","h2",100,-50,50,100,-50,50);  
  TH2F *h3 = new TH2F("h3","h3",100,-50,50,100,-50,50);  
  TH2F *h1c = new TH2F("h1c","h1c",100,-50,50,100,-50,50);  
  TH2F *h2c = new TH2F("h2c","h2c",100,-50,50,100,-50,50);   
  TH2F *h3c = new TH2F("h3c","h3c",100,-50,50,100,-50,50);   
  
  TH1F *hncx1 = new TH1F("hncx1","hncx1",50,0,50); 
  TH1F *hncy1 = new TH1F("hncy1","hncy1",50,0,50);  
  TH1F *hncx2 = new TH1F("hncx2","hncx2",50,0,50);  
  TH1F *hncy2 = new TH1F("hncy2","hncy2",50,0,50);  
  TH1F *hncx3 = new TH1F("hncx3","hncx3",50,0,50);  
  TH1F *hncy3 = new TH1F("hncy3","hncy3",50,0,50);  
  hncx1->SetTitle("x cluster multiplicity, plane 1"); 
  hncx2->SetTitle("x cluster multiplicity, plane 2");  
  hncx3->SetTitle("x cluster multiplicity, plane 3");  
  hncy1->SetTitle("y cluster multiplicity, plane 1");  
  hncy2->SetTitle("y cluster multiplicity, plane 2");   
  hncy3->SetTitle("y cluster multiplicity, plane 3");   
  
  TH1F *hchi2xz = new TH1F("hchi2xz","hchi2xz",200,0,50); 
  TH1F *hchi2yz = new TH1F("hchi2yz","hchi2yz",200,0,50);  
  hchi2xz->SetTitle("#chi^{2} (x-z plane)"); 
  hchi2yz->SetTitle("#chi^{2} (y-z plane)"); 
 
  TH1F *hxquartic1 = new TH1F("hxquartic1","hxquartic1",100,-50,50); 
  TH1F *hyquartic1 = new TH1F("hyquartic1","hyquartic1",100,-50,50);  
  TH2F *hxyquartic1 = new TH2F("hxyquartic1","hxyquartic1",500,-50,50,500,-50,50); 
  // Testing absolute alignment of quartic and tracker w.r.t. beamline 
  //  TH2F *hxyquartic1 = new TH2F("hxyquartic1","hxyquartic1",14,-2.0,40.0,8,-22.0,-4.0); 
  hxquartic1->SetTitle("Track x-position at quartic #1 [mm]"); 
  hyquartic1->SetTitle("Track y-position at quartic #1 [mm]");  
  hxyquartic1->SetTitle("Track x-y position at quartic #1 [mm]"); 

  // Testing absolute alignment of quartic and tracker w.r.t. beamline 
  //  TH2F *hquarticocc = new TH2F("hquarticocc","hquarticocc",14,-2.0,40.0,8,-22.0,-4.0);

  TH2F *hquarticocc = new TH2F("hquarticocc","hquarticocc",14,-2.5,39.5,8,-2.5,21.5);

  /*
   * Tracks and fitting
   */
  Float_t xs[3]; 
  Float_t ys[3]; 
  Float_t zs[3] = {-462, 0, 267};   // From RD51 macros - z-position of the 3 GEM planes 
  
  Float_t xerr[3] = {0.1,0.1,0.1}; 
  Float_t yerr[3] = {0.1,0.1,0.1}; 
  Float_t zerr[3] = {0,0,0}; 

  TGraphErrors *trkxz; 
  TGraphErrors *trkyz; 
  TF1 *f1 = new TF1("f1","pol1",-3000,3000);   
  TF1 *f2 = new TF1("f2","pol1",-3000,3000);    
  TFitResultPtr ftxz; 
  TFitResultPtr ftyz; 

  // Testing absolute alignment of quartic and tracker w.r.t. beamline 
  Float_t beamshiftxlayer3 = 50.0;

  // Simple alignment offset constants - derived from peak of cluster position distribution, relative to middle plane 
  Float_t xoffset1 = 4.5; 
  Float_t xoffset2 = 0;    
  Float_t xoffset3 = 1.2;  
  Float_t yoffset1 = -1.8;   
  Float_t yoffset2 = 0;         
  Float_t yoffset3 = -6.5;     

  // Testing absolute alignment of quartic and tracker w.r.t. beamline 
  //  Float_t xoffset1 = beamshiftxlayer3 + 3.3;
  //  Float_t xoffset2 = beamshiftxlayer3 - 1.2;
  //  Float_t xoffset3 = 0 + beamshiftxlayer3;
  //  Float_t yoffset1 = 4.7;
  //  Float_t yoffset2 = 6.5;
  //  Float_t yoffset3 = 0;

  // Canvas for "event displays" - declare this before the main loop so we can keep updating it with new tracks
  // Define position of GEM planes and quartic
  TCanvas *c1 = new TCanvas("c1","c1"); 
  c1->Divide(2,1); 
  TH2F *hbox = new TH2F("hbox","hbox",2,-1000,1700,2,-100,100);  
  hbox->SetStats(0); 
  TLine *lgem1 = new TLine(-462,-50,-462,50);  
  lgem1->SetLineWidth(3);  
  TLine *lgem2 = new TLine(0,-50,0,50);   
  lgem2->SetLineWidth(3);   
  TLine *lgem3 = new TLine(267,-50,267,50);   
  lgem3->SetLineWidth(3);   
  TLine *lq1 = new TLine(267+1330,-100,267+1330,100);  
  lq1->SetLineColor(4); lq1->SetLineWidth(3);  


  c1->cd(1);    
  hbox->SetTitle("GEM tracks (x-z plane)"); 
  hbox->Draw();  
  lgem1->Draw("same"); lgem2->Draw("same"); lgem3->Draw("same");   

  c1->cd(2); 
  hbox->SetTitle("GEM tracks (y-z plane)"); 
  hbox->Draw();   
  lgem1->Draw("same"); lgem2->Draw("same"); lgem3->Draw("same");  

 
  Int_t x_bin;  
  Int_t y_bin;  
 
  Int_t nclustxlayer1=0; 
  Int_t nclustylayer1=0; 
  Int_t nclustxlayer2=0; 
  Int_t nclustylayer2=0; 
  Int_t nclustxlayer3=0;  
  Int_t nclustylayer3=0;  
 
  Int_t npassing=0; 

  //Define the input and output files
  TFile finput(fileinput);

  // Commented out output for now
  //  TFile foutput(fileoutput,"recreate");
  //  ofstream ofs(fileoutput);
 
  //store the tdc tree in the inputtree variable
  TTree *inputtree = (TTree*)finput.Get("TCluster");
  
  //Number of entries (events) in the input tree
  Int_t nentries = (Int_t)inputtree->GetEntries();

  // Read Branches you are interested in from the input tree
  Int_t evtID; inputtree->SetBranchAddress("evtID",&evtID);
  Int_t nclust; inputtree->SetBranchAddress("nclust",&nclust);
  float clustPos[nentries]; inputtree->SetBranchAddress("clustPos",&clustPos);
  float clustADCs[nentries]; inputtree->SetBranchAddress("clustADCs",&clustADCs);
  Int_t clustSize[nentries]; inputtree->SetBranchAddress("clustSize",&clustSize);
  Int_t clustTimebin[nentries]; inputtree->SetBranchAddress("clustTimebin",&clustTimebin);
  UInt_t detID[nentries]; inputtree->SetBranchAddress("detID",&detID);
  UInt_t planeID[nentries]; inputtree->SetBranchAddress("planeID",&planeID);
  Short_t etaSector[nentries];   inputtree->SetBranchAddress("etaSector",&etaSector);

  // Start main loop over entries in tracking ntuple
  Int_t jentry=0;
  for (jentry=0; jentry<nentries;jentry++)   
    {  
      inputtree->GetEntry(jentry);   

      // Keep reading until we get to an event we're interested in, but not beyond
      if(event > 0 && event > evtID)   
	continue;   
      if(event > 0 && event < evtID)  
	break;  
        
      if(evtID % 5000 == 0)  
	cout << "Evt ID: " << evtID << endl;  

      // Loop over the 3 layers of the GEM tracker
      for(unsigned int layer = 0;layer<3;layer++) 
	{ 
	  Float_t cluster_x    = 0;  
	  Float_t charge_x     = 0;  
	  Float_t size_x       = 0;  
            
	  Float_t cluster_y    = 0;  
	  Float_t charge_y     = 0;  
	  Float_t size_y       = 0;  
            
	  Int_t IdMaxCluster_x=0;  
	  Int_t IdMaxCluster_y=0;  

	  // Find the maximally charged 1D cluster in x and y
	  Float_t IdMaxClusterCharge_x=0;  
	  Float_t IdMaxClusterCharge_y=0;  
	  for (Int_t i = 0; i<nclust; i++){  
	    if (detID[i]==layer && planeID[i]==0)  
	      {  
		if (IdMaxClusterCharge_x<=clustADCs[i])   
		  {  
		    IdMaxCluster_x = i;  
		    IdMaxClusterCharge_x = clustADCs[i];  
		  }  
		if(layer==0) 
		  nclustxlayer1++; 
		if(layer==1)  
		  nclustxlayer2++;  
		if(layer==2)  
		  nclustxlayer3++;  
 
	      }  
	    else if (detID[i]==layer && planeID[i]==1)  
	      {  
		if (IdMaxClusterCharge_y<=clustADCs[i])  
		  {  
		    IdMaxCluster_y = i;  
		    IdMaxClusterCharge_y = clustADCs[i];  
		  }  
		if(layer==0)  
		  nclustylayer1++;  
		if(layer==1)   
		  nclustylayer2++;   
		if(layer==2)   
		  nclustylayer3++;   
	      }  
	  }  
            
	  charge_x     = IdMaxClusterCharge_x;  
	  charge_y     = IdMaxClusterCharge_y;  
            
	  // Require at least 1 cluster with non-zero charge
	  if (charge_x>0 && charge_y>0)  
	    {  
	      cluster_x = clustPos[IdMaxCluster_x];  
	      x_bin = TMath::FloorNint(cluster_x);  
	      size_x = clustSize[IdMaxCluster_x];  
	      cluster_y = clustPos[IdMaxCluster_y];  
	      y_bin = TMath::FloorNint(cluster_y);  
	      size_y = clustSize[IdMaxCluster_y];  

	      // Apply layer-by-layer alignment offsets and fill cluster position maps
	      if(layer == 0)  
		{ 
		  h1->Fill(cluster_x,cluster_y);  
		  xs[0]=cluster_x + xoffset1; 
		  ys[0]=cluster_y + yoffset1; 
		  h1c->Fill(xs[0],ys[0]); 
		} 
	      if(layer == 1)  
		{ 
		  h2->Fill(cluster_x,cluster_y);  
		  xs[1]=cluster_x + xoffset2;  
		  ys[1]=cluster_y + yoffset2;  
		  h2c->Fill(xs[1],ys[1]); 
		} 
	      if(layer == 2)  
		{ 
		  // After the last (3rd) layer, fit tracks, apply event-level cuts, and make plots

		  h3->Fill(cluster_x,cluster_y);  
		  xs[2]=cluster_x + xoffset3; 
		  ys[2]=cluster_y + yoffset3; 
		  h3c->Fill(xs[2],ys[2]); 
 
		  // Fit track in x-z plane 
		  trkxz = new TGraphErrors(3,zs,xs,zerr,xerr); 
		  trkxz->SetMarkerColor(2); 
		  trkxz->SetMarkerStyle(20); 
		  ftxz = trkxz->Fit(f1,"SQ","",-500,1800);  
		  float par0xz = f1->GetParameter(0);  
		  float par1xz = f1->GetParameter(1);  
 
		  // Fit track in y-z plane 
		  trkyz = new TGraphErrors(3,zs,ys,zerr,yerr);  
		  trkyz->SetMarkerColor(2);  
		  trkyz->SetMarkerStyle(20);  
		  ftyz = trkyz->Fit(f2,"SQ","",-500,1800);   
		  float par0yz = f2->GetParameter(0);   
		  float par1yz = f2->GetParameter(1);   
 
		  // Extrapolate to position of quartic#1 (quartz+sapphire) 
		  float zquartic1 = zs[2] + 1330.0;                
		  float xquartic1 = (zquartic1*par1xz) + par0xz; 
		  float yquartic1 = (zquartic1*par1yz) + par0yz; 
 
		  // Apply any track quality cuts passed by the user
		  if((ftxz->Chi2() < chicut) && (ftyz->Chi2() < chicut)) 
		    { 
		      if(nclustxlayer1<=maxclustperlayer && 
                          nclustylayer1<=maxclustperlayer &&  
                          nclustxlayer2<=maxclustperlayer &&  
                          nclustylayer2<=maxclustperlayer &&   
                          nclustxlayer3<=maxclustperlayer &&  
			 nclustylayer3<=maxclustperlayer) 
			{   
			  npassing++; 
 
			  hchi2xz->Fill(ftxz->Chi2()); 
			  hchi2yz->Fill(ftyz->Chi2()); 
 
			  hxquartic1->Fill(xquartic1); 
			  hyquartic1->Fill(yquartic1); 
			  hxyquartic1->Fill(xquartic1,yquartic1); 
 
			  hncx1->Fill(nclustxlayer1); 
			  hncy1->Fill(nclustylayer1);  
			  hncx2->Fill(nclustxlayer2);  
			  hncy2->Fill(nclustylayer2);  
			  hncx3->Fill(nclustxlayer3);  
			  hncy3->Fill(nclustylayer3);  

                          // Now get quartic information for this event if it passes tracking cuts, based on trigger # 
			  //			  ofs << evtID << "\t" << xquartic1 << "\t" << yquartic1 << endl;
			  int nquartichits = runquartic(quarticfileinput, evtID , "num_measurements");
			  for(int x=0; x<nquartichits; x++)
			    {
			      float leadedge = getmeasurements(quarticfileinput, evtID, "leading_edge" ,x); 
                              float tot = getmeasurements(quarticfileinput, evtID, "tot" ,x);  
                              int chid = getmeasurements(quarticfileinput, evtID, "channel_id" ,x);  

			      int qx=0;
			      int qy=0;
                              float qxoffset = 0.0; 
                              float qyoffset = 0.0; 
 			      
			      // Testing absolute alignment of quartic and tracker w.r.t. beamline
			      //			      float qxoffset = 13.0;
			      //			      float qyoffset = -16.0;

			      // Now apply any quartic cuts, mapping, and fill histogram
			      // Cuts here are illustrative, may want to make these configurable
			      if(leadedge > minquarticleadedge && leadedge < maxquarticleadedge && tot > minquartictot && tot < maxquartictot)
				{
				  // Mapping from HPTDC channel to physical quartic bar position
				  switch (chid) 
				    {
				    case 0:  qx = qxoffset + 13.5; qy = qyoffset + 4.5; break;
				    case 1:  qx = qxoffset + 13.5; qy = qyoffset + 10.5; break;
				    case 4:  qx = qxoffset + 13.5; qy = qyoffset + 1.5; break;
				    case 5:  qx = qxoffset + 13.5; qy = qyoffset + 7.5; break;
				    case 8:  qx = qxoffset + 10.5; qy = qyoffset + 4.5; break;
				    case 9:  qx = qxoffset + 10.5; qy = qyoffset + 10.5; break;
				    case 12: qx = qxoffset + 10.5; qy = qyoffset + 1.5; break;
				    case 13: qx = qxoffset + 10.5; qy = qyoffset + 7.5; break;
				    case 16: qx = qxoffset + 7.5; qy = qyoffset + 4.5; break;
				    case 17: qx = qxoffset + 7.5; qy = qyoffset + 10.5; break;
				    case 20: qx = qxoffset + 7.5; qy = qyoffset + 1.5; break;
				    case 21: qx = qxoffset + 7.5; qy = qyoffset +  7.5; break;
				    case 24: qx = qxoffset + 4.5; qy = qyoffset + 4.5; break;
				    case 25: qx = qxoffset + 4.5; qy = qyoffset + 10.5; break;
				    case 26: qx = qxoffset + 4.5; qy = qyoffset + 1.5; break;
				    case 27: qx = qxoffset + 4.5; qy = qyoffset + 7.5; break;
				    case 28: qx = qxoffset + 1.5*1; qy = qyoffset + 4.5; break;
				    case 29: qx = qxoffset + 1.5*1; qy = qyoffset + 10.5; break;
				    case 30: qx = qxoffset + 1.5*1; qy = qyoffset + 1.5; break;
				    case 31: qx = qxoffset + 1.5*1; qy = qyoffset + 7.5; break;
				    }
				  hquarticocc->Fill(qx,qy);
				}
			    }

			  if(plotvar == 1) 
			    { 
			      c1->cd(1); 
			      trkxz->SetMarkerStyle(20);  
			      trkxz->SetMarkerColor(2);  
			      trkxz->Draw("Psame"); 
			      ftxz->Draw("same");  
			      c1->cd(2); 
			      trkyz->SetMarkerStyle(20); 
			      trkyz->SetMarkerColor(2); 
			      trkyz->Draw("Psame"); 
			      ftxz->Draw("same"); 
			    } 
			} 
		    } 
 
		  nclustxlayer1=0; 
		  nclustxlayer2=0; 
		  nclustxlayer3=0; 
		  nclustylayer1=0;  
		  nclustylayer2=0;  
		  nclustylayer3=0;  
		  //              delete trkxz; 
                    
		} 
	    } 
	} 
    }

  /*
   * We're done looping over all events - produce summary plots
   */

  if(plotvar == 1)  // "Event display" 
    { 
      c1->cd(1); 
      lq1->Draw("same");   
      c1->cd(2); 
      lq1->Draw("same"); 
    } 
  if(plotvar == 2)  // Chi2 of track fits 
    { 
      TCanvas *c5 = new TCanvas("c5","c5"); 
      c5->Divide(2,1); 
      c5->cd(1); 
      hchi2xz->Draw("hist");  
      c5->cd(2); 
      hchi2yz->Draw("hist"); 
    } 
  if(plotvar == 3)  // x- and y-position of tracks extrapolated to z-position of quartic 
    { 
      TCanvas *c2 = new TCanvas("c2","c2"); 
      c2->Divide(2,1); 
      c2->cd(1); 
      hxquartic1->Draw("hist");  
      c2->cd(2); 
      hyquartic1->Draw("hist"); 
    } 
  if(plotvar == 4) // 1D cluster multiplicities per GEM layer
    {  
      TCanvas *c3 = new TCanvas("c3","c3");  
      c3->Divide(3,2);  
      c3->cd(1);  
      hncx1->Draw("hist");   
      c3->cd(2);   
      hncy1->Draw("hist");    
      c3->cd(3);   
      hncx2->Draw("hist");    
      c3->cd(4);   
      hncy2->Draw("hist");    
      c3->cd(5);   
      hncx3->Draw("hist");    
      c3->cd(6);   
      hncy3->Draw("hist");    
    }  
  if(plotvar == 5) // Cluster maps per GEM layer
    { 
      TCanvas *c4 = new TCanvas("c4","c4");  
      c4->Divide(3,2);  
      c4->cd(1); 
      h1->Draw("col2z"); 
      c4->cd(2); 
      h2->Draw("col2z"); 
      c4->cd(3); 
      h3->Draw("col2z"); 
      c4->cd(4);  
      h1c->Draw("col2z");  
      c4->cd(5);  
      h2c->Draw("col2z");  
      c4->cd(6);  
      h3c->Draw("col2z");  
    } 
  if(plotvar == 6)  // x-y position scatterplot of tracks extrapolated to z-position of quartic  
    { 
      TCanvas *c6 = new TCanvas("c6","c6"); 
      c6->cd(1);
      hxyquartic1->Draw("col2z"); 
    } 
  if(plotvar == 7) // Map of quartic bars hit in events passing tracking cuts
    {
      TCanvas *c7 = new TCanvas("c7","c7");
      c7->cd(1);
      hquarticocc->Draw("col2z");
    }

  cout << "Events passing all tracking cuts: " << npassing << endl; 

}
