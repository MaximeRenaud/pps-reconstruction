#ifndef GastofGrid_h
#define GastofGrid_h

#include "GastofCanvas.h"
#include "GastofCommons.h"
#include "TH2.h"
#include "TDatime.h"

/**
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \date 25 Jul 2015
 */
class GastofGrid : public GastofCanvas
{
  public:
    inline GastofGrid() :
      GastofCanvas(), fDrawOptions("colz"), fRunId(0), fRunDate(TDatime().AsString()) {;}
    inline GastofGrid(TString name, unsigned int width=500, unsigned int height=500, TString upper_label="") :
      GastofCanvas(name, width, height, upper_label), fDrawOptions("colz"),
      fRunId(0), fRunDate(TDatime().AsString()) { Build(); }
    inline GastofGrid(TString name, TString upper_label) :
      GastofCanvas(name, upper_label), fDrawOptions("colz"),
      fRunId(0), fRunDate(TDatime().AsString()) { Build(); }
    inline virtual ~GastofGrid() {
      if (fHist) delete fHist;
    }

    inline void SetDrawOptions(Option_t* opt) { fDrawOptions = opt; }

    inline void SetLogscale(bool set=true) {
      GastofCanvas::Pad()->SetLogz(set);
      if (set) { if (strstr(GetName(), "_logz")!=0) GastofCanvas::SetName(TString(GetName())+"_logz"); }
      else     { if (strstr(GetName(), "_logz")==0) GastofCanvas::SetName(TString(GetName()).ReplaceAll("_logz", "")); }
    }

    inline void FillChannel(unsigned short nino_id, unsigned short channel_id, double content) {
      FillChannel(GetCoordinates(nino_id, channel_id), content);
    }
    inline void FillChannel(const Coord& c, double content) {
      if (!fHist) return;
      fHist->Fill(c.x, c.y, content);
    }
    inline TH2D* Grid() { return fHist; }

    inline void Save(TString ext="png", TString path=".") {
      //if (!fLabelsDrawn) {
      /*  fLabel3 = new TPaveText(.5, .0, .98, .05, "NDC");
        fLabel3->AddText(Form("Run %d - Spill %d - %s", fRunId, fSpillId, fRunDate.Data()));
        fLabel3->SetMargin(0.);
        fLabel3->SetFillColor(kWhite);
        fLabel3->SetLineColor(kWhite);
        fLabel3->SetLineWidth(0);
        fLabel3->SetShadowColor(kWhite);
        fLabel3->SetTextFont(43);
        fLabel3->SetTextAlign(32);
        fLabel3->SetTextSize(16);
        fLabel3->Draw();*/
      DrawGrid();
      if (!fLabel4) {
        fLabel4 = new TPaveText(.5, .0, .51, .05, "NDC");
        fLabel4->AddText("#downarrow beam #downarrow");
        fLabel4->SetMargin(0.);
        fLabel4->SetFillColor(kWhite);
        fLabel4->SetLineColor(kWhite);
        fLabel4->SetLineWidth(0);
        fLabel4->SetShadowColor(kWhite);
        fLabel4->SetTextFont(43);
        fLabel4->SetTextAlign(22);
        fLabel4->SetTextSize(18);
        fLabel4->Draw("same");
      }
      GastofCanvas::Save(ext, path);
    }

  private:
    inline void Build() {
      fHist = new TH2D(Form("hist_%s", TCanvas::GetName()), "", 8, 1., 9., 8, 1., 9.);
    }
    inline void DrawGrid() {
      GastofCanvas::Pad()->cd();
      GastofCanvas::Pad()->SetRightMargin(0.16);
      fHist->Draw(fDrawOptions);
      fHist->GetXaxis()->SetNdivisions(110);
      fHist->GetXaxis()->CenterLabels();
      fHist->GetYaxis()->SetNdivisions(110);
      fHist->GetYaxis()->CenterLabels();

      //fHist->SetMarkerStyle(20);
      //fHist->SetMarkerSize(.87);
      fHist->SetTitleFont(43, "XYZ");
      fHist->SetTitleSize(22, "XYZ");
      fHist->SetLabelFont(43, "XYZ");
      //fHist->SetTitleOffset(2., "Y");
      fHist->SetLabelSize(24, "XY");
      fHist->SetLabelSize(18, "Z");
      fHist->SetTitleOffset(1.3, "Y");
    }

    TH2D* fHist;
    TPaveText *fLabel3, *fLabel4;
    Option_t* fDrawOptions;
    unsigned int fBoardId, fRunId, fSpillId;
    TString fRunDate;
};

#endif
