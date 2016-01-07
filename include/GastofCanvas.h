#ifndef GastofCanvas_h
#define GastofCanvas_h

#include "TCanvas.h"
#include "TH1.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TStyle.h"

/**
 * \author Laurent Forthomme <laurent.forthomme@cern.ch>
 * \date 25 Jul 2015
 */
class GastofCanvas : public TCanvas
{
  public:
    enum PlotType { TB, Simulation };
    inline GastofCanvas() :
      TCanvas("null"), fBuilt(false), fType(TB),
      fLegend(0), fLegendX(.48), fLegendY(.76), fLegendNumEntries(0),
      fUpperLabel(0), fLabelsDrawn(false) {;}
    inline GastofCanvas(TString name, unsigned int width=500, unsigned int height=500, TString upper_label="") :
      TCanvas(name, "", width, height), fBuilt(false), fWidth(width), fHeight(height), fType(TB),
      fLegend(0), fLegendX(.48), fLegendY(.76), fLegendNumEntries(0),
      fUpperLabelText(upper_label), fUpperLabel(0), fLabelsDrawn(false) { Build(); }
    inline GastofCanvas(TString name, TString upper_label) :
      TCanvas(name, "", 500, 500), fBuilt(false), fWidth(500), fHeight(500), fType(TB),
      fLegend(0), fLegendX(.48), fLegendY(.76), fLegendNumEntries(0),
      fUpperLabelText(upper_label), fUpperLabel(0), fLabelsDrawn(false) { Build(); }
    inline virtual ~GastofCanvas() {
      if (fLegend) delete fLegend;
      if (fUpperLabel) delete fUpperLabel;
    }

    inline void SetPlotType(const PlotType& pt) { fType = pt; }
    inline void SetLegendX(double x) { fLegendX = x; }

    inline void SetUpperLabel(TString text="") {
      if (text=="") return;
      fUpperLabelText = text;
      if (fUpperLabel) return;
      fUpperLabel = new TPaveText(.5, .922, .885, .952, "ndc");
      fUpperLabel->SetMargin(0.);
      fUpperLabel->SetFillColor(kWhite);
      fUpperLabel->SetLineColor(kWhite);
      fUpperLabel->SetLineWidth(0);
      fUpperLabel->SetShadowColor(kWhite);
      fUpperLabel->SetTextFont(43);
      fUpperLabel->SetTextAlign(33);
      fUpperLabel->SetTextSize(18);
      fUpperLabel->AddText(fUpperLabelText);
    }

    inline void Save(TString ext="png", TString path=".") {
      bool valid_ext = true;
      valid_ext |= (strcmp(ext.Data(), "png")!=0);
      valid_ext |= (strcmp(ext.Data(), "pdf")!=0);
      valid_ext |= (strcmp(ext.Data(), "root")!=0);
      if (!valid_ext) return;
      DrawLabels();
      //printf("File saved as %s\n", Form("%s/%s.%s", path.Data(), TCanvas::GetName(), ext.Data()));
      TCanvas::SaveAs(Form("%s/%s.%s", path.Data(), TCanvas::GetName(), ext.Data()));
      //c1->SetLogz();
      //TCanvas::SaveAs(Form("%s/%s_logscale.%s", path.Data(), TCanvas::GetName(), ext.Data()));
    }
    inline TPad* Pad() { return c1; }

    void AddLegendEntry(const TObject* obj_, TString label_, Option_t* option_="lf") {
      if (!fLegend) {
        fLegend = new TLegend(fLegendX, fLegendY, fLegendX+.3, fLegendY+.11);
        fLegend->SetFillColor(kWhite);
        fLegend->SetLineColor(kWhite);
        fLegend->SetLineWidth(0);
        fLegend->SetTextFont(43);
        fLegend->SetTextSize(18);
      }
      fLegend->AddEntry(obj_, label_, option_);
      fLegendNumEntries++;
      if (fLegendNumEntries>3) fLegend->SetY1(fLegend->GetY1()-(fLegendNumEntries-3)*0.01);
    }

    void Prettify(TH1* h) {
      h->GetXaxis()->SetTitleSize(28);
      h->GetXaxis()->SetTitleFont(43);
      h->GetXaxis()->SetTitleOffset(0.75);
      h->GetYaxis()->SetTitleOffset(0.9);
      h->GetYaxis()->SetTitleSize(28);
      h->GetYaxis()->SetTitleFont(43);
      h->GetXaxis()->SetLabelSize(18);
      h->GetXaxis()->SetLabelFont(43);
      h->GetYaxis()->SetLabelSize(18);
      h->GetYaxis()->SetLabelFont(43);
      h->SetLineColor(kBlack);
      h->SetLineWidth(2);
    }

  protected:
    inline void Build() {
      if (fBuilt) return;
      DrawGrid();
      fBuilt = true;
    }
    
  private:
    inline void DrawLabels() {
      if (fLabelsDrawn) {
        fLabel1->Draw("same");
        fLabel2->Draw("same");
        return;
      }
      fLabel1 = new TPaveText(.112, .925, .2, .955, "ndc");
      fLabel1->AddText("GasToF");
      fLabel1->SetMargin(0.);
      fLabel1->SetFillColor(kWhite);
      fLabel1->SetLineColor(kWhite);
      fLabel1->SetLineWidth(0);
      fLabel1->SetShadowColor(kWhite);
      fLabel1->SetTextFont(63);
      fLabel1->SetTextAlign(13);
      fLabel1->SetTextSize(22);
      fLabel1->Draw("same");
      
      fLabel2 = new TPaveText(.28, .925, .36, .955, "NDC");
      fLabel2->AddText(GetTypeLabel(fType));
      fLabel2->SetMargin(0.);
      fLabel2->SetFillColor(kWhite);
      fLabel2->SetLineColor(kWhite);
      fLabel2->SetLineWidth(0);
      fLabel2->SetShadowColor(kWhite);
      fLabel2->SetTextFont(43);
      fLabel2->SetTextAlign(13);
      fLabel2->SetTextSize(22);
      fLabel2->Draw("same");
      
      if (fLegend and fLegend->GetNRows()!=0) fLegend->Draw();
      SetUpperLabel(fUpperLabelText);
      if (fUpperLabel) fUpperLabel->Draw("same");
      fLabelsDrawn = true;

      gStyle->SetMarkerStyle(20);
      gStyle->SetMarkerSize(.87);
      gStyle->SetTitleFont(43, "XYZ");
      gStyle->SetTitleSize(24, "XYZ");
      //gStyle->SetTitleOffset(2., "Y");
      gStyle->SetLabelFont(43, "XYZ");
      gStyle->SetLabelSize(20, "XY");
      gStyle->SetLabelSize(15, "Z");
      gStyle->SetTitleOffset(0.9, "X");
      gStyle->SetTitleOffset(1.1, "Y");
      gStyle->SetHistLineColor(kBlack);
      gStyle->SetHistLineWidth(2);          
    }
    inline void DrawGrid() {
      TCanvas::cd();
      gStyle->SetOptStat(0);

      TCanvas::Divide(1,2);
      c1 = (TPad*)TCanvas::GetPad(1);
      c2 = (TPad*)TCanvas::GetPad(2);
      c1->SetPad(0.,0.,1.,1.);
      c2->SetPad(0.,0.,1.,0.);
      c1->SetBottomMargin(0.1);
      c1->SetLeftMargin(0.12);
      c1->SetRightMargin(0.11);
      c1->SetTopMargin(0.1);
      TCanvas::cd(1);
     
      c1->SetTicks(1, 1);
      //c1->SetGrid(1, 1);
    }

    inline TString GetTypeLabel(const PlotType& pt) {
      switch (pt) {
        case TB: return "TB2015";
        case Simulation: return "simulation";
        default: return "";
      }
    }

    bool fBuilt;
    TPad *c1, *c2;
    double fWidth, fHeight;
    PlotType fType;
    TLegend *fLegend;
    double fLegendX, fLegendY;
    unsigned int fLegendNumEntries;
    TPaveText *fLabel1, *fLabel2;
    TString fUpperLabelText;
    TPaveText *fUpperLabel;
    bool fLabelsDrawn;
};

#endif
