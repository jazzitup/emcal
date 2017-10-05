//Create a 2-D histogram from an image.
//Author: Olivier Couet

#include "commonUtility.h"
#include <TASImage.h>
double zeroInt = 1;

TH2D* rebinhist(TH2D* h, int nRebin);


void checkUniformity()
{
  gStyle->SetPalette(52);
  
  //  TASImage image("inputPics/kodak2_sept11/block2.jpg");
  TString infName = "inputPics/anablesPics/Oct3_2.JPG";
  TASImage image(infName);
  //  TASImage image("/Users/yongsunkim/uiucAnalysis/emcal/inputPics/anablesPics/100_0183_trimmed_small-1.JPG");
  
  UInt_t yPixels = image.GetHeight();
  UInt_t xPixels = image.GetWidth();
  cout << " xPixels = " << xPixels<<endl;
  cout << " yPixels = " << yPixels<<endl;
  
  TH1D* hNfib = new TH1D("hNfib",";number of clusters;Entries",4000,0,4000);


  UInt_t *argb   = image.GetArgbArray();
  
  TH2D* hOriginal = new TH2D("hOriginal","",xPixels,.5,xPixels+1,yPixels+1,.5,yPixels);
  TH2D* h = hOriginal;

  TH1D* h1d = new TH1D("h1d","1D histogram",256,0,256);
  TH1D* henergy = new TH1D("henergy",";energy of clusters",300,0,300000);
  TH1D* henergyNorm = new TH1D("henergyNorm",";Self-normalized intensity",300,0,3);
  TH1D* hMeanE = (TH1D*)henergy->Clone("hMeanE");
  hMeanE->SetXTitle("<cluster intensity>");

  TH1D* hRMSE = new TH1D("hRMSE",";RMS;Entries",1000,0,100000);
  TH1D* hRMSNorm = new TH1D("hRMSNorm",";RMS normalized by mean;Entries",1000,0,1);
  
  h1d->Sumw2();
  for (int row=0; row<xPixels; ++row) {
    for (int col=0; col<yPixels; ++col) {
      int index = col*xPixels+row;
      short grey = argb[index]&0xff ;
      h->SetBinContent(row+1,yPixels-col,grey);
      h1d->Fill(grey);
    }
  }
  
  TCanvas* c0 = new TCanvas("c0","",800,400);
  c0->Divide(2,1);
  c0->cd(1);
  h1d->SetAxisRange(0.5,2e6,"Y");
  h1d->Draw();
  gPad->SetLogy();
  c0->cd(2);
  h->SetAxisRange(0,256,"z");
  h->Draw("colz");
  gPad->SetRightMargin(0.2);
  
  c0->SaveAs("histo-1d.gif");
  
  TH2D* reBinnedH = rebinhist(h,50);
  TCanvas* c1 = new TCanvas("c1","",800,800);
  TH1D* hx = (TH1D*)reBinnedH->ProjectionX();
  TH1D* hy = (TH1D*)reBinnedH->ProjectionY();
  c1->Divide(2,2);
  c1->cd(1);
  reBinnedH->Draw("LEGO");
  c1->cd(3);
  hx->Draw();
  drawText("X projection", 0.5,0.5);
  c1->cd(4);
  hy->Draw();
  drawText("Y projection", 0.5,0.5);

}


TH2D* rebinhist(TH2D* h, int nRebin) {
  int nx = h->GetNbinsX() ;
  int ny = h->GetNbinsY() ;
  int nxNew = nx/nRebin;
  int nyNew = ny/nRebin;
  
  
  TH2D* ret = new TH2D(Form("%s_rebinned",h->GetName()),"",nxNew,0,nxNew,nyNew,0,nyNew);
  for ( int ix = 1 ; ix<= nx ; ix++) {
    for ( int iy = 1 ; iy<= ny ; iy++) {
      double val = h->GetBinContent(ix,iy);
      ret->Fill( ix/nRebin, iy/nRebin, val);
    }}
  return ret;
}
