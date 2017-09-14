//Create a 2-D histogram from an image.
//Author: Olivier Couet

#include <vector>
#include "commonUtility.h"
#include <TASImage.h>
#include <TEllipse.h>
double zeroInt = -1;

bool isDone ( TH2D* hh=0) ;
double getWgtCenter( vector<int>& x, vector<double>& e) ;
void findNeighbors(TH2F* hInd, vector<int>& px, vector<int>& py);
double getSum( vector<double>& e) ;

void fiberCounter()
{
  double seedThr = 0.9;
  double bkgThr = 0.5;
  const int searchRange = 50;
  float rad = 6;

  gStyle->SetPalette(52);


  //  TASImage image("inputPics/picTakenAtNPL_piece1.png");
  //  TASImage image("inputPics/picTakenAtNPL.png");
  //  TASImage image("inputPics/exmample_defacts.png");
  //  TASImage image("inputPics/exmample_goodClusters.png");
  //  TASImage image("inputPics/pic_aug22_samll.png");
  //  TASImage image("inputPics/pic_aug22.png");
  TASImage image("inputPics/kodak1/timmed_100_0015.png");
  //  TASImage image("/Users/yongsunkim/uiucAnalysis/emcal/inputPics/kodak2_sept11/block1.jpg");
  // TASImage image("/Users/yongsunkim/uiucAnalysis/emcal/inputPics/fibersTungstenHindered.png");
  //  TASImage image("/Users/yongsunkim/uiucAnalysis/emcal/inputPics/kodak2_sept11/black1-lowRes2.jpg");

   UInt_t yPixels = image.GetHeight();
   UInt_t xPixels = image.GetWidth();
   UInt_t *argb   = image.GetArgbArray();

   TH2D* h = new TH2D("h","",xPixels,.5,xPixels+1,yPixels+1,.5,yPixels);
   TH1D* h1d = new TH1D("h1d","1D histogram",256,0,1);
   TH1D* henergy = new TH1D("henergy",";energy of clusters",300,0,300);
   

   h1d->Sumw2();
   for (int row=0; row<xPixels; ++row) {
      for (int col=0; col<yPixels; ++col) {
	int index = col*xPixels+row;
	//	 cout << "argb[index] = " << argb[index] << endl;
         float grey = float(argb[index]&0xff)/256;
         h->SetBinContent(row+1,yPixels-col,grey);
	 h1d->Fill(grey);
      }
   }
   



   TCanvas* c0 = new TCanvas("c0","",400,400);
   //   h1d->SetAxisRange(0.5,1,"X");
   h1d->Draw();
   gPad->SetLogy();
   c0->SaveAs("histo-1d.gif");
 
   TCanvas* c1 = new TCanvas("c1","",400,400);
   h->SetAxisRange(0,1,"z");
   //   h->SetAxisRange(300,400,"X");    h->SetAxisRange(300,400,"Y");
   h->Draw("colz");
   gPad->SetRightMargin(0.2);
   


   TH2D* h3 = (TH2D*)h->Clone("h3");
   
   std::vector<int> px;
   std::vector<int> py;
   int kIn=1;  int kOut=-1;   int kUnde=0;
   std::vector<double> inten;
   
   double arrInten[2*searchRange+1,2*searchRange+1];
   int arrFlag[2*searchRange+1,2*searchRange+1];
      
   
   std::vector<int> wgtx;
   std::vector<int> wgty;
   std::vector<double> energy;

   
   
   
   // Clean up backgrounds ;
   for ( int ix0=1 ; ix0<=h3->GetNbinsX() ; ix0++) {
     for ( int iy0=1 ; iy0<=h3->GetNbinsY() ; iy0++) {

       double val0 = h3->GetBinContent(ix0,iy0);
       if (val0 < bkgThr ) 
	 h3->SetBinContent(ix0, iy0, zeroInt);
     }}
   
   
   int nClst = 0;
   
   TCanvas* c2 = new TCanvas("c2","",400,400);
   
   for ( int ix0=1 ; ix0<=h3->GetNbinsX() ; ix0++) {
     for ( int iy0=1 ; iy0<=h3->GetNbinsY() ; iy0++) {
       
       double val0 = h3->GetBinContent(ix0,iy0);

       if ( val0 < seedThr ) 
	 continue;

       nClst++;
       
       TH2F* hInd = new TH2F(Form("ind_%d",nClst),"", 
			     2*searchRange + 1, ix0-searchRange -0.5, ix0+searchRange +0.5,
			     2*searchRange + 1, iy0-searchRange -0.5, iy0+searchRange +0.5);
       // bin ranges from ix0-searchRange -> ix0+searchRange 
       int centerBin = hInd->FindBin(ix0,iy0);
       
       hInd->SetBinContent(centerBin, val0);
       h3->SetBinContent  (centerBin, zeroInt);
       
       
       // Begin clustering
       
       int ix0a =  searchRange;     // index for array.   ii=x0-searchRange => iia = 0 
       int iy0a =  searchRange;     // index for array.   ii=x0-searchRange => iia = 0 
       
       // Fill the array 
       for ( int ii = ix0-searchRange ; ii<= ix0+searchRange ; ii++) {
	 for ( int jj = iy0-searchRange ; jj<= iy0+searchRange ; jj++) {
	   
	   int iia = ii - (ix0-searchRange);    // index for array.   ii=x0-searchRange => iia = 0 
	   int jja = jj - (iy0-searchRange);    // index for array.   ii=x0-searchRange => iia = 0 
	   
	   // initial setup for the array
	   arrFlag[iia,jja] = -1 ;	   
	   arrInten[iia,jja] = 0 ;
	   //////////////////////////////

	   if ( (ii < 1)||(jj<1) ) 
	     continue;
	   if ( (ii > h3->GetNbinsX() ) || (jj > h3->GetNbinsY() ) )
	     continue;
	   if ( valij < bkgThr ) continue;    // Ignore backgrounds
	   
	   double valij = h3->GetBinContent(ii,jj);
	   arrInten[iia,jja,] = valij;
	   arrFlag[iia,jja] = 0 ;	   
	 }}
       
       // Begin Clustering 
       px.clear();
       py.clear();
       arrInten[ix0a,iy0a] = kIn; 
       px.push_back(ix0a); 
       py.push_back(iy0a);
       bool completeFlag = false;
       while ( completeFlag==false)  
	 {
	   completeFlag = true;
	   for ( int vi = 0 ; vi<= px.size() ; vi++) 
	     {
	       int xa = px[vi]  ;
	       int ya = py[vi]  ;
	       if ( xa < 2*searchRange  ) {  // right end
		 if (  (arrFlag[xa+1,ya]==kUnde) && (arrInten[xa+1,ya] > bkgThr) )  {
		   completeFlag=false; arrFlag[xa+1,ya]=kIn;   px.push_back(xa+1); py.push_back(ya); 
		 }}
	       if ( ya < 2*searchRange  ) {  // right end
		 if (  (arrFlag[xa,ya+1]==kUnde) && (arrInten[xa,ya+1] > bkgThr) )  {
		   completeFlag=false; arrFlag[xa,ya+1]=kIn;   px.push_back(xa);   py.push_back(ya+1);
		 }}
	       if ( xa > 0  ) {  // right end
		 if (  (arrFlag[xa-1,ya]==kUnde) && (arrInten[xa-1,ya] > bkgThr) )  {
		   completeFlag=false; arrFlag[xa-1,ya]=kIn;   px.push_back(xa-1); py.push_back(ya);
		 }}
	       if ( ya > 0  ) {  // right end
		 if (  (arrFlag[xa,ya-1]==kUnde) && (arrInten[xa,ya-1] > bkgThr) )  {
		   completeFlag=false; arrFlag[xa,ya-1]=kIn;   px.push_back(xa);   py.push_back(ya-1);
		 }}
	     } // end of px pixels
	 }
       // end of clustering 
       
       
       
       
       
       if ( nClst%100 ==0 ) {
       	 cout << "Iteration " << nClst << "size of the cluster = " << inten.size() << endl;
       }
       wgtx.push_back(getWgtCenter( px, inten) );
       wgty.push_back(getWgtCenter( py, inten) );
       energy.push_back( getSum(inten) );
       henergy->Fill( getSum(inten) ) ;

       //       cout << " x,y,totE = " << wgtx[nClst] << ", " <<  wgty[nClst] <<", " <<  energy[nClst] <<", " << endl;
       

       //       if (nClst == 5)        	 return;

       h3->SetAxisRange(0,2,"z");
       h3->Draw("colz");
       h3->SetXTitle(Form("Fiber counts : %d",nClst));

       //       c2->SaveAs(Form("figs/iter_%d.png",nClst));
       px.clear();   py.clear();   inten.clear();
       //       return;
     }} // end of first iter
   
   
   cout << "number of clusters = " << nClst << endl;
   
   
   TCanvas* c3 = new TCanvas("c3","",400,400);
   h->Draw("colz");
   cout << "n = " << wgtx.size() << endl;
   for ( int ii = 0 ; ii<wgtx.size() ; ii++) { 
     TEllipse *el3 = new TEllipse( wgtx[ii], wgty[ii], rad);
     el3->SetLineColor(kYellow);
     el3->SetLineWidth(2);
     el3->SetFillStyle(0);
     el3->Draw();
   }
   //   c3->Update();
   c3->SaveAs("result.gif");

   TCanvas* c4 = new TCanvas("c4","",400,400);
   handsomeTH1(henergy);
   henergy->Draw();
   c4->SaveAs("energyDist.gif");

   TFile* hout = new TFile("result.root","recreate");
   h->Write();
   h1d->Write();
   henergy->Write();
   hout->Close();
}

bool isDone ( TH2D* hh ) { 
  bool flag = true;
  for ( int ix=1 ; ix<=hh->GetNbinsX() ; ix++) {
    for ( int iy=1 ; iy<=hh->GetNbinsY() ; iy++) {
      if ( hh->GetBinContent(ix,iy) != zeroInt ) {
	flag = false; 
	break;
      }
    }
  }
  return flag;
}


void findNeighbors(TH2F* hInd, vector<int>& px, vector<int>& py) { 
  
  for ( 
  int centerBin = hInd->FindBin(ix0,iy0);  
  
}


double getWgtCenter( vector<int>& x, vector<double>& e) { 
  double sum=0;
  double totEnergy = 0 ;
  for ( int ii = 0 ; ii < x.size() ; ii++) {
    sum = sum + x[ii]*e[ii];
    totEnergy = totEnergy + e[ii] ;
  }
  return sum/totEnergy;
}

double getSum( vector<double>& e) { 
  double sum = 0;
  for ( int ii = 0 ; ii < e.size() ; ii++) {
    sum = sum + e[ii];
  }
  
  return sum;
}
