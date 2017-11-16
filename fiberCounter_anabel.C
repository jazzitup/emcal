//Create a 2-D histogram from an image.
//Author: Olivier Couet

#include <vector>
#include "commonUtility.h"
#include <TASImage.h>
#include <TEllipse.h>
double zeroInt = 1;

double getWgtCenter( vector<int>& x, vector<double>& e) ;
//void findNeighbors(TH2F* hInd, vector<int>& px, vector<int>& py);
double getSum( vector<double>& e);

void fiberCounter(int num=1)
{
  short seedThr = 70;
  short bkgThr = 130;
  const int searchRange = 50;
  //  float rad = 10;
  
  
  const int maxX = 5000;
  const int maxY = 5000;
  gStyle->SetPalette(52);
   
  //TString infName = "test.png";
  TString infName = "diffuser_3/DBN_60-BL_22-N.JPG";
  //TString infName = Form("inputPics/anablesPics/Oct3_%d.JPG",num);
  //TString infName = "inputPics/anablesPics/both_ends/DBN_61-BL_22-NG.JPG";
  TASImage image(infName);
  //TASImage image("/Users/yongsunkim/uiucAnalysis/emcal/inputPics/anablesPics/100_0183_trimmed_small-1.JPG");
  
  UInt_t yPixels = image.GetHeight();
  UInt_t xPixels = image.GetWidth();
  cout << " xPixels = " << xPixels<<endl;
  cout << " yPixels = " << yPixels<<endl;
  
  int xc,yc;
  xc = xPixels/2;
  yc = yPixels/2;
  
  if ( (maxX<xPixels) || ( maxY <yPixels) )   {
    cout <<"  (maxX<xPixels) || ( maxY <yPixels) ! " << endl;
    return;
  }
  
  TH1D* hNfib = new TH1D("hNfib",";number of clusters;Entries",4000,0,4000);

  UInt_t *argb   = image.GetArgbArray();
  
  TH2D* hOriginal = new TH2D("hOriginal","",xPixels,.5,xPixels+1,yPixels+1,.5,yPixels);
  TH2D* h = hOriginal;

  TH1D* h1d = new TH1D("h1d","Intensity of each pixel",256,0,256);
  TH1D* henergy = new TH1D(Form("henergy_%d",num),";energy of clusters",300,0,30000);
  TH1D* henergyNorm = new TH1D(Form("henergyNorm_%d",num),";Self-normalized intensity",300,0,3);
  TH1D* hMeanE = (TH1D*)henergy->Clone("hMeanE");
  hMeanE->SetXTitle("<cluster intensity>");

  TH1D* hRMSE = new TH1D("hRMSE",";RMS;Entries",1000,0,100000);
  TH1D* hRMSNorm = new TH1D("hRMSNorm",";RMS normalized by mean;Entries",1000,0,1);

  TH2D* gFibers = new TH2D("nFibers",";X;Y",20,0,xPixels, 20, 0, yPixels);
  TH2D* gInten = (TH2D*)gFibers->Clone("gInten");
  TH2D* densityInten   = (TH2D*)gFibers->Clone("densityInten");
  TH2D* gNpix  = (TH2D*)gFibers->Clone("gNpix");
  TH2D* gensityNpix   = (TH2D*)gFibers->Clone("densityNpix");

  TH1D* hx = new TH1D("hx","X position",12,0,xPixels);
  TH1D* hy = new TH1D("hy","Y position",14,0,yPixels);
  
  

  h1d->Sumw2();
  for (int row=0; row<xPixels; ++row) {
    for (int col=0; col<yPixels; ++col) {
      int index = col*xPixels+row;
      short grey = argb[index]&0xff ;
      h->SetBinContent(row+1,yPixels-col,grey);
      h1d->Fill(grey);
    }
  }
  
  TCanvas* c0 = new TCanvas("c0","",450,450);
  //c0->Divide(2,1);
  //c0->cd(1);
  h1d->SetAxisRange(0.5,2e6,"Y");
  h1d->Draw();
  jumSun(seedThr,0.5,seedThr,2e6);
  jumSun(bkgThr,0.5,bkgThr,2e6);
  gPad->SetLogy();
  //c0->cd(2);
  //h->SetAxisRange(0,256,"z");
  //h->Draw("colz");
  gPad->SetRightMargin(0.2);
  
  c0->SaveAs("histo-1d.gif");

  
  TH2D* h3 = (TH2D*)h->Clone("h3");
  
  TH2D* h4 = (TH2D*)h->Clone("h4");  
  for ( int ix0=1 ; ix0<=h4->GetNbinsX() ; ix0++) {
    for ( int iy0=1 ; iy0<=h4->GetNbinsY() ; iy0++) {
      h4->SetBinContent(ix0,iy0,0.0001);
    }}
  
  std::vector<int> px;
  std::vector<int> py;
  std::vector<int> pinten;
  
  std::vector<int> vClstE;  
  std::vector<float> vClstX;  
  std::vector<float> vClstY;  


  //  std::vector<short> inten;
     

  
  // Clean up backgrounds ;
  short arrInten[maxX+1][maxY+1];
  short arrFlag[maxX+1][maxY+1];
  int kIn=1;  int kOut=-1;   int kUnde=0;
  int nXbins = h3->GetNbinsX() ; 
  int nYbins = h3->GetNbinsY() ; 

  for ( int ix0=1 ; ix0<=nXbins ; ix0++) {
    for ( int iy0=1 ; iy0<=nYbins ; iy0++) {
      short val0 = h3->GetBinContent(ix0,iy0);
      
      arrInten[ix0][iy0] = val0;
      arrFlag[ix0][iy0] = kUnde;	
      if (val0 < bkgThr )  {
	h3->SetBinContent(ix0, iy0, zeroInt);
	arrInten[ix0][iy0] = 0;
	arrFlag[ix0][iy0] =  kOut;
      }
    }
  }
  
  int nClst = 0;
  
  double Itot=0.,I1=0.,I2=0.,I3=0.,I4=0.;
  double Itotsq=0.,I1sq=0.,I2sq=0.,I3sq=0.,I4sq=0.;
  
  int nClst1=0,nClst2=0,nClst3=0,nClst4=0;
  
  for ( int ix0=1 ; ix0<=nXbins ; ix0++) {
    for ( int iy0=1 ; iy0<=nYbins ; iy0++) {
      short val0 = arrInten[ix0][iy0];
      
      if ( val0 < seedThr ) 
	continue;
      if ( arrFlag[ix0][iy0] != kUnde ) // if it is already included in other clusters.
	continue;

      //cout << " found new seed!   " ;
      //cout <<"(x,y,intensity) = " <<  ix0<<", "<<iy0<<", "<<val0<<endl;


      // Found a seed!     (ix0, iy0, val0) are the seed! 
      nClst++;
      
      // Begin Clustering 
      //      inten.clear();
      px.clear();
      py.clear();
      pinten.clear();

      //      inten.push_back(val0);
      px.push_back(ix0);   
      py.push_back(iy0);
      pinten.push_back( val0 ) ;

      arrFlag[ix0][iy0] = kIn;
      
      bool completeFlag = false;
      int nIter = 0;
      int countedHits=0;
      int itrBegin=0;  // <== this number will evolve over the iteration :  itrBegin = itrEnd;
      int itrEnd;

      while ( completeFlag==false)  
	{
	  //	  cout << " px.size() = " << px.size() << endl;
	  nIter++;
	  completeFlag = true;
	  
	  itrEnd = int(px.size());
	  //	  cout << " itrBegin = " << itrBegin ;
	  //	  cout << ", itrEnd = " << itrEnd<<endl;
	  for ( int vi = itrBegin ; vi < itrEnd ; vi++) {
	    
	    int xCand = px[vi];     int yCand = py[vi];      short intenCand = pinten[vi];
	    
	    if ( (xCand < nXbins) && (arrFlag[xCand+1][yCand]==kUnde) ) {	       // right end 
	      completeFlag=false; 
	      countedHits++;
	      px.push_back( xCand+1 ) ;
	      py.push_back( yCand   ) ;
	      pinten.push_back( arrInten[xCand+1][yCand]   ) ;
	      arrFlag[xCand+1][yCand]=kIn; 
	    }
	    
	    if ( (xCand > 1    ) && (arrFlag[xCand-1][yCand]==kUnde) ) {	       // Left end 
	      completeFlag=false; 
	      countedHits++;
	      px.push_back( xCand-1 ) ;
	      py.push_back( yCand   ) ;
	      pinten.push_back( arrInten[xCand-1][yCand]   ) ;
	      arrFlag[xCand-1][yCand]=kIn; 
	    }
	    if ( (yCand < nYbins) && (arrFlag[xCand][yCand+1]==kUnde) ) { 	       // top end 
	      completeFlag=false; 
	      countedHits++;
	      px.push_back( xCand ) ;
	      py.push_back( yCand+1   ) ;
	      pinten.push_back( arrInten[xCand][yCand+1]   ) ;
	      arrFlag[xCand][yCand+1]=kIn; 
	    }
	    if ( (yCand > 1    ) && (arrFlag[xCand][yCand-1]==kUnde) ) {  	       // bottom end 
	      completeFlag=false; 
	      countedHits++;
	      px.push_back( xCand ) ;
	      py.push_back( yCand-1   ) ;
	      pinten.push_back( arrInten[xCand][yCand-1]  ) ;
	      arrFlag[xCand][yCand-1]=kIn; 
	    }
	  }
	  itrBegin = itrEnd;

	  /*	  if ( nClst==2 ) {
            TCanvas* ctemp = new TCanvas("ctemp","", 400,400);
            TH2D* htemp = new TH2D("htemp","",50, px[0]-15, px[0]+35, 50, py[0]-25, py[0]+25);
            for ( int ii=0 ; ii<px.size(); ii++)   {
              int theBin = htemp->FindBin( px[ii], py[ii]) ;
              htemp->SetBinContent( theBin, pinten[ii] );
            }
	    htemp->SetAxisRange(200,256,"Z");
	    htemp->Draw("colz");
	    drawText(Form("Iter. %d",nIter),0.2,0.8);
	    ctemp->SaveAs(Form("clustering_itr%d.png",nIter));
	      
	    }*/
	}
      //cout << "number of hits in "<<nClst<<"th cluster: " << px.size() << ",  nIteration = "<<nIter<<endl;  
      int sumEnergy = 0;
      for ( int vi=0 ; vi < px.size() ; vi++) {
		sumEnergy = sumEnergy + pinten[vi];
      }
      //cout << " total energy = " << sumEnergy << endl;

      float xmean=0;
      float ymean=0;
      int nPixels= px.size();
      for ( int vi=0 ; vi < nPixels ; vi++) { 
		xmean = xmean + px[vi]*pinten[vi] ;
		ymean = ymean + py[vi]*pinten[vi] ;
      }
      xmean = xmean / sumEnergy ;
      ymean = ymean / sumEnergy ;

      henergy->Fill(sumEnergy);
      vClstE.push_back(sumEnergy);
      vClstX.push_back(xmean);
      vClstY.push_back(ymean);
      
      hx->Fill(xmean);
      hy->Fill(ymean);
      
      if (xmean<xc && ymean>yc) {
		  I1=I1+sumEnergy;
		  I1sq=I1sq+sumEnergy*sumEnergy;
		  nClst1=nClst1+1;}
      if (xmean>=xc && ymean>yc) {
		  I2=I2+sumEnergy;
		  I2sq=I2sq+sumEnergy*sumEnergy;
		  nClst2=nClst2+1;}
      if (xmean<xc && ymean<=yc) {
		  I3=I3+sumEnergy;
		  I3sq=I3sq+sumEnergy*sumEnergy;
		  nClst3=nClst3+1;}
      if (xmean>=xc && ymean<=yc) {
		  I4=I4+sumEnergy;
		  I4sq=I4sq+sumEnergy*sumEnergy;
		  nClst4=nClst4+1;}
      
      Itot=Itot+sumEnergy;
      Itotsq=Itotsq+sumEnergy*sumEnergy;
      
      
      gFibers->Fill(xmean, ymean);
      gInten-> Fill(xmean, ymean,sumEnergy);
      gNpix-> Fill(xmean, ymean,nPixels);
      
      
      //      for ( int vi=0 ; vi < px.size() ; vi++) {
      //	h4->SetBinContent( px[vi],  py[vi], pinten[vi]);
      //      }
      //      TCanvas* c11 = new TCanvas("c11","",400,400);
      //      h4->Draw("colz");
      //      c11->SaveAs(Form("Cluster_%d.png",nClst));
      
    }}
  
  cout << "Picture = " << infName << "\n";
  cout << "seedThr = " << seedThr << "\n";
  cout << "bkgThr = " << bkgThr << "\n";
  
  
  cout<< "Itot = " << Itot << "\n";
  cout<< "I1+I2+I3+I4 = " << I1+I2+I3+I4 << "\n";
  cout<< "I1 = " << I1 << "\n";
  cout<< "I2 = " << I2 << "\n";
  cout<< "I3 = " << I3 << "\n";
  cout<< "I4 = " << I4 << "\n";
  
  cout<< "nClst = " << nClst << "\n";
  cout<< "Sum of nClst =" << nClst1+nClst2+nClst3+nClst4 << "\n";
  cout<< "nClst1 = " << nClst1 << "\n";
  cout<< "nClst2 = " << nClst2 << "\n";
  cout<< "nClst3 = " << nClst3 << "\n";
  cout<< "nClst4 = " << nClst4 << "\n";
  
  double devI= sqrt((Itotsq/nClst)-((Itot*Itot)/(nClst*nClst)));
  double devI1= sqrt((I1sq/nClst1)-((I1*I1)/(nClst1*nClst1)));
  double devI2= sqrt((I2sq/nClst2)-((I2*I2)/(nClst2*nClst2)));
  double devI3= sqrt((I3sq/nClst3)-((I3*I3)/(nClst3*nClst3)));
  double devI4= sqrt((I4sq/nClst4)-((I4*I4)/(nClst4*nClst4)));
  
  cout << "devI = " << devI << "\n";
  cout << "devI1 = " << devI1 << "\n";
  cout << "devI2 = " << devI2 << "\n";
  cout << "devI3 = " << devI3 << "\n";
  cout << "devI4 = " << devI4 << "\n";
  
  double FOMI=(I1/Itot-0.25)*(I1/Itot-0.25)+(I2/Itot-0.25)*(I2/Itot-0.25)+(I3/Itot-0.25)*(I3/Itot-0.25)+(I4/Itot-0.25)*(I4/Itot-0.25);
  cout << "FOM Intensity = " << FOMI << "\n";
  
  double n1= (double)nClst1/(double)nClst, n2= (double)nClst2/(double)nClst, n3= (double)nClst3/(double)nClst, n4= (double)nClst4/(double)nClst;
  
  double FOMN=(n1-0.25)*(n1-0.25)+(n2-0.25)*(n2-0.25)+(n3-0.25)*(n3-0.25)+(n4-0.25)*(n4-0.25);
  cout << "FOM Nfibers = " << FOMN << "\n";
  
  
  double max = gInten->GetMaximum();
  double norm = 1.0/max;
  //cout << "The maximum intensity is: " << max << "\n";
  
  gInten->Scale(norm);
  
  densityInten->Add(gInten);
  densityInten->Divide(gFibers);
  densityNpix->Add(gNpix);
  densityNpix->Divide(gFibers);
  
  TCanvas* c5 = new TCanvas("c5","",600,600);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1);
  gPad->SetRightMargin(0.2);
  //c5->Divide(2,2);
  //c5->cd(1);
  //cleverRangeZ(gInten);
  //gInten->SetMinimum(0);
  //gInten->SetXTitle("Light intensity");
  //gInten->Draw("colz");

  //c5->cd(2);
  //cleverRangeZ(gFibers);
  //gFibers->SetMinimum(0);
  //gFibers->SetMaximum(11.0);
  //gFibers->SetXTitle("Number of fibers");
  //gFibers->Draw("colz");

  //c5->cd(3);
  //cleverRangeZ(densityInten);
  densityInten->SetMinimum(0);
  densityInten->SetMaximum(0.2);
  densityInten->SetXTitle("Intensity/Number of fibers");
  densityInten->Draw("colz");

  //c5->cd(4);
  //cleverRangeZ(densityNpix);
  //densityNpix->SetMinimum(0);
  //densityNpix->SetMaximum(30.0);
  //densityNpix->SetXTitle("Number of Pixel per fiber");
  //densityNpix->Draw("colz");
  
  //cout << "total number of clusters = " << nClst << "\n";
  
  int ncor=0, n2s=0, n3s=0;
  double I2s = 1-2*devI/(henergy->GetMean());
  double I3s = 1-3*devI/(henergy->GetMean());
  
  for ( int ci = 0 ; ci< vClstE.size() ; ci++) {
	  double normint = vClstE[ci] / henergy->GetMean();
	  henergyNorm->Fill(normint) ;
    
    if ((1.7<=normint)&&(normint<2.5)) {
		ncor = ncor +1;}
		
    if ((2.7<=normint)&&(normint<3.5)) {
		ncor = ncor +2;}
		
    if ((3.7<=normint)&&(normint<4.5)) {
		ncor = ncor +3;}
	
	if (normint<I2s) n2s++;
	
	if (normint<I3s) n3s++;
	
	if (normint<0.75) n75++;
  }
  
  cout << "n correction = " << ncor << "\n";
  
  int ntotcor = nClst + ncor;

  cout << "Corrected n fibers = " << ntotcor << "\n";

  cout << "n2s = " << n2s << "\n";
  cout << "n3s = " << n3s << "\n";
  
  cout << "r2s (%) = " << (double)n2s*100/(double)ntotcor << "\n";
  cout << "r3s (%) = " << (double)n3s*100/(double)ntotcor << "\n";
  cout << "r75 (%) = " << (double)n75*100/(double)ntotcor << "\n";

  
  
  hNfib->Fill(nClst);
  hMeanE->Fill( henergy->GetMean() );
  hRMSE->Fill( henergy->GetRMS() );
  hRMSNorm->Fill(  henergy->GetRMS() / henergy->GetMean() ) ;
  
  TCanvas* c15 = new TCanvas("c15","",600,600);
  //  gStyle->SetPalette(52);
  h->SetAxisRange(0,256,"z");
  h->Draw("colz");
  gPad->SetRightMargin(0.2);
  for ( int ci = 0 ; ci< vClstX.size() ; ci++) {
    TEllipse *el3 = new TEllipse( vClstX[ci], vClstY[ci], 5);
    el3->SetLineColor(kBlack);
    el3->SetLineWidth(1.5);
    el3->SetFillStyle(0);
    el3->Draw();
  }

  //double max2 = henergy->GetMaximum();
  //double norm2 = 1.0/max2;
  //henergy->Scale(norm);

  // end of clustering 
  TCanvas* c2 = new TCanvas("c2","",400,400);
  handsomeTH1(henergy,1);
  henergyNorm->SetXTitle("Normalized intensity of each cluster");
  henergyNorm->SetYTitle("Entries");
  henergyNorm->Draw();

  //TCanvas* c6 = new TCanvas("c6","",900,500);
  //c6->Divide(2,1);
  //c6->cd(1);
  //hx->Draw();
  //c6->cd(2);
  //hy->Draw();
  
  TFile* fout = new TFile(Form("%s_outputHistograms.root",infName.Data()),"RECREATE");
  henergy->Write();
  henergyNorm->Write();
  hNfib->Write();
  hMeanE->Write();
  hRMSE->Write();
  hRMSNorm->Write();
  
  //gFibers->Write();
  //gInten->Write();
  //densityInten->Write();
  //gNpix->Write();
  //gensityNpix->Write();
  
  fout->Close();
  
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
