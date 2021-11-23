//plots probability density distributions of fs, a, tstart, tdur for several
   //SMvals.dat files (several distances)
//Aryil Bechtel 2021

void plotFitParams(){

  const Int_t numfiles=2;
  TString dists[numfiles]={"Violet","Lblue"};

  Int_t linecols[numfiles]={2,4};
 
  TString indir="SMvals/SMvals_wc560kt30prct/Violet_vs_Lblue_5_2kpc/";
  TString numruns="100";
  TString modeltype="_m3_100to490ms_10to100dur_wc560kt30prct";

  Double_t alpha=0.1; //threshold pfi value
  Double_t SM_start=-20;
  Double_t SM_end=400;
  Double_t binsize=.05;
  Int_t nbins=(SM_end-SM_start)/binsize;

  TGraph **ROCgraphs=new TGraph*[numfiles];
  TLegend* leg1 = new TLegend(0.7,0.7,.87,.87);

  TH1D **fshists=new TH1D*[numfiles];
  TLegend* leg2=new TLegend(0.6,0.7,0.8,0.75);
  Int_t fbins=20;
  Double_t f_start=57;
  Double_t f_end=228;

  TH1D **ahists=new TH1D*[numfiles];
  TLegend* leg3=new TLegend(0.6,0.7,0.8,0.75);
  Int_t abins=50;
  Double_t a_start=0;
  Double_t a_end=0.6;

  TH1D **tdurHists = new TH1D*[numfiles];
  TLegend* leg4 = new TLegend(0.6,0.7,0.8,0.75);
  Int_t tdurBins = 50;
  Double_t tdurStart = 0;
  Double_t tdurEnd = 100;

  TH1D **tstartHists = new TH1D*[numfiles];
  TLegend* leg5 = new TLegend(0.6,0.7,0.8,0.75);
  Int_t tstartBins = 50;
  Double_t tstartStart = 0;
  Double_t tstartEnd = 545;

  TH1D *lnSMgivenShist=new TH1D("lnSMgivenShist"," ;lnSM;prob",nbins,SM_start,SM_end); //histogram for lnSMgivenS vals
  TH1D *lnSMgivenNShist=new TH1D("lnSMgivenNShist","lnSM",nbins,SM_start,SM_end); //histogram for lnSMgivenNS vals

  Int_t k;
  for(k=0;k<numfiles;k++){

    lnSMgivenShist->Reset();
    lnSMgivenNShist->Reset();

    TString infilename=indir+"SMvals_"+numruns+modeltype+"_"+dists[k]+".dat";
    TString infilename2 = indir+"tdur_"+numruns+modeltype+"_"+dists[k]+".dat";
    TString infilename3 = indir+"tstart_"+numruns+modeltype+"_"+dists[k]+".dat";
    cout<<"opened "<<infilename<<endl;

    ifstream in;
    ifstream in2;
    ifstream in3;
    in.open(infilename);
    in2.open(infilename2);
    in3.open(infilename3);
    Int_t i=0;
    Int_t maxpoints=10000;
    Double_t lnSMgivenS[maxpoints];
    Double_t lnSMgivenNS[maxpoints];

    //params given SASI
    Double_t fsbest[maxpoints];
    Double_t abest[maxpoints];
    Double_t tstartBest[maxpoints];
    Double_t tdurBest[maxpoints];

    //params given no-SASI
    Double_t fsbest2[maxpoints];
    Double_t abest2[maxpoints];
    Double_t tstartBest2[maxpoints];
    Double_t tdurBest2[maxpoints];
  
    while(1){
      in>>lnSMgivenS[i]>>lnSMgivenNS[i]>>fsbest[i]>>abest[i]>>fsbest2[i]>>abest2[i];
      in2>>tdurBest[i]>>tdurBest2[i];
      in3>>tstartBest[i]>>tstartBest2[i];
      
      if(!in.good()) break;
      i++;
    }
    Int_t numpoints=i;
    in.close();
    in2.close();
    in3.close();

    for(i=0;i<numpoints;i++){
      lnSMgivenShist->Fill(lnSMgivenS[i]);
      lnSMgivenNShist->Fill(lnSMgivenNS[i]);
    }

    //normalize the histograms to give probability
    Double_t norm=1;
    Double_t scale = norm/lnSMgivenShist->Integral();
    lnSMgivenShist->Scale(scale);
    scale = norm/lnSMgivenNShist->Integral();
    lnSMgivenNShist->Scale(scale);

    cout<<"lnSMgivenShist integral "<<lnSMgivenShist->Integral()<<endl;
    cout<<"lnSMgivenNShist integral "<<lnSMgivenNShist->Integral()<<endl;

    //calculate prob(false ID) and prob(detection)
    Double_t pfi[nbins];
    Double_t pd[nbins];
    Double_t thresh;
    Int_t ithresh;
    for(i=0;i<nbins;i++){
      pfi[i]=lnSMgivenNShist->Integral(i,nbins);
      pd[i]=lnSMgivenShist->Integral(i,nbins);
      //cout<<"pfi "<<pfi[i]<<endl;
      if(pfi[i]>alpha){
        thresh=lnSMgivenShist->GetXaxis()->GetBinCenter(i+1);
        ithresh=i;
      }
    }
    cout<<"final thresh "<<thresh<<endl;
    cout<<"pfi above thresh "<<lnSMgivenNShist->Integral(ithresh,nbins)<<endl;

  
    ROCgraphs[k]=new TGraph(nbins,pfi,pd);
    ROCgraphs[k]->SetLineColor(linecols[k]);
    ROCgraphs[k]->SetLineWidth(3);
    leg1->AddEntry(ROCgraphs[k],dists[k]+" kpc","l");

    fshists[k]=new TH1D("fs"+dists[k],"fshist",fbins,f_start,f_end);
    ahists[k]=new TH1D("a"+dists[k],"ahist",abins,a_start,a_end);
    tdurHists[k] = new TH1D("tdur"+dists[k],"tdurHist",tdurBins,tdurStart,tdurEnd);
    tstartHists[k] = new TH1D("tstart"+dists[k],"tstartHist",tstartBins,tstartStart,tstartEnd);
    Int_t numpoints2=0;
    Int_t numpoints2p=0;
    for(i=0;i<numpoints;i++){
      
      //add params for all runs with SM>thresh
      if(lnSMgivenS[i]>thresh){
        fshists[k]->Fill(fsbest[i]);
        ahists[k]->Fill(abest[i]);
        tdurHists[k]->Fill(tdurBest2[i]);
        tstartHists[k]->Fill(tstartBest[i]);
        numpoints2+=1;
      }
    }

    cout<<"num runs above thresh "<<numpoints2<<endl;
    //cout<<"numpointsp "<<numpoints2p<<endl;
    fshists[k]->SetLineColor(linecols[k]);
    fshists[k]->SetLineWidth(3);
    leg2->AddEntry(fshists[k],dists[k]+" kpc","l");

    //normalize fshists[k]
    norm=1;
    scale = norm/fshists[k]->Integral();
    fshists[k]->Scale(scale);
    cout<<"fshists integral "<<k<<" "<<fshists[k]->Integral()<<endl;

    ahists[k]->SetLineColor(linecols[k]);
    ahists[k]->SetLineWidth(3);
    //leg3->AddEntry(ahists[k],dists[k]+" kpc","l");
    leg3->AddEntry(ahists[k],dists[k],"l");

    //normalize ahists[k]
    scale = norm/ahists[k]->Integral();
    ahists[k]->Scale(scale);
    cout<<"ahists integral "<<k<<" "<<ahists[k]->Integral()<<endl;

    tdurHists[k]->SetLineColor(linecols[k]);
    tdurHists[k]->SetLineWidth(3);
    //leg4->AddEntry(tdurHists[k],dists[k]+" kpc","l");
    leg4->AddEntry(tdurHists[k],dists[k],"l");

    //normalize tdurHists[k]
    scale = norm/tdurHists[k]->Integral();
    tdurHists[k]->Scale(scale);
    cout<<"tdurHists integral "<<k<<" "<<tdurHists[k]->Integral()<<endl;

    tstartHists[k]->SetLineColor(linecols[k]);
    tstartHists[k]->SetLineWidth(3);
    //leg5->AddEntry(tstartHists[k],dists[k]+" kpc","l");
    leg5->AddEntry(tstartHists[k],dists[k],"l");

    //normalize tstartHists[k]
    scale = norm/tstartHists[k]->Integral();
    tstartHists[k]->Scale(scale);
    cout<<"tstartHists integral "<<k<<" "<<tstartHists[k]->Integral()<<endl;

    //print mean and std of each hist
    //cout<<dists[k]<<endl;
    //cout<<" fs mean and std: "<<fshists[k]->GetMean()<<" "<<fshists[k]->GetStdDev()<<endl;
    //cout<<" a mean and std: "<<ahists[k]->GetMean()<<" "<<ahists[k]->GetStdDev()<<endl;

    cout<<"****************************"<<endl;
  }

  //plot ROC
  gStyle->SetOptStat(0);

  TCanvas *c1=new TCanvas("c1","ROC",800,700);
  ROCgraphs[0]->SetTitle("");
  ROCgraphs[0]->GetXaxis()->SetTitle("P_{FI}");
  ROCgraphs[0]->GetYaxis()->SetTitle("P_{D}");
  ROCgraphs[0]->SetMaximum(1);
  ROCgraphs[0]->Draw("al");
  for(k=1;k<numfiles;k++){
    ROCgraphs[k]->Draw("same");
  }
  leg1->SetBorderSize(0);
  leg1->Draw();

  c1->SetLeftMargin(0.15);

  //plot fshists
  TCanvas *c2=new TCanvas("c2","fshist",800,700);
  fshists[0]->SetTitle("fshist");
  fshists[0]->GetXaxis()->SetTitle("f (Hz)");
  fshists[0]->GetYaxis()->SetTitle("prob/Hz");
  fshists[0]->SetMaximum(1);
  fshists[0]->Draw("hist ");
  
  for(k=1;k<numfiles;k++){
    fshists[k]->Draw("hist  same");
  }

  leg2->SetBorderSize(0);
  leg2->Draw();

  c2->SetLeftMargin(0.15);

  //plot ahists
  TCanvas *c3=new TCanvas("c3","ahist",800,700);
  ahists[0]->SetTitle("relative amplitude");
  ahists[0]->GetXaxis()->SetTitle("a");
  ahists[0]->GetYaxis()->SetTitle("prob");
  ahists[0]->SetMaximum(1);
  ahists[0]->Draw("hist c");
  
  for(k=1;k<numfiles;k++){
    ahists[k]->Draw("hist c same");
  }
  
  leg3->SetBorderSize(0);
  leg3->Draw();

  c3->SetLeftMargin(0.15);

  //plot tdurHists
  TCanvas *c4=new TCanvas("c4","tdurHists",800,700);
  tdurHists[0]->SetTitle("duration");
  tdurHists[0]->GetXaxis()->SetTitle("time (ms)");
  tdurHists[0]->GetYaxis()->SetTitle("prob/ms");
  tdurHists[0]->SetMaximum(1);
  tdurHists[0]->Draw("hist c");
  
  for(k=1;k<numfiles;k++){
    tdurHists[k]->Draw("hist c same");
  }
  
  leg4->SetBorderSize(0);
  leg4->Draw();

  c4->SetLeftMargin(0.15);  

  //plot tstartHists
  TCanvas *c5=new TCanvas("c5","tstartHists",800,700);
  tstartHists[0]->SetTitle("start time");
  tstartHists[0]->GetXaxis()->SetTitle("time (ms)");
  tstartHists[0]->GetYaxis()->SetTitle("prob/ms");
  tstartHists[0]->SetMaximum(1);
  tstartHists[0]->Draw("hist c");
  
  for(k=1;k<numfiles;k++){
    tstartHists[k]->Draw("hist c same");
  }
  
  leg5->SetBorderSize(0);
  leg5->Draw();

  c5->SetLeftMargin(0.15);
}
