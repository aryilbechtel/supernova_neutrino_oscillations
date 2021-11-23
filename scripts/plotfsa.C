//plots probability density distributions of fs and a for several
   //SMvals.dat files (several distances)
//includes option to combine distributions from multiple SMvals.dat 
   //(ex: same distance, distributions from multiple detectors)
//Aryil Bechtel 2020

void plotfsa(){

  const Int_t numfiles=5;

  TString dists[numfiles]={"1","3","4","5","10"};
  //TString dists[numfiles]={"1"};

  Int_t linecols[numfiles]={2,3,4,1,28};
  //Int_t linecols[numfiles]={2};
 
  TString indir="SMvals/SMvals_ar40kt/Tamborra_Lblue/tstart_tdur/";
  TString numruns="100";
  TString modeltype="_m3_100to486ms_10to100dur_ar40kt";

  TString indir2="SMvals/SMvals_scint20kt/";
  TString modeltype2="_scint20kt";

  Double_t alpha=0.1; //threshold pfi value
  Double_t SM_start=-20;
  Double_t SM_end=300;
  Double_t binsize=.05;
  Int_t nbins=(SM_end-SM_start)/binsize;
  //Int_t nbins=1000;

  TGraph **ROCgraphs=new TGraph*[numfiles];
  TLegend* leg1 = new TLegend(0.7,0.7,.87,.87);

  TH1D **fshists=new TH1D*[numfiles];
  TLegend* leg2=new TLegend(0.6,0.7,0.8,0.75);
  Int_t fbins=35;
  Double_t f_start=57;
  Double_t f_end=228;

  TH1D **ahists=new TH1D*[numfiles];
  TLegend* leg3=new TLegend(0.6,0.7,0.8,0.75);
  Int_t abins=30;
  Double_t a_start=0;
  Double_t a_end=0.6;

  TH1D *lnSMgivenShist=new TH1D("lnSMgivenShist"," ;lnSM;prob",nbins,SM_start,SM_end); //histogram for lnSMgivenS vals
  TH1D *lnSMgivenNShist=new TH1D("lnSMgivenNShist","lnSM",nbins,SM_start,SM_end); //histogram for lnSMgivenNS vals

  TH1D *lnSMgivenSphist=new TH1D("lnSMgivenSphist"," ;lnSM;prob",nbins,SM_start,SM_end); //histogram for lnSMgivenSp vals
  TH1D *lnSMgivenNSphist=new TH1D("lnSMgivenNSphist","lnSM",nbins,SM_start,SM_end); //histogram for lnSMgivenNSp vals


  Int_t k;
  for(k=0;k<numfiles;k++){

    lnSMgivenShist->Reset();
    lnSMgivenNShist->Reset();

    TString infilename=indir+"SMvals_"+numruns+modeltype+"_"+dists[k]+"kpc.dat";
    //TString infilename=indir+"chi2ratio_"+numruns+modeltype+"_"+dists[k]+"kpc.dat";
    cout<<"opened "<<infilename<<endl;

  ifstream in;
  in.open(infilename);
  Int_t i=0;
  Int_t maxpoints=10000;
  Double_t lnSMgivenS[maxpoints];
  Double_t lnSMgivenNS[maxpoints];
  Double_t fsbest[maxpoints]; //fs given S
  Double_t abest[maxpoints];
  Double_t fsbest2[maxpoints]; //fs given NS
  Double_t abest2[maxpoints];
  
  TString infilename2=indir2+"SMvals_"+numruns+modeltype2+"_"+dists[k]+"kpc.dat";
  cout<<"opened "<<infilename2<<endl;

  ifstream in2;
  in2.open(infilename2);
  Double_t lnSMgivenSp[maxpoints];
  Double_t lnSMgivenNSp[maxpoints];
  Double_t fsbestp[maxpoints]; //fs given S
  Double_t abestp[maxpoints];
  Double_t fsbest2p[maxpoints]; //fs given NS
  Double_t abest2p[maxpoints];

  Double_t totlnSMgivenS[maxpoints];
  Double_t totlnSMgivenNS[maxpoints];
  
  while(1){
    in>>lnSMgivenS[i]>>lnSMgivenNS[i]>>fsbest[i]>>abest[i]>>fsbest2[i]>>abest2[i];
    
    in2>>lnSMgivenSp[i]>>lnSMgivenNSp[i]>>fsbestp[i]>>abestp[i]>>fsbest2p[i]>>abest2p[i];
    totlnSMgivenS[i]=lnSMgivenS[i]+lnSMgivenSp[i];
    totlnSMgivenNS[i]=lnSMgivenNS[i]+lnSMgivenNSp[i];
    
    if(!in.good()) break;
    i++;
  }
  Int_t numpoints=i;
  in.close();
  in2.close();

  for(i=0;i<numpoints;i++){
    lnSMgivenShist->Fill(lnSMgivenS[i]);
    lnSMgivenNShist->Fill(lnSMgivenNS[i]);
    
    lnSMgivenSphist->Fill(lnSMgivenSp[i]);
    lnSMgivenNSphist->Fill(lnSMgivenNSp[i]);
    
  }
 
  //normalize the histograms to give probability
  Double_t norm=1;
  Double_t scale = norm/lnSMgivenShist->Integral();
  lnSMgivenShist->Scale(scale);
  scale = norm/lnSMgivenNShist->Integral();
  lnSMgivenNShist->Scale(scale);
  
  scale = norm/lnSMgivenSphist->Integral();
  lnSMgivenSphist->Scale(scale);
  scale = norm/lnSMgivenNSphist->Integral();
  lnSMgivenNSphist->Scale(scale);
  
  cout<<"lnSMgivenShist integral "<<lnSMgivenShist->Integral()<<endl;
  cout<<"lnSMgivenNShist integral "<<lnSMgivenNShist->Integral()<<endl;
  
  cout<<"lnSMgivenSphist integral "<<lnSMgivenSphist->Integral()<<endl;
  cout<<"lnSMgivenNSphist integral "<<lnSMgivenNSphist->Integral()<<endl;
  
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
  
  Double_t pfip[nbins];
  Double_t pdp[nbins];
  Double_t threshp;
  Int_t ithreshp;
  for(i=0;i<nbins;i++){
    pfip[i]=lnSMgivenNSphist->Integral(i,nbins);
    pdp[i]=lnSMgivenSphist->Integral(i,nbins);
    //cout<<"pfi "<<pfi[i]<<endl;
    if(pfip[i]>alpha){
      threshp=lnSMgivenSphist->GetXaxis()->GetBinCenter(i+1);
      ithreshp=i;
    }
  }
  cout<<"final threshp "<<threshp<<endl;
  cout<<"pfip above threshp "<<lnSMgivenNSphist->Integral(ithreshp,nbins)<<endl;
  

  ROCgraphs[k]=new TGraph(nbins,pfi,pd);
  ROCgraphs[k]->SetLineColor(linecols[k]);
  ROCgraphs[k]->SetLineWidth(3);
  leg1->AddEntry(ROCgraphs[k],dists[k]+" kpc","l");


  fshists[k]=new TH1D("fs"+dists[k],"fshist",fbins,f_start,f_end);
  ahists[k]=new TH1D("a"+dists[k],"ahist",abins,a_start,a_end);
  Int_t numpoints2=0;
  Int_t numpoints2p=0;
  for(i=0;i<numpoints;i++){
    
    
    if(lnSMgivenS[i]>thresh){
      fshists[k]->Fill(fsbest[i]);
      ahists[k]->Fill(abest[i]);
      numpoints2+=1;
    }
    
    /*
    if(lnSMgivenSp[i]>threshp){
      fshists[k]->Fill(fsbestp[i]);
      ahists[k]->Fill(abestp[i]);
      numpoints2p+=1;
    }
    */
    /*
    fshists[k]->Fill(fsbest[i]);
    ahists[k]->Fill(abest[i]);
    */
  }

  cout<<"numpoints "<<numpoints2<<endl;
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
  leg3->AddEntry(ahists[k],dists[k]+" kpc","l");

  //normalize ahists[k]
  scale = norm/ahists[k]->Integral();
  ahists[k]->Scale(scale);
  cout<<"ahists integral "<<k<<" "<<ahists[k]->Integral()<<endl;

  //print mean and std of each hist
  cout<<dists[k]<<endl;
  cout<<" fs mean and std: "<<fshists[k]->GetMean()<<" "<<fshists[k]->GetStdDev()<<endl;
  cout<<" a mean and std: "<<ahists[k]->GetMean()<<" "<<ahists[k]->GetStdDev()<<endl;

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
  fshists[0]->Draw("hist C");
  
  for(k=1;k<numfiles;k++){
    fshists[k]->Draw("hist C same");
  }

  leg2->SetBorderSize(0);
  leg2->Draw();

  c2->SetLeftMargin(0.15);

  //plot ahists
  TCanvas *c3=new TCanvas("c3","ahist",800,700);
  ahists[0]->SetTitle("ahist");
  ahists[0]->GetXaxis()->SetTitle("a");
  ahists[0]->GetYaxis()->SetTitle("prob");
  ahists[0]->SetMaximum(1);
  ahists[0]->Draw("hist C");
  
  for(k=1;k<numfiles;k++){
    ahists[k]->Draw("hist C same");
  }
  
  leg3->SetBorderSize(0);
  leg3->Draw();

  c3->SetLeftMargin(0.15);

}
