// Reads the digitized files, makes pinched_info.dat for use with Snowglobes pinched.C file
// Also makes key.dat to associate flux number with time slide
// Linearly interpolates flux values at linear timebins within the overlapping time range of 
// all parameters and flavors
// Note: produces luminosities appropriate for *fluences* in Snowglobes-- luminosities are in ergs/s, then integrated over the given time bin.
// Aryil Bechtel 2020, adapted from  make_garching_flux.C by Kate Scholberg

void makeflux(){

  TString indir="Tamborra_Lblue_sasi_model/Tamborra_Lblue_"; //directory of flux files
  //TString indir = "garching_sasi_model/garching_";

  ifstream infile;
  // Two timescales, 6 flavors, three types of quantities (lum, eavg, alpha)
  // One graph array per quantity, indexed by flavor.

  const Int_t numflavor = 3;
  //  TString flavname[6] = {"nue","nuebar","numu","numubar","nutau","nutaubar"};
  TString flavname[numflavor] = {"nue","nuebar","numu"};

  // numu is same as nux

  TGraph** luminosity_graphs = new TGraph*[numflavor];
  TGraph** avgen_graphs = new TGraph*[numflavor];
  TGraph** alpha_graphs = new TGraph*[numflavor];

  TGraph** luminosity_graphs_interp = new TGraph*[numflavor];
  TGraph** avgen_graphs_interp = new TGraph*[numflavor];
  TGraph** alpha_graphs_interp = new TGraph*[numflavor];

  const Int_t maxpoints=10000;
  Double_t time[maxpoints];
  Double_t yval[maxpoints];

  Int_t i=0;
  Int_t j=0;
  Double_t t1, t2, t3; //to get start times for each param and flavor if needed
  Double_t s1, s2, s3; //to get end time for each param and flavor if needed
  Double_t timescale=1.; //to convert time to seconds, if needed

  for (i;i<numflavor;i++){

    // Luminosity
    ifstream in;

    TString infilename2 = indir+"luminosity_"+flavname[i]+".txt";

    cout << "Reading "<< infilename2<<endl;

    in.open(infilename2);
    j=0;
    Int_t numpoints;

    while(1) {
      in >> time[j]>>yval[j];
      time[j]/=timescale;
      //      cout << "time "<<time[j]<<" "<<yval[j]<<endl;
      if (yval[j]<0) {yval[j]=0;}

      if (!in.good()) break;
      j++;

    } // End of file reading loop

    numpoints = j;
    s1=time[0];
    t1=time[j-1];
    cout<<"endtime lum "<<t1<<endl;
    in.close();

    luminosity_graphs[i]= new TGraph(numpoints,time,yval);

    // Average energy

    TString infilename3 = indir+"avgen_"+flavname[i]+".txt";

    cout <<"Reading " << infilename3<<endl;

    in.open(infilename3);
    j=0;

    while(1) {
      in >> time[j]>>yval[j];
      time[j]/=timescale;
      //     cout << "time "<<time[j]<<" "<<yval[j]<<endl;

      if (yval[j]<0) {yval[j]=0;}

      if (!in.good()) break;
      j++;

    } // End of file reading loop

    numpoints = j;
    s2=time[0];
    t2=time[j-1];
    cout<<"endtime lum "<<t2<<endl;
    in.close();

    avgen_graphs[i]= new TGraph(numpoints,time,yval);

    // Alpha

    TString infilename4 = indir+"alpha_"+flavname[i]+".txt";

    cout << "Reading "<<infilename4<<endl;

    in.open(infilename4);
    j=0;

    while(1) {
      in >> time[j]>>yval[j];
      time[j]/=timescale;
      //      cout << "time "<<time[j]<<" "<<yval[j]<<endl;
      if (yval[j]<0) {yval[j]=0;}

      if (!in.good()) break;
      j++;

    } // End of file reading loop

    numpoints = j;
    t3=time[j-1];
    s3=time[0];
    cout<<"endtime lum "<<t3<<endl;
    in.close();

    alpha_graphs[i]= new TGraph(numpoints,time,yval);
      //    graphs[i]->Print();

  } // End of loop over flavors

  //*****************************************
  //starts negative so put an offset
  
  // Double_t t_offset = 0.02;
  //Double_t t_start = -0.02+t_offset+0.001 ;
  //Double_t t_end=9.+t_offset;
  //t_end +=t_offset;
  //  cout<<"t_end "<<t_end<<endl;
  //*****************************************
   
// Get the ranges of each of the graphs
  Int_t k;
  Double_t xx,yy;
  Double_t alpha_first[numflavor];
  Double_t avgen_first[numflavor];
  Double_t luminosity_first[numflavor];
  Double_t alpha_last[numflavor];
  Double_t avgen_last[numflavor];
  Double_t luminosity_last[numflavor];

  Int_t numtimes=numflavor*3;
  Double_t starttimes[numtimes];
  Double_t endtimes[numtimes];

  //Saved for nue flux
  Double_t alpha_val_first[numflavor];
  Int_t ret;

  for (k=0;k<numflavor;k++) {
    Int_t numpt;
    numpt = alpha_graphs[k]->GetN();
    if (numpt>0) {
      ret = alpha_graphs[k]->GetPoint(0,xx,yy);
      alpha_first[k] = xx;
      starttimes[k]=alpha_first[k];

      alpha_val_first[k] = yy;
      ret = alpha_graphs[k]->GetPoint(numpt-1,xx,yy);
      alpha_last[k] = xx;
      endtimes[k]=alpha_last[k];

    }

    numpt = avgen_graphs[k]->GetN();
    if (numpt>0) {
      ret = avgen_graphs[k]->GetPoint(0,xx,yy);
      avgen_first[k] = xx;
      starttimes[k+3]=avgen_first[k];

      ret = avgen_graphs[k]->GetPoint(numpt-1,xx,yy);
      avgen_last[k] = xx;
      endtimes[k+3]=avgen_last[k];

    }

    numpt = luminosity_graphs[k]->GetN();
    if (numpt>0) {
      ret = luminosity_graphs[k]->GetPoint(0,xx,yy);
      luminosity_first[k] = xx;
      starttimes[k+6]=luminosity_first[k];

      ret = luminosity_graphs[k]->GetPoint(numpt-1,xx,yy);
      luminosity_last[k] = xx;
      endtimes[k+6]=luminosity_last[k];

    }

  }

  cout<<"start times:    "<<"end times:"<<endl;
  for(k=0;k<numtimes;k++){
    cout<<starttimes[k]<<" "<<endtimes[k]<<endl;
  }

  //find max start time and min end time
  Double_t t_start=starttimes[0];
  Double_t t_end=endtimes[0];

  for(j=0;j<numtimes;j++){

    if(starttimes[j]>t_start) t_start=starttimes[j];

    if(t_end>endtimes[j]) t_end=endtimes[j];
  }

  cout<<"t_start "<<t_start<<endl;
  cout<<"t_end "<<t_end<<endl;

  //now create linear time intervals
  Double_t timebinsize=1;//.001; //in seconds
  const Int_t ntimebins=(t_end-t_start)/timebinsize;
  cout<<"ntimebins "<<ntimebins<<endl;
  Double_t timebins[ntimebins+1];

  Double_t t_step=(t_end-t_start)/double(ntimebins);
  cout<<"t_step "<<t_step<<endl;
  
  Double_t tot_dt=0;

    // Output files
    ofstream outfile;
    outfile.open("Tamborra_Lblue_pinched_info_1msTEST.dat");

    ofstream outfile2;
    outfile2.open("Tamborra_Lblue_pinched_info_key_1msTEST.dat");
    
   
    Double_t alpha_nue[ntimebins], avgen_nue[ntimebins], luminosity_nue[ntimebins];
    Double_t alpha_nuebar[ntimebins], avgen_nuebar[ntimebins], luminosity_nuebar[ntimebins];
    Double_t alpha_nux[ntimebins], avgen_nux[ntimebins], luminosity_nux[ntimebins];


    for (i=0;i<ntimebins;i++){
      timebins[i]=i*t_step+t_start;
      //      cout<<"timebin "<<i<<" "<<timebins[i]<<endl;
      Double_t dt=t_step;
      tot_dt += dt;  // for check

      // Get the values from the plots.  But first check we are within the meaningful range for each graph for each flavor.

      /*
      //linear interpolate points from graph, no fine-tuning
      alpha_nue = alpha_graphs[0]->Eval(timebins[i],0,"");
      avgen_nue = avgen_graphs[0]->Eval(timebins[i],0,"");
      luminosity_nue = luminosity_graphs[0]->Eval(timebins[i],0,"");

      alpha_nuebar = alpha_graphs[1]->Eval(timebins[i],0,"");
      avgen_nuebar = avgen_graphs[1]->Eval(timebins[i],0,"");
      luminosity_nuebar = luminosity_graphs[1]->Eval(timebins[i],0,"");

      alpha_nux = alpha_graphs[2]->Eval(timebins[i],0,"");
      avgen_nux = avgen_graphs[2]->Eval(timebins[i],0,"");
      luminosity_nux = luminosity_graphs[2]->Eval(timebins[i],0,"");
      */
            
// Nue
if (timebins[i]>=alpha_first[0] && timebins[i]<=alpha_last[0]
         && timebins[i]>=avgen_first[0] && timebins[i]<=avgen_last[0]
    && timebins[i]>=luminosity_first[0] && timebins[i]<=luminosity_last[0]) {

   alpha_nue[i] = alpha_graphs[0]->Eval(timebins[i],0,"");
   avgen_nue[i] = avgen_graphs[0]->Eval(timebins[i],0,"");
   luminosity_nue[i] = luminosity_graphs[0]->Eval(timebins[i],0,"");

 } else {
   alpha_nue[i]=0.;
   avgen_nue[i] = 0.;
   luminosity_nue[i]=0.;
  }


// Nuebar
 if (timebins[i]>=alpha_first[1] && timebins[i]<=alpha_last[1]
          && timebins[i]>=avgen_first[1] && timebins[i]<=avgen_last[1]
     && timebins[i]>=luminosity_first[1] && timebins[i]<=luminosity_last[1]) {
	
   alpha_nuebar[i] = alpha_graphs[1]->Eval(timebins[i],0,"");
   avgen_nuebar[i] = avgen_graphs[1]->Eval(timebins[i],0,"");
   luminosity_nuebar[i] = luminosity_graphs[1]->Eval(timebins[i],0,"");
   
 } else {
   alpha_nuebar[i]=0.;
   avgen_nuebar[i] = 0.;
   luminosity_nuebar[i]=0.;
 }

 //Numu

 if (timebins[i]>=alpha_first[2] && timebins[i]<=alpha_last[2]
          && timebins[i]>=avgen_first[2] && timebins[i]<=avgen_last[2]
     && timebins[i]>=luminosity_first[2] && timebins[i]<=luminosity_last[2]) {
   
   alpha_nux[i] = alpha_graphs[2]->Eval(timebins[i],0,"");
   avgen_nux[i] = avgen_graphs[2]->Eval(timebins[i],0,"");
   luminosity_nux[i] = luminosity_graphs[2]->Eval(timebins[i],0,"");
   
 } else {
   alpha_nux[i]=0.;
   avgen_nux[i] = 0.;
   luminosity_nux[i]=0.;
 }
       
 // Save the time and dt for each flux file number in the key file
 outfile2 << i<< " "<<timebins[i]<<" "<<dt<<endl;

 // Correct luminosity units

 //Double_t lumscale=1.e52;
 Double_t lumscale=1;
 outfile << i<< " "
	 <<alpha_nue[i]<<" "
	 <<alpha_nuebar[i]<<" "
	 <<alpha_nux[i]<<" "
	 <<avgen_nue[i]<<" "
	 <<avgen_nuebar[i]<<" "
	 <<avgen_nux[i]<<" "
	 <<luminosity_nue[i]*lumscale<<" "
	 <<luminosity_nuebar[i]*lumscale<<" "
	 <<luminosity_nux[i]*lumscale<<endl;
    }

    //    cout << "tot_dt "<< tot_dt <<" "<<t_end-t_start<<endl;

   outfile.close();


   TLegend* leg1 = new TLegend(0.6,0.7,0.8,0.75);

   TCanvas * c1 = new TCanvas("c1","luminosity_interp",800,700);

   Int_t lineWidth = 2;
   
   luminosity_graphs_interp[0]=new TGraph(ntimebins,timebins,luminosity_nue);
   luminosity_graphs_interp[0]->GetXaxis()->SetTitle("time (ms)");
   luminosity_graphs_interp[0]->GetYaxis()->SetTitle("Luminosity (10^{52} ergs/s)");
   luminosity_graphs_interp[0]->SetTitle("");
   luminosity_graphs_interp[0]->SetLineColor(2);
   luminosity_graphs_interp[0]->SetMinimum(0);
   luminosity_graphs_interp[0]->GetXaxis()->SetLimits(t_start,t_end);
   luminosity_graphs_interp[0]->SetLineWidth(lineWidth);
   luminosity_graphs_interp[0]->Draw("al");
   leg1->AddEntry(luminosity_graphs_interp[0],"#nu_{e}","l");
   
   luminosity_graphs_interp[1]=new TGraph(ntimebins,timebins,luminosity_nuebar);
   luminosity_graphs_interp[1]->SetLineColor(3);
   luminosity_graphs_interp[1]->SetLineWidth(lineWidth);
   luminosity_graphs_interp[1]->Draw("same");
   leg1->AddEntry(luminosity_graphs_interp[1],"#bar{#nu}_{e}","l");

   luminosity_graphs_interp[2]=new TGraph(ntimebins,timebins,luminosity_nux);
   luminosity_graphs_interp[2]->SetLineColor(4);
   luminosity_graphs_interp[2]->SetLineWidth(lineWidth);
   luminosity_graphs_interp[2]->Draw("same");
   leg1->AddEntry(luminosity_graphs_interp[2],"#nu_{x}","l");
   leg1->SetBorderSize(0);
   leg1->Draw();


   TLegend* leg2 = new TLegend(0.6,0.7,0.8,0.75);

   TCanvas * c2 = new TCanvas("c2","avgen_interp",800,700);

   avgen_graphs_interp[0]=new TGraph(ntimebins,timebins,avgen_nue);
   avgen_graphs_interp[0]->GetXaxis()->SetTitle("time (ms)");
   avgen_graphs_interp[0]->GetYaxis()->SetTitle("Average energy (MeV)");
   avgen_graphs_interp[0]->SetTitle("");
   avgen_graphs_interp[0]->SetLineColor(2);
   avgen_graphs_interp[0]->SetMinimum(0);
   avgen_graphs_interp[0]->SetMaximum(19.);
   avgen_graphs_interp[0]->GetXaxis()->SetLimits(t_start,t_end);
   avgen_graphs_interp[0]->SetLineWidth(lineWidth);
   avgen_graphs_interp[0]->Draw("al");
   leg2->AddEntry(avgen_graphs_interp[0],"#nu_{e}","l");

   avgen_graphs_interp[1]=new TGraph(ntimebins,timebins,avgen_nuebar);
   avgen_graphs_interp[1]->SetLineColor(3);
   avgen_graphs_interp[1]->SetLineWidth(lineWidth);
   avgen_graphs_interp[1]->Draw("same");
   leg2->AddEntry(avgen_graphs_interp[1],"#bar{#nu}_{e}","l");

   avgen_graphs_interp[2]=new TGraph(ntimebins,timebins,avgen_nux);
   avgen_graphs_interp[2]->SetLineColor(4);
   avgen_graphs_interp[2]->SetLineWidth(lineWidth);
   avgen_graphs_interp[2]->Draw("same");
   leg2->AddEntry(avgen_graphs_interp[2],"#nu_{x}","l");
   leg2->SetBorderSize(0);
   //leg2->Draw();

   TLegend* leg3 = new TLegend(0.6,0.7,0.8,0.75);

   TCanvas * c3 = new TCanvas("c3","alpha_interp",800,700);

   alpha_graphs_interp[0]=new TGraph(ntimebins,timebins,alpha_nue);
   alpha_graphs_interp[0]->GetXaxis()->SetTitle("time (ms)");
   alpha_graphs_interp[0]->GetYaxis()->SetTitle("#alpha");
   alpha_graphs_interp[0]->SetTitle("");
   alpha_graphs_interp[0]->SetLineColor(2);
   alpha_graphs_interp[0]->SetMinimum(0);
   alpha_graphs_interp[0]->GetXaxis()->SetLimits(t_start,t_end);
   alpha_graphs_interp[0]->SetLineWidth(lineWidth);
   alpha_graphs_interp[0]->Draw("al");
   leg3->AddEntry(alpha_graphs_interp[0],"#nu_{e}","l");


   alpha_graphs_interp[1]=new TGraph(ntimebins,timebins,alpha_nuebar);
   alpha_graphs_interp[1]->SetLineColor(3);
   alpha_graphs_interp[1]->SetLineWidth(lineWidth);
   alpha_graphs_interp[1]->Draw("same");
   leg3->AddEntry(alpha_graphs_interp[1],"#bar{#nu}_{e}","l");

   alpha_graphs_interp[2]=new TGraph(ntimebins,timebins,alpha_nux);
   alpha_graphs_interp[2]->SetLineColor(4);
   alpha_graphs_interp[2]->SetLineWidth(lineWidth);
   alpha_graphs_interp[2]->Draw("same");
   leg3->AddEntry(alpha_graphs_interp[2],"#nu_{x}","l");
   leg3->SetBorderSize(0);
   //leg3->Draw();

   //COMMENT OUT SECTION BELOW FOR NO SMOOTHING
   //fit polynomial to parameters to create smooth signal
   /*
   Double_t alpha_nue_smooth[ntimebins], avgen_nue_smooth[ntimebins], luminosity_nue_smooth[ntimebins];
   Double_t alpha_nuebar_smooth[ntimebins], avgen_nuebar_smooth[ntimebins], luminosity_nuebar_smooth[ntimebins];
   Double_t alpha_nux_smooth[ntimebins], avgen_nux_smooth[ntimebins], luminosity_nux_smooth[ntimebins];

   TGraph** luminosity_graphs_smooth = new TGraph*[numflavor];
   TGraph** avgen_graphs_smooth = new TGraph*[numflavor];
   TGraph** alpha_graphs_smooth = new TGraph*[numflavor];

   for(k=0;k<numflavor;k++){

   TF1* f1=new TF1("f1","pol7"); //fit a 5-parameter polynomial to smoothed data
   luminosity_graphs_interp[k]->Fit(f1);
   //get parameters
   Double_t param0=f1->GetParameter(0);
   Double_t param1=f1->GetParameter(1);
   Double_t param2=f1->GetParameter(2);
   Double_t param3=f1->GetParameter(3);
   Double_t param4=f1->GetParameter(4);
   Double_t param5=f1->GetParameter(5);
   Double_t param6=f1->GetParameter(6);
   Double_t param7=f1->GetParameter(7);
   //calculate fit value at each time point in range

   //nue
   if(k==0){
   for(i=0;i<ntimebins;i++){
     luminosity_nue_smooth[i]=param0+param1*timebins[i]+param2*TMath::Power(timebins[i],2)+param3*TMath::Power(timebins[i],3)+param4*TMath::Power(timebins[i],4)+param5*TMath::Power(timebins[i],5)+param6*TMath::Power(timebins[i],6)+param7*TMath::Power(timebins[i],7);
   }
   }

   //nuebar
   if(k==1){
     for(i=0;i<ntimebins;i++){
       luminosity_nuebar_smooth[i]=param0+param1*timebins[i]+param2*TMath::Power(timebins[i],2)+param3*TMath::Power(timebins[i],3)+param4*TMath::Power(timebins[i],4)+param5*TMath::Power\
	 (timebins[i],5)+param6*TMath::Power(timebins[i],6)+param7*TMath::Power(timebins[i],7);
     }
   }

   //nux
   if(k==2){
     for(i=0;i<ntimebins;i++){
       luminosity_nux_smooth[i]=param0+param1*timebins[i]+param2*TMath::Power(timebins[i],2)+param3*TMath::Power(timebins[i],3)+param4*TMath::Power(timebins[i],4)+param5*TMath::Power\
	 (timebins[i],5)+param6*TMath::Power(timebins[i],6)+param7*TMath::Power(timebins[i],7);
     }
   }
   }


   for(k=0;k<numflavor;k++){

     TF1* f1=new TF1("f1","pol7"); //fit a 5-parameter polynomial to smoothed data
     avgen_graphs_interp[k]->Fit(f1);
     //get parameters
     Double_t param0=f1->GetParameter(0);
     Double_t param1=f1->GetParameter(1);
     Double_t param2=f1->GetParameter(2);
     Double_t param3=f1->GetParameter(3);
     Double_t param4=f1->GetParameter(4);
     Double_t param5=f1->GetParameter(5);
     Double_t param6=f1->GetParameter(6);
     Double_t param7=f1->GetParameter(7);
     //calculate fit value at each time point in range

     //nue
     if(k==0){
       for(i=0;i<ntimebins;i++){
	 avgen_nue_smooth[i]=param0+param1*timebins[i]+param2*TMath::Power(timebins[i],2)+param3*TMath::Power(timebins[i],3)+param4*TMath::Power(timebins[i],4)+param5\
	   *TMath::Power(timebins[i],5)+param6*TMath::Power(timebins[i],6)+param7*TMath::Power(timebins[i],7);
       }
     }

     //nuebar
     if(k==1){
       for(i=0;i<ntimebins;i++){
	 avgen_nuebar_smooth[i]=param0+param1*timebins[i]+param2*TMath::Power(timebins[i],2)+param3*TMath::Power(timebins[i],3)+param4*TMath::Power(timebins[i],4)+\
	   param5*TMath::Power\
	   (timebins[i],5)+param6*TMath::Power(timebins[i],6)+param7*TMath::Power(timebins[i],7);
       }
     }

     //nux
     if(k==2){
       for(i=0;i<ntimebins;i++){
	 avgen_nux_smooth[i]=param0+param1*timebins[i]+param2*TMath::Power(timebins[i],2)+param3*TMath::Power(timebins[i],3)+param4*TMath::Power(timebins[i],4)+\
	   param5*TMath::Power\
	   (timebins[i],5)+param6*TMath::Power(timebins[i],6)+param7*TMath::Power(timebins[i],7);
       }
     }
   }
   


   for(k=0;k<numflavor;k++){

     TF1* f1=new TF1("f1","pol7"); //fit a 5-parameter polynomial to smoothed data
     alpha_graphs_interp[k]->Fit(f1);
     //get parameters
     Double_t param0=f1->GetParameter(0);
     Double_t param1=f1->GetParameter(1);
     Double_t param2=f1->GetParameter(2);
     Double_t param3=f1->GetParameter(3);
     Double_t param4=f1->GetParameter(4);
     Double_t param5=f1->GetParameter(5);
     Double_t param6=f1->GetParameter(6);
     Double_t param7=f1->GetParameter(7);
     //calculate fit value at each time point in range

     //nue
     if(k==0){
       for(i=0;i<ntimebins;i++){
         alpha_nue_smooth[i]=param0+param1*timebins[i]+param2*TMath::Power(timebins[i],2)+param3*TMath::Power(timebins[i],3)+param4*TMath::Power(timebins[i],4)+param5\
           *TMath::Power(timebins[i],5)+param6*TMath::Power(timebins[i],6)+param7*TMath::Power(timebins[i],7);
       }
     }
       //nuebar
       if(k==1){
	 for(i=0;i<ntimebins;i++){
	   alpha_nuebar_smooth[i]=param0+param1*timebins[i]+param2*TMath::Power(timebins[i],2)+param3*TMath::Power(timebins[i],3)+param4*TMath::Power(timebins[i],4)+\
	     param5*TMath::Power\
	     (timebins[i],5)+param6*TMath::Power(timebins[i],6)+param7*TMath::Power(timebins[i],7);
	 }
       }

       //nux
       if(k==2){
	 for(i=0;i<ntimebins;i++){
	   alpha_nux_smooth[i]=param0+param1*timebins[i]+param2*TMath::Power(timebins[i],2)+param3*TMath::Power(timebins[i],3)+param4*TMath::Power(timebins[i],4)+\
	     param5*TMath::Power\
	     (timebins[i],5)+param6*TMath::Power(timebins[i],6)+param7*TMath::Power(timebins[i],7);
	 }
       }
   }

   //save smoothed parameters in outfile
   outfile.open("sasi_pinched_info_1ms_smooth.dat");
   for(i=0;i<ntimebins;i++){ 
     Double_t lumscale=1.e52;
     outfile << i<< " "
	   <<alpha_nue_smooth[i]<<" "
	   <<alpha_nuebar_smooth[i]<<" "
	   <<alpha_nux_smooth[i]<<" "
	   <<avgen_nue_smooth[i]<<" "
	   <<avgen_nuebar_smooth[i]<<" "
	   <<avgen_nux_smooth[i]<<" "
	   <<luminosity_nue_smooth[i]*lumscale<<" "
	   <<luminosity_nuebar_smooth[i]*lumscale<<" "
	   <<luminosity_nux_smooth[i]*lumscale<<endl;
   }
   outfile.close();

   //plot smoothed parameters
   TCanvas *c4=new TCanvas("c4","smooth luminosity",800,700);
   luminosity_graphs_smooth[0]=new TGraph(ntimebins,timebins,luminosity_nue_smooth);
   luminosity_graphs_smooth[0]->SetMinimum(0);
   luminosity_graphs_smooth[0]->SetLineColor(2);
   luminosity_graphs_smooth[0]->Draw("al");

   luminosity_graphs_smooth[1]=new TGraph(ntimebins,timebins,luminosity_nuebar_smooth);
   luminosity_graphs_smooth[1]->SetLineColor(3);
   luminosity_graphs_smooth[1]->Draw("same");

   luminosity_graphs_smooth[2]=new TGraph(ntimebins,timebins,luminosity_nux_smooth);
   luminosity_graphs_smooth[2]->SetLineColor(4);
   luminosity_graphs_smooth[2]->Draw("same");

   TCanvas *c5=new TCanvas("c5","smooth avgen",800,700);
   avgen_graphs_smooth[0]=new TGraph(ntimebins,timebins,avgen_nue_smooth);
   avgen_graphs_smooth[0]->SetMinimum(0);
   avgen_graphs_smooth[0]->SetLineColor(2);
   avgen_graphs_smooth[0]->Draw("al");

   avgen_graphs_smooth[1]=new TGraph(ntimebins,timebins,avgen_nuebar_smooth);
   avgen_graphs_smooth[1]->SetLineColor(3);
   avgen_graphs_smooth[1]->Draw("same");

   avgen_graphs_smooth[2]=new TGraph(ntimebins,timebins,avgen_nux_smooth);
   avgen_graphs_smooth[2]->SetLineColor(4);
   avgen_graphs_smooth[2]->Draw("same");

   TCanvas *c6=new TCanvas("c6","smooth alpha",800,700);
   alpha_graphs_smooth[0]=new TGraph(ntimebins,timebins,alpha_nue_smooth);
   alpha_graphs_smooth[0]->SetMinimum(0);
   alpha_graphs_smooth[0]->SetLineColor(2);
   alpha_graphs_smooth[0]->Draw("al");

   alpha_graphs_smooth[1]=new TGraph(ntimebins,timebins,alpha_nuebar_smooth);
   alpha_graphs_smooth[1]->SetLineColor(3);
   alpha_graphs_smooth[1]->Draw("same");

   alpha_graphs_smooth[2]=new TGraph(ntimebins,timebins,alpha_nux_smooth);
   alpha_graphs_smooth[2]->SetLineColor(4);
   alpha_graphs_smooth[2]->Draw("same");
   */
}

