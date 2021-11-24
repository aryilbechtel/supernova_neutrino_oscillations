/*
  creates #runs events per bin time series from random Poisson process
    from eventspertime.dat file of event info

  PARAMETERS: frequency, relative amplitude, start time, and duration
  
  calculates #runs "SASI meter" (likelihood ratio=L(alternative)/L(null)) values from series

  TERMINOLOGY:
  givenS: given SASI oscillations
  givenNS: given no SASI oscillations
  
  saves SM values in SMvals.dat with format:
    lnSMgivenS lnSMgivenNS fsbestgivenS abestgivenS fsbestgivenNS abestgivenNS
  saves start times in tstart.dat with format:
    tstartgivenS tstartgivenNS
  saves duration times in tdur.dat with format:
    tdurgivenS tdurgivenNS
  
  generated no-SASI model is big-binned spline of moving avg'd SASI model
  template method:
    null-template: smoothed big-binning of moving-avg'd data
    2p-template: injected sinusoid over smoothed signal

V4: modified to allow for larger event rates and oscillation amplitudes:
    when Bessel fn input is too high, asymptotic form is used.
    logs also used to prevent inf values. 

*/

//Aryil Bechtel 2021

//set boundaries of algorithm
Double_t tdstep = 10; //duration step in ms
Double_t tdstart = 10;
Double_t tdend = 110;
Int_t numtd = (tdend - tdstart + 1)/tdstep;

Double_t timestep=10; //start time step in ms
Double_t t_start=900; //start time range in ms
Double_t t_end=1200;
Double_t f_start=57; //fq range in Hz
Double_t f_end=228; 

Int_t wsize=9; //moving window size (units of bins) for moving average of shot-noised data
Int_t nsbinsize=9; //big bin size: sample interval to get smoothed template
Int_t wedge=wsize/2;

Int_t runs = 35;

Int_t maxpoints=2000;
Double_t dist2 = 1.25;
Double_t dist_scale = 100/(TMath::Power(dist2,2));
Double_t scale=1000; //scaling factor to convert to ms and events/ms

#include <TVirtualFFT.h>

Double_t *tdurfn(Int_t fftnum, Double_t eventsrp[maxpoints], Double_t nseventsrp[maxpoints],
                Double_t timep[maxpoints],Double_t dtscaled, Int_t input_array[3]);

void sasimeterruns_tstart_tdur_v4(){

    TString outdir="/Users/aryilbechtel/Documents/Duke_18_22/neutrino_lab/sasi/sasistudybundle/plotssnowglobes/snewpy_models/Zha_2021/detectors/ar40kt/s23/";
    TString outfilename=outdir+"SMvals_35_m3_900to1200ms_10to100dur_ar40kt_1.25kpc.dat";
    TString outfilename2=outdir+"tstart_35_m3_900to1200ms_10to100dur_ar40kt_1.25kpc.dat";
    TString outfilename3=outdir+"tdur_35_m3_900to1200ms_10to100dur_ar40kt_1.25kpc.dat";

    ofstream out;
    ofstream out2;
    ofstream out3;
    out.open(outfilename);
    out2.open(outfilename2);
    out3.open(outfilename3);

    //get event info from eventspertime.dat file
    TString infilename="/Users/aryilbechtel/Documents/Duke_18_22/neutrino_lab/sasi/sasistudybundle/plotssnowglobes/snewpy_models/Zha_2021/detectors/ar40kt/s23/eventspertime_Zha_2021_s23_1ms.dat";
    Double_t time[maxpoints];
    Double_t events[maxpoints];
    Double_t rate[maxpoints];
    Double_t nseventsp[maxpoints];
    Double_t nsratep[maxpoints];
    //scale event rate and events per bin to mimic changes in fiducial distance 
    Double_t fdist=TMath::Sqrt(100/dist_scale); //calculate fiducial distance from dist_scale
    //cout<<"fiducial distance "<<fdist<<endl;

    ifstream in;
    in.open(infilename);  
    Int_t i=0;
    while(1){
        in>>time[i]>>events[i]>>rate[i];
        time[i]*=scale;
        rate[i]/=scale;
        events[i]*=dist_scale;
        rate[i]*=dist_scale;
        if(!in.good()) break;
        i++;
    }
    in.close();
    Int_t numpoints=i;
    Double_t dt=(time[numpoints-1]-time[numpoints-2]);
    Double_t dtscaled = dt/scale;
    //cout<<"dt "<<dt<<endl; 
    cout<<numpoints<<endl;

    Double_t t_final = numpoints - (nsbinsize + wedge + tdstart + 1);
    if(t_final<t_end){
        cout<<"ERROR: END START TIME IS OUT OF BOUNDS"<<endl;
        cout<<"MAX POSSIBLE END START TIME IS "<<t_final<<
            " MS FOR INPUT START DURATION"<<endl;
        exit(1);
    }

    if(t_start<(t_start-wedge)){
        cout<<"ERROR: INITIAL START TIME IS OUT OF BOUNDS"<<endl;
        cout<<"MIN POSSIBLE INITIAL START TIME IS "<<t_start-wedge<<
            " MS FOR INPUT MOVING AVG WINDOW SIZE"<<endl;
        exit(1);
    }

    //create no-SASI model with spline of big-binning of moving avg of SASI model
    Int_t i_start = 1+ (t_start-time[0])/dt;
    Int_t i_end = 1+ (t_end-time[0])/dt;
    Int_t i_startp = wedge;
    Int_t i_endp = (numpoints-1)-nsbinsize-wedge;
    Int_t numpointskp=i_endp-i_startp+1;
    //cout<<"numpointskp "<<numpointskp<<endl;

    Double_t avgevents[numpointskp+nsbinsize]; //array for moving avg'd SASI model
    for(i=0;i<numpointskp+nsbinsize;i++){
        avgevents[i]=0;
        for(Int_t c=0;c<wsize;c++){
            avgevents[i]+=events[i+i_startp-wedge+c];
        }
        avgevents[i]/=wsize;
    }
    Int_t numbigbins_ns=1+numpointskp/nsbinsize;
    Double_t timebigbin_ns[numbigbins_ns]; //array for big bin time array
    for(i=0;i<numbigbins_ns;i++){
        timebigbin_ns[i]=time[i_startp+i*nsbinsize];
    }
    TGraph *bigbin_ns=new TGraph(numbigbins_ns,timebigbin_ns,timebigbin_ns);
    //sample moving avg'd SASI in big bins and spline interpolate between
    Double_t events_nsbb[numbigbins_ns];
    for(i=0;i<numbigbins_ns;i++){
        events_nsbb[i]=avgevents[i*nsbinsize];
        bigbin_ns->SetPoint(i,timebigbin_ns[i],events_nsbb[i]);
    }
    for(i=0;i<numpointskp;i++){
        nseventsp[i]=bigbin_ns->Eval(time[i+i_startp],0,"S");
    }
    
    //crop sasi and time array
    Double_t timep[numpointskp]; //cropped time array
    Double_t eventsp[numpointskp];
    Int_t numbigbins=numbigbins_ns;

    for(i=0;i<numpointskp;i++){
        timep[i]=time[i+i_startp];
        eventsp[i]=events[i+i_startp];
    }

    //start run loop here
    Int_t m;
    //these vvv must be initialized within tdurfn
    //TGraph *sevbigbin=new TGraph(numbigbins);
    //TGraph *nsevbigbin=new TGraph(numbigbins);
    TRandom3 *pois1 = new TRandom3();
    TRandom3 *pois2 = new TRandom3();
    for(m=0;m<runs;m++){

        //Initialize final best fit parameters
        Double_t maxlnSMgivenSp=-9e7;
        Double_t maxlnSMgivenNSp=-9e7;
        Int_t s_startSp, s_startNSp;
        Double_t tstartSp, tstartNSp;
        Double_t abestpp, fsbestpp, abestpp2, fsbestpp2;
        Int_t tdur, tdurS, tdurNS;

        //shot-noise sasi and no-sasi models
        Double_t eventsrp[numpointskp]; 
        Double_t nseventsrp[numpointskp];
        for(i=0;i<numpointskp;i++){
            pois1->SetSeed(0);
            pois2->SetSeed(0);    
            eventsrp[i]=pois1->Poisson(eventsp[i]); //sasi run
            nseventsrp[i]=pois2->Poisson(nseventsp[i]); //no sasi run       
            //uncomment following 2 lines to remove shot-noise
            //eventsrp[i]=eventsp[i];
            //nseventsrp[i]=nseventsp[i];
            
        }

        //loop over time durations here
        for(i=0;i<numtd;i++){
            tdur = i*tdstep + tdstart;
            //create input array with ints needed to run tdurfn
            Int_t numtimes;
            Int_t i_last;
            Double_t t_last = time[numpoints - (nsbinsize + wedge + tdur + 1)];
            cout<<"t_last "<<t_last<<endl;
            cout<<"t_end "<<t_end<<endl;
            if(t_last<t_end){
                i_last = 1+ (t_last-time[0])/dt;
                numtimes = -1+(i_last-i_start+1)/timestep;
            }
            else numtimes = (i_end-i_start+1)/timestep;
            cout<<"numtimes for "<<tdur<<"ms: "<<numtimes<<endl;
            Int_t input_array[3] = {i_start,i_startp,numtimes};

            Double_t *tdurfn_return=tdurfn(tdur,eventsrp,nseventsrp,timep,dtscaled,input_array);
            
            Double_t maxlnSMgivenS=tdurfn_return[0];               
            Double_t maxlnSMgivenNS=tdurfn_return[1];

            Double_t fsbestp=tdurfn_return[2];
            Double_t abestp=tdurfn_return[3];
            Double_t fsbestp2=tdurfn_return[4];
            Double_t abestp2=tdurfn_return[5];

            Double_t tstartS=tdurfn_return[6];
            Double_t tstartNS=tdurfn_return[7];

            if(maxlnSMgivenS>maxlnSMgivenSp){
                maxlnSMgivenSp=maxlnSMgivenS;
                fsbestpp=fsbestp;
                abestpp=abestp;
                tstartSp=tstartS;
                tdurS=tdur;
            }
            if(maxlnSMgivenNS>maxlnSMgivenNSp){
                maxlnSMgivenNSp=maxlnSMgivenS;
                fsbestpp2=fsbestp2;
                abestpp2=abestp2;
                tstartNSp=tstartNS;
                tdurNS=tdur;
            }
            
            cout<<"run "<<m<<" duration "<<tdur<<" complete"<<endl;
            cout<<"----------------------"<<endl;
        }

        //save SMvals to outfile
        out<<maxlnSMgivenSp<<" "<<maxlnSMgivenNSp<<" "<<fsbestpp<<" "<<abestpp<<" "<<fsbestpp2<<" "<<abestpp2<<endl;

        //save start times to outfile2
        out2<<tstartSp<<" "<<tstartNSp<<endl;

        //save tdurs to outfile3
        out3<<tdurS<<" "<<tdurNS<<endl;


        cout<<"lnSMgivenS: "<<maxlnSMgivenSp
        <<" lnSMgivenNS: "<<maxlnSMgivenNSp
        <<" fsbest: "<<fsbestpp
        <<" abset: "<<abestpp
        <<" fsbest2: "<<fsbestpp2
        <<" abest2 "<<abestpp2
        <<endl;
        
        cout<<"best SASI start time and duration: "<<tstartSp<<" "<<tdurS<<endl;
        cout<<"best No-SASI start time and duration: "<<tstartNSp<<" "<<tdurNS<<endl;;
        cout<<"run "<<m<<" complete"<<endl;
        cout<<"**********************"<<endl;
        cout<<"**********************"<<endl;

    }


    out.close();
    out2.close();
    out3.close();

}

/*
    tdurfn (function):
        inputs: SASI duration (# pts in fft), shot-noised sasi, no-sasi models, 
                time array, input_array[3]
                    input_array[0]: index of first start time
                    input_array[1]: index of earliest possible start time
                    input_array[2]: number of start times

        outputs: lnSMgivenS and lnSMgivenNS maximized with respect
                to frequency, rel amp, start time 
                
                stored in array in following order:
                (maxlnSMgivenS, maxlnSMgivenNS, 
                fsbestp, abestp, tstartS, tdurS, 
                fsbestp2, abestp2, tstartNS,tdurNS)

*/
Double_t *tdurfn(Int_t fftnum, Double_t eventsrp[maxpoints], Double_t nseventsrp[maxpoints], Double_t timep[maxpoints],
                Double_t dtscaled, Int_t input_array[3]){
                    
    //Unpack input array
    Int_t i_start = input_array[0];
    Int_t i_startp = input_array[1];
    Int_t numtimes = input_array[2];

    //Initialize best fit parameters for start time, amp, freq
    Double_t maxlnSMgivenS=-999;
    Double_t maxlnSMgivenNS=-999;
    Int_t s_startS, s_startNS;
    Double_t tstartS, tstartNS;
    Double_t abestp, fsbestp, abestp2, fsbestp2;

    //Initialize fft pointer
    TVirtualFFT *trafo = TVirtualFFT::FFT(1, &fftnum, "R2HC");
    Int_t numpointsk = fftnum;

    //Initialize big binned graphs
    Int_t numbigbins = 1 + fftnum/nsbinsize;
    TGraph *sevbigbin=new TGraph(numbigbins);
    TGraph *nsevbigbin=new TGraph(numbigbins);

    //create fftpower lambda
    auto fftpower=[=](Double_t numpoints2,Double_t points2[maxpoints],Double_t dt2){
        trafo->SetPoints(points2);
        trafo->Transform();
        Double_t re, im;
        Double_t *power=new Double_t[maxpoints];
        Int_t fnyq=1/(2*dt2); //max resolvable freq, nyquist freq
        for(Int_t i=0; i<numpoints2+1; i++){
            trafo->GetPointComplex(i, re, im);
            //calculate the powerspectrum
            Double_t absoutsqr=(re*re)+(im*im);
            //cout<<"absoutsqr "<<i<<" "<<absoutsqr<<endl;
            Int_t xint=(i)/(dt2*numpoints2);
            if(xint==0 || xint==fnyq){
                power[i]=absoutsqr/(TMath::Power(numpoints2,2));
            }
            else if(i==numpoints2) power[i]=0; //for ease of plotting
            else{
                power[i]=2*(absoutsqr)/(TMath::Power(numpoints2,2)); //for frequencies 0<f<nyquist freq (SASI freq range is within that)
            }
        }
        return power;
        delete [] trafo;
    };

    //create twopfit lambda
    auto twopfit=[=](Int_t numpoints,Double_t time[maxpoints],Double_t power[maxpoints],Double_t A[maxpoints],Int_t if_start,Int_t if_end){
        Double_t n1=0.; //mean value of background rate
        Double_t phase=0; //phase is set constant before parameter fitting loops
        Double_t scale=1000;
        Double_t dtscaled=(time[1]-time[0])/scale;

        Int_t j;
        Double_t sasiLbest=0;
        Double_t abest=0;
        Double_t a_start=0.001;
        Double_t a_end=0.6;
        Double_t astep=.001; //determines precision of amplitude parameter
        Int_t anum=1+(a_end-a_start)/astep;

        Int_t l;
        Double_t fsbest=0;
        Double_t fs_start=f_start;
        Double_t fs_end=f_end;
        Double_t fsstep=1; //determines precision of freq parameter
        Int_t fsnum=1+(fs_end-fs_start)/fsstep;

        //loop over freq values
        for(l=0;l<fsnum;l++){
            Double_t fs=l*fsstep+fs_start;
            //loop over amplitude values
            for(j=0;j<anum;j++){
                Double_t a=j*astep+a_start;
                Double_t sasiparam[numpoints];
                Double_t sparamtotevents=0;
                for(Int_t i=0;i<numpoints;i++){
                    //built-in conversion from hz to rad/ms vvv
                    sasiparam[i]=(A[i]-n1)*(1.+a*TMath::Sin((2/scale)*TMath::Pi()*fs*time[i]+phase))+n1;
                    sparamtotevents+=sasiparam[i];
                }
                //calculate likelihood that sasi power series vector is realization of sasi hypothesis

                //fft and find power spectrum of 2-param model for given params
                Double_t *sasiparampower=fftpower(numpointsk,sasiparam,dtscaled);

                Double_t sigmasqr=sparamtotevents/2;
                Int_t k;
                Double_t sasiL2;
                Double_t logsasiL2;

                //approximate Bessel fn output when Bessel input > bessThresh
                //approximate using first two terms of asymptotic form of Bessel fn
                Double_t bessThresh = 700;
                
                //loop over frequencies
                for(k=if_start;k<if_end;k++){
                    //calculate prob distr of power spectrum at given freq
                    Double_t b=(TMath::Power(numpointsk,2))/(4*sigmasqr);
                    Double_t bessIn = 2*b*TMath::Sqrt(power[k]*sasiparampower[k]);
                    Double_t logbessApprox = bessIn - TMath::Log(TMath::Sqrt(2*TMath::Pi()*bessIn)) + 
                                        TMath::Log(1+1/(8*bessIn));
                    if(bessIn<bessThresh){
                        if(k==if_start){
                            logsasiL2=TMath::Log(b) + (-b*(power[k]+sasiparampower[k]))+ TMath::Log(TMath::BesselI0(bessIn));
                            sasiL2 = exp(logsasiL2);
                        }
                        else{
                            logsasiL2=TMath::Log(b) + (-b*(power[k]+sasiparampower[k]))+ TMath::Log(TMath::BesselI0(bessIn));
                            sasiL2 *= exp(logsasiL2);
                        }
                    }
                    else{
                        if(k==if_start){
                            logsasiL2=TMath::Log(b) + (-b*(power[k]+sasiparampower[k])) + logbessApprox;
                            sasiL2 = exp(logsasiL2);
                        }
                        else{
                            logsasiL2=TMath::Log(b) + (-b*(power[k]+sasiparampower[k])) + logbessApprox;
                            sasiL2 *= exp(logsasiL2);
                        }

                    }
                }
                //cout<<sasiL2<<endl;
                if(sasiL2>sasiLbest){
                    sasiLbest=sasiL2;
                //cout<<sasiLbest<<endl;
                    abest=a;
                    fsbest=fs;
                }
                delete [] sasiparampower;
                } //close amplitude loop
            } //close freq loop
        Double_t *fitreturn= new Double_t[3];
        fitreturn[0]=sasiLbest;
        fitreturn[1]=abest;
        fitreturn[2]=fsbest;
        return fitreturn;
    };
  
    //loop over start times
    Int_t i_startT;
    //Int_t numpointsk=fftnum;
    Int_t tdur_int = fftnum;
    Int_t q;
    Int_t i;
    for(q=0;q<numtimes;q++){
        i_startT=i_start+q*timestep;
        //cout<<"testing start time: "<<time[i_startT]<<endl;

        Double_t timek[numpointsk];
        Double_t eventsr[numpointsk];
        Double_t toteventsr=0;
        Double_t nseventsr[numpointsk];
        Double_t tot_nseventsr=0;

        //apply final crop to random runs
        for(i=0;i<numpointsk;i++){
            timek[i]=timep[i+i_startT-i_startp];
            eventsr[i]=eventsrp[i+i_startT-i_startp];
            toteventsr+=eventsr[i];
            nseventsr[i]=nseventsrp[i+i_startT-i_startp];
            tot_nseventsr+=nseventsr[i];
        }
    
        //take moving avg of random runs
        Double_t avgsevr[numpointsk+nsbinsize]; //array for moving avg'd SASI run
        Double_t avgnsevr[numpointsk+nsbinsize];
        for(i=0;i<numpointsk+nsbinsize;i++){
            avgsevr[i]=0;
            avgnsevr[i]=0;
                for(Int_t c=0;c<wsize;c++){
                avgsevr[i]+=eventsrp[i+i_startT-i_startp-wedge+c];
                avgnsevr[i]+=nseventsrp[i+i_startT-i_startp-wedge+c];
            }
            avgsevr[i]/=wsize;
            avgnsevr[i]/=wsize;
        }
        //sample moving avg'd runs in big bins
        Double_t timebigbin[numbigbins];
        for(i=0;i<numbigbins;i++){
            timebigbin[i]=timek[i*nsbinsize];
        }
        Double_t sA[numpointsk];
        Double_t nsA[numpointsk];
        Double_t sbigbin[numbigbins];
        Double_t nsbigbin[numbigbins];
        for(i=0;i<numbigbins;i++){
            sbigbin[i]=avgsevr[i*nsbinsize];
            sevbigbin->SetPoint(i,timebigbin[i],sbigbin[i]);
            nsbigbin[i]=avgnsevr[i*nsbinsize];
            nsevbigbin->SetPoint(i,timebigbin[i],nsbigbin[i]);
        }
        for(i=0;i<numpointsk;i++){
            sA[i]=sevbigbin->Eval(timek[i],0,"S");
            nsA[i]=nsevbigbin->Eval(timek[i],0,"S");
        }
        //CALCULATE POWER SPECTRUM OF SASI RUN
        Double_t x[numpointsk];
        Double_t *sasipower=fftpower(numpointsk,eventsr,dtscaled);
        for(i=0;i<numpointsk+1;i++) x[i]=(i)/(dtscaled*numpointsk);
        
        //find indices of frequency boundaries
        Double_t dx=x[1]; //freq step size in Hz
        //cout<<"dx "<<dx<<endl;
        Int_t if_start=f_start/dx; //may need to manually adjust these indices by one bc of rounding
        Int_t if_end=(f_end/dx)+1;
        Int_t if_interval=if_end-if_start;
        //cout<<"if_start "<<if_start<<" if_end "<<if_end<<endl;

        //FIT 2 PARAM TEMPLATE TO SASI RUN
        Double_t t_interval=timek[numpointsk-1]-timek[0];
        Double_t n1=0.; //mean value of background rate
        Double_t phase=0; //phase is set constant before parameter fitting loops

        Double_t *twopfitreturn=twopfit(numpointsk,timek,sasipower,sA,if_start,if_end);
        Double_t sasiLbest=twopfitreturn[0];
        Double_t abest=twopfitreturn[1];
        Double_t fsbest=twopfitreturn[2];

        //SAVE SASIPARAM ARRAY AND POWERSPECTRUM FOR SASI MODEL FOR PLOTTING
        Double_t sasiparamgivenS[numpointsk];
        for(i=0;i<numpointsk;i++){
            sasiparamgivenS[i]=(sA[i]-n1)*(1.+abest*TMath::Sin((2/scale)*TMath::Pi()*fsbest*timek[i]+phase))+n1;
        }
        Double_t *sasiparamgivenSpower=fftpower(numpointsk,sasiparamgivenS,dtscaled);

        //CREATE NULL-PARAM TEMPLATE FOR SASI RUN
        Double_t nosasiparam[numpointsk];
        Double_t nsparamtotevents=0;
        for(i=0;i<numpointsk;i++){
            nosasiparam[i]=sA[i]; //time-average event rate over time interval
            nsparamtotevents+=nosasiparam[i];
        }
        //cout<<"nsparamtotevents "<<nsparamtotevents<<endl;

        //CALCULATE POWER SPECTRUM OF NULL-PARAM TEMPLATE FOR SASI RUN
        Double_t *nosasiparampower=fftpower(numpointsk,nosasiparam,dtscaled);

        //CALCULATE LIKELIHOOD FOR SASIPOWER GIVEN NOSASI HYPOTHESIS
        Double_t sigmasqr=nsparamtotevents/2;
        Int_t k;
        Double_t nosasiL2;
        Double_t lognosasiL2;
        //approximate Bessel fn output when Bessel input > bessThresh 
        //approximate using first two terms of asymptotic form of Bessel fn
        Double_t bessThresh = 700;
        //loop over frequencies
        for(k=if_start;k<if_end;k++){

            Double_t b=(TMath::Power(numpointsk,2))/(4*sigmasqr);
            Double_t bessIn = 2*b*TMath::Sqrt(sasipower[k]*nosasiparampower[k]);
            Double_t logbessApprox = bessIn - TMath::Log(TMath::Sqrt(2*TMath::Pi()*bessIn)) + 
                                TMath::Log(1+1/(8*bessIn));
            if(bessIn<bessThresh){
                if(k==if_start){
                    lognosasiL2=TMath::Log(b)+(-b*(sasipower[k]+nosasiparampower[k]))+log(TMath::BesselI0(bessIn));
                    nosasiL2=exp(lognosasiL2);
                }
                else{
                    lognosasiL2=TMath::Log(b)+(-b*(sasipower[k]+nosasiparampower[k]))+log(TMath::BesselI0(bessIn));
                    nosasiL2*=exp(lognosasiL2);
                }
            }
            else{
                if(k==if_start){
                    lognosasiL2=TMath::Log(b)+(-b*(sasipower[k]+nosasiparampower[k]))+logbessApprox;
                    nosasiL2=exp(lognosasiL2);
                }
                else{
                    lognosasiL2=TMath::Log(b)+(-b*(sasipower[k]+nosasiparampower[k]))+logbessApprox;
                    nosasiL2*=exp(lognosasiL2);
                }
            }
            //if(lognosasiL2>750) cout<<"lognosasiL2 "<<lognosasiL2<<endl;

        /*
            //calculate prob distr of power spectrum at given freq
            Double_t b=(TMath::Power(numpointsk,2))/(4*sigmasqr);
            if(k==if_start){
                nosasiL2=b*exp(-b*(sasipower[k]+nosasiparampower[k]))*TMath::BesselI0(2*b*TMath::Sqrt(sasipower[k]*nosasiparampower[k]));
            }
            else{
                nosasiL2*=b*exp(-b*(sasipower[k]+nosasiparampower[k]))*TMath::BesselI0(2*b*TMath::Sqrt(sasipower[k]*nosasiparampower[k]));
            }
            */
        }
    
        //cout<<"nosasiL2 given SASI "<<nosasiL2<<endl;

        //CALCULATE LN(SASI METER) (SM) GIVEN SASI FEATURE IS ACTUALLY PRESENT
        Double_t SMgivenS=sasiLbest/nosasiL2; //sasi meter given sasi feature is actually present
        Double_t lnSMgivenS=log(SMgivenS);
        //cout<<"lnSMgivenS "<<lnSMgivenS<<endl;
        
        //CALCULATE POWER SPECTRUM OF NO SASI RUN
        Double_t *nosasipower=fftpower(numpointsk,nseventsr,dtscaled);
        //for(i=0;i<numpointsk;i++) cout<<nosasipower[i]<<endl;

        //FIT 2 PARAM TEMPLATE TO NO SASI RUN
        Double_t *twopfitreturn2=twopfit(numpointsk,timek,nosasipower,nsA,if_start,if_end);
        Double_t sasiLbest2=twopfitreturn2[0];
        Double_t abest2=twopfitreturn2[1];
        Double_t fsbest2=twopfitreturn2[2];
        //cout<<"sasiLbest2 given NO SASI "<<sasiLbest2<<" fsbest2 "<<fsbest2<<" abest2 "<<abest2<<endl;

        //SAVE SASIPARAM ARRAY AND POWERSPECTRUM FOR NO SASI MODEL FOR PLOTTING
        Double_t sasiparamgivenNS[numpointsk];
        for(i=0;i<numpointsk;i++){
            //built-in conversion from hz to rad/ms vvv
            sasiparamgivenNS[i]=(nsA[i]-n1)*(1.+abest2*TMath::Sin((2/scale)*TMath::Pi()*fsbest2*timek[i]+phase))+n1;
        }
        Double_t *sasiparamgivenNSpower=fftpower(numpointsk,sasiparamgivenNS,dtscaled);

        //CREATE NULL PARAM TEMPLATE FOR NO SASI RUN
        nsparamtotevents=0;
        for(i=0;i<numpointsk;i++){
            nosasiparam[i]=nsA[i]; //time-average event rate over time interval
            nsparamtotevents+=nosasiparam[i];
        }
        //cout<<"nsparamtotevents "<<nsparamtotevents<<endl;

        //CALCULATE POWER SPECTRUM OF NULL-PARAM TEMPLATE FOR SASI RUN
        Double_t *nosasiparampower2=fftpower(numpointsk,nosasiparam,dtscaled);

        //CALCULATE LIKELIHOOD FOR NOSASIPOWER GIVEN NOSASI HYPOTHESIS
        sigmasqr=nsparamtotevents/2;
        //loop over frequencies
        for(k=if_start;k<if_end;k++){

            //calculate prob distr of power spectrum at given freq
            Double_t b=(TMath::Power(numpointsk,2))/(4*sigmasqr);
            Double_t bessIn = 2*b*TMath::Sqrt(nosasipower[k]*nosasiparampower2[k]);
            //if(nosasipower[k]>1e10) cout<<"nosasipower["<<k<<"] "<<nosasipower[k]<<endl;
            //if(bessIn>1e4) cout<<"bessIn "<<bessIn<<endl;
            Double_t logbessApprox = bessIn - TMath::Log(TMath::Sqrt(2*TMath::Pi()*bessIn)) + 
                                TMath::Log(1+1/(8*bessIn));
            if(bessIn<bessThresh){
                if(k==if_start){
                    lognosasiL2=TMath::Log(b)+(-b*(nosasipower[k]+nosasiparampower2[k]))+log(TMath::BesselI0(bessIn));
                    nosasiL2=exp(lognosasiL2);
                }
                else{
                    lognosasiL2=TMath::Log(b)+(-b*(nosasipower[k]+nosasiparampower2[k]))+log(TMath::BesselI0(bessIn));
                    nosasiL2*=exp(lognosasiL2);
                }
            }
            else{
                if(k==if_start){
                    lognosasiL2=TMath::Log(b)+(-b*(nosasipower[k]+nosasiparampower2[k]))+logbessApprox;
                    nosasiL2=exp(lognosasiL2);
                }
                else{
                    lognosasiL2=TMath::Log(b)+(-b*(nosasipower[k]+nosasiparampower2[k]))+logbessApprox;
                    nosasiL2*=exp(lognosasiL2);
                }
            }
        }
        //cout<<"nosasiL2 given NO SASI "<<nosasiL2<<endl;
        //CALCULATE LN(SASI METER) (SM) GIVEN SASI FEATURE IS NOT PRESENT
        Double_t SMgivenNS=sasiLbest2/nosasiL2; //sasi meter given sasi feature is actually present
        Double_t lnSMgivenNS=log(SMgivenNS);
        //cout<<"lnSMgivenNS "<<lnSMgivenNS<<endl;

        if(lnSMgivenS>maxlnSMgivenS){
            maxlnSMgivenS=lnSMgivenS;
            fsbestp=fsbest;
            abestp=abest;
            s_startS=i_startT;
            //convert time from index to ms
            tstartS=timep[s_startS-i_startp];
        }
        if(lnSMgivenNS>maxlnSMgivenNS){
            maxlnSMgivenNS=lnSMgivenNS;
            fsbestp2=fsbest2;
            abestp2=abest2;
            s_startNS=i_startT;
            tstartNS=timep[s_startNS-i_startp];
        }
    
        delete [] sasipower;
        delete [] sasiparamgivenSpower;
        delete [] nosasiparampower2;
        delete [] sasiparamgivenNSpower;
        delete [] nosasipower;
        delete [] nosasiparampower;
        delete [] twopfitreturn;
        delete [] twopfitreturn2;

        if(q==0) cout<<"0 %"<<endl;
        if(q==-1+numtimes/4) cout<<"25 %"<<endl;
        if(q==-1+numtimes/2) cout<<"50 %"<<endl;
        if(q==-1+3*numtimes/4) cout<<"75 %"<<endl;
        if(q==numtimes-1) cout<<"100 %"<<endl;  
    }//close loop over start times

    /*
    cout<<"lnSMgivenS: "<<maxlnSMgivenS
    <<" lnSMgivenNS: "<<maxlnSMgivenNS
    <<" fsbest: "<<fsbestp
    <<" abset: "<<abestp
    <<" fsbest2: "<<fsbestp2
    <<" abest2 "<<abestp2
    <<endl;
    
    cout<<"best SASI start time: "<<tstartS<<endl;
    cout<<"best SASI start index: "<<s_startS<<endl;
    cout<<"best No-SASI start time: "<<tstartNS<<endl;;
    cout<<"run "<<m<<" of "<<fftnum<<" ms duration complete"<<endl;
    cout<<"**********************"<<endl;
    */
    
    Double_t *tdurfn_return = new Double_t[8];

    tdurfn_return[0]=maxlnSMgivenS;
    tdurfn_return[1]=maxlnSMgivenNS;

    tdurfn_return[2]=fsbestp;
    tdurfn_return[3]=abestp;
    tdurfn_return[4]=fsbestp2;
    tdurfn_return[5]=abestp2;

    tdurfn_return[6]=tstartS;
    tdurfn_return[7]=tstartNS;

    return tdurfn_return;

} //close tdurfn function