#ifndef FITFIN_C
#define FITFIN_C


#include <TFile.h>
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TArrayD.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include <TGraph.h>
#include "TROOT.h"
#include "TVirtualFitter.h"
#include <fstream>
#include <iostream>
#include <array>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <tuple>


// intial peak search param
Int_t npeaks = 30;   // # of peaks to fit

Double_t fpeaks(Double_t *x, Double_t *par)
    {
    Double_t result = par[0] + par[1] * x[0];
    for (Int_t p = 0; p < npeaks; p++) {
        Double_t norm = par[3 * p + 2]; // "height" or "area"
        Double_t mean = par[3 * p + 3];
        Double_t sigma = par[3 * p + 4];
    #if defined(__PEAKS_C_FIT_AREAS__)
        norm /= sigma * (TMath::Sqrt(TMath::TwoPi())); // "area"
    #endif                                              /* defined(__PEAKS_C_FIT_AREAS__) */
        result += norm * TMath::Gaus(x[0], mean, sigma);
    }
    return result;
    }


// to scale bins by energy values
 TF1* Poly(Double_t *x, Int_t n, Double_t *E){ //x -> peak positions, n-> #peaks, E-> known energies of peaks
    TGraph *g = new TGraph(n, x, E);

    g->SetTitle("Energy calibration;Channel;Energy");

    TF1 *cal = new TF1("cal","pol1");   // start linear

    //fitting graph
    g->Fit(cal,"Q");

    //new canvas
    TCanvas *c4 = new TCanvas("c4","Calibration");
    g->Draw("AP");

    //returning cal for later use
    return cal;
 }

 
//to read in and bound initial file
void read(string filepath, int l_cut,int u_cut) {

    string line;
    ifstream f(filepath);

    if (!f.is_open()) {
        cerr << "Error opening file: " << filepath << endl;
        exit(EXIT_FAILURE);
        return;
    }

    // Find marker & check if existent 
    bool foundData = false;
    while(getline(f,line)){
        if(line.find("$DATA") != string::npos){
            foundData = true;
            break;
        }
    }
    
    if (!foundData) {
        cerr << "Error: $DATA marker not found" << endl;
        return;
    }

    int low, high;
    f >> low >> high;

    //////////////////////////TEMP DEBUGG///////////////////////
    // Validate bin count
    if (high <= low) {
        cerr << "Error: high (" << high << ") <= low (" << low << ")" << endl;
        return;
    }
    /////////////////////////////////////////////////////////////


    cout<< "CHECKPOINT 1" << endl;
    double counts;
    
    //for initialization & bounds protection
    if (l_cut< 0 || l_cut < low) l_cut = low;
    if (u_cut < 0 || u_cut > high) u_cut = high;

    int ch = 0;
    int nbins = u_cut - l_cut + 1; 

    TH1F *h = new TH1F("h","MAESTRO Spectrum",nbins,low,high);
    
    while(f >> counts){
        if(ch >= l_cut && ch <= u_cut){
            h->SetBinContent(ch-l_cut+1, counts);
        }
        ch++;
    }
    cout<< "CHECKPOINT 2" << endl;

    h->Draw();
    TFile c("spectra.root","RECREATE");
    h->Write();
    c.Close();
    cout<< "CHECKPOINT 3" << endl;

}



//to calibrate the histogram as required
TH1F* calibration(TH1F* h, Double_t offset,Double_t gain, Double_t xmin,Double_t xmax) { //-> original histogram, new offset, new gain (per bin)
    //changing axis limits and getting bin # same
    TH1F* hcal = new TH1F("hcal","calibrated",
                    h->GetNbinsX(),
                    xmin*gain+offset,
                    xmax*gain+offset);

    //filling new bins with counts as prior stored 
    for(int b=1; b<=h->GetNbinsX(); ++b){
        hcal->SetBinContent(b, h->GetBinContent(b));
        hcal->SetBinError(b, h->GetBinError(b));
    }

    return hcal;

}

 //to allow for mult returns from fitter func
 struct Fitresults{
    TH1F* hist;
    Double_t *xp;
    Double_t sig;
    Double_t thsh;
    vector<tuple<double,double,double>> glob;
 };


Fitresults Fitter(TH1F* h, Double_t sigma, Double_t thresh){
    //cleaning canvases
    gROOT->GetListOfCanvases(); 
    delete gROOT->FindObject("c2"); 

    // canvas to see :)
    TCanvas *c2 = new TCanvas("c2", "Calibration", 1200, 800);
    c2->Divide(1, 2);
    c2->cd(1);
    h->Draw();

    struct Fitresults r;

    vector<tuple<double,double,double>> glob;

    // Get histogram range
    Double_t xmin = h->GetXaxis()->GetXmin(); //temp static 
    Double_t xmax = h->GetXaxis()->GetXmax();

    ///////////////////////////////////////////////////////////////////
    // Initialize parameters array (adjust size as needed)
    const int max_peaks = 3; //physical peak limit
    Double_t par[2 + 3*max_peaks];
    /////////////////////////////////////////////////////////////////////

   // Use TSpectrum to find the peak candidates
   TSpectrum *s = new TSpectrum(2 * npeaks); //# of peaks to search for
   //important to properly set sigma/thresh here
   Int_t nfound = s->Search(h, sigma, "", thresh); // const TH1 *hist, sigma of peaks searched for, Option_t *option="", Double_t threshold=0.05
   cout << "Found " << nfound << " candidate peaks to fit" << endl;

   //Estimate background using TSpectrum::Background
   TH1 *hb = s->Background(h, 20); //h is og histo, smoothing, name 

   //clone histogram to fit peaks above background
   TH1 *h2 = (TH1*)h->Clone("h2");
   h2->SetDirectory(0);
 
   // estimate linear background using a fitting method
   TF1 *fline = new TF1("fline", "pol1", xmin, xmax);
   h->Fit("fline", "qn");


   par[0] = fline->GetParameter(0);
   par[1] = fline->GetParameter(1);

   //resetting npeaks
   npeaks = 0;

   //loop on all found peaks Eliminate peaks at the background level
   Int_t p;
   Double_t *xpeaks;
   xpeaks = s->GetPositionX();
   for (p = 0; p < nfound; p++) {
        Double_t xp = xpeaks[p];
        Int_t bin = h->GetXaxis()->FindBin(xp);
        Double_t yp = h->GetBinContent(bin);
        if (yp - TMath::Sqrt(yp) < fline->Eval(xp))
            continue;

        // forced break for exceeding array bounds
        if (npeaks >= max_peaks) {
            cout << "Warning: Maximum number of peaks (" << max_peaks << ") reached!" << endl;
            break;
        }

        par[3 * npeaks + 2] = yp; // "height"
        par[3 * npeaks + 3] = xp; // "mean"
        //par[3 * npeaks + 4] = 5; 
        par[3*npeaks+4] = 3. * h->GetXaxis()->GetBinWidth(bin); // "new sigma"

    #if defined(__PEAKS_C_FIT_AREAS__)
        par[3 * npeaks + 2] *= par[3 * npeaks + 4] * (TMath::Sqrt(TMath::TwoPi())); // "area"
    #endif                                                                            /* defined(__PEAKS_C_FIT_AREAS__) */
        npeaks++;
    }

    cout << "Found " << npeaks << " useful peaks to fit" <<endl;


    if(npeaks == 0) {
        cout << "No peaks found above background!" << endl;
        r.hist = h;
        r.xp = nullptr;
        r.sig = sigma;
        r.thsh = thresh;
        r.glob.clear();
        return r;
        }

    TF1 *fit = new TF1("fit", fpeaks, xmin, xmax, 2 + 3 * npeaks);

    printf("Now fitting: Be patient\n");

    // We may have more than the default 25 parameters
    TVirtualFitter::Fitter(h2, 10 + 3 * npeaks);

    // Set initial parameters
    fit->SetParameters(par);
    fit->SetNpx(1000);

    h2->Fit(fit,"RQM0");

    //extracting sigma FWHM & mean for each peak
    glob.clear();

    for(int p = 0; p < npeaks; ++p){

        double mean  = fit->GetParameter(3*p + 3);
        double sigma = fit->GetParameter(3*p + 4);
        double fwhm  = 2.35482 * sigma;

        glob.emplace_back(mean, sigma, fwhm);
    }


    //Drawing results
    c2->cd(2);
    h2->Draw();
    fit->Draw("same");

    c2->cd(1);
    h->Draw();
    c2->Update();

    //returning via struct
    r.hist = h;
    r.xp = xpeaks;
    r.sig = sigma;
    r.thsh = thresh;
    r.glob = glob;
    return r;
}


void FitFin(){

    //giving filepath to desired spectra
    string filepath;

    cout << "enter the full filepath: " << endl;

    // Clearing any potential newline
    //cin.ignore(numeric_limits<streamsize>::max(), '\n');
    getline(cin, filepath);

    char again = 'n';
    int l_cut = -1 , u_cut = -1;

    do{

    //read func
    read(filepath,l_cut,u_cut);

    //input 2
    cout << "please input a lower bin threshold: ";
    cin >> l_cut;
    cout << "please input an upper bin threshold: ";
    cin >> u_cut;

    cout << "Do you want to run the intital thresholding again? (y/n): ";
    cin >> again;

    }while(again == 'y');

    // pulling histo from previously saved root file
    TFile i("spectra.root");
    if (!i.IsOpen()) {
    cout << "Error: Could not open root file!" << endl;
    return;
    }

    TH1F *h = (TH1F*)i.Get("h");

    h->SetDirectory(0); // Detach histogram from file

    if (!h) {
    cout << "Error: Histogram 'h' not found in spectra.root!" << endl;
    return;
    }

    // get # peaks to search for
    int peaks_num;
    cout << "please input the number of peaks to search for: ";
    cin >> peaks_num;

    //user inputted search parameters
    double inpsigma = 1.0;
    double inpthresh = 1.0;

    // loop to ensure correct parameterisation
    do
    {
    
    //Input
    cout << "please input the sigma for peak searching: ";
    cin >> inpsigma;
    cout << "please input the threshold for peak searching: ";
    cin >> inpthresh;

    //running peak&gaussian fitting
    Fitter(h, inpsigma, inpthresh);

    cout << "Do you want to run the fitting again? (y/n): ";
    cin >> again;

    } while (again == 'y');

    //calling results
    Fitresults r = Fitter(h, inpsigma, inpthresh);
    
    TH1F* hfit = r.hist;
    auto xp2 = r.xp;
    auto sig =r.sig;
    auto thresh = r.thsh;

    // Extract peak positions and known energies for calibration
    TArrayD E_array(peaks_num); // intiliasing energy array

    //filling with user input
    cout << "There are : " << peaks_num << " peaks, " << "please provide their energies in keV sequentially:" << endl;
    for (Int_t i = 0; i < peaks_num; i++) {
        cout << "Energy for peak " << i + 1 << ": ";
        cin >> E_array[i];
    }

    //poly func
    TF1 *cal = Poly(xp2, peaks_num, E_array.GetArray());

    //applying calibration to histogram
    Double_t offset = cal->GetParameter(0);  //fit using a + b*x -> a
    Double_t gain   = cal->GetParameter(1); // b

    //getting axis limits
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    //visual check of calibration
    TCanvas* c14 = new TCanvas("c14","calib histo");
    h->GetXaxis()->SetLimits(xmin*gain+offset, xmax*gain+offset);
    h->Draw();
    c14->Update();

    //calibration function
    TH1F* hcal = calibration(h, offset, gain, xmin, xmax);

    //2nd fitting loop to ensure correct parameterisation
    do
    {
    //input 2
    cout << "please input a new sigma for peak searching: ";
    cin >> inpsigma;
    cout << "please input a new threshold for peak searching: ";
    cin >> inpthresh;

    //running peak&gaussian fitting
    Fitter(hcal, inpsigma, inpthresh);

    cout << "Do you want to run the fitting again? (y/n): ";
    cin >> again;

    } while (again == 'y');

    Fitresults r2 = Fitter(hcal, inpsigma, inpthresh);

    auto glob = r2.glob;

    size_t c;

    for(c = 0; c < glob.size(); c++){
        cout << "peak # " << c 
        << ". mean: " << get<0>(glob[c])
        << ". sigma: " <<get<1>(glob[c])
        << ". FWHM: " <<get<2>(glob[c]) <<endl;
    }

    cout << "\nFinito" << endl;

    return;
} 

#endif