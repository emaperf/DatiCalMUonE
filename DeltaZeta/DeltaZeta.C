#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"

#include "TRandom3.h"
//noise parameters
const double noise=25.;
const double fixedTerm=0.001;

TRandom *tGen;

//Theory and experiment parameters
const double r_ = 0.999996345272132;
const double eMass_ = 0.510999;
const double muMass_ = 105.658;
const double beamEnergy_ = 150000.;

const double minTheta = 0.;

// Crystal parameters
const double innerSize_ = 25. + 0.5*2; //lato+separazione (mm)
const double outerSize_ = 97.8 + 0.5*2;
// Position measurement parameters
const int showerSize_ = 2;
const double w0_ = 4.;

//Fit parameters
const double a=-0.0002244; //GeV
const double b=0.04133;//MeV
const double c=-1.33;

double DelZeta (double r, double Theta, double deltaTh){
		double DelZetaVal=(-r/pow(TMath::Sin(Theta),2))*deltaTh;
		return DelZetaVal;
		}

void DeltaZeta(TString file="/lustre/cmswork/mpresill/MUonE/Fall2019/depth23/ntuD23T20E130.root", double maxTheta=5., int nEvents=-1)
{

TString suffix=file(file.Index ("ntu")+3, file.Index(".root")-file.Index("ntu")-3);

	cout<<"---You choose:"<<endl;
	cout<<"---- file = "<<file<<endl;
	cout<<"---- nEvents = "<<nEvents<<endl;
	cout<<"---- suffix = "<<suffix<<endl;

	gStyle->SetOptStat(0);
    	gStyle->SetOptFit(1);

	TFile *f=new TFile (file); //apertura file

	cout<<"--- File opened"<<endl;

	TTree *trTracks = (TTree*)f->Get("tracks"); //cerca il tree tracks
	TTree *trHits   = (TTree*)f->Get("hits");
	TTree *trCalo   = (TTree*)f->Get("calo");
	TTree *trConfig = (TTree*)f->Get("config");

	cout<<"--- Tree read"<<endl;

	//Read config
	double cCalZentry;
	trConfig->SetBranchAddress("calZentry", &cCalZentry); //punto fisso impatto calorimetro

	// Get n events, max tracks/hits/calo
	int nTracks;
	trTracks->SetBranchAddress("nTracks", &nTracks);

	if (nEvents<0) nEvents = trTracks->GetEntries();
	if (nEvents>trTracks->GetEntries()) nEvents = trTracks->GetEntries();

	//this numebrs are hadcoded in the G4 simulation
	int maxTracks=500;
	int maxcrystals=5000;

	// Tracks
    	vector<int> t_PDG, t_parentPDG, t_trackID;
    	vector<double> t_kinEnergy, t_xVertex, t_yVertex, t_zVertex, t_theta;

    	t_PDG.resize(maxTracks, 0);
    	t_parentPDG.resize(maxTracks, 0);
    	t_trackID.resize(maxTracks, 0);
    	t_kinEnergy.resize(maxTracks, 0);
    	t_xVertex.resize(maxTracks, 0);
    	t_yVertex.resize(maxTracks, 0);
    	t_zVertex.resize(maxTracks, 0);
    	t_theta.resize(maxTracks, 0);

    	trTracks->SetBranchAddress("PDG", &(t_PDG[0]));
    	trTracks->SetBranchAddress("parentPDG", &(t_parentPDG[0]));
    	trTracks->SetBranchAddress("trackID", &(t_trackID[0]));
    	trTracks->SetBranchAddress("kinEnergy", &(t_kinEnergy[0]));
    	trTracks->SetBranchAddress("xVertex", &(t_xVertex[0]));
    	trTracks->SetBranchAddress("yVertex", &(t_yVertex[0]));
    	trTracks->SetBranchAddress("zVertex", &(t_zVertex[0]));//z impatto target
    	trTracks->SetBranchAddress("theta", &(t_theta[0]));


	// Calo
    	float c_totalEnergy;
    	int c_ncrystals;
    	vector<double> c_deposit, c_calposx, c_calposy;
    	vector<int> c_crystalType;

    	c_deposit.resize(maxcrystals, 0);
    	c_calposx.resize(maxcrystals, 0);
    	c_calposy.resize(maxcrystals, 0);
    	c_crystalType.resize(maxcrystals, 0);

    	trCalo->SetBranchAddress("totalEnergy", &c_totalEnergy);//somma energia tutti i cristalli
    	trCalo->SetBranchAddress("ncrystals", &c_ncrystals);
    	trCalo->SetBranchAddress("crystalType", &(c_crystalType[0]));//interno o esterno
    	trCalo->SetBranchAddress("deposit", &(c_deposit[0]));//energia depositata su ciascun cristallo
    	trCalo->SetBranchAddress("calposx", &(c_calposx[0]));
    	trCalo->SetBranchAddress("calposy", &(c_calposy[0]));


	cout<<"--- Branches setted"<<endl;

    	trConfig->GetEvent(0);

	tGen=new TRandom3();
	 // Define result histograms
    	TH1D *measuredAngle = new TH1D("measuredAngle", "#theta_{calo} "+suffix+";#theta [mrad]", 50, minTheta, maxTheta);
    	TH1D *angleResolution = new TH1D("angleResolution", "#delta(#theta_{calo}) "+suffix+";#delta(#theta_{calo}) [mrad]", 50, -(maxTheta-minTheta)/2, +(maxTheta-minTheta)/2);
    	TH1D *measuredEnergy = new TH1D("measuredEnergy", "E_{loss} "+suffix+";E_{loss}", 150, 0, 150);
    	TH1D *measuredEnFraction = new TH1D("measuredEnFraction", "E_{loss}/E_{gen} "+suffix+";E_{loss}/E_{gen}", 100, 0.5, 1.2);
    	TH1D *esumCalo=new TH1D ("esumCalo", "esumCalo", 150, 0, 150); //isto noise matrice
    	TH1D *esumCaloFraction=new TH1D ("esumCaloFraction", "esumCaloFraction", 100, 0.5, 1.2);//isto noise matrice
    	TH1D *cDeposit=new TH1D("cDeposit", "E_noise", 150, 0, 150); //isto noise cristalli interni
    	TH1D *cDepositFraction=new TH1D("cDepositFraction", "E_noise", 100, 0.5, 1.2); //isto noise cristalli interni

    	//istogrammi cristalli
    	TH1D *nCryst=new TH1D("nCryst", "Cristalli illuminati (noise)", 50, 0, 300);




    	// Define utility variables
    	int nShowerCrystals = pow(2*showerSize_ + 1, 2);
    
    	// Events loop begins here
    	cout<<"--- Event loop begins"<<endl;

	int ngood = 0;
	int ninner = 0;
	int nouter = 0;
	int ninneredge = 0;
	int nouteredge = 0;


	for(int i=0; i<nEvents; ++i){

	trTracks->GetEvent(i);
	trCalo->GetEvent(i);

	if(i%1000 == 0) cout<<"event "<<i<<endl;

	// Event Selection and Primary Tracks
        int eleIndex = -1;

	for(int j=0; j<nTracks; ++j){
		if(t_trackID[j]==2){  //ID ele
		   eleIndex = j;
		   break;
	    }
	}
	if(eleIndex == -1) continue;
        if(t_PDG[eleIndex] != 11) continue; //PDG ele=11

	ngood++;
	
	//smearing
	double nTotcrys=0;
	double Esum=0;
	double EsumFrac=0;
	for(int j=0; j<c_ncrystals; j++)
	{
		double s1=tGen->Gaus()*noise;
		double s2=tGen->Gaus()*fixedTerm*c_deposit[j];

		if(c_crystalType[j]!=1) continue;
		c_deposit[j]+=(s1+s2);

		if(c_deposit[j]>10) nTotcrys+=1;
	
		Esum+=c_deposit[j];
		EsumFrac+=(c_deposit[j]/t_kinEnergy[eleIndex]);
		
	}

	cDeposit->Fill(Esum/1000);
	cDepositFraction->Fill(EsumFrac);
	nCryst->Fill(nTotcrys);

	// --- ELOSS ANALYSIS
        measuredEnergy->Fill(c_totalEnergy/1000); // convert MeV in GeV
        measuredEnFraction->Fill(c_totalEnergy/t_kinEnergy[eleIndex]);

	 // --- CRYSTAL ANALYSIS
        double maxDeposit = 0;
        int maxCrystal = -1;

	// Search for crystal with the largest energy deposit
	for (int j=0; j<c_ncrystals; ++j){
		if(c_crystalType[j] == 0) continue;
		double en = c_deposit[j];
		if(en>maxDeposit){
		     maxDeposit = en;
		     maxCrystal = j;
		}
	}

	
	double maxCalX = c_calposx[maxCrystal];
        double maxCalY = c_calposy[maxCrystal];
        int    maxType = c_crystalType[maxCrystal]; // 1 == inner calo, 2 == outer calo
	
	if(maxType == 1) ninner++;
        if(maxType == 2) nouter++;

	// Define shower
        double size = (maxType == 1 ? innerSize_ : outerSize_);
        double showerSizeMM = (float)showerSize_*size + size*0.05; // add 5% tolerance

        // Measure total shower energy
        double xMean = 0.;
        double yMean = 0.;
        double norm  = 0.;
        double etot = 0.;
        int    nshower = 0; //cristalli matrice
        vector<int> showerCrystals; //vettore indici cristalli matrice

	
        for(int j=0; j<c_ncrystals; ++j){
            
	    if(c_deposit[j] == 0) continue;
            if(c_crystalType[j] == 0) continue;
            
	    if(abs(c_calposx[j] - maxCalX) > showerSizeMM ) continue;
            if(abs(c_calposy[j] - maxCalY) > showerSizeMM ) continue;
            nshower++;
            etot += c_deposit[j];
            showerCrystals.push_back(j); // store shower crystals indices
        }
	
	esumCalo->Fill(etot/1000);
        esumCaloFraction->Fill(etot/t_kinEnergy[eleIndex]);

	 if(nshower < nShowerCrystals){
            
            ninneredge++;
        }

        if(nshower > nShowerCrystals){
           
            nouteredge++;
        }

	for(auto it:showerCrystals){   //scorre gli indici crist matrice
            double en = c_deposit[it]; //en depositata sul crist con indice it
            double w = TMath::Max(0., w0_ + log(en/etot)); //peso
            norm  += w;
            xMean += w*c_calposx[it];
            yMean += w*c_calposy[it];
				
	}

	xMean /= norm;
        yMean /= norm;

	// Compute theta calo
        double r = sqrt(pow(xMean,2)+pow(yMean,2));
        double dist = cCalZentry - t_zVertex[eleIndex]; //zimpatto ele su target
        double theta = TMath::ATan(r/dist);
        double reso = (theta-t_theta[eleIndex]);
	
	//Compute delta zeta        
	double theta_v=t_theta[eleIndex]*1000;
	double deltaTheta=a*pow((etot/1000),2)+b*(etot/1000)+c;
	double deltaZeta=DelZeta (r, theta_v, deltaTheta);
	
	
	cout<<"norm: "<<norm<<endl;
	cout<<"xMean: "<<xMean<<" "<<"yMean: "<<yMean<<endl;
	cout<<"r "<<r<<endl;
	cout<<"Delta Zeta: "<<deltaZeta<<endl;
	cout<<endl;
	
	double theta_del = TMath::ATan(r/(dist-deltaZeta));
        double reso_del = (theta_del-t_theta[eleIndex]);

	measuredAngle->Fill(theta_del*1000);
        angleResolution->Fill(reso_del*1000);

	}

	cout<<endl;
    	cout<<"ngood "<<ngood<<endl;
    	cout<<"ninner "<<ninner<<endl;
    	cout<<"ninneredge "<<ninneredge<<endl;
    	cout<<"nouteredge "<<nouteredge<<endl;
    	cout<<"nouter "<<nouter<<endl;


	TCanvas *c1 = new TCanvas("c1","c1",1000,600);
    	c1->Divide(2,2);
    	c1->cd(1);
    	gStyle->SetOptStat("emr");
    	gPad->SetGrid();
    	measuredAngle->SetMarkerStyle(20);
    	measuredAngle->SetMarkerSize(.5);
    	measuredAngle->Draw("PE");
    	c1->cd(2);
    	gPad->SetGrid();
    	angleResolution->SetMarkerStyle(20);
    	angleResolution->SetMarkerSize(.5);
    	angleResolution->Draw("PE");


	measuredAngle->Fit("gaus", "Q");
    	angleResolution->Fit("gaus", "Q");

	
	c1->cd(3);
    	gPad->SetGrid();
    	measuredEnergy->SetMarkerStyle(20);
    	measuredEnergy->SetMarkerSize(.5);
    	
	double meanE=measuredEnergy->GetMean();
    	double sigmaE=measuredEnergy->GetStdDev();
    	cout<<"MeanE= "<<meanE<<" "<<"StdE= "<<sigmaE<<endl;
    	
	measuredEnergy->GetXaxis()->SetRangeUser(meanE-3.5*sigmaE, meanE+2*sigmaE);
    	measuredEnergy->Draw("PE");
    	
	esumCalo->SetMarkerStyle(20);
    	esumCalo->SetMarkerSize(.5);
    	esumCalo->SetMarkerColor(2);
    	esumCalo->SetLineColor(2);
    	esumCalo->Draw("sames");
    	
	cDeposit->SetMarkerStyle(20);
    	cDeposit->SetMarkerSize(.5);
    	cDeposit->SetMarkerColor(3);
    	cDeposit->SetLineColor(3);
    	cDeposit->Draw("sames");


	c1->cd(4);
    	gPad->SetGrid();
    	measuredEnFraction->SetMarkerStyle(20);
    	measuredEnFraction->SetMarkerSize(.5);
    	
	double meanErel=measuredEnFraction->GetMean();
    	double sigmaErel=measuredEnFraction->GetStdDev();
    	cout<<"MeanErel= "<<meanErel<<" "<<"StdErel= "<<sigmaErel<<endl;
    	
	measuredEnFraction->GetXaxis()->SetRangeUser(meanErel-3.5*sigmaErel, meanErel+2*sigmaErel);
    	measuredEnFraction->Draw("PE");
    	esumCaloFraction->SetMarkerStyle(20);
    	
	esumCaloFraction->SetMarkerSize(.5);
    	esumCaloFraction->SetMarkerColor(2);
    	esumCaloFraction->SetLineColor(2);
    	esumCaloFraction->Draw("sames");

	cDepositFraction->SetMarkerStyle(20);
	cDepositFraction->SetMarkerSize(.5);
	cDepositFraction->SetMarkerColor(3);
	cDepositFraction->SetLineColor(3);
	cDepositFraction->Draw("sames");


	TCanvas *c2= new TCanvas("c2", "c2", 400, 300);
    	gPad->SetGrid();
    	
	double meanCrys=nCryst->GetMean();
    	double sigmaCrys=nCryst->GetStdDev();
    	
	nCryst->GetXaxis()->SetRangeUser(meanCrys-3.5*sigmaCrys, meanCrys+3.5*sigmaCrys);

    	nCryst->Draw();



	return;

}
















