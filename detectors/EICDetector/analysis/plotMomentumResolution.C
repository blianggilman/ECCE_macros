
#include <iostream>
#include <string>
#include <vector>
using namespace std;


//Global Variables
TFile *file = new TFile("/project/projectdirs/m3763/blianggi/out_ECCE/60218804/tracking_output/ALL_G4EICDetector_g4tracking_eval.root");
// TFile *file = new TFile("/project/projectdirs/m3763/blianggi/out_ECCE/60218804/tracking_output/G4EICDetector_1_g4tracking_eval.root");
TTree *tree = (TTree*) file->Get("tracks");

float gpx, gpy, gpz, px, py, pz;

double etavals[15] = {-3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5};
char etachars[15][25] = {"n35", "n30", "n25", "n20", "n15", "n10", "n05", "00", "05", "10", "15", "20", "25", "30", "35"};
char etabins[14][25] = {"-3.5 ;< #eta < -3.0", "-3.0 < #eta < -2.5", "-2.5 < #eta < -2.0", 
    "-2.0 < #eta < -1.5", "-1.5 < #eta < -1.0", "-1.0 < #eta < -0.5", "-0.5 < #eta < 0.0", 
    "0.0 < #eta < 0.5", "0.5 < #eta < 1.0", "1.0 < #eta < 1.5", "1.5 < #eta < 2.0", 
    "2.0 < #eta < 2.5", "2.5 < #eta < 3.0", "3.0 < #eta < 3.5"};

char pbins[10][20] = {"0 < p < 2", "2 < p < 4", "4 < p < 6", "6 < p < 8", "8 < p < 10", "10 < p < 12",
    "12 < p < 14", "14 < p < 16", "16 < p < 18", "18 < p < 20"};

//first need to initialize histogram
vector<vector<TH1F*>> mom_res_preliminary_hists(14);
vector<vector<TH1F*>> mom_res_hists(14);




//calculate the three vector momentum
float calculateP(float px, float py, float pz){
    float sum = pow(px,2) + pow(py,2) + pow(pz,2);
    return sqrt(sum);
}


//calculate eta from theta
float calculateEta(float px, float py, float pz){
    //to find theta: theta is the COM scattering angle = the angle between the z-axis and 2D momentum
    //cos theta = pz/p
    float p = calculateP(px, py, pz);
    float theta = acos(pz/p);

    //given theta, calculate eta = -ln(tan(theta/2))
    return -log(tan(theta/2.));
}

//calculate the momentum resolution
float calculateMomRes(float gpx, float gpy, float gpz, float px, float py, float pz){
    // //subtract the vectors, then take the magnitude
    // float gp = calculateP(gpx, gpy, gpz); 
    // // float p = calculateP(px, py, pz); //originally tried to divide by this
    // float num = calculateP(px-gpx, py-gpy, pz-gpz);
    // return abs(num)/abs(gp);

    //take the magnitude of each momentum vector, then subtract and take the absolute value
    float p = calculateP(px, py, pz);
    float gp = calculateP(gpx, gpy, gpz);
    float num = p-gp;
    return (num/abs(gp));
}


//make gaussian fits
Double_t** makeFits(){

    TF1* g1 = new TF1("g1", "gaus",-0.1,0.1);
    Double_t** st_dev = 0;
    st_dev = new double*[14]; //[eta range][p range] in order from min to max (neg->pos)
    for (int i=0; i<14; i++){
        st_dev[i] = new double[10];
    }
    // int etactr = 0;
    // int pctr = 0;
    
    for (int i=0; i<14; i++){
        // if (i==1 || i==2 || i== 3) continue; //CHANGE THIS AFTER TEST
        for (int j=0; j<10; j++){
            mom_res_preliminary_hists[i][j]->Fit(g1, "R");
            st_dev[i][j] = g1->GetParameter(2);
        }
        // if (i==4) break; //CHANGE THIS AFTER TEST
        
    }


    // cout << "param fits" << st_dev[0] << " " << st_dev[1] << " " << st_dev[2] << " " << st_dev[3] << " " << st_dev[4] << " " << st_dev[5]<<endl;

    for (int i=0; i<10; i++){
        cout << st_dev[0][i] << endl;
    }

    // return st_dev;
    return st_dev;


}

//make gaussian fits
Double_t** makeFits2(){

    TF1* g1 = new TF1("g1", "gaus",-0.1,0.1);
    Double_t** st_dev = 0;
    st_dev = new double*[14]; //[eta range][p range] in order from min to max (neg->pos)
    for (int i=0; i<14; i++){
        st_dev[i] = new double[10];
    }
    
    for (int i=0; i<14; i++){
        for (int j=0; j<10; j++){
            mom_res_hists[i][j]->Fit(g1, "R");
            st_dev[i][j] = g1->GetParameter(2);
        }
    }
   
    for (int i=0; i<10; i++){
        cout << st_dev[0][i] << endl;
    }

    // return st_dev;
    return st_dev;
        
}


//calculate PWG requirements
// Double_t** makeFits(){

//     TF1* g1 = new TF1("g1", "gaus",-0.1,0.1);
//     Double_t** st_dev = 0;
//     st_dev = new double*[14];
double* calculatePWGreqs(int i, double p[]){
    double A[] = {0.1, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1};
    double B[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1., 1., 1., 2., 2.};
    double* dp_p = 0;
    dp_p = new double[14];
    
    for (int j=0; j<14; j++){
        double sum = pow(A[i]*p[j],2) + pow(B[i],2);
        dp_p[j] = sqrt(sum);
    }

    return dp_p;

}


//plot SD results
void plotSD(Double_t** st_dev){
    
    //initialize array of hists
    // TGraph* momPlots_by_eta[14];
    std::vector<TGraph*> momPlots_by_eta(14);
    std::vector<TGraph*> pwg_req_eqs(14);
    const Int_t n = 10;
    double p[n] = {2.,4.,6.,8.,10.,12.,14.,16.,18.,20.};

    // for (int i=0; i<10; i++){
    //     cout << st_dev[0][i] << endl;
    // }
    

    //make eta labels
    std::vector<TLatex*> etalabels(14);
    cout << "length of eta " << sizeof(etabins) << " " << sizeof(etalabels) << endl;
    for (int i=0; i<14; i++) {
        double* dp_p = calculatePWGreqs(i, p);
        double* max;
        max = max_element(&dp_p[i], &dp_p[i]+10);

        double width = 16.;
        double height = *max*0.5;
        // cout << "MAXIMUM!!!! " << *max << ", " << height << endl;
        etalabels[i] = new TLatex(width,height,etabins[i]);
        // cout << "eta label " << i << ": " << etalabels[i] << endl;

    }


    //fill graphs and plot, also get PWG requirements
    TCanvas *c15 = new TCanvas("c15","c15",1200,900);
    c15->Divide(3,5);
    
    for(int i=0; i<14; i++){
        // cout << "iteration: " << i << endl;

        cout << st_dev[1][1] << ", " << st_dev[5][5] << endl;
        //convert st dev's into percentages
        for (int j=0; j<10; j++){
            st_dev[i][j] *= 100.0;
        }
        
        cout << st_dev[1][1] << ", " << st_dev[5][5] << endl;
        momPlots_by_eta[i] = new TGraph(n,p,st_dev[i]);

        double* dp_p = calculatePWGreqs(i, p);
        // TGraph *g2 = new TGraph(n, p, dp_p);
        pwg_req_eqs[i] = new TGraph(n,p,dp_p);

        c15->cd(i+1);
        TMultiGraph *mg = new TMultiGraph();
        mg->Add(momPlots_by_eta[i]);
        mg->Add(pwg_req_eqs[i]);
        mg->Draw("ALP");
        // momPlots_by_eta[i]->Draw("ACP");
        momPlots_by_eta[i]->SetLineColor(4);
        momPlots_by_eta[i]->SetLineWidth(4);
        momPlots_by_eta[i]->SetMarkerColor(4);
        momPlots_by_eta[i]->SetMarkerSize(1.5);
        momPlots_by_eta[i]->SetMarkerStyle(21);

        mg->GetXaxis()->SetTitle("Track p [GeV/c]");
        mg->GetYaxis()->SetTitle("#Deltap/p (%)");
        mg->GetXaxis()->CenterTitle();
        mg->GetYaxis()->CenterTitle();
        // mg->GetXaxis()->SetLabelFont(36);
        // mg->GetYaxis()->SetLabelFont(36);
        // momPlots_by_eta[i]->GetXaxis()->SetTitle("Track p [GeV/c]");
        // momPlots_by_eta[i]->GetYaxis()->SetTitle("#Deltap/p");
        // // momPlots_by_eta[i]->GetHistogram()->SetMinimum(0.);
        // // momPlots_by_eta[i]->GetHistogram()->SetMinimum(0.2);

        etalabels[i]->SetTextFont(43); etalabels[i]->SetTextSize(12);
        etalabels[i]->Draw("same");

        // pwg_req_eqs[i]->Draw("ACP");

    }

    c15->Print("plots/mom_res_SD_by_eta.pdf");


}



//plot histogram of dp/p in multiple bins of eta and p (see graph for specifics)
//eta: [-3.5,3.5] in bins of 0.5
//p: [0,10] in bins of 1 or 2? //try 1 first
void plotMomRes(int nEntries){    

    int ctr[14];

    int max[14];
    int savei[14];


    //fill histograms based on the bins
    for (int i=0; i<nEntries; i++) {
        if (i%1000000==0) cout << "jet: " << i << " out of: " << nEntries << endl;
        
        tree->GetEntry(i);
        float eta = calculateEta(px,py,pz);
        float p = calculateP(px,py,pz);
        float momRes = calculateMomRes(gpx,gpy,gpz,px,py,pz);
        

        if (i<10) cout << "p for " <<i << " is " << px << ", " << py << ", " << pz << endl;
        if (i<10) cout << "eta for " <<i << " is " << eta << endl;
        if (i<10) cout << "p for " <<i << " is " << p << endl;
        if (i<10) cout << "theta for " <<i << " is " << acos(pz/p) << endl;


        for (int i=0; i<14; i++){
            for (int j=0; j<10; j++){
                if (eta > etavals[i] && eta <= etavals[i+1]) {
                    ctr[i]++;
                    if (momRes > max[i]) {
                        max[i] = momRes;
                        savei[i] = i;
                    }
                    if (p > j*2 && p <= j*2+2) mom_res_preliminary_hists[i][j]->Fill(momRes); //does comparison have to be a double?
                }
            }
        }
        


    }

    cout << "COUNTERS" << endl;
    cout << ctr[0] << " " << ctr[1] << " " << ctr[2] << " " << ctr[3]  << " " << ctr[4] << " " << ctr[5] << " " << ctr[6] << " " << ctr[7] << endl;
    cout << ctr[8] << " " << ctr[9] << " " << ctr[10] << " " << ctr[11] << " " << ctr[12] << " " << ctr[13] << endl;
    cout << "MAX VALS PER ETA (event#: max resolution)" << endl;
    cout << savei[0] << ": " << max[0] << " " << savei[1] << ": " << max[1] << " " << savei[2] << ": " << max[2] << " " << savei[3] << ": " << max[3] << endl;
    cout << savei[4] << ": " << max[4] << " " << savei[5] << ": " << max[5] << " " << savei[6] << ": " << max[6] << " " << savei[7] << ": " << max[7] << endl;
    cout << savei[8] << ": " << max[8] << " " << savei[9] << ": " << max[9] << " " << savei[10] << ": " << max[10] << " " << savei[11] << ": " << max[11] << endl;
    cout << savei[12] << ": " << max[12] << " " << savei[13] << ": " << max[13] << endl;


    //calculate standard deviations
    Double_t** st_dev_preliminary = 0;
    st_dev_preliminary = new double*[14]; //[eta range][p range] in order from min to max (neg->pos)
    for (int i=0; i<14; i++){
        st_dev_preliminary[i] = new double[10];
        for (int j=0; j<10; j++){
            st_dev_preliminary[i][j] = mom_res_preliminary_hists[i][j]->GetStdDev();
        }
    }


    //make fits
    Double_t** st_dev = makeFits();

    cout << "compare st devs!!!!" << endl;
    for (int i=0; i<14; i++){
        for (int j=0; j<10; j++){
            cout << st_dev_preliminary[i][j] << " ?= " << st_dev[i][j] << endl;
        }
    }
    

    //redefine histograms with new binning
    for (int i=0; i<14; i++){
        for(int j=0; j<10; j++){
            char title[1024];
            sprintf(title, "mom_res_eta%s-%s_p%d-%d_rebinned", etachars[i], etachars[i+1], j*2, j*2+2);
            char label[1024];
            sprintf(label, ";dp/p, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            mom_res_hists[i][j] = new TH1F(title, label, 100, -3*st_dev[i][j], 3*st_dev[i][j]);
        }
    }

    //fill histograms with new binning
    for (int i=0; i<nEntries; i++) {
        if (i%1000000==0) cout << "jet: " << i << " out of: " << nEntries << endl;
        
        tree->GetEntry(i);
        float eta = calculateEta(px,py,pz);
        float p = calculateP(px,py,pz);
        float momRes = calculateMomRes(gpx,gpy,gpz,px,py,pz);

        for (int i=0; i<14; i++){
            for (int j=0; j<10; j++){
                if (eta > etavals[i] && eta <= etavals[i+1]) {
                    if (p > j*2 && p <= j*2+2) mom_res_hists[i][j]->Fill(momRes); //does comparison have to be a double?
                }
            }
        }
    }

    Double_t** st_dev_2 = makeFits2();




    //Study this more later - why isn't it showing up??
    // TF1* g1 = new TF1("g1", "gaus",-0.1,0.1);
    // mom_res_eta0005_p0708->Fit(g1, "R"); 
    // cout << "PARAM FOR ONE THATS NOT SHOWING UP" << g1->GetParameter(2) << endl;


    //make canvas to draw
    vector<TCanvas*> canvas(14);
    for (int i=0; i<14; i++){
        char canvas_title[1024];
        sprintf(canvas_title, "c%d", i);
        canvas[i] = new TCanvas(canvas_title,canvas_title,1200,900);
        canvas[i]->Divide(5,2);

    }


    // // TCanvas c1 = makeCanvas();
    // TCanvas *c1 = new TCanvas("c1","c1",1200,900);
    // c1->Divide(5,2);
    // // c1->SetTitle("Jet Energy Drop by p");
    // // gPad->SetLogy();

    // // String outname = "mom_res_eta"

    for (int i=0; i<14; i++){
        for (int j=0; j<10; j++){
            canvas[i]->cd(j+1);
            // mom_res_preliminary_hists[i][j]->Draw();
            mom_res_hists[i][j]->Draw();

        }
        char outname[1024];
        sprintf(outname, "plots/mom_res_etan%s-%s_bins100_rebinned.pdf", etachars[i], etachars[i+1]);
        canvas[i]->Print(outname);
    }


    cout << "line 1109" << endl;

    //plot SD
    // for (int i=0; i<10; i++){
    //     cout << st_dev[0][i] << endl;
    // }
    plotSD(st_dev_2);

}



//THIS DOESN'T WORK
TCanvas* makeCanvas(){
    //p canvas
    TCanvas *c1 = new TCanvas("c1","c1",1200,900);
    c1->Divide(5,2);
    // c1->SetTitle("Jet Energy Drop by p");
    gPad->SetLogy();

    return c1;
}



//combine
void plotMomentumResolution(){
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    for (int i=0; i<14; i++){
        mom_res_preliminary_hists[i] = vector<TH1F*>(10);
        mom_res_hists[i] = vector<TH1F*>(10);

        for(int j=0; j<10; j++){
            // cout << "iteration: " << i << ", " << j << endl;

            // char title[] = ("mom_res_eta%c-%c_p%d-%d", etachars[i], etachars[i+1], j*2, j*2+2);
            // char label[] = (";dp/p, %d < #eta < %d, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            char title[1024];
            sprintf(title, "mom_res_eta%s-%s_p%d-%d", etachars[i], etachars[i+1], j*2, j*2+2);
            char label[1024];
            sprintf(label, ";dp/p, %f < #eta < %f, %d < p < %d", etavals[i], etavals[i+1], j*2, j*2+2);
            mom_res_preliminary_hists[i][j] = new TH1F(title,label,100,-1,1);

        }
    }


    tree->SetBranchAddress("gpx", &gpx);
    tree->SetBranchAddress("gpy", &gpy);
    tree->SetBranchAddress("gpz", &gpz);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);

    
    const int nEntries = tree->GetEntries();

    plotMomRes(nEntries); 

}