


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
    return -log(tan(theta/2));
}

//calculate the momentum resolution
float calculateMomRes(float gpx, float gpy, float gpz, float px, float py, float pz){
    //subtract the vectors, then take the magnitude
    float gp = calculateP(gpx, gpy, gpz);
    float num = calculateP(px-gpx, py-gpy, pz-gpz);
    return abs(num)/abs(gp);

    // //take the magnitude of each momentum vector, then subtract and take the absolute value
    // float p = calculateP(px, py, pz);
    // float gp = calculateP(gpx, gpy, gpz);
    // float num = abs(p-gp);
    // return num/abs(gp);
}

//plot histogram of dp/p in multiple bins of eta and p (see graph for specifics)
//eta: [-3.5,3.5] in bins of 0.5
//pt: [0,10] in bins of 1 or 2? //try 1 first
void plotMomRes(TTree* tree, int nEntries, float gpx, float gpy, float gpz, float px, float py, float pz){
        
    //first need to initialize histogram
    TH1F *mom_res_etan3530_pt0001 = new TH1F("mom_res_etan3530_pt0001",";Momentum Resolution, -3.5 < #eta < -3, 0 < pT < 1",25,0,1);
    TH1F *mom_res_etan3530_pt0102 = new TH1F("mom_res_etan3530_pt0102",";Momentum Resolution, -3.5 < #eta < -3, 1 < pT < 2",25,0,1);
    TH1F *mom_res_etan3530_pt0203 = new TH1F("mom_res_etan3530_pt0203",";Momentum Resolution, -3.5 < #eta < -3, 2 < pT < 3",25,0,1);
    TH1F *mom_res_etan3530_pt0304 = new TH1F("mom_res_etan3530_pt0304",";Momentum Resolution, -3.5 < #eta < -3, 3 < pT < 4",25,0,1);
    TH1F *mom_res_etan3530_pt0405 = new TH1F("mom_res_etan3530_pt0405",";Momentum Resolution, -3.5 < #eta < -3, 4 < pT < 5",25,0,1);
    TH1F *mom_res_etan3530_pt0506 = new TH1F("mom_res_etan3530_pt0506",";Momentum Resolution, -3.5 < #eta < -3, 5 < pT < 6",25,0,1);
    TH1F *mom_res_etan3530_pt0607 = new TH1F("mom_res_etan3530_pt0607",";Momentum Resolution, -3.5 < #eta < -3, 6 < pT < 7",25,0,1);
    TH1F *mom_res_etan3530_pt0708 = new TH1F("mom_res_etan3530_pt0708",";Momentum Resolution, -3.5 < #eta < -3, 7 < pT < 8",25,0,1);
    TH1F *mom_res_etan3530_pt0809 = new TH1F("mom_res_etan3530_pt0809",";Momentum Resolution, -3.5 < #eta < -3, 8 < pT < 9",25,0,1);
    TH1F *mom_res_etan3530_pt0910 = new TH1F("mom_res_etan3530_pt0910",";Momentum Resolution, -3.5 < #eta < -3, 9 < pT < 10",25,0,1);
    
    TH1F *mom_res_etan3025_pt0001 = new TH1F("mom_res_etan3025_pt0001",";Momentum Resolution, -3 < #eta < -2.5, 0 < pT < 1",25,0,1);
    TH1F *mom_res_etan3025_pt0102 = new TH1F("mom_res_etan3025_pt0102",";Momentum Resolution, -3 < #eta < -2.5, 1 < pT < 2",25,0,1);
    TH1F *mom_res_etan3025_pt0203 = new TH1F("mom_res_etan3025_pt0203",";Momentum Resolution, -3 < #eta < -2.5, 2 < pT < 3",25,0,1);
    TH1F *mom_res_etan3025_pt0304 = new TH1F("mom_res_etan3025_pt0304",";Momentum Resolution, -3 < #eta < -2.5, 3 < pT < 4",25,0,1);
    TH1F *mom_res_etan3025_pt0405 = new TH1F("mom_res_etan3025_pt0405",";Momentum Resolution, -3 < #eta < -2.5, 4 < pT < 5",25,0,1);
    TH1F *mom_res_etan3025_pt0506 = new TH1F("mom_res_etan3025_pt0506",";Momentum Resolution, -3 < #eta < -2.5, 5 < pT < 6",25,0,1);
    TH1F *mom_res_etan3025_pt0607 = new TH1F("mom_res_etan3025_pt0607",";Momentum Resolution, -3 < #eta < -2.5, 6 < pT < 7",25,0,1);
    TH1F *mom_res_etan3025_pt0708 = new TH1F("mom_res_etan3025_pt0708",";Momentum Resolution, -3 < #eta < -2.5, 7 < pT < 8",25,0,1);
    TH1F *mom_res_etan3025_pt0809 = new TH1F("mom_res_etan3025_pt0809",";Momentum Resolution, -3 < #eta < -2.5, 8 < pT < 9",25,0,1);
    TH1F *mom_res_etan3025_pt0910 = new TH1F("mom_res_etan3025_pt0910",";Momentum Resolution, -3 < #eta < -2.5, 9 < pT < 10",25,0,1);
    
    TH1F *mom_res_etan2520_pt0001 = new TH1F("mom_res_etan2520_pt0001",";Momentum Resolution, -2.5 < #eta < -2, 0 < pT < 1",25,0,1);
    TH1F *mom_res_etan2520_pt0102 = new TH1F("mom_res_etan2520_pt0102",";Momentum Resolution, -2.5 < #eta < -2, 1 < pT < 2",25,0,1);
    TH1F *mom_res_etan2520_pt0203 = new TH1F("mom_res_etan2520_pt0203",";Momentum Resolution, -2.5 < #eta < -2, 2 < pT < 3",25,0,1);
    TH1F *mom_res_etan2520_pt0304 = new TH1F("mom_res_etan2520_pt0304",";Momentum Resolution, -2.5 < #eta < -2, 3 < pT < 4",25,0,1);
    TH1F *mom_res_etan2520_pt0405 = new TH1F("mom_res_etan2520_pt0405",";Momentum Resolution, -2.5 < #eta < -2, 4 < pT < 5",25,0,1);
    TH1F *mom_res_etan2520_pt0506 = new TH1F("mom_res_etan2520_pt0506",";Momentum Resolution, -2.5 < #eta < -2, 5 < pT < 6",25,0,1);
    TH1F *mom_res_etan2520_pt0607 = new TH1F("mom_res_etan2520_pt0607",";Momentum Resolution, -2.5 < #eta < -2, 6 < pT < 7",25,0,1);
    TH1F *mom_res_etan2520_pt0708 = new TH1F("mom_res_etan2520_pt0708",";Momentum Resolution, -2.5 < #eta < -2, 7 < pT < 8",25,0,1);
    TH1F *mom_res_etan2520_pt0809 = new TH1F("mom_res_etan2520_pt0809",";Momentum Resolution, -2.5 < #eta < -2, 8 < pT < 9",25,0,1);
    TH1F *mom_res_etan2520_pt0910 = new TH1F("mom_res_etan2520_pt0910",";Momentum Resolution, -2.5 < #eta < -2, 9 < pT < 10",25,0,1);
    
    TH1F *mom_res_etan2015_pt0001 = new TH1F("mom_res_etan2015_pt0001",";Momentum Resolution, -2 < #eta < -1.5, 0 < pT < 1",25,0,1);
    TH1F *mom_res_etan2015_pt0102 = new TH1F("mom_res_etan2015_pt0102",";Momentum Resolution, -2 < #eta < -1.5, 1 < pT < 2",25,0,1);
    TH1F *mom_res_etan2015_pt0203 = new TH1F("mom_res_etan2015_pt0203",";Momentum Resolution, -2 < #eta < -1.5, 2 < pT < 3",25,0,1);
    TH1F *mom_res_etan2015_pt0304 = new TH1F("mom_res_etan2015_pt0304",";Momentum Resolution, -2 < #eta < -1.5, 3 < pT < 4",25,0,1);
    TH1F *mom_res_etan2015_pt0405 = new TH1F("mom_res_etan2015_pt0405",";Momentum Resolution, -2 < #eta < -1.5, 4 < pT < 5",25,0,1);
    TH1F *mom_res_etan2015_pt0506 = new TH1F("mom_res_etan2015_pt0506",";Momentum Resolution, -2 < #eta < -1.5, 5 < pT < 6",25,0,1);
    TH1F *mom_res_etan2015_pt0607 = new TH1F("mom_res_etan2015_pt0607",";Momentum Resolution, -2 < #eta < -1.5, 6 < pT < 7",25,0,1);
    TH1F *mom_res_etan2015_pt0708 = new TH1F("mom_res_etan2015_pt0708",";Momentum Resolution, -2 < #eta < -1.5, 7 < pT < 8",25,0,1);
    TH1F *mom_res_etan2015_pt0809 = new TH1F("mom_res_etan2015_pt0809",";Momentum Resolution, -2 < #eta < -1.5, 8 < pT < 9",25,0,1);
    TH1F *mom_res_etan2015_pt0910 = new TH1F("mom_res_etan2015_pt0910",";Momentum Resolution, -2 < #eta < -1.5, 9 < pT < 10",25,0,1);
    
    TH1F *mom_res_etan1510_pt0001 = new TH1F("mom_res_etan1510_pt0001",";Momentum Resolution, -1.5 < #eta < -1, 0 < pT < 1",25,0,1);
    TH1F *mom_res_etan1510_pt0102 = new TH1F("mom_res_etan1510_pt0102",";Momentum Resolution, -1.5 < #eta < -1, 1 < pT < 2",25,0,1);
    TH1F *mom_res_etan1510_pt0203 = new TH1F("mom_res_etan1510_pt0203",";Momentum Resolution, -1.5 < #eta < -1, 2 < pT < 3",25,0,1);
    TH1F *mom_res_etan1510_pt0304 = new TH1F("mom_res_etan1510_pt0304",";Momentum Resolution, -1.5 < #eta < -1, 3 < pT < 4",25,0,1);
    TH1F *mom_res_etan1510_pt0405 = new TH1F("mom_res_etan1510_pt0405",";Momentum Resolution, -1.5 < #eta < -1, 4 < pT < 5",25,0,1);
    TH1F *mom_res_etan1510_pt0506 = new TH1F("mom_res_etan1510_pt0506",";Momentum Resolution, -1.5 < #eta < -1, 5 < pT < 6",25,0,1);
    TH1F *mom_res_etan1510_pt0607 = new TH1F("mom_res_etan1510_pt0607",";Momentum Resolution, -1.5 < #eta < -1, 6 < pT < 7",25,0,1);
    TH1F *mom_res_etan1510_pt0708 = new TH1F("mom_res_etan1510_pt0708",";Momentum Resolution, -1.5 < #eta < -1, 7 < pT < 8",25,0,1);
    TH1F *mom_res_etan1510_pt0809 = new TH1F("mom_res_etan1510_pt0809",";Momentum Resolution, -1.5 < #eta < -1, 8 < pT < 9",25,0,1);
    TH1F *mom_res_etan1510_pt0910 = new TH1F("mom_res_etan1510_pt0910",";Momentum Resolution, -1.5 < #eta < -1, 9 < pT < 10",25,0,1);
    
    TH1F *mom_res_etan1005_pt0001 = new TH1F("mom_res_etan1005_pt0001",";Momentum Resolution, -1 < #eta < -0.5, 0 < pT < 1",25,0,1);
    TH1F *mom_res_etan1005_pt0102 = new TH1F("mom_res_etan1005_pt0102",";Momentum Resolution, -1 < #eta < -0.5, 1 < pT < 2",25,0,1);
    TH1F *mom_res_etan1005_pt0203 = new TH1F("mom_res_etan1005_pt0203",";Momentum Resolution, -1 < #eta < -0.5, 2 < pT < 3",25,0,1);
    TH1F *mom_res_etan1005_pt0304 = new TH1F("mom_res_etan1005_pt0304",";Momentum Resolution, -1 < #eta < -0.5, 3 < pT < 4",25,0,1);
    TH1F *mom_res_etan1005_pt0405 = new TH1F("mom_res_etan1005_pt0405",";Momentum Resolution, -1 < #eta < -0.5, 4 < pT < 5",25,0,1);
    TH1F *mom_res_etan1005_pt0506 = new TH1F("mom_res_etan1005_pt0506",";Momentum Resolution, -1 < #eta < -0.5, 5 < pT < 6",25,0,1);
    TH1F *mom_res_etan1005_pt0607 = new TH1F("mom_res_etan1005_pt0607",";Momentum Resolution, -1 < #eta < -0.5, 6 < pT < 7",25,0,1);
    TH1F *mom_res_etan1005_pt0708 = new TH1F("mom_res_etan1005_pt0708",";Momentum Resolution, -1 < #eta < -0.5, 7 < pT < 8",25,0,1);
    TH1F *mom_res_etan1005_pt0809 = new TH1F("mom_res_etan1005_pt0809",";Momentum Resolution, -1 < #eta < -0.5, 8 < pT < 9",25,0,1);
    TH1F *mom_res_etan1005_pt0910 = new TH1F("mom_res_etan1005_pt0910",";Momentum Resolution, -1 < #eta < -0.5, 9 < pT < 10",25,0,1);
    
    TH1F *mom_res_etan0500_pt0001 = new TH1F("mom_res_etan0500_pt0001",";Momentum Resolution, -0.5 < #eta < 0, 0 < pT < 1",25,0,1);
    TH1F *mom_res_etan0500_pt0102 = new TH1F("mom_res_etan0500_pt0102",";Momentum Resolution, -0.5 < #eta < 0, 1 < pT < 2",25,0,1);
    TH1F *mom_res_etan0500_pt0203 = new TH1F("mom_res_etan0500_pt0203",";Momentum Resolution, -0.5 < #eta < 0, 2 < pT < 3",25,0,1);
    TH1F *mom_res_etan0500_pt0304 = new TH1F("mom_res_etan0500_pt0304",";Momentum Resolution, -0.5 < #eta < 0, 3 < pT < 4",25,0,1);
    TH1F *mom_res_etan0500_pt0405 = new TH1F("mom_res_etan0500_pt0405",";Momentum Resolution, -0.5 < #eta < 0, 4 < pT < 5",25,0,1);
    TH1F *mom_res_etan0500_pt0506 = new TH1F("mom_res_etan0500_pt0506",";Momentum Resolution, -0.5 < #eta < 0, 5 < pT < 6",25,0,1);
    TH1F *mom_res_etan0500_pt0607 = new TH1F("mom_res_etan0500_pt0607",";Momentum Resolution, -0.5 < #eta < 0, 6 < pT < 7",25,0,1);
    TH1F *mom_res_etan0500_pt0708 = new TH1F("mom_res_etan0500_pt0708",";Momentum Resolution, -0.5 < #eta < 0, 7 < pT < 8",25,0,1);
    TH1F *mom_res_etan0500_pt0809 = new TH1F("mom_res_etan0500_pt0809",";Momentum Resolution, -0.5 < #eta < 0, 8 < pT < 9",25,0,1);
    TH1F *mom_res_etan0500_pt0910 = new TH1F("mom_res_etan0500_pt0910",";Momentum Resolution, -0.5 < #eta < 0, 9 < pT < 10",25,0,1);
    
    
    TH1F *mom_res_eta0005_pt0001 = new TH1F("mom_res_eta0005_pt0001",";Momentum Resolution, 0 < #eta < 0.5, 0 < pT < 1",25,0,1);
    TH1F *mom_res_eta0005_pt0102 = new TH1F("mom_res_eta0005_pt0102",";Momentum Resolution, 0 < #eta < 0.5, 1 < pT < 2",25,0,1);
    TH1F *mom_res_eta0005_pt0203 = new TH1F("mom_res_eta0005_pt0203",";Momentum Resolution, 0 < #eta < 0.5, 2 < pT < 3",25,0,1);
    TH1F *mom_res_eta0005_pt0304 = new TH1F("mom_res_eta0005_pt0304",";Momentum Resolution, 0 < #eta < 0.5, 3 < pT < 4",25,0,1);
    TH1F *mom_res_eta0005_pt0405 = new TH1F("mom_res_eta0005_pt0405",";Momentum Resolution, 0 < #eta < 0.5, 4 < pT < 5",25,0,1);
    TH1F *mom_res_eta0005_pt0506 = new TH1F("mom_res_eta0005_pt0506",";Momentum Resolution, 0 < #eta < 0.5, 5 < pT < 6",25,0,1);
    TH1F *mom_res_eta0005_pt0607 = new TH1F("mom_res_eta0005_pt0607",";Momentum Resolution, 0 < #eta < 0.5, 6 < pT < 7",25,0,1);
    TH1F *mom_res_eta0005_pt0708 = new TH1F("mom_res_eta0005_pt0708",";Momentum Resolution, 0 < #eta < 0.5, 7 < pT < 8",25,0,1);
    TH1F *mom_res_eta0005_pt0809 = new TH1F("mom_res_eta0005_pt0809",";Momentum Resolution, 0 < #eta < 0.5, 8 < pT < 9",25,0,1);
    TH1F *mom_res_eta0005_pt0910 = new TH1F("mom_res_eta0005_pt0910",";Momentum Resolution, 0 < #eta < 0.5, 9 < pT < 10",25,0,1);
    
    TH1F *mom_res_eta0510_pt0001 = new TH1F("mom_res_eta0510_pt0001",";Momentum Resolution, 0.5 < #eta < 1, 0 < pT < 1",25,0,1);
    TH1F *mom_res_eta0510_pt0102 = new TH1F("mom_res_eta0510_pt0102",";Momentum Resolution, 0.5 < #eta < 1, 1 < pT < 2",25,0,1);
    TH1F *mom_res_eta0510_pt0203 = new TH1F("mom_res_eta0510_pt0203",";Momentum Resolution, 0.5 < #eta < 1, 2 < pT < 3",25,0,1);
    TH1F *mom_res_eta0510_pt0304 = new TH1F("mom_res_eta0510_pt0304",";Momentum Resolution, 0.5 < #eta < 1, 3 < pT < 4",25,0,1);
    TH1F *mom_res_eta0510_pt0405 = new TH1F("mom_res_eta0510_pt0405",";Momentum Resolution, 0.5 < #eta < 1, 4 < pT < 5",25,0,1);
    TH1F *mom_res_eta0510_pt0506 = new TH1F("mom_res_eta0510_pt0506",";Momentum Resolution, 0.5 < #eta < 1, 5 < pT < 6",25,0,1);
    TH1F *mom_res_eta0510_pt0607 = new TH1F("mom_res_eta0510_pt0607",";Momentum Resolution, 0.5 < #eta < 1, 6 < pT < 7",25,0,1);
    TH1F *mom_res_eta0510_pt0708 = new TH1F("mom_res_eta0510_pt0708",";Momentum Resolution, 0.5 < #eta < 1, 7 < pT < 8",25,0,1);
    TH1F *mom_res_eta0510_pt0809 = new TH1F("mom_res_eta0510_pt0809",";Momentum Resolution, 0.5 < #eta < 1, 8 < pT < 9",25,0,1);
    TH1F *mom_res_eta0510_pt0910 = new TH1F("mom_res_eta0510_pt0910",";Momentum Resolution, 0.5 < #eta < 1, 9 < pT < 10",25,0,1);
    
    TH1F *mom_res_eta1015_pt0001 = new TH1F("mom_res_eta1015_pt0001",";Momentum Resolution, 1 < #eta < 1.5, 0 < pT < 1",25,0,1);
    TH1F *mom_res_eta1015_pt0102 = new TH1F("mom_res_eta1015_pt0102",";Momentum Resolution, 1 < #eta < 1.5, 1 < pT < 2",25,0,1);
    TH1F *mom_res_eta1015_pt0203 = new TH1F("mom_res_eta1015_pt0203",";Momentum Resolution, 1 < #eta < 1.5, 2 < pT < 3",25,0,1);
    TH1F *mom_res_eta1015_pt0304 = new TH1F("mom_res_eta1015_pt0304",";Momentum Resolution, 1 < #eta < 1.5, 3 < pT < 4",25,0,1);
    TH1F *mom_res_eta1015_pt0405 = new TH1F("mom_res_eta1015_pt0405",";Momentum Resolution, 1 < #eta < 1.5, 4 < pT < 5",25,0,1);
    TH1F *mom_res_eta1015_pt0506 = new TH1F("mom_res_eta1015_pt0506",";Momentum Resolution, 1 < #eta < 1.5, 5 < pT < 6",25,0,1);
    TH1F *mom_res_eta1015_pt0607 = new TH1F("mom_res_eta1015_pt0607",";Momentum Resolution, 1 < #eta < 1.5, 6 < pT < 7",25,0,1);
    TH1F *mom_res_eta1015_pt0708 = new TH1F("mom_res_eta1015_pt0708",";Momentum Resolution, 1 < #eta < 1.5, 7 < pT < 8",25,0,1);
    TH1F *mom_res_eta1015_pt0809 = new TH1F("mom_res_eta1015_pt0809",";Momentum Resolution, 1 < #eta < 1.5, 8 < pT < 9",25,0,1);
    TH1F *mom_res_eta1015_pt0910 = new TH1F("mom_res_eta1015_pt0910",";Momentum Resolution, 1 < #eta < 1.5, 9 < pT < 10",25,0,1);
    
    TH1F *mom_res_eta1520_pt0001 = new TH1F("mom_res_eta1520_pt0001",";Momentum Resolution, 1.5 < #eta < 2, 0 < pT < 1",25,0,1);
    TH1F *mom_res_eta1520_pt0102 = new TH1F("mom_res_eta1520_pt0102",";Momentum Resolution, 1.5 < #eta < 2, 1 < pT < 2",25,0,1);
    TH1F *mom_res_eta1520_pt0203 = new TH1F("mom_res_eta1520_pt0203",";Momentum Resolution, 1.5 < #eta < 2, 2 < pT < 3",25,0,1);
    TH1F *mom_res_eta1520_pt0304 = new TH1F("mom_res_eta1520_pt0304",";Momentum Resolution, 1.5 < #eta < 2, 3 < pT < 4",25,0,1);
    TH1F *mom_res_eta1520_pt0405 = new TH1F("mom_res_eta1520_pt0405",";Momentum Resolution, 1.5 < #eta < 2, 4 < pT < 5",25,0,1);
    TH1F *mom_res_eta1520_pt0506 = new TH1F("mom_res_eta1520_pt0506",";Momentum Resolution, 1.5 < #eta < 2, 5 < pT < 6",25,0,1);
    TH1F *mom_res_eta1520_pt0607 = new TH1F("mom_res_eta1520_pt0607",";Momentum Resolution, 1.5 < #eta < 2, 6 < pT < 7",25,0,1);
    TH1F *mom_res_eta1520_pt0708 = new TH1F("mom_res_eta1520_pt0708",";Momentum Resolution, 1.5 < #eta < 2, 7 < pT < 8",25,0,1);
    TH1F *mom_res_eta1520_pt0809 = new TH1F("mom_res_eta1520_pt0809",";Momentum Resolution, 1.5 < #eta < 2, 8 < pT < 9",25,0,1);
    TH1F *mom_res_eta1520_pt0910 = new TH1F("mom_res_eta1520_pt0910",";Momentum Resolution, 1.5 < #eta < 2, 9 < pT < 10",25,0,1);
    
    TH1F *mom_res_eta2025_pt0001 = new TH1F("mom_res_eta2025_pt0001",";Momentum Resolution, 2 < #eta < 2.5, 0 < pT < 1",25,0,1);
    TH1F *mom_res_eta2025_pt0102 = new TH1F("mom_res_eta2025_pt0102",";Momentum Resolution, 2 < #eta < 2.5, 1 < pT < 2",25,0,1);
    TH1F *mom_res_eta2025_pt0203 = new TH1F("mom_res_eta2025_pt0203",";Momentum Resolution, 2 < #eta < 2.5, 2 < pT < 3",25,0,1);
    TH1F *mom_res_eta2025_pt0304 = new TH1F("mom_res_eta2025_pt0304",";Momentum Resolution, 2 < #eta < 2.5, 3 < pT < 4",25,0,1);
    TH1F *mom_res_eta2025_pt0405 = new TH1F("mom_res_eta2025_pt0405",";Momentum Resolution, 2 < #eta < 2.5, 4 < pT < 5",25,0,1);
    TH1F *mom_res_eta2025_pt0506 = new TH1F("mom_res_eta2025_pt0506",";Momentum Resolution, 2 < #eta < 2.5, 5 < pT < 6",25,0,1);
    TH1F *mom_res_eta2025_pt0607 = new TH1F("mom_res_eta2025_pt0607",";Momentum Resolution, 2 < #eta < 2.5, 6 < pT < 7",25,0,1);
    TH1F *mom_res_eta2025_pt0708 = new TH1F("mom_res_eta2025_pt0708",";Momentum Resolution, 2 < #eta < 2.5, 7 < pT < 8",25,0,1);
    TH1F *mom_res_eta2025_pt0809 = new TH1F("mom_res_eta2025_pt0809",";Momentum Resolution, 2 < #eta < 2.5, 8 < pT < 9",25,0,1);
    TH1F *mom_res_eta2025_pt0910 = new TH1F("mom_res_eta2025_pt0910",";Momentum Resolution, 2 < #eta < 2.5, 9 < pT < 10",25,0,1);
    
    TH1F *mom_res_eta2530_pt0001 = new TH1F("mom_res_eta2530_pt0001",";Momentum Resolution, 2.5 < #eta < 3, 0 < pT < 1",25,0,1);
    TH1F *mom_res_eta2530_pt0102 = new TH1F("mom_res_eta2530_pt0102",";Momentum Resolution, 2.5 < #eta < 3, 1 < pT < 2",25,0,1);
    TH1F *mom_res_eta2530_pt0203 = new TH1F("mom_res_eta2530_pt0203",";Momentum Resolution, 2.5 < #eta < 3, 2 < pT < 3",25,0,1);
    TH1F *mom_res_eta2530_pt0304 = new TH1F("mom_res_eta2530_pt0304",";Momentum Resolution, 2.5 < #eta < 3, 3 < pT < 4",25,0,1);
    TH1F *mom_res_eta2530_pt0405 = new TH1F("mom_res_eta2530_pt0405",";Momentum Resolution, 2.5 < #eta < 3, 4 < pT < 5",25,0,1);
    TH1F *mom_res_eta2530_pt0506 = new TH1F("mom_res_eta2530_pt0506",";Momentum Resolution, 2.5 < #eta < 3, 5 < pT < 6",25,0,1);
    TH1F *mom_res_eta2530_pt0607 = new TH1F("mom_res_eta2530_pt0607",";Momentum Resolution, 2.5 < #eta < 3, 6 < pT < 7",25,0,1);
    TH1F *mom_res_eta2530_pt0708 = new TH1F("mom_res_eta2530_pt0708",";Momentum Resolution, 2.5 < #eta < 3, 7 < pT < 8",25,0,1);
    TH1F *mom_res_eta2530_pt0809 = new TH1F("mom_res_eta2530_pt0809",";Momentum Resolution, 2.5 < #eta < 3, 8 < pT < 9",25,0,1);
    TH1F *mom_res_eta2530_pt0910 = new TH1F("mom_res_eta2530_pt0910",";Momentum Resolution, 2.5 < #eta < 3, 9 < pT < 10",25,0,1);
    
    TH1F *mom_res_eta3035_pt0001 = new TH1F("mom_res_eta3035_pt0001",";Momentum Resolution, 3 < #eta < 3.5, 0 < pT < 1",25,0,1);
    TH1F *mom_res_eta3035_pt0102 = new TH1F("mom_res_eta3035_pt0102",";Momentum Resolution, 3 < #eta < 3.5, 1 < pT < 2",25,0,1);
    TH1F *mom_res_eta3035_pt0203 = new TH1F("mom_res_eta3035_pt0203",";Momentum Resolution, 3 < #eta < 3.5, 2 < pT < 3",25,0,1);
    TH1F *mom_res_eta3035_pt0304 = new TH1F("mom_res_eta3035_pt0304",";Momentum Resolution, 3 < #eta < 3.5, 3 < pT < 4",25,0,1);
    TH1F *mom_res_eta3035_pt0405 = new TH1F("mom_res_eta3035_pt0405",";Momentum Resolution, 3 < #eta < 3.5, 4 < pT < 5",25,0,1);
    TH1F *mom_res_eta3035_pt0506 = new TH1F("mom_res_eta3035_pt0506",";Momentum Resolution, 3 < #eta < 3.5, 5 < pT < 6",25,0,1);
    TH1F *mom_res_eta3035_pt0607 = new TH1F("mom_res_eta3035_pt0607",";Momentum Resolution, 3 < #eta < 3.5, 6 < pT < 7",25,0,1);
    TH1F *mom_res_eta3035_pt0708 = new TH1F("mom_res_eta3035_pt0708",";Momentum Resolution, 3 < #eta < 3.5, 7 < pT < 8",25,0,1);
    TH1F *mom_res_eta3035_pt0809 = new TH1F("mom_res_eta3035_pt0809",";Momentum Resolution, 3 < #eta < 3.5, 8 < pT < 9",25,0,1);
    TH1F *mom_res_eta3035_pt0910 = new TH1F("mom_res_eta3035_pt0910",";Momentum Resolution, 3 < #eta < 3.5, 9 < pT < 10",25,0,1);
    
    

    int ctr1 = 0;
    int ctr2 = 0;
    int ctr3 = 0;
    int ctr4 = 0;
    int ctr5 = 0;
    int ctr6 = 0;
    int ctr7 = 0;
    int ctr8 = 0;
    int ctr9 = 0;
    int ctr10 = 0;
    int ctr11 = 0;
    int ctr12 = 0;
    int ctr13 = 0;
    int ctr14 = 0;

    //fill histograms based on the bins
    for (int i=0; i<nEntries; i++) {
        if (i%1000000==0) cout << "jet: " << i << " out of: " << nEntries << endl;
        
        tree->GetEntry(i);
        float eta = calculateEta(px,py,pz);
        float pT = calculateP(px,py,pz);
        float momRes = calculateMomRes(gpx,gpy,gpz,px,py,pz);

        if (i<10) cout << "p for " <<i << " is " << px << ", " << py << ", " << pz << endl;
        if (i<10) cout << "eta for " <<i << " is " << eta << endl;
        if (i<10) cout << "pT for " <<i << " is " << pT << endl;
        if (i<10) cout << "theta for " <<i << " is " << acos(pz/pT) << endl;


        if (eta <= -3.5) cout << "jet < -3.5!" << endl;
        else if (eta <= -3) {
            ctr1++;
            if (pT <= 1) mom_res_etan3530_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_etan3530_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_etan3530_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_etan3530_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_etan3530_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_etan3530_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_etan3530_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_etan3530_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_etan3530_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_etan3530_pt0910->Fill(momRes);
        } else if (eta <= -2.5) {
            ctr2++;
            if (pT <= 1) mom_res_etan3025_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_etan3025_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_etan3025_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_etan3025_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_etan3025_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_etan3025_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_etan3025_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_etan3025_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_etan3025_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_etan3025_pt0910->Fill(momRes);
        } else if (eta <= -2) {
            ctr3++;
            if (pT <= 1) mom_res_etan2520_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_etan2520_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_etan2520_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_etan2520_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_etan2520_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_etan2520_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_etan2520_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_etan2520_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_etan2520_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_etan2520_pt0910->Fill(momRes);
        } else if (eta <= -1.5) {
            ctr4++;
            if (pT <= 1) mom_res_etan2015_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_etan2015_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_etan2015_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_etan2015_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_etan2015_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_etan2015_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_etan2015_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_etan2015_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_etan2015_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_etan2015_pt0910->Fill(momRes);
        } else if (eta <= -1) {
            ctr5++;
            if (pT <= 1) mom_res_etan1510_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_etan1510_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_etan1510_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_etan1510_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_etan1510_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_etan1510_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_etan1510_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_etan1510_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_etan1510_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_etan1510_pt0910->Fill(momRes);
        } else if (eta <= -0.5) {
            ctr6++;
            if (pT <= 1) mom_res_etan1005_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_etan1005_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_etan1005_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_etan1005_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_etan1005_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_etan1005_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_etan1005_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_etan1005_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_etan1005_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_etan1005_pt0910->Fill(momRes);
        } else if (eta <= 0) {
            ctr7++;
            if (pT <= 1) mom_res_etan0500_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_etan0500_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_etan0500_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_etan0500_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_etan0500_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_etan0500_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_etan0500_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_etan0500_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_etan0500_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_etan0500_pt0910->Fill(momRes);
        } else if (eta <= 0.5) {
            ctr8++;
            if (pT <= 1) mom_res_eta0005_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_eta0005_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_eta0005_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_eta0005_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_eta0005_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_eta0005_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_eta0005_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_eta0005_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_eta0005_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_eta0005_pt0910->Fill(momRes);
        } else if (eta <= 1) {
            ctr9++;
            if (pT <= 1) mom_res_eta0510_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_eta0510_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_eta0510_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_eta0510_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_eta0510_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_eta0510_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_eta0510_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_eta0510_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_eta0510_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_eta0510_pt0910->Fill(momRes);
        } else if (eta <= 1.5) {
            ctr10++;
            if (pT <= 1) mom_res_eta1015_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_eta1015_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_eta1015_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_eta1015_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_eta1015_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_eta1015_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_eta1015_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_eta1015_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_eta1015_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_eta1015_pt0910->Fill(momRes);
        } else if (eta <= 2) {
            ctr11++;
            if (pT <= 1) mom_res_eta1520_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_eta1520_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_eta1520_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_eta1520_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_eta1520_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_eta1520_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_eta1520_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_eta1520_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_eta1520_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_eta1520_pt0910->Fill(momRes);
        } else if (eta <= 2.5) {
            ctr12++;
            if (pT <= 1) mom_res_eta2025_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_eta2025_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_eta2025_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_eta2025_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_eta2025_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_eta2025_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_eta2025_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_eta2025_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_eta2025_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_eta2025_pt0910->Fill(momRes);
        } else if (eta <= 3) {
            ctr13++;
            if (pT <= 1) mom_res_eta2530_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_eta2530_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_eta2530_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_eta2530_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_eta2530_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_eta2530_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_eta2530_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_eta2530_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_eta2530_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_eta2530_pt0910->Fill(momRes);
        } else if (eta <= 3.5) {
            ctr14++;
            if (pT <= 1) mom_res_eta3035_pt0001->Fill(momRes);
            else if (pT <= 2) mom_res_eta3035_pt0102->Fill(momRes);
            else if (pT <= 3) mom_res_eta3035_pt0203->Fill(momRes);
            else if (pT <= 4) mom_res_eta3035_pt0304->Fill(momRes);
            else if (pT <= 5) mom_res_eta3035_pt0405->Fill(momRes);
            else if (pT <= 6) mom_res_eta3035_pt0506->Fill(momRes);
            else if (pT <= 7) mom_res_eta3035_pt0607->Fill(momRes);
            else if (pT <= 8) mom_res_eta3035_pt0708->Fill(momRes);
            else if (pT <= 9) mom_res_eta3035_pt0809->Fill(momRes);
            else if (pT <= 10) mom_res_eta3035_pt0910->Fill(momRes);
        } 
        // else cout << "jet > 3.5!" << endl;

    }


    cout << ctr1 << " " << ctr2 << " " << ctr3 <<" " << ctr3 <<" " << ctr4 << " " << ctr5 << " " << ctr6 << " " << ctr7 << " " << ctr8 <<endl;
    cout << ctr9 << " " << ctr10 << " " << ctr11 << " " << ctr12 << " " << ctr13 << " " << ctr14 << endl;

    //make canvas to draw
    // TCanvas c1 = makeCanvas();
    TCanvas *c1 = new TCanvas("c1","c1",1200,900);
    c1->Divide(5,2);
    // c1->SetTitle("Jet Energy Drop by pt");
    // gPad->SetLogy();

    c1->cd(1); mom_res_etan3530_pt0001->Draw();
    c1->cd(2); mom_res_etan3530_pt0102->Draw();
    c1->cd(3); mom_res_etan3530_pt0203->Draw();
    c1->cd(4); mom_res_etan3530_pt0304->Draw();
    c1->cd(5); mom_res_etan3530_pt0405->Draw();
    c1->cd(6); mom_res_etan3530_pt0506->Draw();
    c1->cd(7); mom_res_etan3530_pt0607->Draw();
    c1->cd(8); mom_res_etan3530_pt0708->Draw();
    c1->cd(9); mom_res_etan3530_pt0809->Draw();
    c1->cd(10); mom_res_etan3530_pt0910->Draw();
    c1->Print("plots/mom_res_etan3530.png");

    c1->cd(1); mom_res_etan3025_pt0001->Draw();
    c1->cd(2); mom_res_etan3025_pt0102->Draw();
    c1->cd(3); mom_res_etan3025_pt0203->Draw();
    c1->cd(4); mom_res_etan3025_pt0304->Draw();
    c1->cd(5); mom_res_etan3025_pt0405->Draw();
    c1->cd(6); mom_res_etan3025_pt0506->Draw();
    c1->cd(7); mom_res_etan3025_pt0607->Draw();
    c1->cd(8); mom_res_etan3025_pt0708->Draw();
    c1->cd(9); mom_res_etan3025_pt0809->Draw();
    c1->cd(10); mom_res_etan3025_pt0910->Draw();
    c1->Print("plots/mom_res_etan3025.png");

    // c1->cd(1); mom_res_etan2520_pt0001->Draw();
    // c1->cd(2); mom_res_etan2520_pt0102->Draw();
    // c1->cd(3); mom_res_etan2520_pt0203->Draw();
    // c1->cd(4); mom_res_etan2520_pt0304->Draw();
    // c1->cd(5); mom_res_etan2520_pt0405->Draw();
    // c1->cd(6); mom_res_etan2520_pt0506->Draw();
    // c1->cd(7); mom_res_etan2520_pt0607->Draw();
    // c1->cd(8); mom_res_etan2520_pt0708->Draw();
    // c1->cd(9); mom_res_etan2520_pt0809->Draw();
    // c1->cd(10); mom_res_etan2520_pt0910->Draw();
    // c1->Print("plots/mom_res_etan2520.pdf");

    // c1->cd(1); mom_res_etan2015_pt0001->Draw();
    // c1->cd(2); mom_res_etan2015_pt0102->Draw();
    // c1->cd(3); mom_res_etan2015_pt0203->Draw();
    // c1->cd(4); mom_res_etan2015_pt0304->Draw();
    // c1->cd(5); mom_res_etan2015_pt0405->Draw();
    // c1->cd(6); mom_res_etan2015_pt0506->Draw();
    // c1->cd(7); mom_res_etan2015_pt0607->Draw();
    // c1->cd(8); mom_res_etan2015_pt0708->Draw();
    // c1->cd(9); mom_res_etan2015_pt0809->Draw();
    // c1->cd(10); mom_res_etan2015_pt0910->Draw();
    // c1->Print("plots/mom_res_etan2015.pdf");

    // c1->cd(1); mom_res_etan1510_pt0001->Draw();
    // c1->cd(2); mom_res_etan1510_pt0102->Draw();
    // c1->cd(3); mom_res_etan1510_pt0203->Draw();
    // c1->cd(4); mom_res_etan1510_pt0304->Draw();
    // c1->cd(5); mom_res_etan1510_pt0405->Draw();
    // c1->cd(6); mom_res_etan1510_pt0506->Draw();
    // c1->cd(7); mom_res_etan1510_pt0607->Draw();
    // c1->cd(8); mom_res_etan1510_pt0708->Draw();
    // c1->cd(9); mom_res_etan1510_pt0809->Draw();
    // c1->cd(10); mom_res_etan1510_pt0910->Draw();
    // c1->Print("plots/mom_res_etan1510.pdf");

    // c1->cd(1); mom_res_etan1005_pt0001->Draw();
    // c1->cd(2); mom_res_etan1005_pt0102->Draw();
    // c1->cd(3); mom_res_etan1005_pt0203->Draw();
    // c1->cd(4); mom_res_etan1005_pt0304->Draw();
    // c1->cd(5); mom_res_etan1005_pt0405->Draw();
    // c1->cd(6); mom_res_etan1005_pt0506->Draw();
    // c1->cd(7); mom_res_etan1005_pt0607->Draw();
    // c1->cd(8); mom_res_etan1005_pt0708->Draw();
    // c1->cd(9); mom_res_etan1005_pt0809->Draw();
    // c1->cd(10); mom_res_etan1005_pt0910->Draw();
    // c1->Print("plots/mom_res_etan1005.pdf");

    // c1->cd(1); mom_res_etan0500_pt0001->Draw();
    // c1->cd(2); mom_res_etan0500_pt0102->Draw();
    // c1->cd(3); mom_res_etan0500_pt0203->Draw();
    // c1->cd(4); mom_res_etan0500_pt0304->Draw();
    // c1->cd(5); mom_res_etan0500_pt0405->Draw();
    // c1->cd(6); mom_res_etan0500_pt0506->Draw();
    // c1->cd(7); mom_res_etan0500_pt0607->Draw();
    // c1->cd(8); mom_res_etan0500_pt0708->Draw();
    // c1->cd(9); mom_res_etan0500_pt0809->Draw();
    // c1->cd(10); mom_res_etan0500_pt0910->Draw();
    // c1->Print("plots/mom_res_etan0500.pdf");

    // c1->cd(1); mom_res_eta0005_pt0001->Draw();
    // c1->cd(2); mom_res_eta0005_pt0102->Draw();
    // c1->cd(3); mom_res_eta0005_pt0203->Draw();
    // c1->cd(4); mom_res_eta0005_pt0304->Draw();
    // c1->cd(5); mom_res_eta0005_pt0405->Draw();
    // c1->cd(6); mom_res_eta0005_pt0506->Draw();
    // c1->cd(7); mom_res_eta0005_pt0607->Draw();
    // c1->cd(8); mom_res_eta0005_pt0708->Draw();
    // c1->cd(9); mom_res_eta0005_pt0809->Draw();
    // c1->cd(10); mom_res_eta0005_pt0910->Draw();
    // c1->Print("plots/mom_res_eta0005.pdf");

    // c1->cd(1); mom_res_eta0510_pt0001->Draw();
    // c1->cd(2); mom_res_eta0510_pt0102->Draw();
    // c1->cd(3); mom_res_eta0510_pt0203->Draw();
    // c1->cd(4); mom_res_eta0510_pt0304->Draw();
    // c1->cd(5); mom_res_eta0510_pt0405->Draw();
    // c1->cd(6); mom_res_eta0510_pt0506->Draw();
    // c1->cd(7); mom_res_eta0510_pt0607->Draw();
    // c1->cd(8); mom_res_eta0510_pt0708->Draw();
    // c1->cd(9); mom_res_eta0510_pt0809->Draw();
    // c1->cd(10); mom_res_eta0510_pt0910->Draw();
    // c1->Print("plots/mom_res_eta0510.pdf");

    // c1->cd(1); mom_res_eta1015_pt0001->Draw();
    // c1->cd(2); mom_res_eta1015_pt0102->Draw();
    // c1->cd(3); mom_res_eta1015_pt0203->Draw();
    // c1->cd(4); mom_res_eta1015_pt0304->Draw();
    // c1->cd(5); mom_res_eta1015_pt0405->Draw();
    // c1->cd(6); mom_res_eta1015_pt0506->Draw();
    // c1->cd(7); mom_res_eta1015_pt0607->Draw();
    // c1->cd(8); mom_res_eta1015_pt0708->Draw();
    // c1->cd(9); mom_res_eta1015_pt0809->Draw();
    // c1->cd(10); mom_res_eta1015_pt0910->Draw();
    // c1->Print("plots/mom_res_eta1015.pdf");

    // c1->cd(1); mom_res_eta1520_pt0001->Draw();
    // c1->cd(2); mom_res_eta1520_pt0102->Draw();
    // c1->cd(3); mom_res_eta1520_pt0203->Draw();
    // c1->cd(4); mom_res_eta1520_pt0304->Draw();
    // c1->cd(5); mom_res_eta1520_pt0405->Draw();
    // c1->cd(6); mom_res_eta1520_pt0506->Draw();
    // c1->cd(7); mom_res_eta1520_pt0607->Draw();
    // c1->cd(8); mom_res_eta1520_pt0708->Draw();
    // c1->cd(9); mom_res_eta1520_pt0809->Draw();
    // c1->cd(10); mom_res_eta1520_pt0910->Draw();
    // c1->Print("plots/mom_res_eta1520.pdf");

    // c1->cd(1); mom_res_eta2025_pt0001->Draw();
    // c1->cd(2); mom_res_eta2025_pt0102->Draw();
    // c1->cd(3); mom_res_eta2025_pt0203->Draw();
    // c1->cd(4); mom_res_eta2025_pt0304->Draw();
    // c1->cd(5); mom_res_eta2025_pt0405->Draw();
    // c1->cd(6); mom_res_eta2025_pt0506->Draw();
    // c1->cd(7); mom_res_eta2025_pt0607->Draw();
    // c1->cd(8); mom_res_eta2025_pt0708->Draw();
    // c1->cd(9); mom_res_eta2025_pt0809->Draw();
    // c1->cd(10); mom_res_eta2025_pt0910->Draw();
    // c1->Print("plots/mom_res_eta2025.pdf");

    // c1->cd(1); mom_res_eta2530_pt0001->Draw();
    // c1->cd(2); mom_res_eta2530_pt0102->Draw();
    // c1->cd(3); mom_res_eta2530_pt0203->Draw();
    // c1->cd(4); mom_res_eta2530_pt0304->Draw();
    // c1->cd(5); mom_res_eta2530_pt0405->Draw();
    // c1->cd(6); mom_res_eta2530_pt0506->Draw();
    // c1->cd(7); mom_res_eta2530_pt0607->Draw();
    // c1->cd(8); mom_res_eta2530_pt0708->Draw();
    // c1->cd(9); mom_res_eta2530_pt0809->Draw();
    // c1->cd(10); mom_res_eta2530_pt0910->Draw();
    // c1->Print("plots/mom_res_eta2530.pdf");

    // c1->cd(1); mom_res_eta3035_pt0001->Draw();
    // c1->cd(2); mom_res_eta3035_pt0102->Draw();
    // c1->cd(3); mom_res_eta3035_pt0203->Draw();
    // c1->cd(4); mom_res_eta3035_pt0304->Draw();
    // c1->cd(5); mom_res_eta3035_pt0405->Draw();
    // c1->cd(6); mom_res_eta3035_pt0506->Draw();
    // c1->cd(7); mom_res_eta3035_pt0607->Draw();
    // c1->cd(8); mom_res_eta3035_pt0708->Draw();
    // c1->cd(9); mom_res_eta3035_pt0809->Draw();
    // c1->cd(10); mom_res_eta3035_pt0910->Draw();
    // c1->Print("plots/mom_res_eta3035.pdf");
    


}


//THIS DOESN'T WORK
TCanvas* makeCanvas(){
    //pt canvas
    TCanvas *c1 = new TCanvas("c1","c1",1200,900);
    c1->Divide(5,2);
    // c1->SetTitle("Jet Energy Drop by pt");
    gPad->SetLogy();

    return c1;
}

//fit plot with a gaussian and find SD


//combine
void plotMomentumResolution(){
    //read in file
    // TFile *file = new TFile("/project/projectdirs/m3763/blianggi/out_ECCE/60218804/tracking_output/ALL_G4EICDetector_g4tracking_eval.root");
    TFile *file = new TFile("/project/projectdirs/m3763/blianggi/out_ECCE/60218804/tracking_output/G4EICDetector_18_g4tracking_eval.root");
    TTree *tree = (TTree*) file->Get("tracks");

    float gpx, gpy, gpz, px, py, pz;
    const int nEntries = tree->GetEntries();

    tree->SetBranchAddress("gpx", &gpx);
    tree->SetBranchAddress("gpy", &gpy);
    tree->SetBranchAddress("gpz", &gpz);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);

    plotMomRes(tree, nEntries, gpx, gpy, gpz, px, py, pz); 

}