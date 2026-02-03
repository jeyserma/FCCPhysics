

namespace FCCAnalyses {

bool is_ww_hadronic(Vec_mc mc, Vec_i ind) {
   int l1 = 0;
   int l2 = 0;
   //cout << "*********" << endl;
   for(size_t i = 0; i < mc.size(); ++i) {
        auto & p = mc[i];
        if(std::abs(p.PDG) == 24) {
            int ds = p.daughters_begin;
            int de = p.daughters_end;
            for(int k=ds; k<de; k++) {
                int pdg = abs(mc[ind[k]].PDG);
                if(pdg == 24) continue;
                //std::cout << "W " << pdg << endl;
                if(l1 == 0) l1 = pdg;
                else l2 = pdg;
            }
        }
   }
   if((l1 < 6 && l2 < 6)) {
       //std::cout << "HADRONIC-----------" << l1 << " " << l2 << endl;
       return true;
   }
   return false;
}

bool is_zz_hadronic(Vec_mc mc, Vec_i ind) {
   int l1 = 0;
   int l2 = 0;
   //cout << "*********" << endl;
   for(size_t i = 0; i < mc.size(); ++i) {
        auto & p = mc[i];
        if(std::abs(p.PDG) == 23) {
            int ds = p.daughters_begin;
            int de = p.daughters_end;
            for(int k=ds; k<de; k++) {
                int pdg = abs(mc[ind[k]].PDG);
                if(pdg == 24) continue;
                //std::cout << "W " << pdg << endl;
                if(l1 == 0) l1 = pdg;
                else l2 = pdg;
            }
        }
   }
   if((l1 < 6 && l2 < 6)) {
       //std::cout << "HADRONIC-----------" << l1 << " " << l2 << endl;
       return true;
   }
   return false;
}


// Original chi2 for WW and ZZ
// Slightly problematic for ZZ: peaks at 80 GeV
Vec_i pairing_WW_ZZ(Vec_tlv J, float target) {
    int nJets6 = 6;
    float d1W[6][6], dW1[6][6], d2W[6][6],dW2[6][6], d1Z[6][6],d2Z[6][6], d1H[6][6][6][6],d2H[6][6][6][6];

    for(int i=0;i<6;i++) {
        for(int j=0;j<6;j++) {
            d1W[i][j]=-1000.;
            d2W[i][j]=-1000.;
            d1Z[i][j]=-1000.;
            d2Z[i][j]=-1000.;
            dW1[i][j]=-1000.;
            dW2[i][j]=-1000.;
        }
    }

    for(int i=0;i<6;i++){
        for(int j=0;j<6;j++){
            for(int k=0;k<6;k++){
                for(int s=0;s<6;s++) {
                    d1H[i][j][k][s]=-1000.;
                    d2H[i][j][k][s]=-1000.;
                }
            }
        }
    }

    double DWmin_1=1000.,DWmin_2=1000.,chi2_1=10000000000., chi2_2=10000000000.,chi2M=10000000000.,chi2min=100000000.,chi2First=10000000.,chi2Second=100000000;float dmin=1000000.,DZmin=1000000.;

    int iW1=-1.,jW1=-1.,iW2=-1.,jW2=-1.,iZ=-1.,jZ=-1.;
    int iW1_Z=-1.,jW1_Z=-1.,iW2_Z=-1.,jW2_Z=-1.,iZ_Z=-1.,jZ_Z=-1.;
    int iW1_1=-1.,jW1_1=-1.,iW2_1=-1.,jW2_1=-1.,iZ_1=-1.,jZ_1=-1.;
    int iW1_2=-1.,jW1_2=-1.,iW2_2=-1.,jW2_2=-1.,iZ_2=-1.,jZ_2=-1.;
    int iW1_3=-1.,jW1_3=-1.,iW2_3=-1.,jW2_3=-1.,iZ_3=-1.,jZ_3=-1.;

    for(int i=0;i<nJets6;i++) {
        for(int j=i+1;j<nJets6;j++) {
            for(int k=i+1;k<nJets6;k++) {
                if(!(k==i)&&!(k==j)) {
                    for(int s=k+1;s<nJets6;s++) {
                        if(!(s==i)&&!(s==j)&&!(s==k)) {
                            for(int l=0;l<nJets6;l++) {
                                if(!(l==i)&&!(l==j)&&!(l==k)&&!(l==s)) {
                                    for(int m=l+1;m<nJets6;m++) {
                                        if(!(m==i)&&!(m==j)&&!(m==k)&&!(m==s)&&!(m==l)) {

                                            d1W[i][j]      =fabs( target -((J[i]+J[j]).M()));
                                            dW1[l][m]      =fabs( target -((J[l]+J[m]).M()));
                                            d1Z[k][s]      =fabs( 91.19  -((J[k]+J[s]).M()));
                                            d1H[i][j][l][m]=fabs(125.0   -((J[i]+J[j]+J[l]+J[m]).M()));
                                      
                                            chi2First=(d1W[i][j])*(d1W[i][j])+(d1Z[k][s])*(d1Z[k][s]); // chi2 with only V
                                            chi2_1   =(d1W[i][j])*(d1W[i][j])+(d1Z[k][s])*(d1Z[k][s])+d1H[i][j][l][m]*d1H[i][j][l][m]; // chi2 with V and Higgs

                                            if(d1Z[k][s]<DZmin) {
                                                DZmin=d1Z[k][s];
                                                iZ_Z=k;
                                                jZ_Z=s;
                                                if(d1W[i][j]<dW1[l][m]) { iW1_Z=i; jW1_Z=j; iW2_Z=l;jW2_Z=m;}
                                                else                    { iW1_Z=l; jW1_Z=m; iW2_Z=i;jW2_Z=j;}
                                            }
                                            if(chi2First<chi2min) {chi2min=chi2First;iW1_1=i;jW1_1=j; iZ_1=k; jZ_1=s;iW2_1=l;jW2_1=m;}
                                            if(chi2_1<chi2M)      {chi2M=chi2_1;     iW1_2=i;jW1_2=j; iZ_2=k; jZ_2=s;iW2_2=l;jW2_2=m;}
                                            if( DWmin_1<dmin)     {dmin= DWmin_1;    iW1_3=i;jW1_3=j; iZ_3=k; jZ_3=s;iW2_3=l;jW2_3=m;}

                                            d2W[k][s]      =fabs( target -((J[k]+J[s]).M()));
                                            dW2[l][m]      =fabs( target -((J[l]+J[m]).M()));
                                            d2Z[i][j]      =fabs( 91.19  -((J[i]+J[j]).M()));
                                            d2H[k][s][l][m]=fabs(125.0   -((J[k]+J[s]+J[l]+J[m]).M()));

                                            chi2Second=(d2W[k][s])*(d2W[k][s])+(d2Z[i][j])*(d2Z[i][j]); 
                                            chi2_2=(d2W[k][s])*(d2W[k][s])+(d2Z[i][j])*(d2Z[i][j])+(d2H[k][s][l][m]*d2H[k][s][l][m]); 

                                            if(d2Z[i][j]<DZmin) {
                                                DZmin=d2Z[i][j];
                                                iZ_Z=i;
                                                jZ_Z=j; 

                                                if(d2W[k][s]<dW2[l][m]) {iW1_Z=k; jW1_Z=s; iW2_Z=l; jW2_Z=m;}
                                                else {iW1_Z=l; jW1_Z=m; iW2_Z=k; jW2_Z=s;}
                                            }

                                            if(chi2Second<chi2min) {chi2min=chi2Second;iW1_1=k;jW1_1=s;iZ_1=i; jZ_1=j;iW2_1=l;jW2_1=m;}
                                            if(chi2_2<chi2M)       {chi2M  =chi2_2;    iW1_2=k;jW1_2=s;iZ_2=i; jZ_2=j;iW2_2=l;jW2_2=m;}
                                            if( DWmin_2<dmin)      {dmin   =DWmin_2;   iW1_3=k;jW1_3=s; iZ_3=i;jZ_3=j;iW2_3=l;jW2_3=m;}
                                        }
                                    } // m
                                }
                            }// l
                        }
                    } // s
                }
            }// k
        } // j
    }// i

    Vec_i ret;
    ret.push_back(iW1_2); // jet1 onshell V
    ret.push_back(jW1_2); // jet2 onshell V
    ret.push_back(iW2_2); // jet1 offshell V
    ret.push_back(jW2_2); // jet2 offshell V
    ret.push_back(iZ_2); // jet1 production Z
    ret.push_back(jZ_2); // jet2 production Z
    ret.push_back(chi2M); // jet2 production Z
    return ret;
}


//Vec_tlv pairing_WW(Vec_tlv J, float target=80.385) {
Vec_i pairing_test(Vec_tlv J, float target) {
    // chi2 from Mila
    int nJets6=6;
    float d1W[6][6], dW1[6][6], d2W[6][6],dW2[6][6], d1Z[6][6],d2Z[6][6], d1H[6][6][6][6],d2H[6][6][6][6];

    for(int i=0;i<6;i++) {
        for(int j=0;j<6;j++) {
            d1W[i][j]=-1000.;
            d2W[i][j]=-1000.;
            d1Z[i][j]=-1000.;
            d2Z[i][j]=-1000.;
            dW1[i][j]=-1000.;
            dW2[i][j]=-1000.;
        }
    }

    for(int i=0;i<6;i++){
        for(int j=0;j<6;j++){
            for(int k=0;k<6;k++){
                for(int s=0;s<6;s++) {
                    d1H[i][j][k][s]=-1000.;
                    d2H[i][j][k][s]=-1000.;
                }
            }
        }
    }

    double DWmin_1=1000.,DWmin_2=1000.,chi2_1=10000000000., chi2_2=10000000000.,chi2M=10000000000.,chi2min=100000000.,chi2First=10000000.,chi2Second=100000000;float dmin=1000000.,DZmin=1000000.;

    int iW1=-1.,jW1=-1.,iW2=-1.,jW2=-1.,iZ=-1.,jZ=-1.;
    int iW1_Z=-1.,jW1_Z=-1.,iW2_Z=-1.,jW2_Z=-1.,iZ_Z=-1.,jZ_Z=-1.;
    int iW1_1=-1.,jW1_1=-1.,iW2_1=-1.,jW2_1=-1.,iZ_1=-1.,jZ_1=-1.;
    int iW1_2=-1.,jW1_2=-1.,iW2_2=-1.,jW2_2=-1.,iZ_2=-1.,jZ_2=-1.;
    int iW1_3=-1.,jW1_3=-1.,iW2_3=-1.,jW2_3=-1.,iZ_3=-1.,jZ_3=-1.;


    for(int i=0;i<nJets6;i++) {
        for(int j=i+1;j<nJets6;j++) {
            for(int k=i+1;k<nJets6;k++) {
                if(!(k==i)&&!(k==j)) {
                    for(int s=k+1;s<nJets6;s++) {
                        if(!(s==i)&&!(s==j)&&!(s==k)) {
                            for(int l=0;l<nJets6;l++) {
                                if(!(l==i)&&!(l==j)&&!(l==k)&&!(l==s)) {
                                    for(int m=l+1;m<nJets6;m++) {
                                        if(!(m==i)&&!(m==j)&&!(m==k)&&!(m==s)&&!(m==l)) {

                                            d1W[i][j]      =fabs( target -((J[i]+J[j]).M()));
                                            dW1[l][m]      =fabs( target -((J[l]+J[m]).M()));
                                            d1Z[k][s]      =fabs( 91.19  -((J[k]+J[s]).M()));
                                            d1H[i][j][l][m]=fabs(125.0   -((J[i]+J[j]+J[l]+J[m]).M()));
                                      
                                            chi2First=(d1W[i][j])*(d1W[i][j])+(d1Z[k][s])*(d1Z[k][s]); // chi2 with only V
                                            chi2_1   =(d1W[i][j])*(d1W[i][j])+(d1Z[k][s])*(d1Z[k][s])+d1H[i][j][l][m]*d1H[i][j][l][m]; // chi2 with V and Higgs
                                            chi2_1 = (d1Z[k][s])*(d1Z[k][s]) + d1H[i][j][l][m]*d1H[i][j][l][m];

                                            if(d1Z[k][s]<DZmin) {
                                                DZmin=d1Z[k][s];
                                                iZ_Z=k;
                                                jZ_Z=s;
                                                if(d1W[i][j]<dW1[l][m]) { iW1_Z=i; jW1_Z=j; iW2_Z=l;jW2_Z=m;}
                                                else                    { iW1_Z=l; jW1_Z=m; iW2_Z=i;jW2_Z=j;}
                                            }
                                            if(chi2First<chi2min) {chi2min=chi2First;iW1_1=i;jW1_1=j; iZ_1=k; jZ_1=s;iW2_1=l;jW2_1=m;}
                                            if(chi2_1<chi2M)      {chi2M=chi2_1;     iW1_2=i;jW1_2=j; iZ_2=k; jZ_2=s;iW2_2=l;jW2_2=m;}
                                            if( DWmin_1<dmin)     {dmin= DWmin_1;    iW1_3=i;jW1_3=j; iZ_3=k; jZ_3=s;iW2_3=l;jW2_3=m;}

                                            d2W[k][s]      =fabs( target -((J[k]+J[s]).M()));
                                            dW2[l][m]      =fabs( target -((J[l]+J[m]).M()));
                                            d2Z[i][j]      =fabs( 91.19  -((J[i]+J[j]).M()));
                                            d2H[k][s][l][m]=fabs(125.0   -((J[k]+J[s]+J[l]+J[m]).M()));

                                            chi2Second=(d2W[k][s])*(d2W[k][s])+(d2Z[i][j])*(d2Z[i][j]); 
                                            chi2_2=(d2W[k][s])*(d2W[k][s])+(d2Z[i][j])*(d2Z[i][j])+(d2H[k][s][l][m]*d2H[k][s][l][m]);
                                            chi2_2 = (d2Z[i][j])*(d2Z[i][j]) + (d2H[k][s][l][m]*d2H[k][s][l][m]);

                                            if(d2Z[i][j]<DZmin) {
                                                DZmin=d2Z[i][j];
                                                iZ_Z=i;
                                                jZ_Z=j; 

                                                if(d2W[k][s]<dW2[l][m]) {iW1_Z=k; jW1_Z=s; iW2_Z=l; jW2_Z=m;}
                                                else {iW1_Z=l; jW1_Z=m; iW2_Z=k; jW2_Z=s;}
                                            }

                                            if(chi2Second<chi2min) {chi2min=chi2Second;iW1_1=k;jW1_1=s;iZ_1=i; jZ_1=j;iW2_1=l;jW2_1=m;}
                                            if(chi2_2<chi2M)       {chi2M  =chi2_2;    iW1_2=k;jW1_2=s;iZ_2=i; jZ_2=j;iW2_2=l;jW2_2=m;}
                                            if( DWmin_2<dmin)      {dmin   =DWmin_2;   iW1_3=k;jW1_3=s; iZ_3=i;jZ_3=j;iW2_3=l;jW2_3=m;}
                                        }
                                    } // m
                                }
                            }// l
                        }
                    } // s
                }
            }// k
        } // j
    
    }// i


    int iW1_2_n = iW1_2, jW1_2_n = jW1_2, iW2_2_n = iW2_2, jW2_2_n = jW2_2; // Higgs jets
    float chi12 = fabs(target - (J[iW1_2]+J[jW1_2]).M());
    float chi13 = fabs(target - (J[iW1_2]+J[iW2_2]).M());
    float chi14 = fabs(target - (J[iW1_2]+J[jW2_2]).M());
    
    float chi23 = fabs(target - (J[jW1_2]+J[iW2_2]).M());
    float chi24 = fabs(target - (J[jW1_2]+J[jW2_2]).M());

    float chi34 = fabs(target - (J[iW2_2]+J[jW2_2]).M());

    //if(chi12<chi13 && chi12<chi14 && chi12<chi23 && chi12<chi24 && chi12<chi34) {
    //    cout << chi12 << " " << chi13 << " " << chi14 << " " << chi23 << " " << chi24 << " " << chi34 << endl;
    //}

    if(chi12<chi13 && chi12<chi14 && chi12<chi23 && chi12<chi24 && chi12<chi34) {iW1_2_n=iW1_2; jW1_2_n=jW1_2;   iW2_2_n=iW2_2; jW2_2_n=jW2_2; };
    if(chi13<chi12 && chi13<chi14 && chi13<chi23 && chi13<chi24 && chi13<chi34) {iW1_2_n=iW1_2; jW1_2_n=iW2_2;   iW2_2_n=jW1_2; jW2_2_n=jW2_2; };
    if(chi14<chi12 && chi14<chi13 && chi14<chi23 && chi14<chi24 && chi14<chi34) {iW1_2_n=iW1_2; jW1_2_n=jW2_2;   iW2_2_n=jW1_2; jW2_2_n=iW2_2; };
    
    if(chi23<chi13 && chi23<chi14 && chi23<chi23 && chi23<chi24 && chi23<chi34) {iW1_2_n=jW1_2; jW1_2_n=iW2_2;   iW2_2_n=iW1_2; jW2_2_n=jW2_2; };
    if(chi24<chi12 && chi24<chi13 && chi24<chi14 && chi24<chi23 && chi24<chi34) {iW1_2_n=jW1_2; jW1_2_n=jW2_2;   iW2_2_n=iW1_2; jW2_2_n=iW2_2; };
    
    if(chi34<chi12 && chi34<chi13 && chi34<chi14 && chi34<chi23 && chi34<chi24) {iW1_2_n=iW2_2; jW1_2_n=jW2_2;   iW2_2_n=iW1_2; jW2_2_n=jW1_2; };

    //if(chi2 < chi1 && chi2 < chi3) {iW1_2_n=iW1_2; jW1_2_n=iW2_2; iW2_2_n=jW1_2; jW2_2_n=jW2_2; cout << "CHI2" << endl; };
    //if(chi3 < chi1 && chi3 < chi2) {iW1_2_n=jW1_2; jW1_2_n=jW2_2; iW2_2_n=iW1_2; jW2_2_n=iW2_2; cout << "CHI3" << endl; };

    Vec_i ret;
    ret.push_back(iW1_2_n); // jet1 onshell V
    ret.push_back(jW1_2_n); // jet2 onshell V
    ret.push_back(iW2_2_n); // jet1 offshell V
    ret.push_back(jW2_2_n); // jet2 offshell V
    ret.push_back(iZ_2); // jet1 production Z
    ret.push_back(jZ_2); // jet2 production Z
    ret.push_back(chi2M); // jet2 production Z
    return ret;
}

//Vec_tlv pairing_WW(Vec_tlv J, float target=80.385) {
Vec_i pairing_WW(Vec_tlv J, float target=80.385) {
    // chi2 from Mila
    int nJets6=6;
    /*
    float d1W[6][6], dW1[6][6], d2W[6][6],dW2[6][6], d1Z[6][6],d2Z[6][6], d1H[6][6][6][6],d2H[6][6][6][6];

    for(int i=0;i<6;i++) {
        for(int j=0;j<6;j++) {
            d1W[i][j]=-1000.;
            d2W[i][j]=-1000.;
            d1Z[i][j]=-1000.;
            d2Z[i][j]=-1000.;
            dW1[i][j]=-1000.;
            dW2[i][j]=-1000.;
        }
    }

    for(int i=0;i<6;i++){
        for(int j=0;j<6;j++){
            for(int k=0;k<6;k++){
                for(int s=0;s<6;s++) {
                    d1H[i][j][k][s]=-1000.;
                    d2H[i][j][k][s]=-1000.;
                }
            }
        }
    }
*/
    double DWmin_1=1000.,DWmin_2=1000.,chi2_1=10000000000., chi2_2=10000000000.,chi2M=9e99,chi2min=100000000.,chi2First=10000000.,chi2Second=100000000;float dmin=1000000.,DZmin=1000000.;

    int iW1=-1.,jW1=-1.,iW2=-1.,jW2=-1.,iZ=-1.,jZ=-1.;
    int iW1_Z=-1.,jW1_Z=-1.,iW2_Z=-1.,jW2_Z=-1.,iZ_Z=-1.,jZ_Z=-1.;
    int iW1_1=-1.,jW1_1=-1.,iW2_1=-1.,jW2_1=-1.,iZ_1=-1.,jZ_1=-1.;
    int iW1_2=-1.,jW1_2=-1.,iW2_2=-1.,jW2_2=-1.,iZ_2=-1.,jZ_2=-1.;
    int iW1_3=-1.,jW1_3=-1.,iW2_3=-1.,jW2_3=-1.,iZ_3=-1.,jZ_3=-1.;


    TLorentzVector init, recoil;
    init.SetPxPyPzE(0, 0, 0, 240);

    float C_Von, C_Voff, C_Z, C_H, C_H_rec, chi2;
    for(int i=0;i<nJets6;i++) { // i 1
        for(int j=0;j<nJets6;j++) { // j 2
            if(j==i) continue;
            for(int k=0;k<nJets6;k++) { // k 3
                if(k==i or k==j) continue;
                for(int l=0;l<nJets6;l++) { // l 4
                    if(l==i or l==j or l==k) continue;
                    for(int m=0;m<nJets6;m++) { // m 5
                        if(m==i or m==j or m==k or m==l) continue;
                        for(int n=0;n<nJets6;n++) { // n 6
                            if(n==i or n==j or n==k or n==l or n==m) continue;

                            // we have i,j,k,l,m,n which are all different
                            // form all pairs

                            C_Von = fabs(target -((J[i]+J[j]).M())); // onshell V
                            C_Voff = fabs(target -((J[k]+J[l]).M())); // offshel V
                            C_Z = fabs(91.19  -((J[m]+J[n]).M())); // Z
                            C_H = fabs(125.0  -((J[i]+J[j]+J[k]+J[l]).M())); // Higgs
                            C_H_rec = fabs(125.0  -((init-J[m]-J[n]).M())); // Higgs recoil
                            chi2 = C_Von*C_Von + C_Z*C_Z + 0.0000001*C_H*C_H + 0.*C_H_rec*C_H_rec;
                            if(chi2 < chi2M) { chi2M=chi2; iW1_2=i; jW1_2=j; iW2_2=k; jW2_2=l; iZ_2=m; jZ_2=n; }


                            C_Von = fabs(target -((J[m]+J[n]).M())); // onshell V
                            C_Voff = fabs(target -((J[i]+J[j]).M())); // offshel V
                            C_Z = fabs(91.19  -((J[k]+J[l]).M())); // Z
                            C_H = fabs(125.0  -((J[m]+J[n]+J[i]+J[j]).M())); // Higgs
                            C_H_rec = fabs(125.0  -((init-J[k]-J[l]).M())); // Higgs recoil
                            chi2 = C_Von*C_Von + C_Z*C_Z + 0.0000001*C_H*C_H + 0.*C_H_rec*C_H_rec;
                            if(chi2 < chi2M) { chi2M=chi2; iW1_2=m; jW1_2=n; iW2_2=i; jW2_2=j; iZ_2=k; jZ_2=l; }


                            C_Von = fabs(target -((J[k]+J[l]).M())); // onshell V
                            C_Voff = fabs(target -((J[m]+J[n]).M())); // offshel V
                            C_Z = fabs(91.19  -((J[i]+J[j]).M())); // Z
                            C_H = fabs(125.0  -((J[k]+J[l]+J[m]+J[n]).M())); // Higgs
                            C_H_rec = fabs(125.0  -((init-J[i]-J[j]).M())); // Higgs recoil
                            chi2 = C_Von*C_Von + C_Z*C_Z + 0.0000001*C_H*C_H + 0.*C_H_rec*C_H_rec;
                            if(chi2 < chi2M) { chi2M=chi2; iW1_2=k; jW1_2=l; iW2_2=m; jW2_2=n; iZ_2=i; jZ_2=j; }

                        } // n
                    } // m
                } // l
            }// k
        } // j
    }// i


    // now we have 4 jets from Higgs iW1_2  jW1_2 iW2_2 jW2_2
    // repair them to from good V
    //cout <<  chi2 << endl;
    Vec_i ret;

    ret.push_back(iW1_2); // jet1 onshell V
    ret.push_back(jW1_2); // jet2 onshell V
    ret.push_back(iW2_2); // jet1 offshell V
    ret.push_back(jW2_2); // jet2 offshell V
    ret.push_back(iZ_2); // jet1 production Z
    ret.push_back(jZ_2); // jet2 production Z
    ret.push_back(chi2M);
/*
    ret.push_back(0); // jet1 onshell V
    ret.push_back(1); // jet2 onshell V
    ret.push_back(2); // jet1 offshell V
    ret.push_back(3); // jet2 offshell V
    ret.push_back(4); // jet1 production Z
    ret.push_back(5); // jet2 production Z
    */
    return ret;
}


Vec_i pairing_ZZ(Vec_tlv J, float target=91.2) {

    int Njets=6;
    float d1W[6][6], dW1[6][6], d2W[6][6],dW2[6][6], d1Z[6][6],d2Z[6][6], d1H[6][6][6][6],d2H[6][6][6][6];

    for(int i=0;i<6;i++) {
        for(int j=0;j<6;j++) {
            d1W[i][j]=-1000.;
            d2W[i][j]=-1000.;
            d1Z[i][j]=-1000.;
            d2Z[i][j]=-1000.;
            dW1[i][j]=-1000.;
            dW2[i][j]=-1000.;
        }
    }

    for(int i=0;i<6;i++){
        for(int j=0;j<6;j++){
            for(int k=0;k<6;k++){
                for(int s=0;s<6;s++) {
                    d1H[i][j][k][s]=-1000.;
                    d2H[i][j][k][s]=-1000.;
                }
            }
        }
    }

    double DWmin_1=1000.,DWmin_2=1000.,chi2_1=10000000000., chi2_2=10000000000.,chi2M=10000000000.,chi2min=100000000.,chi2First=10000000.,chi2Second=100000000;float dmin=1000000.,DZmin=1000000.;

    int iW1=-1.,jW1=-1.,iW2=-1.,jW2=-1.,iZ=-1.,jZ=-1.;
    int iW1_Z=-1.,jW1_Z=-1.,iW2_Z=-1.,jW2_Z=-1.,iZ_Z=-1.,jZ_Z=-1.;
    int iW1_1=-1.,jW1_1=-1.,iW2_1=-1.,jW2_1=-1.,iZ_1=-1.,jZ_1=-1.;
    int iW1_2=-1.,jW1_2=-1.,iW2_2=-1.,jW2_2=-1.,iZ_2=-1.,jZ_2=-1.;
    int iW1_3=-1.,jW1_3=-1.,iW2_3=-1.,jW2_3=-1.,iZ_3=-1.,jZ_3=-1.;


   	for(int i=0;i<Njets;i++){
	for(int j=i+1;j<Njets;j++){
	for(int k=i+1;k<Njets;k++){
	  if(!(k==i)&&!(k==j)){
	for(int s=k+1;s<Njets;s++){
       if(!(s==i)&&!(s==j)&&!(s==k)){
	for(int l=0;l<Njets;l++){
	  if(!(l==i)&&!(l==j)&&!(l==k)&&!(l==s)){
	for(int m=l+1;m<Njets;m++){
	  if(!(m==i)&&!(m==j)&&!(m==k)&&!(m==s)&&!(m==l)){
	

 	d1W[i][j]      =fabs( target -((J[i]+J[j]).M()));
	dW1[l][m]      =fabs( target -((J[l]+J[m]).M()));
 	d1Z[k][s]      =fabs( 91.19  -((J[k]+J[s]).M()));
	d1H[i][j][l][m]=fabs(125.0   -((J[i]+J[j]+J[l]+J[m]).M()));
	  
        chi2First=(d1W[i][j])*(d1W[i][j])+(d1Z[k][s])*(d1Z[k][s]); 
	chi2_1=(d1W[i][j])*(d1W[i][j])+(d1Z[k][s])*(d1Z[k][s])+d1H[i][j][l][m]*d1H[i][j][l][m]; 
	 
	 if(d1Z[k][s]<DZmin)
	  {
		DZmin=d1Z[k][s];

		 iZ_Z=k;
		 jZ_Z=s;
		if(d1W[i][j]<dW1[l][m]){ iW1_Z=i; jW1_Z=j; iW2_Z=l;jW2_Z=m;}   
		else{ iW1_Z=l; jW1_Z=m;  iW2_Z=i; jW2_Z=j;}
//		cout<<"1: "<<chi2min<<endl;
	  }
	  
        if(chi2First<chi2min) {chi2min=chi2First;iW1_1=i;jW1_1=j; iZ_1=k; jZ_1=s;iW2_1=l;jW2_1=m;  }
	if(chi2_1<chi2M)      {chi2M=chi2_1;     iW1_2=i;jW1_2=j; iZ_2=k; jZ_2=s;iW2_2=l;jW2_2=m;  }
	  if( DWmin_1<dmin)   {dmin= DWmin_1;    iW1_3=i;jW1_3=j; iZ_3=k; jZ_3=s;iW2_3=l;jW2_3=m; }
	  
/*	  cout<<"after chiFirst: "<<endl;
        cout<< "chi2_1"<<" "<< "chi2_2"<<" "<<"chi2M"<<" "<<"chi2min"<<" "<<"chi2First"<<" "<<"chi2Second"<<" "<<" "<<dmin<<" "<<DZmin<<endl;
        cout<< chi2_1<<" "<< chi2_2<<" "<<chi2M<<" "<<chi2min<<" "<<chi2First<<" "<<chi2Second<<" "<<" "<<dmin<<" "<<DZmin<<endl;
      	 cout<< iW1_2<<jW1_2<<iZ_2<<jZ_2<<iW2_2<<jW2_2<<endl;*/
  //   
          d2W[k][s]      =fabs( target -((J[k]+J[s]).M()));//pogledamo obe kombinacije jer setamo samo 12 34 ali ne i 34 12
	  dW2[l][m]      =fabs( target -((J[l]+J[m]).M()));//pogledamo obe kombinacije jer amo samo 12 34 ali ne i 34 12
	  d2Z[i][j]      =fabs( 91.19  -((J[i]+J[j]).M()));
          d2H[k][s][l][m]=fabs(125.0   -((J[k]+J[s]+J[l]+J[m]).M()));
	  
          
/*_1*/    chi2Second=(d2W[k][s])*(d2W[k][s])+(d2Z[i][j])*(d2Z[i][j]); 
/*_2*/	  chi2_2=(d2W[k][s])*(d2W[k][s])+(d2Z[i][j])*(d2Z[i][j])+  d2H[k][s][l][m]*d2H[k][s][l][m]; 
/*_3*/	  DWmin_2=(fabs(d2W[k][s])+fabs(d2Z[i][j]));
	  
	  if(d2Z[i][j]<DZmin)
	  {
		DZmin=d2Z[i][j];

		 iZ_Z=i;
		 jZ_Z=j;
		 
	    if(d2W[k][s]<dW2[l][m]) 
		{   
		  iW1_Z=k;
		  jW1_Z=s; 
		  iW2_Z=l;
		  jW2_Z=m;
		}   
		else
		{
		 iW1_Z=l;
		 jW1_Z=m; 
		 iW2_Z=k;
		 jW2_Z=s;
		}
//		cout<<"1: "<<chi2min<<endl;
	  }

	  if(chi2Second<chi2min) {chi2min=chi2Second;iW1_1=k;jW1_1=s; iZ_1=i; jZ_1=j; iW2_1=l;jW2_1=m; }
  	  if(chi2_2<chi2M)	 {chi2M=chi2_2;      iW1_2=k;jW1_2=s; iZ_2=i; jZ_2=j; iW2_2=l;jW2_2=m;  }
          if( DWmin_2<dmin)      {dmin= DWmin_2;     iW1_3=k;jW1_3=s; iZ_3=i; jZ_3=j; iW2_3=l;jW2_3=m;	}		

 /*	  cout<<"after chiSecond "<<endl;
        cout<< "chi2_1"<<" "<< "chi2_2"<<" "<<"chi2M"<<" "<<"chi2min"<<" "<<"chi2First"<<" "<<"chi2Second"<<" "<<" "<<dmin<<" "<<DZmin<<endl;
        cout<< chi2_1<<" "<< chi2_2<<" "<<chi2M<<" "<<chi2min<<" "<<chi2First<<" "<<chi2Second<<" "<<" "<<dmin<<" "<<DZmin<<endl;
      	 cout<< iW1_2<<jW1_2<<iZ_2<<jZ_2<<iW2_2<<jW2_2<<endl;*/
	 }}}}//l,m
	}}}}//k.s
	}}//i,j
//	cout<<"chi2min: "<<chi2min/1000000<<endl;





    Vec_i ret;
    ret.push_back(iW1_2); // jet1 onshell V        iW1_2 or  iW1_1
    ret.push_back(jW1_2); // jet2 onshell V        jW1_2 or  jW1_1
    ret.push_back(iW2_2); // jet1 offshell V       iW2_2 or  iW2_1
    ret.push_back(jW2_2); // jet2 offshell V       jW2_2 or  jW2_1
    ret.push_back(iZ_2); // jet1 production Z      iZ_2 or  iZ_1
    ret.push_back(jZ_2); // jet2 production Z      jZ_2 or  jZ_1
    return ret;
}


Vec_i pairing_ZZ_old(Vec_tlv J, float target=91.2) {
    // chi2 from Mila
    // is identical to WW it seems
    // energy resolution
    //double sigmaOverE = 0.05;
    //int count=0;
    //random_device rd;
    //mt19937 gen(rd());

    Float_t jet_p_1,jet_p_2,jet_p_3,jet_p_4,jet_p_5,jet_p_6;
    Float_t jet_e_1,jet_e_2,jet_e_3,jet_e_4,jet_e_5,jet_e_6;
    Float_t InvMassZ1_1,   InvMassZ2_1,   InvMassZ_1,InvMassHiggs_1;
    Float_t InvMassZ1_2,   InvMassZ2_2,   InvMassZ_2,InvMassHiggs_2;
    Float_t InvMassZ1_3,   InvMassZ2_3,   InvMassZ_3,InvMassHiggs_3;

    Float_t jet_theta_1,jet_theta_2,jet_theta_3,jet_theta_4,jet_theta_5,jet_theta_6;
    Float_t jet_phi_1,jet_phi_2,jet_phi_3,jet_phi_4,jet_phi_5,jet_phi_6;
    TLorentzVector JHiggs,JW1,JW2,JZ,JTotal;
    JHiggs*=0.,JW1*=0.,	JW2*=0.,JZ*=0.;
    JTotal*=0;
    TLorentzVector Jets2[2],Jets4[4],JW1_Jets4,JW2_Jets4 ;
    TVector3       Jets2_3[2],Jets4_3[4],JW1_Jets4_3,JW2_Jets4_3 ;
    float d1W[6][6],dW1[6][6], d2W[6][6],dW2[6][6], d1Z[6][6], d2Z[6][6], d1H[6][6][6][6],d2H[6][6][6][6], m_gen[6][6];
    for(int i=0;i<6;i++){for(int j=0;j<6;j++){
        d1W[i][j]=-1000.;
        d2W[i][j]=-1000.;
        d1Z[i][j]=-1000.;
        d2Z[i][j]=-1000.;
        dW1[i][j]=-1000.;
        dW2[i][j]=-1000.;
        m_gen[i][j]=-1000;
    }}
    for(int i=0;i<6;i++){for(int j=0;j<6;j++){for(int k=0;k<6;k++){for(int s=0;s<6;s++){
        d1H[i][j][k][s]=-1000.;
        d2H[i][j][k][s]=-1000.;
    }}}}
     
    double DWmin_1=1000.,DWmin_2=1000.,chi2_1=10000000000., chi2_2=10000000000.,chi2M=10000000000.,chi2min=100000000.,chi2First=10000000.,chi2Second=100000000;float dmin=1000000.,DZmin=1000000.;

    int iW1=-1.,jW1=-1.,iW2=-1.,jW2=-1.,iZ=-1.,jZ=-1.;
    int iW1_Z=-1.,jW1_Z=-1.,iW2_Z=-1.,jW2_Z=-1.,iZ_Z=-1.,jZ_Z=-1.;
    int iW1_1=-1.,jW1_1=-1.,iW2_1=-1.,jW2_1=-1.,iZ_1=-1.,jZ_1=-1.;
    int iW1_2=-1.,jW1_2=-1.,iW2_2=-1.,jW2_2=-1.,iZ_2=-1.,jZ_2=-1.;
    int iW1_3=-1.,jW1_3=-1.,iW2_3=-1.,jW2_3=-1.,iZ_3=-1.,jZ_3=-1.;

    int Njets=6;

    for(int i=0;i<Njets;i++) {
        for(int j=i+1;j<Njets;j++) {
            for(int k=i+1;k<Njets;k++) {
                if(!(k==i)&&!(k==j)) {
                    for(int s=k+1;s<Njets;s++) {
                        if(!(s==i)&&!(s==j)&&!(s==k)) {
                            for(int l=0;l<Njets;l++) {
                                if(!(l==i)&&!(l==j)&&!(l==k)&&!(l==s)) {
                                    for(int m=l+1;m<Njets;m++) {
                                        if(!(m==i)&&!(m==j)&&!(m==k)&&!(m==s)&&!(m==l)) {

                                            d1W[i][j]      =fabs(target -((J[i]+J[j]).M()));
                                            dW1[l][m]      =fabs(target -((J[l]+J[m]).M()));
                                            d1Z[k][s]      =fabs(91.19  -((J[k]+J[s]).M()));
                                            d1H[i][j][l][m]=fabs(125.0   -((J[i]+J[j]+J[l]+J[m]).M()));

                                            chi2First=(d1W[i][j])*(d1W[i][j])+(d1Z[k][s])*(d1Z[k][s]); 
                                            chi2_1=(d1W[i][j])*(d1W[i][j])+(d1Z[k][s])*(d1Z[k][s])+d1H[i][j][l][m]*d1H[i][j][l][m]; 
                                            DWmin_1=(fabs(d1W[i][j])+fabs(d1Z[k][s]));

                                            if(d1Z[k][s]<DZmin) {
                                                DZmin=d1Z[k][s];

                                                iZ_Z=k;
                                                jZ_Z=s;
                                                if(d1W[i][j]<dW1[l][m]) {
                                                    iW1_Z=i;
                                                    jW1_Z=j; 
                                                    iW2_Z=l;
                                                    jW2_Z=m;
                                                }   
                                                else {
                                                    iW1_Z=l;
                                                    jW1_Z=m; 
                                                    iW2_Z=i;
                                                    jW2_Z=j;
                                                }
                                            }

                                            if(chi2First<chi2min) { chi2min=chi2First; iW1_1=i; jW1_1=j; iZ_1=k; jZ_1=s; iW2_1=l; jW2_1=m; }
                                            if(chi2_1<chi2M) { chi2M=chi2_1; iW1_2=i; jW1_2=j; iZ_2=k; jZ_2=s; iW2_2=l; jW2_2=m; }
                                            if( DWmin_1<dmin) { dmin= DWmin_1; iW1_3=i; jW1_3=j; iZ_3=k; jZ_3=s; iW2_3=l; jW2_3=m; }

                                            d2W[k][s]      =fabs(target -((J[k]+J[s]).M()));
                                            dW2[l][m]      =fabs(target -((J[l]+J[m]).M()));
                                            d2Z[i][j]      =fabs(91.19  -((J[i]+J[j]).M()));
                                            d2H[k][s][l][m]=fabs(125.0   -((J[k]+J[s]+J[l]+J[m]).M()));

                                            chi2Second=(d2W[k][s])*(d2W[k][s])+(d2Z[i][j])*(d2Z[i][j]); 
                                            chi2_2=(d2W[k][s])*(d2W[k][s])+(d2Z[i][j])*(d2Z[i][j])+  d2H[k][s][l][m]*d2H[k][s][l][m]; 
                                            DWmin_2=(fabs(d2W[k][s])+fabs(d2Z[i][j]));

                                            if(d2Z[i][j]<DZmin) {
                                                DZmin=d2Z[i][j];
                                                iZ_Z=i;
                                                jZ_Z=j;

                                                if(d2W[k][s]<dW2[l][m]) {
                                                    iW1_Z=k;
                                                    jW1_Z=s; 
                                                    iW2_Z=l;
                                                    jW2_Z=m;
                                                }
                                                else {
                                                    iW1_Z=l;
                                                    jW1_Z=m; 
                                                    iW2_Z=k;
                                                    jW2_Z=s;
                                                }

                                            }

                                            if(chi2Second<chi2min) { chi2min=chi2Second; iW1_1=k; jW1_1=s; iZ_1=i; jZ_1=j; iW2_1=l; jW2_1=m; }
                                            if(chi2_2<chi2M) { chi2M=chi2_2; iW1_2=k; jW1_2=s; iZ_2=i; jZ_2=j; iW2_2=l; jW2_2=m; }
                                            if( DWmin_2<dmin) { dmin= DWmin_2; iW1_3=k; jW1_3=s; iZ_3=i; jZ_3=j; iW2_3=l; jW2_3=m; }
                                        }
                                    }
                                }
                            }//l,m
                        }
                    }
                }
            }//k.s
        }
    }//i,j



    Vec_i ret;
    ret.push_back(iW1_2); // jet1 onshell V        iW1_2 or  iW1_1
    ret.push_back(jW1_2); // jet2 onshell V        jW1_2 or  jW1_1
    ret.push_back(iW2_2); // jet1 offshell V       iW2_2 or  iW2_1
    ret.push_back(jW2_2); // jet2 offshell V       jW2_2 or  jW2_1
    ret.push_back(iZ_2); // jet1 production Z      iZ_2 or  iZ_1
    ret.push_back(jZ_2); // jet2 production Z      jZ_2 or  jZ_1
    return ret;
}


//Vec_tlv pairing_WW(Vec_tlv J, float target=80.385) {
Vec_i pairing_ZZ_flavor_orig(Vec_tlv J, Vec_f recojet_isB, Vec_f recojet_isC, Vec_f recojet_isS) {

    int nJets6 = 6;
    std::vector<std::pair<float, size_t>> chi2_1;
    Vec_i idxs_1_1, idxs_2_1;
    Vec_f fs_b_1, fs_c_1, fs_s_1; // flavor sum

    // first select the 2 jets that form the Z with combined score > 1
    int k = 0;
    for(int i = 0; i < nJets6; i++) { // i
        for(int j=i+1; j < nJets6; j++) { // j
            auto Z = J[i] + J[j];
            float chi2_t = std::pow(Z.M() - 91.2, 2);
            float flavor_sum_B = recojet_isB[i] + recojet_isB[j];
            float flavor_sum_C = recojet_isC[i] + recojet_isC[j];
            float flavor_sum_S = recojet_isS[i] + recojet_isS[j];
            chi2_1.emplace_back(chi2_t, k);
            fs_b_1.push_back(flavor_sum_B);
            fs_c_1.push_back(flavor_sum_C);
            fs_s_1.push_back(flavor_sum_S);
            idxs_1_1.push_back(i);
            idxs_2_1.push_back(j);
            k += 1;
        }
    }
    


    // check if Z candidate found: sort on chi2
    std::sort(chi2_1.begin(), chi2_1.end());
    int idx_1_1 = idxs_1_1[chi2_1.at(0).second]; // by default take lowest chi2
    int idx_2_1 = idxs_2_1[chi2_1.at(0).second]; // by default take lowest chi2
    for(const auto& pair : chi2_1) { // override when flavor score > 1
        if(fs_b_1[pair.second] > 1.2 or fs_c_1[pair.second] > 1 or fs_s_1[pair.second] > 0.8) {
            idx_1_1 = idxs_1_1[pair.second];
            idx_2_1 = idxs_2_1[pair.second];
            break;
        }
    }



    // now we have one good Z, find the second one 
    // same procedure, but neglect first 2 found jets
    std::vector<std::pair<float, size_t>> chi2_2;
    Vec_i idxs_1_2, idxs_2_2;
    Vec_f fs_b_2, fs_c_2, fs_s_2; // flavor sum

    // first select the 2 jets that form the Z with combined score > 1
    k = 0;
    for(int i = 0; i < nJets6; i++) { // i
        if(i == idx_1_1 or i == idx_2_1) continue;
        for(int j=i+1; j < nJets6; j++) { // j
            if(j == idx_1_1 or j == idx_2_1) continue;
            auto Z = J[i] + J[j];
            float chi2_t = std::pow(Z.M() - 91.2, 2);
            float flavor_sum_B = recojet_isB[i] + recojet_isB[j];
            float flavor_sum_C = recojet_isC[i] + recojet_isC[j];
            float flavor_sum_S = recojet_isS[i] + recojet_isS[j];
            chi2_2.emplace_back(chi2_t, k);
            fs_b_2.push_back(flavor_sum_B);
            fs_c_2.push_back(flavor_sum_C);
            fs_s_2.push_back(flavor_sum_S);
            idxs_1_2.push_back(i);
            idxs_2_2.push_back(j);
            k += 1;
        }
    }

    // check if Z candidate found: sort on chi2
    std::sort(chi2_2.begin(), chi2_2.end());
    int idx_1_2 = idxs_1_2[chi2_2.at(0).second]; // by default take lowest chi2
    int idx_2_2 = idxs_2_2[chi2_2.at(0).second]; // by default take lowest chi2
    for(const auto& pair : chi2_2) { // override when flavor score > 1
        if(fs_b_2[pair.second] > 1.2 or fs_c_2[pair.second] > 1 or fs_s_2[pair.second] > 0.8) {
            idx_1_2 = idxs_1_2[pair.second];
            idx_2_2 = idxs_2_2[pair.second];
            break;
        }
    }


    // remaining jets --> Z*
    int idx_1_3, idx_2_3;
    for(int i = 0; i < nJets6; i++) { // i
        if(i == idx_1_1 or i == idx_2_1 or i == idx_1_2 or i == idx_2_2) continue;
        for(int j=i+1; j < nJets6; j++) { // j
            if(j == idx_1_1 or j == idx_2_1 or j == idx_1_2 or j == idx_2_2) continue;
            idx_1_3 = i;
            idx_2_3 = j;
        }
    }

    

    Vec_i ret;
    // pair Z1 or Z2 with Z* to form H
    //auto c1 = J[idx_1_1] + J[idx_2_1] + J[0] + J[0];
    auto c1 = J[idx_1_1] + J[idx_2_1] + J[idx_1_3] + J[idx_2_3];
    auto c2 = J[idx_1_2] + J[idx_2_2] + J[idx_1_3] + J[idx_2_3];
    
    //std::cout << idx_1_1 << " " << idx_2_1 << " " << idx_1_3 << " " << idx_2_3 << " " << c1.M() << std::endl;

    if(abs(c1.M() - 125) < abs(c2.M() - 125)) { // Z1 + Z* = Higgs
        ret.push_back(idx_1_1); // jet1 onshell V
        ret.push_back(idx_2_1); // jet2 onshell V
        ret.push_back(idx_1_3); // jet1 offshell V
        ret.push_back(idx_2_3); // jet2 offshell V
        ret.push_back(idx_1_2); // jet1 production Z
        ret.push_back(idx_2_2); // jet2 production Z
    }
    else { // Z2 + Z* = Higgs
        ret.push_back(idx_1_2); // jet1 onshell V
        ret.push_back(idx_2_2); // jet2 onshell V
        ret.push_back(idx_1_3); // jet1 offshell V
        ret.push_back(idx_2_3); // jet2 offshell V
        ret.push_back(idx_1_1); // jet1 production Z
        ret.push_back(idx_2_1); // jet2 production Z
    }

    

    return ret;
}


//Vec_tlv pairing_WW(Vec_tlv J, float target=80.385) {
Vec_i pairing_ZZ_flavor(Vec_tlv J, Vec_f recojet_isB, Vec_f recojet_isC, Vec_f recojet_isS) {

    TLorentzVector init, recoil;
    init.SetPxPyPzE(0, 0, 0, 240);

    int nJets6 = 6;
    std::vector<std::pair<float, size_t>> chi2_1;
    Vec_i idxs_1_1, idxs_2_1;
    Vec_f fs_b_1, fs_c_1, fs_s_1; // flavor sum

    // first select the 2 jets that form the Z with combined score > 1
    int k = 0;
    for(int i = 0; i < nJets6; i++) { // i
        for(int j=i+1; j < nJets6; j++) { // j
            auto Z = J[i] + J[j];
            auto R = init - J[i] - J[j];
            float chi2_t = std::pow(Z.M() - 91.2, 2) + std::pow(R.M() - 125.0, 2);
            chi2_1.emplace_back(chi2_t, k);
            idxs_1_1.push_back(i);
            idxs_2_1.push_back(j);
            k += 1;
        }
    }
    


    // check if Z candidate found: sort on chi2
    std::sort(chi2_1.begin(), chi2_1.end());
    int idx_1_1 = idxs_1_1[chi2_1.at(0).second]; // by default take lowest chi2
    int idx_2_1 = idxs_2_1[chi2_1.at(0).second]; // by default take lowest chi2
    /*
    for(const auto& pair : chi2_1) { // override when flavor score > 1
        if(fs_b_1[pair.second] > 1.2 or fs_c_1[pair.second] > 1 or fs_s_1[pair.second] > 0.8) {
            //idx_1_1 = idxs_1_1[pair.second];
            //idx_2_1 = idxs_2_1[pair.second];
            break;
        }
    }
    */



    // now we have one good Z, find the second one 
    // same procedure, but neglect first 2 found jets
    std::vector<std::pair<float, size_t>> chi2_2;
    Vec_i idxs_1_2, idxs_2_2;
    Vec_f fs_b_2, fs_c_2, fs_s_2; // flavor sum

    // first select the 2 jets that form the Z with combined score > 1
    k = 0;
    for(int i = 0; i < nJets6; i++) { // i
        if(i == idx_1_1 or i == idx_2_1) continue;
        for(int j=i+1; j < nJets6; j++) { // j
            if(j == idx_1_1 or j == idx_2_1) continue;
            auto Z = J[i] + J[j];
            auto R = init - J[i] - J[j];
            float chi2_t = std::pow(Z.M() - 91.2, 2) + 0.0*std::pow(R.M() - 125.0, 2);
            chi2_2.emplace_back(chi2_t, k);
            idxs_1_2.push_back(i);
            idxs_2_2.push_back(j);
            k += 1;
        }
    }

    // check if Z candidate found: sort on chi2
    std::sort(chi2_2.begin(), chi2_2.end());
    int idx_1_2 = idxs_1_2[chi2_2.at(0).second]; // by default take lowest chi2
    int idx_2_2 = idxs_2_2[chi2_2.at(0).second]; // by default take lowest chi2
    /*
    for(const auto& pair : chi2_2) { // override when flavor score > 1
        if(fs_b_2[pair.second] > 1.2 or fs_c_2[pair.second] > 1 or fs_s_2[pair.second] > 0.8) {
            //idx_1_2 = idxs_1_2[pair.second];
            //idx_2_2 = idxs_2_2[pair.second];
            break;
        }
    }
    */


    // remaining jets --> Z*
    int idx_1_3, idx_2_3;
    for(int i = 0; i < nJets6; i++) { // i
        if(i == idx_1_1 or i == idx_2_1 or i == idx_1_2 or i == idx_2_2) continue;
        for(int j=i+1; j < nJets6; j++) { // j
            if(j == idx_1_1 or j == idx_2_1 or j == idx_1_2 or j == idx_2_2) continue;
            idx_1_3 = i;
            idx_2_3 = j;
        }
    }

    

    Vec_i ret;
    // pair Z1 or Z2 with Z* to form H
    //auto c1 = J[idx_1_1] + J[idx_2_1] + J[0] + J[0];
    auto c1 = J[idx_1_1] + J[idx_2_1] + J[idx_1_3] + J[idx_2_3];
    auto c2 = J[idx_1_2] + J[idx_2_2] + J[idx_1_3] + J[idx_2_3];
    
    //std::cout << idx_1_1 << " " << idx_2_1 << " " << idx_1_3 << " " << idx_2_3 << " " << c1.M() << std::endl;

    if(abs(c1.M() - 125) < abs(c2.M() - 125)) { // Z1 + Z* = Higgs
        ret.push_back(idx_1_1); // jet1 onshell V
        ret.push_back(idx_2_1); // jet2 onshell V
        ret.push_back(idx_1_3); // jet1 offshell V
        ret.push_back(idx_2_3); // jet2 offshell V
        ret.push_back(idx_1_2); // jet1 production Z
        ret.push_back(idx_2_2); // jet2 production Z
    }
    else { // Z2 + Z* = Higgs
        ret.push_back(idx_1_2); // jet1 onshell V
        ret.push_back(idx_2_2); // jet2 onshell V
        ret.push_back(idx_1_3); // jet1 offshell V
        ret.push_back(idx_2_3); // jet2 offshell V
        ret.push_back(idx_1_1); // jet1 production Z
        ret.push_back(idx_2_1); // jet2 production Z
    }

    

    return ret;
}


//Vec_tlv pairing_WW(Vec_tlv J, float target=80.385) {
Vec_i pairing_ZZd(Vec_tlv J, Vec_f recojet_isB, Vec_f recojet_isC, Vec_f recojet_isS) {

    int nJets6 = 6;
    float chi2_min = 9e99;

    // first select the 2 jets that form the Z with combined score > 1
    int idx_1_1 = -1, idx_2_1 = -1;
    for(int i = 0; i < nJets6; i++) { // i
        for(int j=i+1; j < nJets6; j++) { // j
            auto Z = J[i] + J[j];
            float chi2_t = std::pow(Z.M() - 91.2, 2);
            float flavor_sum_B = recojet_isB[i] + recojet_isB[j];
            float flavor_sum_C = recojet_isC[i] + recojet_isC[j];
            float flavor_sum_S = recojet_isS[i] + recojet_isS[j];
            if(flavor_sum_B > 1 or flavor_sum_C > 1 or flavor_sum_S > 1) {
                if(chi2_t < chi2_min) {
                    chi2_min = chi2_t;
                    idx_1_1 = i;
                    idx_2_1 = j;
                }
            }
        }
    }
    if(idx_1_1 == -1) {
        std::cout << "JET PAIRING 1 NOT FOUND" << std::endl;
    }





    // now we have one good Z, find the second one 
    // same procedure, but neglect first 2 found jets
    int idx_1_2 = -1, idx_2_2 = -1;
    chi2_min = 9e99;
    for(int i = 0; i < nJets6; i++) { // i
        if(i == idx_1_1 or i == idx_2_1) continue;
        for(int j=i+1; j < nJets6; j++) { // j
            if(j == idx_1_1 or j == idx_2_1) continue;
            auto Z = J[i] + J[j];
            float chi2_t = std::pow(Z.M() - 91.2, 2);
            float flavor_sum_B = recojet_isB[i] + recojet_isB[j];
            float flavor_sum_C = recojet_isC[i] + recojet_isC[j];
            float flavor_sum_S = recojet_isS[i] + recojet_isS[j];
            if(flavor_sum_B > 1 or flavor_sum_C > 1 or flavor_sum_S > 1) {
                if(chi2_t < chi2_min) {
                    chi2_min = chi2_t;
                    idx_1_2 = i;
                    idx_2_2 = j;
                }
            }
        }
    }

    if(idx_1_1 == -1) {
        std::cout << "JET PAIRING 2 NOT FOUND" << std::endl;
    }

    // remaining jets --> Z*
    int idx_1_3, idx_2_3;
    for(int i = 0; i < nJets6; i++) { // i
        if(i == idx_1_1 or i == idx_2_1 or i == idx_1_2 or i == idx_2_2) continue;
        for(int j=i+1; j < nJets6; j++) { // j
            if(j == idx_1_1 or j == idx_2_1 or j == idx_1_2 or j == idx_2_2) continue;
            idx_1_3 = i;
            idx_2_3 = j;
        }
    }

    

    Vec_i ret;
    // pair Z1 or Z2 with Z* to form H
    //auto c1 = J[idx_1_1] + J[idx_2_1] + J[0] + J[0];
    auto c1 = J[idx_1_1] + J[idx_2_1] + J[idx_1_3] + J[idx_2_3];
    auto c2 = J[idx_1_2] + J[idx_2_2] + J[idx_1_3] + J[idx_2_3];
    
    //std::cout << idx_1_1 << " " << idx_2_1 << " " << idx_1_3 << " " << idx_2_3 << " " << c1.M() << std::endl;

    if(abs(c1.M() - 125) < abs(c2.M() - 125)) { // Z1 + Z* = Higgs
        ret.push_back(idx_1_1); // jet1 onshell V
        ret.push_back(idx_2_1); // jet2 onshell V
        ret.push_back(idx_1_3); // jet1 offshell V
        ret.push_back(idx_2_3); // jet2 offshell V
        ret.push_back(idx_1_2); // jet1 production Z
        ret.push_back(idx_2_2); // jet2 production Z
    }
    else { // Z2 + Z* = Higgs
        ret.push_back(idx_1_2); // jet1 onshell V
        ret.push_back(idx_2_2); // jet2 onshell V
        ret.push_back(idx_1_3); // jet1 offshell V
        ret.push_back(idx_2_3); // jet2 offshell V
        ret.push_back(idx_1_1); // jet1 production Z
        ret.push_back(idx_2_1); // jet2 production Z
    }


    return ret;
    //Vec_i ret = {0,1,2,3,4,5};
    //return ret;

}

}