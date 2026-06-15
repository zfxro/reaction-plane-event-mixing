enum RPO {IN,OUT,ALL};

const float PI = TMath::Pi();
const int nadj = 2;
const int nmix = 2;
const int ntype = 3;

const int ncent = 8;
const int nsumcent = 2;

const int markstyle[] = {20,24,21,25,22,26,23,32,34,47};

namespace GAT {
    #include "gamAnaTree.h"
}

//sys err x
double sex_xi = 0.;
double sex_ptt = 0.;
//bin widths
double xibw = 0.;

//x positions. will include offset if flag set
double xpos_xi[GAT::nptt][GAT::nxi];
double xpos_ptt[GAT::nptt][GAT::nxi];

//TGraphs:
//-------------------- yields --------------------
TGraphErrors* tgrixiyield[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::ndpsi+1] = {};
TGraphErrors* tgrixiyield_syserr[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::ndpsi+1] = {};

TGraphErrors* tgripttyield[nadj][ntype+1][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};
TGraphErrors* tgripttyield_syserr[nadj][ntype+1][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};
//-------------------- centrality ratio --------------------
//central to midcentral ratio as function of xi
TGraphErrors* tgrixicentrat[nadj][ntype+1][GAT::nptt] = {};
TGraphErrors* tgrixicentrat_syserr[nadj][ntype+1][GAT::nptt] = {};

//central to midcentral ratio as function of ptt
TGraphErrors* tgripttcentrat[nadj][ntype+1][GAT::nxi] = {};
TGraphErrors* tgripttcentrat_syserr[nadj][ntype+1][GAT::nxi] = {};

//-------------------- jet v2 --------------------
//v2jet vs xi
TGraphErrors* tgrixiv2j[nadj][ntype+1][GAT::nptt][nsumcent+1] = {};
TGraphErrors* tgrixiv2j_syserr[nadj][ntype+1][GAT::nptt][nsumcent+1] = {};

//v2jet vs ptt
TGraphErrors* tgripttv2j[nadj][ntype+1][nsumcent+1][GAT::nxi] = {};
TGraphErrors* tgripttv2j_syserr[nadj][ntype+1][nsumcent+1][GAT::nxi] = {};

//v2jet vs cent
TGraphErrors* tgricentv2j[nadj][ntype+1][GAT::nptt][GAT::nxi] = {};

//-------------------- jet v2 (last bin)--------------------
//v2jet vs xi
TGraphErrors* tgrixiv2jlastbin[nadj][ntype+1][GAT::nptt][nsumcent+1] = {};
TGraphErrors* tgrixiv2jlastbin_syserr[nadj][ntype+1][GAT::nptt][nsumcent+1] = {};

//v2jet vs ptt
TGraphErrors* tgripttv2jlastbin[nadj][ntype+1][nsumcent+1][GAT::nxi] = {};
TGraphErrors* tgripttv2jlastbin_syserr[nadj][ntype+1][nsumcent+1][GAT::nxi] = {};

//-------------------- end tgraph declarations --------------------

double yield[nadj][ntype+1][GAT::nptt][ncent][GAT::ndpsi+1][GAT::nxi] = {};
double yieldsum[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};

double yielderr[nadj][ntype+1][GAT::nptt][ncent][GAT::ndpsi+1][GAT::nxi] = {};
double yielderrsum[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};

double yieldsum_syserr[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};

double yieldlastbinsum[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};
double yieldlastbinerrsum[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};
double yieldlastbinsum_syserr[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};

double v2j[nadj][ntype+1][GAT::nptt][ncent][GAT::nxi] = {};
double v2jsum[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::nxi] = {};

double v2jerr[nadj][ntype+1][GAT::nptt][ncent][GAT::nxi] = {};
double v2jerrsum[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::nxi] = {};

double v2j_syserr[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::nxi] = {};

double v2jlastbinsum[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::nxi] = {};
double v2jlastbinerrsum[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::nxi] = {};
double v2jlastbin_syserr[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::nxi] = {};

TGraphErrors* tgrsig = nullptr;
TGraphErrors* tgrexpsig = nullptr;
TGraphErrors* tgrksig[nadj][ntype+1][GAT::nptt] = {};

TFile* input = nullptr;
TFile* inputsig = nullptr;

struct errorinfo_t {
    Int_t iadj, itype, iptt, ixi, icent, idpsi;
    Double_t ntrg;
    Double_t flate[2] = {};//flat term error
    Double_t cose[2] = {};//cosine term error
};
    
errorinfo_t errinfo[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};

int ReadFile(const char* name, TFile*& file);
void GetYield(int iadj, int itype, int iptt, int idpsi, int ixi, int icent, bool issum = false);
void FillSysError();

void integAwaySide(
    const bool offsetx = true,
    const char* inname = "weighted_bg_subtracted.root",
    const char* innamesig = "kComponents/sigmak.root",
    const char* outname = "integawaysideyield.root"
    )
{
    cerr<<"WARNING! Smearing is not yet taken into account. Please include smearing. Search for smear in drawDphiPlots.C"<<endl;
    cerr<<"WARNING! Efficiency error not taken into account. Search for efferr in drawDphiPlots.C"<<endl;
    
    if (ncent != GAT::ncent)
    {
        cerr<<"Local value of ncent, "<<ncent<<", does not match that from gamAnaTree.h "<<GAT::ncent<<endl;
        return;
    }
    
    TFile* output = new TFile(outname,"RECREATE");
    
    if (ReadFile(inname,input)<0) return;
    input->GetObject("npartsigma",tgrsig);

    FillSysError();

    //---------- set x positions and sys err x width ----------
    {
        double xioff[GAT::nptt] = {-.12,-.04,.04,.12};
        sex_xi = .035;

        if (!offsetx)
            fill_n(xioff,GAT::nptt,0.);
    
        xibw = GAT::limxi[1]-GAT::limxi[0];
        sex_xi *= xibw;

        double ptoff[GAT::nxi] = {-.09,.0,.09};
        sex_ptt = .04;

        if (!offsetx)
            fill_n(ptoff,GAT::nxi,0.);

        double pttbw = GAT::limptt[1]-GAT::limptt[0];
        sex_ptt *= pttbw;
    
        for (int iptt = 0; iptt < GAT::nptt; iptt++)
        {
            for (int ixi = 0; ixi < GAT::nxi; ixi++)
            {
                double ximid = (GAT::limxi[ixi+1]+GAT::limxi[ixi])/2.;
                double ptmid = (GAT::limptt[iptt+1]+GAT::limptt[iptt])/2.;

                xpos_xi[iptt][ixi] = ximid + xioff[iptt]*xibw;
                xpos_ptt[iptt][ixi] = ptmid + ptoff[ixi]*pttbw;
            }
        }
    }
    //---------- xi and ptt x positions have been set ----------
    
    if (ReadFile(innamesig,inputsig)<0) return;
    inputsig->GetObject("expsig",tgrexpsig);

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        ostringstream osssuf;
        osssuf<<"dphi_";
        if (iadj)
            osssuf<<"adj_";
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {

                inputsig->GetObject(Form("ksigma_Np_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt),tgrksig[iadj][itype][iptt]);
                
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {                   
                    for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                    {                   
                        for (int icent = 0; icent < ncent; icent++)
                        {
                            GetYield(iadj,itype,iptt,idpsi,ixi,icent);

                            if (idpsi!=ALL) continue;

                            double iny = yield[iadj][itype][iptt][icent][IN][ixi];
                            double inye = yielderr[iadj][itype][iptt][icent][IN][ixi];

                            double outy = yield[iadj][itype][iptt][icent][OUT][ixi];
                            double outye = yielderr[iadj][itype][iptt][icent][OUT][ixi];

                            double diff = iny - outy;
                            double error = sqrt(inye*inye + outye*outye);

                            double sigma = tgrsig->GetY()[icent];

                            diff/=sigma;
                            error/=sigma;

                            double ally = yield[iadj][itype][iptt][icent][ALL][ixi];
                            double allye = yielderr[iadj][itype][iptt][icent][ALL][ixi];
                            
                            diff = diff/ally;
                            error = error/ally;
                            error = sqrt( pow(error,2.) + pow(allye*diff/ally,2.) );

                            v2j[iadj][itype][iptt][icent][ixi] = diff;
                            v2jerr[iadj][itype][iptt][icent][ixi] = error;
                        }

                        for (int icent = 0; icent < nsumcent+1; icent++)
                        {
                            GetYield(iadj,itype,iptt,idpsi,ixi,icent,true);

                            if (idpsi!=ALL) continue;

                            double iny = yieldsum[iadj][itype][iptt][icent][IN][ixi];
                            double inye = yielderrsum[iadj][itype][iptt][icent][IN][ixi];

                            double outy = yieldsum[iadj][itype][iptt][icent][OUT][ixi];
                            double outye = yielderrsum[iadj][itype][iptt][icent][OUT][ixi];

                            double diff = iny - outy;
                            double error = sqrt(inye*inye + outye*outye);

                            double sigma = tgrexpsig->GetY()[icent];

                            sigma *= tgrksig[iadj][itype][iptt]->GetY()[icent];

                            diff/=sigma;
                            error/=sigma;

                            double ally = yieldsum[iadj][itype][iptt][icent][ALL][ixi];
                            double allye = yielderrsum[iadj][itype][iptt][icent][ALL][ixi];
                            
                            diff = diff/ally;
                            error = error/ally;
                            error = sqrt( pow(error,2.) + pow(allye*diff/ally,2.) );

                            v2jsum[iadj][itype][iptt][icent][ixi] = diff;
                            v2jerrsum[iadj][itype][iptt][icent][ixi] = error;

                            //-------------------- last bin --------------------
                            double inylastbin = yieldlastbinsum[iadj][itype][iptt][icent][IN][ixi];
                            double inylastbine = yieldlastbinerrsum[iadj][itype][iptt][icent][IN][ixi];

                            double outylastbin = yieldlastbinsum[iadj][itype][iptt][icent][OUT][ixi];
                            double outylastbine = yieldlastbinerrsum[iadj][itype][iptt][icent][OUT][ixi];

                            double difflastbin = inylastbin - outylastbin;
                            double errorlastbin = sqrt(inylastbine*inylastbine + outylastbine*outylastbine);

                            difflastbin/=sigma;
                            errorlastbin/=sigma;

                            double allylastbin = yieldlastbinsum[iadj][itype][iptt][icent][ALL][ixi];
                            double allylastbine = yieldlastbinerrsum[iadj][itype][iptt][icent][ALL][ixi];
                            
                            difflastbin = difflastbin/allylastbin;
                            errorlastbin = errorlastbin/allylastbin;
                            errorlastbin = sqrt( pow(errorlastbin,2.) + pow(allylastbine*difflastbin/allylastbin,2.) );

                            v2jlastbinsum[iadj][itype][iptt][icent][ixi] = difflastbin;
                            v2jlastbinerrsum[iadj][itype][iptt][icent][ixi] = errorlastbin;

                            //---------- systematic error calculation ----------
                            {                           
                                double ntrgin = errinfo[iadj][itype][iptt][icent][IN][ixi].ntrg;
                                double e1in = errinfo[iadj][itype][iptt][icent][IN][ixi].flate[0] * ntrgin;
                                double e2in = errinfo[iadj][itype][iptt][icent][IN][ixi].flate[1] * ntrgin;
                                double e3in = errinfo[iadj][itype][iptt][icent][IN][ixi].cose[0] * ntrgin;
                                double e4in = errinfo[iadj][itype][iptt][icent][IN][ixi].cose[1] * ntrgin;

                                double ntrgout = errinfo[iadj][itype][iptt][icent][OUT][ixi].ntrg;
                                double e1out = errinfo[iadj][itype][iptt][icent][OUT][ixi].flate[0] * ntrgout;
                                double e2out = errinfo[iadj][itype][iptt][icent][OUT][ixi].flate[1] * ntrgout;
                                double e3out = errinfo[iadj][itype][iptt][icent][OUT][ixi].cose[0] * ntrgout;
                                double e4out = errinfo[iadj][itype][iptt][icent][OUT][ixi].cose[1] * ntrgout;

                                double ntrgall = errinfo[iadj][itype][iptt][icent][ALL][ixi].ntrg;
                                double e1all = errinfo[iadj][itype][iptt][icent][ALL][ixi].flate[0] * ntrgall / 2.;
                                double e3all = errinfo[iadj][itype][iptt][icent][ALL][ixi].cose[0] * ntrgall / 2.;

                                bool agree = true;
                                if (fabs(e1in-e1all) > .01*e1all) agree = false;
                                if (fabs(e1out-e1all) > .01*e1all) agree = false;
                                if (fabs(e3in-e3all) > .01*e3all) agree = false;
                                if (fabs(e3out-e3all) > .01*e3all) agree = false;
                                if (fabs(e2in-e2out) > .01*(e2in+e2out)/2.) agree = false;
                                if (fabs(e4in-e4out) > .01*(e4in+e4out)/2.) agree = false;

                                if (!agree)
                                {
                                    printf("adj %d; type %d; ptt %d; cent %d; xi %d;\endl",iadj,itype,iptt,icent,ixi);
                                    cout<<"e1in = "<<e1in<<" e2in = "<<e2in<<" e3in = "<<e3in<<" e4in = "<<e4in<<endl;
                                    cout<<"e1out = "<<e1out<<" e2out = "<<e2out<<" e3out = "<<e3out<<" e4out = "<<e4out<<endl;
                                    cout<<"e1all = "<<e1all<<" e3all = "<<e3all<<endl;
                                    return;
                                }

                                double e1 = e1in;
                                double e2 = e2in;
                                double e3 = e3in;
                                double e4 = e4in;                               

                                double Bio = diff*ally*sigma;
                                double Ball = ally;

                                double et1 = -(1./ntrgin - 1./ntrgout)*(PI/3.) + (1./ntrgall)*(2.*PI/3.)*Bio/Ball;
                                et1 *= (e1)/(Ball*sigma);

                                double et3 = -(1./ntrgin - 1./ntrgout)*(sqrt(3.)/4.) + (1./ntrgall)*(sqrt(3.)/2.)*Bio/Ball;
                                et3 *= e3/(Ball*sigma);

                                double et2 = -(1./ntrgin + 1./ntrgout)*(PI/3.);
                                et2 *= e2/(Ball*sigma);

                                double et4 = -(1./ntrgin + 1./ntrgout)*(sqrt(3.)/4.);
                                et4 *= e4/(Ball*sigma);

                                double syserr = sqrt(et1*et1 + et2*et2 + et3*et3 + et4*et4);
                                v2j_syserr[iadj][itype][iptt][icent][ixi] = syserr;

                                //-------------------- last bin --------------------
                                {
                                    double Bio = difflastbin*allylastbin*sigma;
                                    double Ball = allylastbin;

                                    double et1 = -(1./ntrgin - 1./ntrgout)*(PI/12.) + (1./ntrgall)*(2.*PI/12.)*Bio/Ball;
                                    et1 *= (e1)/(Ball*sigma);

                                    double et3 = -(1./ntrgin - 1./ntrgout)*(1./4.) + (1./ntrgall)*(1./2.)*Bio/Ball;
                                    et3 *= e3/(Ball*sigma);

                                    double et2 = -(1./ntrgin + 1./ntrgout)*(PI/12.);
                                    et2 *= e2/(Ball*sigma);

                                    double et4 = -(1./ntrgin + 1./ntrgout)*(1./4.);
                                    et4 *= e4/(Ball*sigma);

                                    double syserr = sqrt(et1*et1 + et2*et2 + et3*et3 + et4*et4);
                                    v2jlastbin_syserr[iadj][itype][iptt][icent][ixi] = syserr;
                                }
                                
                            }                           
                        }//icent < nsumcent+1
                    }//idpsi
                }//ixi
            }//iptt
        }//itype
    }//iadj

    //Create TGraphs of interest
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    double centyield = yieldsum[iadj][itype][iptt][0][ALL][ixi];
                    double midcentyield = yieldsum[iadj][itype][iptt][1][ALL][ixi];
                    double ratio = centyield/midcentyield;

                    double e0 = yielderrsum[iadj][itype][iptt][0][ALL][ixi];
                    double e1 = yielderrsum[iadj][itype][iptt][1][ALL][ixi];
                    e0 /= midcentyield;
                    e1 = -e1 * ratio / midcentyield;

                    double error = sqrt(e0*e0 + e1*e1);

                    double syserror = -9999.;
                    {
                        double a0 = PI/3.;
                        double a1 = sqrt(3.)/4.;
                        double fe0c = errinfo[iadj][itype][iptt][0][ALL][ixi].flate[0];                
                        fe0c *= a0;
                        
                        double ce0c = errinfo[iadj][itype][iptt][0][ALL][ixi].cose[0];
                        ce0c *= a1;

                        double ntrgp = errinfo[iadj][itype][iptt][1][ALL][ixi].ntrg;
                        double fe0p = errinfo[iadj][itype][iptt][1][ALL][ixi].flate[0];
                        fe0p *= a0;
                        double ce0p = errinfo[iadj][itype][iptt][1][ALL][ixi].cose[0];
                        ce0p *= a1;

                        double EC = (fe0c*fe0c + ce0c*ce0c)/(centyield*centyield);
                        double EP = (fe0p*fe0p + ce0p*ce0p)/(midcentyield*midcentyield);
                        syserror = sqrt(ratio*ratio*(EC+EP));
                    }
                        
                    //---------- cent ratio vs xi ----------
                    if (!tgrixicentrat[iadj][itype][iptt])
                    {
                        tgrixicentrat[iadj][itype][iptt] = new TGraphErrors();
                        tgrixicentrat[iadj][itype][iptt]->SetName(Form("ixicentrat_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt));
                        tgrixicentrat[iadj][itype][iptt]->SetTitle(";#xi;Y_{0-20%}/Y_{20-60%}");

                        tgrixicentrat_syserr[iadj][itype][iptt] = new TGraphErrors();
                        tgrixicentrat_syserr[iadj][itype][iptt]->SetName(Form("ixicentrat_syserr_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt));
                        tgrixicentrat_syserr[iadj][itype][iptt]->SetTitle(";#xi;Y_{0-20%}/Y_{20-60%}");
                    }
                    
                    tgrixicentrat[iadj][itype][iptt]->SetPoint(ixi,xpos_xi[iptt][ixi],ratio);
                    tgrixicentrat[iadj][itype][iptt]->SetPointError(ixi,0,error);

                    tgrixicentrat_syserr[iadj][itype][iptt]->SetPoint(ixi,xpos_xi[iptt][ixi],ratio);
                    tgrixicentrat_syserr[iadj][itype][iptt]->SetPointError(ixi,sex_xi,syserror);

                    //---------- cent ratio vs ptt ----------
                    if (!tgripttcentrat[iadj][itype][ixi])
                    {
                        tgripttcentrat[iadj][itype][ixi] = new TGraphErrors();
                        tgripttcentrat[iadj][itype][ixi]->SetName(Form("ipttcentrat_adj%i_%s_xi%i",iadj,GAT::spec[itype],ixi));
                        tgripttcentrat[iadj][itype][ixi]->SetTitle(";p_{Tt};Y_{0-20%}/Y_{20-60%}");

                        tgripttcentrat_syserr[iadj][itype][ixi] = new TGraphErrors();
                        tgripttcentrat_syserr[iadj][itype][ixi]->SetName(Form("ipttcentrat_syserr_adj%i_%s_xi%i",iadj,GAT::spec[itype],ixi));
                        tgripttcentrat_syserr[iadj][itype][ixi]->SetTitle(";p_{Tt};Y_{0-20%}/Y_{20-60%}");
                    }
                    
                    tgripttcentrat[iadj][itype][ixi]->SetPoint(iptt,xpos_ptt[iptt][ixi],ratio);
                    tgripttcentrat[iadj][itype][ixi]->SetPointError(iptt,0,error);

                    tgripttcentrat_syserr[iadj][itype][ixi]->SetPoint(iptt,xpos_ptt[iptt][ixi],ratio);
                    tgripttcentrat_syserr[iadj][itype][ixi]->SetPointError(iptt,sex_ptt,syserror);
                    
                                                       
                    for (int icent = 0; icent < nsumcent+1; icent++)
                    {
                        for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                        {
                            //---------- yield vs xi ----------
                            if (!tgrixiyield[iadj][itype][iptt][icent][idpsi])
                            {
                                tgrixiyield[iadj][itype][iptt][icent][idpsi] = new TGraphErrors();
                                tgrixiyield[iadj][itype][iptt][icent][idpsi]->SetName(Form("ixiyield_adj%i_%s_ptt%i_cent%i_dpsi%i",iadj,GAT::spec[itype],iptt,icent,idpsi));
                                tgrixiyield[iadj][itype][iptt][icent][idpsi]->SetTitle(";#xi;1/N_{t} dN/d#xi");

                                tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi] = new TGraphErrors();
                                tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi]->SetName(Form("ixiyield_syserr_adj%i_%s_ptt%i_cent%i_dpsi%i",iadj,GAT::spec[itype],iptt,icent,idpsi));
                                tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi]->SetTitle(";#xi;1/N_{t} dN/d#xi");
                            }
                            tgrixiyield[iadj][itype][iptt][icent][idpsi]->SetPoint(ixi,xpos_xi[iptt][ixi],yieldsum[iadj][itype][iptt][icent][idpsi][ixi]/xibw);
                            tgrixiyield[iadj][itype][iptt][icent][idpsi]->SetPointError(ixi,0,yielderrsum[iadj][itype][iptt][icent][idpsi][ixi]/xibw);

                            tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi]->SetPoint(ixi,xpos_xi[iptt][ixi],yieldsum[iadj][itype][iptt][icent][idpsi][ixi]/xibw);
                            tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi]->SetPointError(ixi,sex_xi,fabs(yieldsum_syserr[iadj][itype][iptt][icent][idpsi][ixi])/xibw);

                            //---------- yield vs pt ----------
                            if (!tgripttyield[iadj][itype][icent][idpsi][ixi])
                            {
                                tgripttyield[iadj][itype][icent][idpsi][ixi] = new TGraphErrors();
                                tgripttyield[iadj][itype][icent][idpsi][ixi]->SetName(Form("ipttyield_adj%i_%s_cent%i_dpsi%i_xi%i",iadj,GAT::spec[itype],icent,idpsi,ixi));
                                tgripttyield[iadj][itype][icent][idpsi][ixi]->SetTitle(";p_{Tt};1/N_{t} dN/d#xi");

                                tgripttyield_syserr[iadj][itype][icent][idpsi][ixi] = new TGraphErrors();
                                tgripttyield_syserr[iadj][itype][icent][idpsi][ixi]->SetName(Form("ipttyield_syserr_adj%i_%s_cent%i_dpsi%i_xi%i",iadj,GAT::spec[itype],icent,idpsi,ixi));
                                tgripttyield_syserr[iadj][itype][icent][idpsi][ixi]->SetTitle(";p_{Tt};1/N_{t} dN/d#xi");
                            }
                        
                            
                            tgripttyield[iadj][itype][icent][idpsi][ixi]->SetPoint(iptt,xpos_ptt[iptt][ixi],yieldsum[iadj][itype][iptt][icent][idpsi][ixi]/xibw);
                            tgripttyield[iadj][itype][icent][idpsi][ixi]->SetPointError(iptt,0,yielderrsum[iadj][itype][iptt][icent][idpsi][ixi]/xibw);

                            tgripttyield_syserr[iadj][itype][icent][idpsi][ixi]->SetPoint(iptt,xpos_ptt[iptt][ixi],yieldsum[iadj][itype][iptt][icent][idpsi][ixi]/xibw);
                            tgripttyield_syserr[iadj][itype][icent][idpsi][ixi]->SetPointError(iptt,sex_ptt,fabs(yieldsum_syserr[iadj][itype][iptt][icent][idpsi][ixi])/xibw);

                        }

                        //--------------- v2j vs xi ---------------
                        if (!tgrixiv2j[iadj][itype][iptt][icent])
                        {
                            tgrixiv2j[iadj][itype][iptt][icent] = new TGraphErrors();
                            tgrixiv2j[iadj][itype][iptt][icent]->SetName(Form("ixiv2j_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent));
                            tgrixiv2j[iadj][itype][iptt][icent]->SetTitle(";#xi;v_{2,away-side} (incorrect)");
                        }
                        cout<<"Additional factor of PI/8 being applied to v2j here"<<endl;
                        tgrixiv2j[iadj][itype][iptt][icent]->SetPoint(ixi,xpos_xi[iptt][ixi],v2jsum[iadj][itype][iptt][icent][ixi]*(PI/8.));
                        tgrixiv2j[iadj][itype][iptt][icent]->SetPointError(ixi,0,v2jerrsum[iadj][itype][iptt][icent][ixi]*(PI/8.));

                        if (!tgrixiv2j_syserr[iadj][itype][iptt][icent])
                        {
                            tgrixiv2j_syserr[iadj][itype][iptt][icent] = new TGraphErrors();
                            tgrixiv2j_syserr[iadj][itype][iptt][icent]->SetName(Form("ixiv2j_syserr_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent));
                            tgrixiv2j_syserr[iadj][itype][iptt][icent]->SetTitle(";#xi;v_{2,away-side} (incorrect)");
                        }
                        
                        tgrixiv2j_syserr[iadj][itype][iptt][icent]->SetPoint(ixi,xpos_xi[iptt][ixi],v2jsum[iadj][itype][iptt][icent][ixi]*(PI/8.));
                        tgrixiv2j_syserr[iadj][itype][iptt][icent]->SetPointError(ixi,sex_xi,v2j_syserr[iadj][itype][iptt][icent][ixi]*(PI/8.));

                        //--------------- v2j vs ptt ---------------
                        if (!tgripttv2j[iadj][itype][icent][ixi])
                        {
                            tgripttv2j[iadj][itype][icent][ixi] = new TGraphErrors();
                            tgripttv2j[iadj][itype][icent][ixi]->SetName(Form("ipttv2j_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi));
                            tgripttv2j[iadj][itype][icent][ixi]->SetTitle(";p_{Tt};v_{2,away-side} (incorrect)");
                        }
                        cout<<"Additional factor of PI/8. being applied to v2j here"<<endl;
                        tgripttv2j[iadj][itype][icent][ixi]->SetPoint(iptt,xpos_ptt[iptt][ixi],v2jsum[iadj][itype][iptt][icent][ixi]*(PI/8.));
                        tgripttv2j[iadj][itype][icent][ixi]->SetPointError(iptt,0,v2jerrsum[iadj][itype][iptt][icent][ixi]*(PI/8.));

                        if (!tgripttv2j_syserr[iadj][itype][icent][ixi])
                        {
                            tgripttv2j_syserr[iadj][itype][icent][ixi] = new TGraphErrors();
                            tgripttv2j_syserr[iadj][itype][icent][ixi]->SetName(Form("ipttv2j_syserr_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi));
                            tgripttv2j_syserr[iadj][itype][icent][ixi]->SetTitle(";p_{Tt};v_{2,away-side} (incorrect)");
                        }
                        
                        tgripttv2j_syserr[iadj][itype][icent][ixi]->SetPoint(iptt,xpos_ptt[iptt][ixi],v2jsum[iadj][itype][iptt][icent][ixi]*(PI/8.));
                        tgripttv2j_syserr[iadj][itype][icent][ixi]->SetPointError(iptt,sex_ptt,v2j_syserr[iadj][itype][iptt][icent][ixi]*(PI/8.));

                        //--------------- v2j vs cent ---------------
                        if (!tgricentv2j[iadj][itype][iptt][ixi])
                        {
                            tgricentv2j[iadj][itype][iptt][ixi] = new TGraphErrors();
                            tgricentv2j[iadj][itype][iptt][ixi]->SetName(Form("icentv2j_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));
                            tgricentv2j[iadj][itype][iptt][ixi]->SetTitle(";icent;v_{2,away-side} (incorrect)");
                        }
                        tgricentv2j[iadj][itype][iptt][ixi]->SetPoint(icent,icent,v2jsum[iadj][itype][iptt][icent][ixi]*(PI/8.));
                        tgricentv2j[iadj][itype][iptt][ixi]->SetPointError(icent,0,v2jerrsum[iadj][itype][iptt][icent][ixi]*(PI/8.));

                        //========================= last bin =========================
                        //--------------- v2jlastbin vs xi ---------------
                        if (!tgrixiv2jlastbin[iadj][itype][iptt][icent])
                        {
                            tgrixiv2jlastbin[iadj][itype][iptt][icent] = new TGraphErrors();
                            tgrixiv2jlastbin[iadj][itype][iptt][icent]->SetName(Form("ixiv2jlastbin_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent));
                            tgrixiv2jlastbin[iadj][itype][iptt][icent]->SetTitle(";#xi;v_{2,away-side}");
                        }
                        cout<<"Additional factor of PI/8 being applied to v2jlastbin here"<<endl;
                        tgrixiv2jlastbin[iadj][itype][iptt][icent]->SetPoint(ixi,xpos_xi[iptt][ixi],v2jlastbinsum[iadj][itype][iptt][icent][ixi]*(PI/8.));
                        tgrixiv2jlastbin[iadj][itype][iptt][icent]->SetPointError(ixi,0,v2jlastbinerrsum[iadj][itype][iptt][icent][ixi]*(PI/8.));

                        if (!tgrixiv2jlastbin_syserr[iadj][itype][iptt][icent])
                        {
                            tgrixiv2jlastbin_syserr[iadj][itype][iptt][icent] = new TGraphErrors();
                            tgrixiv2jlastbin_syserr[iadj][itype][iptt][icent]->SetName(Form("ixiv2jlastbin_syserr_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent));
                            tgrixiv2jlastbin_syserr[iadj][itype][iptt][icent]->SetTitle(";#xi;v_{2,away-side}");
                        }
                        
                        tgrixiv2jlastbin_syserr[iadj][itype][iptt][icent]->SetPoint(ixi,xpos_xi[iptt][ixi],v2jlastbinsum[iadj][itype][iptt][icent][ixi]*(PI/8.));
                        tgrixiv2jlastbin_syserr[iadj][itype][iptt][icent]->SetPointError(ixi,sex_xi,v2jlastbin_syserr[iadj][itype][iptt][icent][ixi]*(PI/8.));

                        //--------------- v2jlastbin vs ptt ---------------
                        if (!tgripttv2jlastbin[iadj][itype][icent][ixi])
                        {
                            tgripttv2jlastbin[iadj][itype][icent][ixi] = new TGraphErrors();
                            tgripttv2jlastbin[iadj][itype][icent][ixi]->SetName(Form("ipttv2jlastbin_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi));
                            tgripttv2jlastbin[iadj][itype][icent][ixi]->SetTitle(";p_{Tt};v_{2,away-side}");
                        }
                        cout<<"Additional factor of PI/8. being applied to v2jlastbin here"<<endl;
                        tgripttv2jlastbin[iadj][itype][icent][ixi]->SetPoint(iptt,xpos_ptt[iptt][ixi],v2jlastbinsum[iadj][itype][iptt][icent][ixi]*(PI/8.));
                        tgripttv2jlastbin[iadj][itype][icent][ixi]->SetPointError(iptt,0,v2jlastbinerrsum[iadj][itype][iptt][icent][ixi]*(PI/8.));

                        if (!tgripttv2jlastbin_syserr[iadj][itype][icent][ixi])
                        {
                            tgripttv2jlastbin_syserr[iadj][itype][icent][ixi] = new TGraphErrors();
                            tgripttv2jlastbin_syserr[iadj][itype][icent][ixi]->SetName(Form("ipttv2jlastbin_syserr_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi));
                            tgripttv2jlastbin_syserr[iadj][itype][icent][ixi]->SetTitle(";p_{Tt};v_{2,away-side}");
                        }
                        
                        tgripttv2jlastbin_syserr[iadj][itype][icent][ixi]->SetPoint(iptt,xpos_ptt[iptt][ixi],v2jlastbinsum[iadj][itype][iptt][icent][ixi]*(PI/8.));
                        tgripttv2jlastbin_syserr[iadj][itype][icent][ixi]->SetPointError(iptt,sex_ptt,v2jlastbin_syserr[iadj][itype][iptt][icent][ixi]*(PI/8.));
                        
                    }//icent
                }//ixi
            }//iptt
        }//itype
    }//iadj

    //write tgraphs
    output->cd();
    for (int iadj = 0; iadj < nadj; iadj++)
    {   
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                tgrixicentrat[iadj][itype][iptt]->Write();
                tgrixicentrat_syserr[iadj][itype][iptt]->Write();
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    tgricentv2j[iadj][itype][iptt][ixi]->Write();
                }
                for (int icent = 0; icent < nsumcent+1; icent++)
                {
                    tgrixiv2j[iadj][itype][iptt][icent]->Write();
                    tgrixiv2j_syserr[iadj][itype][iptt][icent]->Write();
                    
                    tgrixiv2jlastbin[iadj][itype][iptt][icent]->Write();
                    tgrixiv2jlastbin_syserr[iadj][itype][iptt][icent]->Write();
                    
                    for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                    {
                        tgrixiyield[iadj][itype][iptt][icent][idpsi]->Write();
                        tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi]->Write();
                    }
                }
            }

            for (int ixi = 0; ixi < GAT::nxi; ixi++)
            {
                tgripttcentrat[iadj][itype][ixi]->Write();
                tgripttcentrat_syserr[iadj][itype][ixi]->Write();

                for (int icent = 0; icent < nsumcent+1; icent++)
                {
                    tgripttv2j[iadj][itype][icent][ixi]->Write();
                    tgripttv2j_syserr[iadj][itype][icent][ixi]->Write();

                    tgripttv2jlastbin[iadj][itype][icent][ixi]->Write();
                    tgripttv2jlastbin_syserr[iadj][itype][icent][ixi]->Write();
                    
                    for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                    {
                        tgripttyield[iadj][itype][icent][idpsi][ixi]->Write();
                        tgripttyield_syserr[iadj][itype][icent][idpsi][ixi]->Write();
                    }
                }
            }
        }
    }

    cout<<"Output written to "<<output->GetName()<<endl;
    output->Close();
    cout<<"Done!"<<endl;
    return;
}

int ReadFile(const char* name, TFile*& file)
{
    file = new TFile(name,"READ");
    if (!file)
    {
        cerr<<"Could not open "<<name<<endl;
        return -1;
    }
    if (file->IsZombie())
    {
        cerr<<"File "<<name<<" is zombie!"<<endl;
        return -2;
    }
    cout<<"Reading file "<<name<<endl;
    return 0;
}

void GetYield(int iadj, int itype, int iptt, int idpsi, int ixi, int icent, bool issum)
{
            
    ostringstream oss;
    oss<<"dphi";
    if (iadj)
        oss<<"_adj";
    oss<<"_"<<GAT::spec[itype];
    if (!issum)
        oss<<"_mix0";
    oss<<"_ptt"<<iptt<<"_";
    if (issum)
        oss<<"sum";
    oss<<"cent"<<icent<<"_dpsi"<<idpsi<<"_xi"<<ixi;

    TH1F* h1dphi = nullptr;
    input->GetObject(oss.str().c_str(),h1dphi);
    if (!h1dphi)
    {
        cerr<<"Could not find "<<oss.str()<<endl;
    }

    const int bin0 = 9;
    const double expec = .67;

    const int binlast = 12;
    const double expeclast = .92;

    double binlowedge = h1dphi->GetBinLowEdge(bin0);
    if (fabs(binlowedge - expec) > .01)
    {
        cerr<<"ERROR: binlowedge = "<<binlowedge<<endl;
        return;
    }

    int nbin = h1dphi->GetNbinsX();
    
    double error = -9999.;
    double y = h1dphi->IntegralAndError(bin0,nbin,error,"width");
    //width option multiplies by bin width.  the pi has been dropped and must be included manually...
    y*=PI;
    error*=PI;

    double lastbincont = h1dphi->GetBinContent(h1dphi->GetNbinsX());
    double lastbinerr = h1dphi->GetBinError(h1dphi->GetNbinsX());
    lastbincont *= PI;
    lastbinerr *= PI;

    if (issum)
    {
        yieldsum[iadj][itype][iptt][icent][idpsi][ixi] = y;
        yielderrsum[iadj][itype][iptt][icent][idpsi][ixi] = error;

        double fe0 = errinfo[iadj][itype][iptt][icent][idpsi][ixi].flate[0];
        double fe1 = errinfo[iadj][itype][iptt][icent][idpsi][ixi].flate[1];
        
        double ce0 = errinfo[iadj][itype][iptt][icent][idpsi][ixi].cose[0];
        double ce1 = errinfo[iadj][itype][iptt][icent][idpsi][ixi].cose[1];

        double af = PI/3.;
        double ac = sqrt(3.)/4.;

        fe0 *= af;
        fe1 *= af;

        ce0 *= ac;
        ce1 *= ac;

        double syserr = sqrt( fe0*fe0 + fe1*fe1 + ce0*ce0 + ce1*ce1 );
        
        //note: bin width scaling done at set point time
        
        yieldsum_syserr[iadj][itype][iptt][icent][idpsi][ixi] = syserr;

        //==================== last bin ====================

        yieldlastbinsum[iadj][itype][iptt][icent][idpsi][ixi] = lastbincont;
        yieldlastbinerrsum[iadj][itype][iptt][icent][idpsi][ixi] = lastbinerr;

        TF1* fsyserr = nullptr;
        input->GetObject(Form("%s_syserr0",h1dphi->GetName()),fsyserr);
        double lastbinsyserr = fsyserr->Integral(double(binlast-1.)/double(nbin),1.0);  
        if (fabs(double(binlast-1.)/double(nbin) - expeclast) > .01)
        {
            cerr<<"ERROR: last lower bin edge incorrectly calculated = "<<double(binlast-1.)/double(nbin)<<endl;
            return;
        }

        lastbinsyserr *= PI;
        yieldlastbinsum_syserr[iadj][itype][iptt][icent][idpsi][ixi] = lastbinsyserr;

    }
    else
    {
        yield[iadj][itype][iptt][icent][idpsi][ixi] = y;
        yielderr[iadj][itype][iptt][icent][idpsi][ixi] = error;
    }
}

void FillSysError()
{

    TTree* tree;
    input->GetObject("Terrorinfo",tree);
    errorinfo_t info;
    tree->SetBranchAddress("info",&info);

    int nentries = tree->GetEntries();

    for (int ientry = 0; ientry < nentries; ientry++)
    {
        tree->GetEntry(ientry);
        errinfo[info.iadj][info.itype][info.iptt][info.icent][info.idpsi][info.ixi] = info;
    }                       
}
