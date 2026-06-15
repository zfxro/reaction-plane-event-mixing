//draw results from integAwaySide.C

enum RPO {IN,OUT,ALL};
string strplane[3] = {"in-plane","out-of-plane","all-plane"};

const float PI = TMath::Pi();
const int nadj = 3;
const int ntype = 3;

const int ncent = 8;
const int nsumcent = 2;

const int limsumcent[nsumcent+1] = {0,20,60};

const double efferr = .076;//from pgg023

const double v2min = -3.0;
const double v2max = 3.0;

const double ratmin = -3.0;
const double ratmax = 8.0;


const int markstyle[] = {20,21,33,47,34};

const int markcol[] = {kBlue,kRed,kGreen,kMagenta,kCyan,kYellow};

double EX;


namespace GAT {
#include "gamAnaTree.h"
}

//averages from parentptana
const double pt_missed[GAT::nptt] = {6.11, 8.37, 10.87, 14.21};
const double pt_tagged[GAT::nptt] = {6.99, 9.61, 12.20, 15.44};
                                     
//-------------------- yields --------------------
TGraphErrors* tgrixiyield[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::ndpsi+1] = {};
TGraphErrors* tgrixiyield_syserr[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::ndpsi+1] = {};
TGraphErrors* tgrixiyield_fullsyserr[nadj][ntype+1][GAT::nptt][nsumcent+1][GAT::ndpsi+1] = {};

TGraphErrors* tgripttyield[nadj][ntype+1][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};
TGraphErrors* tgripttyield_syserr[nadj][ntype+1][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};
TGraphErrors* tgripttyield_fullsyserr[nadj][ntype+1][nsumcent+1][GAT::ndpsi+1][GAT::nxi] = {};

TGraph* tgrsmear_ratio[nadj][ntype+1][GAT::nptt] = {};

//manually setting these limits
double biglims[nadj][ntype+1][nsumcent+1][GAT::nxi][2] = {};
void SetBigLims();

//confidence interval
double ci = .68;

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

//-------------------- jet v2 (from fit)--------------------
//v2jet vs xi
TGraphErrors* tgrixiv2jfromfit[nadj][ntype+1][GAT::nptt][nsumcent+1] = {};
TGraphErrors* tgrixiv2jfromfit_syserr[nadj][ntype+1][GAT::nptt][nsumcent+1] = {};

//v2jet vs ptt
TGraphErrors* tgripttv2jfromfit[nadj][ntype+1][nsumcent+1][GAT::nxi] = {};
TGraphErrors* tgripttv2jfromfit_syserr[nadj][ntype+1][nsumcent+1][GAT::nxi] = {};

//tried using TObject*& but would not work. so overloading functions instead
bool FindObject(TFile* input, TString name, TGraphErrors*& obj);
bool FindObject(TFile* input, TString name, TGraph*& obj);

void FillFullSysErr(bool);

void FillPrevResult(TGraphErrors*& tgr_prev,TGraphErrors*& tgr_prev_sys);

void drawIntegAwaySide_prevResultCompare(
    const char* inname = "iasy_dir.root",
    bool islateadj = false,
    bool isdiradj2 = true,
    const char* inname_smear = "smeared.root",
    const char* inname_drawdphi = "dphiplotsdrawn.root"
    )
{
    if (islateadj)
    {
        cout<<"Drawing graphs for adj1 dir only since islateadj flag set true"<<endl;
        cout<<"Skipping the drawing of any in or out of plane yields. If interested in these lateAdjust.C must be adjusted to include them in its output."<<endl;
        cout<<"Actually FillFullSysErr() has not been configured for the case of islateadj = true. Exiting"<<endl;
        return;
    }

    TGraphErrors* tgr_prev;
    TGraphErrors* tgr_prev_sys;
    FillPrevResult(tgr_prev,tgr_prev_sys);

    gStyle->SetOptFit();

    SetBigLims();


    TFile* input_smear = new TFile(inname_smear,"READ");
    cout<<"Reading file: "<<inname_smear<<endl;
    if (!input_smear || input_smear->IsZombie())
    {
        cerr<<"Error opening "<<inname_smear<<endl;
        return;
    }

    TFile* input_drawdphi;
    if (!islateadj && !isdiradj2)
    {
        input_drawdphi = new TFile(inname_drawdphi,"READ");
        cout<<"Reading file: "<<inname_drawdphi<<endl;
        if (!input_drawdphi || input_drawdphi->IsZombie())
        {
            cerr<<"Error opening "<<inname_drawdphi<<endl;
            return;
        }
    }
    
    TFile* input = new TFile(inname,"READ");
    cout<<"Reading file: "<<inname<<endl;
    if (!input || input->IsZombie())
    {
        cerr<<"Error opening "<<inname<<endl;
        return;
    }

    for (int iadj = 0; iadj < nadj; iadj++)     
    {
        if (!isdiradj2 && iadj==2) continue;
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            if (isdiradj2)
            {
                if (itype != GAT::DIR || iadj != 2)
                    continue;
            }

            if (islateadj && !(itype == GAT::DIR && iadj == 1)) continue;

            cout<<" GAT::spec[itype] = "<<GAT::spec[itype]<<endl;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                if (!FindObject(input,Form("ixicentrat_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt),tgrixicentrat[iadj][itype][iptt])) return;

                if (!FindObject(input,Form("ixicentrat_syserr_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt),tgrixicentrat_syserr[iadj][itype][iptt])) return;
                
                if (isdiradj2)
                {
                    if (!FindObject(input_smear,Form("ratio_adj%i_%s_ptt%i",0,GAT::spec[GAT::INC],iptt),tgrsmear_ratio[0][GAT::INC][iptt])) return;
                }
                else
                {
                    if (!FindObject(input_smear,Form("ratio_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt),tgrsmear_ratio[iadj][itype][iptt])) return;
                }
                
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    if (!isdiradj2)
                    {
                        if (!FindObject(input,Form("icentv2j_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi),tgricentv2j[iadj][itype][iptt][ixi])) return;
                    }

                    if (iptt == 0)//only get once
                    {
                        if (!FindObject(input,Form("ipttcentrat_adj%i_%s_xi%i",iadj,GAT::spec[itype],ixi),tgripttcentrat[iadj][itype][ixi])) return;                    

                        if (!FindObject(input,Form("ipttcentrat_syserr_adj%i_%s_xi%i",iadj,GAT::spec[itype],ixi),tgripttcentrat_syserr[iadj][itype][ixi])) return;
                    }
                }
                for (int icent = 0; icent < nsumcent+1; icent++)
                {
                    if (!islateadj && !isdiradj2)
                    {

                        if (!FindObject(input,Form("ixiv2j_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent),tgrixiv2j[iadj][itype][iptt][icent])) return;
                        
                        if (!FindObject(input,Form("ixiv2j_syserr_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent),tgrixiv2j_syserr[iadj][itype][iptt][icent])) return;

                        if (!FindObject(input,Form("ixiv2jlastbin_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent),tgrixiv2jlastbin[iadj][itype][iptt][icent])) return;
                        
                        if (!FindObject(input,Form("ixiv2jlastbin_syserr_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent),tgrixiv2jlastbin_syserr[iadj][itype][iptt][icent])) return;

                        if (!FindObject(input_drawdphi,Form("ixiv2jfromfit_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent),tgrixiv2jfromfit[iadj][itype][iptt][icent])) return;
                        
                        if (!FindObject(input_drawdphi,Form("ixiv2jfromfit_syserr_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent),tgrixiv2jfromfit_syserr[iadj][itype][iptt][icent])) return;
                    }
                    else
                    {

                        if (!FindObject(input,Form("ixiv2jfromfit_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent),tgrixiv2jfromfit[iadj][itype][iptt][icent])) return;
                        
                        if (!FindObject(input,Form("ixiv2jfromfit_syserr_adj%i_%s_ptt%i_cent%i",iadj,GAT::spec[itype],iptt,icent),tgrixiv2jfromfit_syserr[iadj][itype][iptt][icent])) return;
                    }

                    if (iptt == 0)//only get once
                    {
                        for (int ixi = 0; ixi < GAT::nxi; ixi++)
                        {
                            if (!islateadj && !isdiradj2)
                            {

                                if (!FindObject(input,Form("ipttv2j_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi),tgripttv2j[iadj][itype][icent][ixi])) return;
                                if (!FindObject(input,Form("ipttv2j_syserr_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi),tgripttv2j_syserr[iadj][itype][icent][ixi])) return;

                                if (!FindObject(input,Form("ipttv2jlastbin_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi),tgripttv2jlastbin[iadj][itype][icent][ixi])) return;
                                if (!FindObject(input,Form("ipttv2jlastbin_syserr_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi),tgripttv2jlastbin_syserr[iadj][itype][icent][ixi])) return;                            

                                if (!FindObject(input_drawdphi,Form("ipttv2jfromfit_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi),tgripttv2jfromfit[iadj][itype][icent][ixi])) return;
                                if (!FindObject(input_drawdphi,Form("ipttv2jfromfit_syserr_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi),tgripttv2jfromfit_syserr[iadj][itype][icent][ixi])) return;
                                
                            }
                            else
                            {

                                printf("Will try to find ipttv2jfromfit_adj%i_%s_cent%i_xi%i\n",iadj,GAT::spec[itype],icent,ixi);
                                if (!FindObject(input,Form("ipttv2jfromfit_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi),tgripttv2jfromfit[iadj][itype][icent][ixi])) return;

                                printf("Will try to find ipttv2jfromfit_syserr_adj%i_%s_cent%i_xi%i\n",iadj,GAT::spec[itype],icent,ixi);
                                if (!FindObject(input,Form("ipttv2jfromfit_syserr_adj%i_%s_cent%i_xi%i",iadj,GAT::spec[itype],icent,ixi),tgripttv2jfromfit_syserr[iadj][itype][icent][ixi])) return;
                            }
                        }
                    }
                    
                    for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                    {
                        if (islateadj && idpsi != GAT::ndpsi) continue;

                        if (!FindObject(input,Form("ixiyield_adj%i_%s_ptt%i_cent%i_dpsi%i",iadj,GAT::spec[itype],iptt,icent,idpsi),tgrixiyield[iadj][itype][iptt][icent][idpsi])) return;
                        if (!FindObject(input,Form("ixiyield_syserr_adj%i_%s_ptt%i_cent%i_dpsi%i",iadj,GAT::spec[itype],iptt,icent,idpsi),tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi])) return;

                        for (int ixi = 0; ixi < tgrixiyield[iadj][itype][iptt][icent][idpsi]->GetN(); ixi++)
                        {
                            double xioff = .04;
                            if (iptt == 0) xioff = -.04;
                            double xibw = GAT::limxi[ixi+1]-GAT::limxi[ixi];
                            double ximid = (GAT::limxi[ixi+1]+GAT::limxi[ixi])/2.;
                            double xpos = ximid + xioff*xibw;
                            tgrixiyield[iadj][itype][iptt][icent][idpsi]->SetPoint(ixi,xpos,tgrixiyield[iadj][itype][iptt][icent][idpsi]->GetY()[ixi]);
                            tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi]->SetPoint(ixi,xpos,tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi]->GetY()[ixi]);
                        }
                        
                        for (int ixi = 0; ixi < GAT::nxi; ixi++)
                        {
                            //only get once
                            if (iptt == 0)
                            {
                                if (!FindObject(input,Form("ipttyield_adj%i_%s_cent%i_dpsi%i_xi%i",iadj,GAT::spec[itype],icent,idpsi,ixi),tgripttyield[iadj][itype][icent][idpsi][ixi])) return;
                                if (!FindObject(input,Form("ipttyield_syserr_adj%i_%s_cent%i_dpsi%i_xi%i",iadj,GAT::spec[itype],icent,idpsi,ixi),tgripttyield_syserr[iadj][itype][icent][idpsi][ixi])) return;
                            }
                            
                        }
                    }
                }
            }
        }
    }

    //add smearing, hadeff, and bg sub sys err in quadrature. save as error in *_fullsyserr[...]
    FillFullSysErr(isdiradj2);

    //Draw integrated away side yield
    //function of xi
    gStyle->SetEndErrorSize(8);
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        if (!isdiradj2 && iadj==2) continue;
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            if (islateadj && !(itype == GAT::DIR && iadj == 1)) continue;

            if (isdiradj2)
            {
                if (itype != GAT::DIR || iadj != 2)
                    continue;
            }

            
            for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
            {
                if (islateadj && idpsi != GAT::ndpsi) continue;
                
                TCanvas* can[nsumcent+1];
                TMultiGraph* mgr[nsumcent+1];
                double ymaxer = numeric_limits<double>::min();
                double yminer = numeric_limits<double>::max();
                for (int icent = 0; icent < nsumcent+1; icent++)
                {
                    can[icent] = new TCanvas(Form("can%i",icent),"",800,900);
                    // can[icent]->SetGrid(0,1);
                    can[icent]->SetLogy();
                    mgr[icent] = new TMultiGraph();
                    TLegend* leg = new TLegend(0.25,0.2,0.89,0.35);
                    leg->SetTextSize(.035);
                    leg->SetBorderSize(0);
                    leg->SetMargin(.1);
                    for (int iptt = 1; iptt >= 0; iptt--)               
                    {
                        TGraphErrors* tgr = tgrixiyield[iadj][itype][iptt][icent][idpsi];

                        tgr->SetMarkerColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                        tgr->SetLineColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                
                        tgr->SetMarkerSize(2.);
                        tgr->SetLineWidth(3);
                        tgr->SetMarkerStyle(markstyle[iptt]);

                        TGraphErrors* tgr_syserr = tgrixiyield_fullsyserr[iadj][itype][iptt][icent][idpsi];
                        tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)),.50);

                        mgr[icent]->Add(tgr_syserr,"2");
                        mgr[icent]->Add(tgr,"P");

                        mgr[icent]->SetTitle(tgr->GetTitle());
                    }

                    tgr_prev->SetMarkerColor(TColor::GetColorPalette(220*2/(GAT::nptt-1)));
                    tgr_prev->SetLineColor(TColor::GetColorPalette(220*2/(GAT::nptt-1)));
                
                    tgr_prev->SetMarkerSize(2.);
                    tgr_prev->SetLineWidth(3);
                    tgr_prev->SetMarkerStyle(kCircle);

                    tgr_prev_sys->SetFillColorAlpha(TColor::GetColorPalette(220*2/(GAT::nptt-1)),.50);

                    mgr[icent]->Add(tgr_prev_sys,"2");
                    mgr[icent]->Add(tgr_prev,"P");

                    for (int iptt = 0; iptt < 2; iptt++)
                    {
                        TGraphErrors* tgr = tgrixiyield[iadj][itype][iptt][icent][idpsi];
                        ostringstream oss;
                        leg->AddEntry(tgr,Form("%.0f #leq  p_{T} < %.0f (2014) 0-60%%  |#Delta#phi - #pi| < #pi/3",GAT::limptt[iptt],GAT::limptt[iptt+1]),"P");
                    }

                    leg->AddEntry(tgr_prev,Form("%.0f #leq  p_{T} < %.0f (2007 + 2010 + 2011) 0-40%%  |#Delta#phi - #pi| < #pi/2",5.,9.),"P");
            
                    mgr[icent]->Draw("A");
                

                    double ymin = mgr[icent]->GetHistogram()->GetYaxis()->GetXmin();
                    double ymax = mgr[icent]->GetHistogram()->GetYaxis()->GetXmax();
                    // mgr[icent]->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*5.);

                    ymaxer = max(ymax,ymaxer);
                    yminer = min(ymin,yminer);
            
                    mgr[icent]->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                    mgr[icent]->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                    mgr[icent]->GetHistogram()->GetXaxis()->SetLimits(0.,2.4);
                    mgr[icent]->GetHistogram()->GetXaxis()->SetNdivisions(-3);
                
                    leg->Draw();

                    double xmin = mgr[icent]->GetHistogram()->GetXaxis()->GetXmin();
                    double xmax = mgr[icent]->GetHistogram()->GetXaxis()->GetXmax();
                    TLine* line = new TLine(xmin,0.0,xmax,0.0);
                    line->SetLineStyle(7);
                    line->Draw();

                    // {
                    // }
                    // {
                    // }
                }

                for (int icent = 0; icent < nsumcent+1; icent++)
                {
                    if (icent != nsumcent) continue;
                    // mgr[icent]->GetHistogram()->GetYaxis()->SetRangeUser(yminer,ymaxer+(ymaxer-yminer)*5.);

                    //make wider version too. maybe this will work.
                    can[icent]->SetWindowSize(1600,900);
                    can[icent]->Modified();
                    can[icent]->SaveAs(Form("integpics_prevResultCompare/y_wide_adj%i_%s_cent%i_dpsi%i.pdf",iadj,GAT::spec[itype],icent,idpsi));
                }
            }
        }
    }

    cout<<"Drawing Complete!"<<endl;
}

void SetBigLims()
{
    //---------- PI0 ----------
    biglims[0][GAT::PI0][nsumcent][0][0] = .016;
    biglims[0][GAT::PI0][nsumcent][0][1] = .057;
    biglims[0][GAT::PI0][nsumcent][1][0] = .125;
    biglims[0][GAT::PI0][nsumcent][1][1] = .245;
    biglims[0][GAT::PI0][nsumcent][2][0] = .26;
    biglims[0][GAT::PI0][nsumcent][2][1] = 1.02;

    biglims[0][GAT::PI0][1][0][0] = .028;
    biglims[0][GAT::PI0][1][0][1] = .078;
    biglims[0][GAT::PI0][1][1][0] = .15;
    biglims[0][GAT::PI0][1][1][1] = .295;
    biglims[0][GAT::PI0][1][2][0] = .32;
    biglims[0][GAT::PI0][1][2][1] = .92;

    biglims[0][GAT::PI0][0][0][0] = .006;
    biglims[0][GAT::PI0][0][0][1] = .054;
    biglims[0][GAT::PI0][0][1][0] = .065;
    biglims[0][GAT::PI0][0][1][1] = .21;
    biglims[0][GAT::PI0][0][2][0] = .15;
    biglims[0][GAT::PI0][0][2][1] = 1.35;

    //---------- DIR ----------
    biglims[2][GAT::DIR][nsumcent][0][0] = -.002;
    biglims[2][GAT::DIR][nsumcent][0][1] = .046;
    biglims[2][GAT::DIR][nsumcent][1][0] = .02;
    biglims[2][GAT::DIR][nsumcent][1][1] = .27;
    biglims[2][GAT::DIR][nsumcent][2][0] = -.02;
    biglims[2][GAT::DIR][nsumcent][2][1] = .66;

    biglims[2][GAT::DIR][1][0][0] = -.015;
    biglims[2][GAT::DIR][1][0][1] = .07;
    biglims[2][GAT::DIR][1][1][0] = .01;
    biglims[2][GAT::DIR][1][1][1] = .29;
    biglims[2][GAT::DIR][1][2][0] = -.15;
    biglims[2][GAT::DIR][1][2][1] = .95;

    biglims[2][GAT::DIR][0][0][0] = -.008;
    biglims[2][GAT::DIR][0][0][1] = .044;
    biglims[2][GAT::DIR][0][1][0] = -.03;
    biglims[2][GAT::DIR][0][1][1] = .29;
    biglims[2][GAT::DIR][0][2][0] = -.35;
    biglims[2][GAT::DIR][0][2][1] = .7;


    // old before proper bin scaling
    //---------- PI0 ----------
    // biglims[0][GAT::PI0][nsumcent][0][0] = .00675;
    // biglims[0][GAT::PI0][nsumcent][0][1] = .02175;
    // biglims[0][GAT::PI0][nsumcent][1][0] = .052;
    // biglims[0][GAT::PI0][nsumcent][1][1] = .094;
    // biglims[0][GAT::PI0][nsumcent][2][0] = .11;
    // biglims[0][GAT::PI0][nsumcent][2][1] = .39;

    // biglims[0][GAT::PI0][1][0][0] = .011;
    // biglims[0][GAT::PI0][1][0][1] = .0295;
    // biglims[0][GAT::PI0][1][1][0] = .062;
    // biglims[0][GAT::PI0][1][1][1] = .118;
    // biglims[0][GAT::PI0][1][2][0] = .13;
    // biglims[0][GAT::PI0][1][2][1] = .36;

    // biglims[0][GAT::PI0][0][0][0] = .0025;
    // biglims[0][GAT::PI0][0][0][1] = .019;
    // biglims[0][GAT::PI0][0][1][0] = .028;
    // biglims[0][GAT::PI0][0][1][1] = .084;
    // biglims[0][GAT::PI0][0][2][0] = .06;
    // biglims[0][GAT::PI0][0][2][1] = .54;
}

bool FindObject(TFile* input, TString name, TGraphErrors*& obj)
{
    
    input->GetObject(name,obj);

    if (!obj)
    {
        cerr<<"Could not find "<<name<<". Exiting..."<<endl;
        return false;
    }

    return true;
}

bool FindObject(TFile* input, TString name, TGraph*& obj)
{
    
    input->GetObject(name,obj);

    if (!obj)
    {
        cerr<<"Could not find "<<name<<". Exiting..."<<endl;
        return false;
    }

    return true;
}

void FillFullSysErr(bool isdiradj2)
{
    for (int iadj = 0; iadj < nadj; iadj++)     
    {
        if (!isdiradj2 && iadj==2) continue;
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (isdiradj2)
            {
                if (itype != GAT::DIR || iadj != 2)
                    continue;
            }

            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int icent = 0; icent < nsumcent+1; icent++)
            {
                for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                {
                    //create the tgraphs and set name
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        tgripttyield_fullsyserr[iadj][itype][icent][idpsi][ixi] = new TGraphErrors(GAT::nptt);
                        tgripttyield_fullsyserr[iadj][itype][icent][idpsi][ixi]->SetTitle(";p_{Tt};1/N_{#gamma} dN/d#xi");
                        tgripttyield_fullsyserr[iadj][itype][icent][idpsi][ixi]->SetName(Form("ipttyield_fullsyserr_adj%i_%s_cent%i_dpsi%i_xi%i",iadj,GAT::spec[itype],icent,idpsi,ixi));
                    }

                    for (int iptt = 0; iptt < GAT::nptt; iptt++)
                    {
                        tgrixiyield_fullsyserr[iadj][itype][iptt][icent][idpsi] = new TGraphErrors(GAT::nxi);
                        tgrixiyield_fullsyserr[iadj][itype][iptt][icent][idpsi]->SetTitle(";#xi;1/N_{#gamma} dN/d#xi");
                        tgrixiyield_fullsyserr[iadj][itype][iptt][icent][idpsi]->SetName(Form("ixiyield_fullsyserr_adj%i_%s_ptt%i_cent%i_dpsi%i",iadj,GAT::spec[itype],iptt,icent,idpsi));
                    }
                    
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        for (int iptt = 0; iptt < GAT::nptt; iptt++)
                        {
                            {
                                TGraphErrors* tgr = tgripttyield_syserr[iadj][itype][icent][idpsi][ixi];

                                double efer = efferr*tgr->GetY()[iptt];

                                //take care: smearing ratio binned in xi
                                double smr = -9999.;
                                if (isdiradj2)
                                    smr = fabs(tgrsmear_ratio[0][GAT::INC][iptt]->GetY()[ixi]-1.)*tgr->GetY()[iptt];
                                else
                                    smr = fabs(tgrsmear_ratio[iadj][itype][iptt]->GetY()[ixi]-1.)*tgr->GetY()[iptt];

                                double se = tgr->GetEY()[iptt];
                            
                                double tot = sqrt(se*se + smr*smr + efer*efer);
                                tgripttyield_fullsyserr[iadj][itype][icent][idpsi][ixi]->SetPoint(iptt,tgr->GetX()[iptt],tgr->GetY()[iptt]);
                                tgripttyield_fullsyserr[iadj][itype][icent][idpsi][ixi]->SetPointError(iptt,tgr->GetEX()[iptt],tot);
                            }

                            {
                                TGraphErrors* tgr = tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi];

                                double efer = efferr*tgr->GetY()[ixi];

                                //take care: smearing ratio binned in xi
                                double smr = -9999.;
                                if (isdiradj2)
                                    smr = fabs(tgrsmear_ratio[0][GAT::INC][iptt]->GetY()[ixi]-1.)*tgr->GetY()[ixi];
                                else
                                    smr = fabs(tgrsmear_ratio[iadj][itype][iptt]->GetY()[ixi]-1.)*tgr->GetY()[ixi];

                                double se = tgr->GetEY()[ixi];
                            
                                double tot = sqrt(se*se + smr*smr + efer*efer);

                                tgrixiyield_fullsyserr[iadj][itype][iptt][icent][idpsi]->SetPoint(ixi,tgr->GetX()[ixi],tgr->GetY()[ixi]);
                                tgrixiyield_fullsyserr[iadj][itype][iptt][icent][idpsi]->SetPointError(ixi,tgr->GetEX()[ixi],tot);
                            }

                        }
                    }
                }
            }
        }
    }
}

void FillPrevResult(TGraphErrors*& tgr_prev, TGraphErrors*& tgr_prev_sys)
{
    // x,AuAu,Astat,Asys,dAu,dstat,dsys
    // 0.2,0.003704,0.002088,0.003312,0.011974,0.001375,0.001428
    // 0.6,0.02772,0.005758,0.007051,0.066163,0.005219,0.006703
    // 1.0,0.13896,0.00819,0.018936,0.13464,0.01444,0.01592
    // 1.4,0.24692,0.03816,0.052152,0.24492,0.02228,0.02483
    // 1.8,0.5665,0.02696,0.04904,0.3351,0.0310,0.0389
    // 2.2,0.5166,0.09055,0.17955,0.3892,0.0289,0.1563

    vector<double> vx = {0.2,0.6,1.0,1.4,1.8,2.2};
    vector<double> vy = {0.003704,0.02772,0.13896,0.24692,0.5665,0.5166};
    vector<double> vye = {0.002088,0.005758,0.00819,0.03816,0.02696,0.09055};
    vector<double> vysys = {0.003312,0.007051,0.018936,0.052152,0.04904,0.17955};
    vector<double> vxsys (vx.size(),0.03);

    tgr_prev = new TGraphErrors(vx.size(),vx.data(),vy.data(),nullptr,vye.data());
    tgr_prev_sys = new TGraphErrors(vx.size(),vx.data(),vy.data(),vxsys.data(),vysys.data());
    
}
