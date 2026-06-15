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

// // const int markstyle[8] = {20,21,22,23,29,33,34,47};
// const int markstyle[] = {20,24,21,25,22,26,23,32,34,47};
// // const int markstyle[] = {20,24,21,25,22,26,23,32,34,28,47,46};

const int markstyle[] = {20,21,33,47,34};

const int markcol[] = {kBlue,kRed,kGreen,kMagenta,kCyan,kYellow};


namespace GAT {
#include "../src/gamAnaTree/gamAnaTree.h"
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
// double ci = .95;
double ci = .68;

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

void drawIntegAwaySide(
    const char* inname = "integawaysideyield.root",
    bool islateadj = false,
    bool isdiradj2 = false,
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

    gStyle->SetOptFit();

    SetBigLims();

    // gStyle->SetTitleSize(0.04,"X");
    // gStyle->SetTitleSize(0.04,"Y");

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
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            // if (itype != GAT::PI0) continue;

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
                                // cout<<"GAT::spec["<<itype<<"] = "<<GAT::spec[itype]<<endl;

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

    //set palette
    // gStyle->SetPalette(kBird);
    
    //function of xi
    //same plot for pt
    // TCanvas* cancentrat = new TCanvas("cancentrat","",1000,1000);
    gStyle->SetEndErrorSize(10);
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        if (!isdiradj2 && iadj==2) continue;
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            // if (itype != GAT::PI0) continue;
            if (islateadj && !(itype == GAT::DIR && iadj == 1)) continue;

            if (isdiradj2)
            {
                if (itype != GAT::DIR || iadj != 2)
                    continue;
            }
            
            // TCanvas cancentrat ("cancentrat","",1000,1000);
            // TCanvas cancentrat ("cancentrat","",2000,1000);
            TCanvas cancentrat ("cancentrat","",1500,800);          
            TMultiGraph* multi = new TMultiGraph();
            TMultiGraph* multi_nosys = new TMultiGraph();
            // TLegend leg (0.15,0.7,.45,.875);
            // TLegend leg (0.3,0.775,0.9,0.9);
            TLegend leg (0.6,0.775,0.9,0.9);
            leg.SetTextSize(.035);
            leg.SetTextAlign(22);
            leg.SetNColumns(2);

            TGraphErrors* tgr_all {};

            for (int iptt = GAT::nptt-1; iptt >= 0; iptt--)
            {
                cout<<"Drawing tgrixicentrat["<<iadj<<"]["<<itype<<"]["<<iptt<<"]"<<endl;
                
                TGraphErrors* tgr = tgrixicentrat[iadj][itype][iptt];
                
                // tgr->SetMarkerColor(TColor::GetColorPalette(254*iptt/(GAT::nptt-1)));
                // tgr->SetLineColor(TColor::GetColorPalette(254*iptt/(GAT::nptt-1)));

                tgr->SetMarkerColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                tgr->SetLineColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));

                // tgr->SetMarkerColor(markcol[iptt]);
                // tgr->SetLineColor(markcol[iptt]);
                
                // tgr->SetMarkerSize(2.);
                // tgr->SetLineWidth(3);
                tgr->SetMarkerSize(2);
                tgr->SetLineWidth(2);
                tgr->SetMarkerStyle(markstyle[iptt]);

                TGraphErrors* tgr_syserr = tgrixicentrat_syserr[iadj][itype][iptt];
                tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)),.50);
                
                multi->Add(tgr_syserr,"2");//add sys error box first
                multi->Add(tgr,"P");
                multi->SetTitle(tgr->GetTitle());

                // copy of tgraph with x at actual position (not shifted for illustrative purposes)
                TGraphErrors* tgr_actualx = new TGraphErrors(*tgr);
                if (!tgr_all) tgr_all = new TGraphErrors(*tgr);
                for (int i = 0; i < GAT::nxi; i++)
                {
                    tgr_actualx->SetPoint(i,(GAT::limxi[i]+GAT::limxi[i+1])/2.,tgr->GetY()[i]);
                    tgr_actualx->SetPointError(i,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));

                    tgr_all->SetPoint(i+iptt*GAT::nxi,(GAT::limxi[i]+GAT::limxi[i+1])/2.,tgr->GetY()[i]);
                    tgr_all->SetPointError(i+iptt*GAT::nxi,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));
                }

                multi_nosys->Add(tgr_actualx,"P");
                multi_nosys->SetTitle(tgr_actualx->GetTitle());
            }

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                TGraphErrors* tgr = tgrixicentrat[iadj][itype][iptt];
                ostringstream oss;
                leg.AddEntry(tgr,Form("%.0f #leq  p_{T} < %.0f",GAT::limptt[iptt],GAT::limptt[iptt+1]),"P");
            }
                
            multi->Draw("A");

            // multi->GetHistogram()->GetYaxis()->SetRangeUser(0.25,1.75);
            double ymin = multi->GetHistogram()->GetYaxis()->GetXmin();
            double ymax = multi->GetHistogram()->GetYaxis()->GetXmax();

            ymin = (ymin < ratmin) ? ratmin : ymin;
            ymax = (ymax > ratmax) ? ratmax : ymax;         
            
            multi->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);
            
            multi->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
            multi->GetHistogram()->GetXaxis()->SetLimits(0.,2.4);
            multi->GetHistogram()->GetXaxis()->SetNdivisions(-3);
            leg.Draw();
            double xmin = multi->GetHistogram()->GetXaxis()->GetXmin();
            double xmax = multi->GetHistogram()->GetXaxis()->GetXmax();
            TLine line (xmin,1.0,xmax,1.0);
            line.SetLineStyle(7);
            line.Draw();
            cancentrat.SaveAs(Form("integpics/centrat_adj%i_%s.pdf",iadj,GAT::spec[itype]));

            //lets draw with just stat errors and fit
            //have to put all data into single tgraph and fit in order to find stats
            tgr_all->Draw("A");
            tgr_all->Fit("pol1","0");
            TF1* fit = (TF1*)tgr_all->GetFunction("pol1");
            fit->SetRange(GAT::limxi[0],GAT::limxi[GAT::nxi]);
            // fit->SetRange((GAT::limxi[0]+GAT::limxi[1])/2.,(GAT::limxi[GAT::nxi-1]+GAT::limxi[GAT::nxi])/2.);
            fit->SetLineColor(kRed-7);
            cancentrat.Modified();
            cancentrat.Update();
            TPaveStats* ps = (TPaveStats *)tgr_all->GetListOfFunctions()->FindObject("stats");
            
            multi_nosys->Draw("A");

            fit->Draw("same");

            ps->SetX1NDC(.3);
            ps->SetX2NDC(.6);
            ps->SetY1NDC(.775);
            ps->SetY2NDC(.9);
            ps->Draw();
            
            multi_nosys->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);
            
            multi_nosys->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
            multi_nosys->GetHistogram()->GetXaxis()->SetLimits(0.,2.4);
            multi_nosys->GetHistogram()->GetXaxis()->SetNdivisions(-3);
            leg.Draw();
            line.Draw();

            /*Create a TGraphErrors to hold the confidence intervals*/
            TGraphErrors* grint = new TGraphErrors(GAT::nxi);
            for (int i=0; i<GAT::nxi; i++)
                grint->SetPoint(i, (GAT::limxi[i+1]+GAT::limxi[i])/2., 0);
            (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,ci);
            // grint->Draw("psame");

            TGraph* grmin = new TGraphErrors(GAT::nxi);
            TGraph* grmax = new TGraphErrors(GAT::nxi);
            for (int i=0; i<GAT::nxi; i++) {
                grmin->SetPoint(i,grint->GetX()[i],grint->GetY()[i]-grint->GetEY()[i]);
                grmax->SetPoint(i,grint->GetX()[i],grint->GetY()[i]+grint->GetEY()[i]);
            }
            
            TGraph* grshade = new TGraph(7*2);
            for (int i=0; i<7; i++)
            {
                grshade->SetPoint(i,i*.4,grmax->Eval(i*.4));
                grshade->SetPoint(7+i,2.4-i*.4,grmin->Eval(2.4-i*.4));
            }

            grshade->SetFillStyle(1001);
            grshade->SetFillColor(kRed-7);
            grshade->SetFillColorAlpha(kRed-7,.2);
            grshade->Draw("fsame");

            multi_nosys->Draw();//bring points to front

            cancentrat.SaveAs(Form("integpics/centrat_fit_adj%i_%s.pdf",iadj,GAT::spec[itype]));

        }
    }

    //---------- copy above centrat vs xi for centrat vs ptt ----------
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        if (!isdiradj2 && iadj==2) continue;
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            // if (itype != GAT::PI0) continue;
            if (islateadj && !(itype == GAT::DIR && iadj == 1)) continue;

            if (isdiradj2)
            {
                if (itype != GAT::DIR || iadj != 2)
                    continue;
            }
            
            // cout<<"iadj "<<iadj<<" itype "<<itype<<" GAT::spec[itype] = "<<GAT::spec[itype]<<endl;
            
            
            // TCanvas cancentrat ("cancentrat","",1000,1000);
            // TCanvas cancentrat ("cancentrat","",2000,1000);
            TCanvas cancentrat ("cancentrat","",1500,800);          
            TMultiGraph* multi = new TMultiGraph();
            TMultiGraph* multi_nosys = new TMultiGraph();
            // TLegend leg (0.15,0.7,.45,.875);
            // TLegend leg (0.3,0.775,0.9,0.9);
            TLegend leg (0.6,0.775,0.9,0.9);
            leg.SetTextSize(.035);
            leg.SetTextAlign(22);
            leg.SetNColumns(2);

            TGraphErrors* tgr_all {};

            for (int ixi = GAT::nxi-1; ixi >= 0; ixi--)
            {
                TGraphErrors* tgr = tgripttcentrat[iadj][itype][ixi];

                tgr->SetMarkerColor(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)));
                tgr->SetLineColor(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)));

                tgr->SetMarkerSize(2);
                tgr->SetLineWidth(2);
                tgr->SetMarkerStyle(markstyle[ixi]);

                TGraphErrors* tgr_syserr = tgripttcentrat_syserr[iadj][itype][ixi];

                tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)),.50);               
                multi->Add(tgr_syserr,"2");//add sys error box first
                
                multi->Add(tgr,"P");

                multi->SetTitle(tgr->GetTitle());

                TGraphErrors* tgr_actualx = new TGraphErrors(*tgr);
                if (!tgr_all) tgr_all = new TGraphErrors(*tgr);
                for (int i = 0; i < GAT::nptt; i++)
                {
                    tgr_actualx->SetPoint(i,(GAT::limptt[i]+GAT::limptt[i+1])/2.,tgr->GetY()[i]);
                    tgr_all->SetPoint(i+ixi*GAT::nptt,(GAT::limptt[i]+GAT::limptt[i+1])/2.,tgr->GetY()[i]);

                    tgr_actualx->SetPointError(i,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));
                    tgr_all->SetPointError(i+ixi*GAT::nptt,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));
                }

                multi_nosys->Add(tgr_actualx,"P");
                multi_nosys->SetTitle(tgr_actualx->GetTitle());
            }

            for (int ixi = 0; ixi < GAT::nxi; ixi++)
            {
                TGraphErrors* tgr = tgripttcentrat[iadj][itype][ixi];
                ostringstream oss;
                if (ixi == GAT::nxi-1)
                    leg.AddEntry((TObject*)0, "", "");
                leg.AddEntry(tgr,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]),"P");
            }
                
            multi->Draw("A");

            // multi->GetHistogram()->GetYaxis()->SetRangeUser(0.25,1.75);
            double ymin = multi->GetHistogram()->GetYaxis()->GetXmin();
            double ymax = multi->GetHistogram()->GetYaxis()->GetXmax();

            ymin = (ymin < ratmin) ? ratmin : ymin;
            ymax = (ymax > ratmax) ? ratmax : ymax;         

            multi->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);
            
            multi->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
            multi->GetHistogram()->GetXaxis()->SetLimits(GAT::limptt[0],GAT::limptt[GAT::nptt]);
            multi->GetHistogram()->GetXaxis()->SetNdivisions(-GAT::limptt[0]+GAT::limptt[GAT::nptt]);
            int xlab[] = {2,4,6,7,9,10};
            for (auto ilab : xlab)
            {
                multi->GetHistogram()->GetXaxis()->ChangeLabel(ilab,-1,-1,-1,-1,-1," ");
            }
                    
            leg.Draw();
            double xmin = multi->GetHistogram()->GetXaxis()->GetXmin();
            double xmax = multi->GetHistogram()->GetXaxis()->GetXmax();
            TLine line (xmin,1.0,xmax,1.0);
            line.SetLineStyle(7);
            line.Draw();
            cancentrat.SaveAs(Form("integpics/centratptt_adj%i_%s.pdf",iadj,GAT::spec[itype]));

            //lets draw with just stat errors and fit
            //have to put all data into single tgraph and fit in order to find stats
            tgr_all->Draw("A");
            tgr_all->Fit("pol1","0");
            TF1* fit = (TF1*)tgr_all->GetFunction("pol1");
            fit->SetRange(GAT::limptt[0],GAT::limptt[GAT::nptt]);
            // fit->SetRange((GAT::limptt[0]+GAT::limptt[1])/2.,(GAT::limptt[GAT::nptt-1]+GAT::limptt[GAT::nptt])/2.);
            fit->SetLineColor(kRed-7);
            cancentrat.Modified();
            cancentrat.Update();
            TPaveStats* ps = (TPaveStats *)tgr_all->GetListOfFunctions()->FindObject("stats");
            
            multi_nosys->Draw("A");

            fit->Draw("same");

            ps->SetX1NDC(.3);
            ps->SetX2NDC(.6);
            ps->SetY1NDC(.775);
            ps->SetY2NDC(.9);
            ps->Draw();
            
            multi_nosys->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);
            
            multi_nosys->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
            multi_nosys->GetHistogram()->GetXaxis()->SetLimits(GAT::limptt[0],GAT::limptt[GAT::nptt]);
            multi_nosys->GetHistogram()->GetXaxis()->SetNdivisions(-GAT::limptt[0]+GAT::limptt[GAT::nptt]);
            for (auto ilab : xlab)
            {
                multi_nosys->GetHistogram()->GetXaxis()->ChangeLabel(ilab,-1,-1,-1,-1,-1," ");
            }
            leg.Draw();
            line.Draw();

            TGraphErrors* grint = new TGraphErrors(GAT::nptt);
            for (int i=0; i<GAT::nptt; i++)
                grint->SetPoint(i, (GAT::limptt[i+1]+GAT::limptt[i])/2., 0);
            (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,ci);
            // grint->Draw("psame");

            TGraph* grmin = new TGraphErrors(GAT::nptt);
            TGraph* grmax = new TGraphErrors(GAT::nptt);
            for (int i=0; i<GAT::nptt; i++) {
                grmin->SetPoint(i,grint->GetX()[i],grint->GetY()[i]-grint->GetEY()[i]);
                grmax->SetPoint(i,grint->GetX()[i],grint->GetY()[i]+grint->GetEY()[i]);
            }
            
            TGraph* grshade = new TGraph(2*(2*GAT::nptt+1));
            for (int i=0; i<2*GAT::nptt+1; i++)
            {
                double x = GAT::limptt[i/2]+(i%2)*(GAT::limptt[i/2+1]-GAT::limptt[i/2])/2.;
                cout<<"Setting point "<<i<<" to ("<<x<<","<<grmax->Eval(x)<<")"<<endl;
                grshade->SetPoint(i,x,grmax->Eval(x));
                cout<<"Setting point "<<2*GAT::nptt+1+i<<" to ("<<GAT::limptt[GAT::nptt]-(x-GAT::limptt[0])<<","<<grmin->Eval(GAT::limptt[GAT::nptt]-(x-GAT::limptt[0]))<<")"<<endl;
                grshade->SetPoint(2*GAT::nptt+1+i,GAT::limptt[GAT::nptt]-(x-GAT::limptt[0]),grmin->Eval(GAT::limptt[GAT::nptt]-(x-GAT::limptt[0])));
                
            }

            grshade->SetFillStyle(1001);
            grshade->SetFillColor(kRed-7);
            grshade->SetFillColorAlpha(kRed-7,.2);
            grshade->Draw("fsame");

            multi_nosys->Draw();//bring points to front

            cancentrat.SaveAs(Form("integpics/centratptt_fit_adj%i_%s.pdf",iadj,GAT::spec[itype]));
        }
    }
    //---------- end centrat vs ptt copied from xi ----------

    //Draw v2 
    //function of xi
    //new plot for each cent
    //same plot for pt
    gStyle->SetEndErrorSize(10);
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        if (!isdiradj2 && iadj==2) continue;
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            // if (itype != GAT::PI0) continue;
            // if (islateadj && !(itype == GAT::DIR && iadj == 1)) continue;
            if (islateadj) continue;
            if (isdiradj2) continue;

            // TCanvas canxiv2 ("canxiv2","",1000,1000);
            TCanvas canxiv2 ("canxiv2","",1500,800);

            for (int icent = 0; icent < nsumcent+1; icent++)
            {

                TMultiGraph* multi = new TMultiGraph();
                TMultiGraph* multi_nosys = new TMultiGraph();
                // TLegend leg (0.3,0.775,0.9,0.9);
                TLegend leg (0.6,0.775,0.9,0.9);
                // TLegend leg (0.3,0.8375,0.9,0.9);
                leg.SetTextSize(.035);
                leg.SetTextAlign(22);
                leg.SetNColumns(2);

                TGraphErrors* tgr_all {};

                for (int iptt = GAT::nptt-1; iptt >= 0; iptt--)
                    // for (int iptt = 1; iptt >= 0; iptt--)
                {
                    TGraphErrors* tgr = tgrixiv2j[iadj][itype][iptt][icent];
                
                    // tgr->SetMarkerColor(TColor::GetColorPalette(254*iptt/(GAT::nptt-1)));
                    // tgr->SetLineColor(TColor::GetColorPalette(254*iptt/(GAT::nptt-1)));

                    tgr->SetMarkerColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                    tgr->SetLineColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));

                    // tgr->SetMarkerColor(markcol[iptt]);
                    // tgr->SetLineColor(markcol[iptt]);
                
                    tgr->SetMarkerSize(2.);
                    tgr->SetLineWidth(2);
                    tgr->SetMarkerStyle(markstyle[iptt]);

                    TGraphErrors* tgr_syserr = tgrixiv2j_syserr[iadj][itype][iptt][icent];
                    tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)),.50);

                    multi->Add(tgr_syserr,"2");
                    multi->Add(tgr,"P");

                    multi->SetTitle(tgr->GetTitle());

                    TGraphErrors* tgr_actualx = new TGraphErrors(*tgr);
                    if (!tgr_all) tgr_all = new TGraphErrors(*tgr);
                    for (int i = 0; i < GAT::nxi; i++)
                    {
                        tgr_actualx->SetPoint(i,(GAT::limxi[i]+GAT::limxi[i+1])/2.,tgr->GetY()[i]);
                        tgr_actualx->SetPointError(i,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));

                        tgr_all->SetPoint(i+iptt*GAT::nxi,(GAT::limxi[i]+GAT::limxi[i+1])/2.,tgr->GetY()[i]);
                        tgr_all->SetPointError(i+iptt*GAT::nxi,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));
                    }

                    multi_nosys->Add(tgr_actualx,"P");
                    multi_nosys->SetTitle(tgr_actualx->GetTitle());

                }

                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                    // for (int iptt = 0; iptt < 2; iptt++)
                {
                    TGraphErrors* tgr = tgrixiv2j[iadj][itype][iptt][icent];
                    ostringstream oss;
                    leg.AddEntry(tgr,Form("%.0f #leq  p_{T} < %.0f",GAT::limptt[iptt],GAT::limptt[iptt+1]),"P");
                }

            
                multi->Draw("A");
                double ymin = multi->GetHistogram()->GetYaxis()->GetXmin();
                double ymax = multi->GetHistogram()->GetYaxis()->GetXmax();

                ymin = (ymin < v2min) ? v2min : ymin;
                ymax = (ymax > v2max) ? v2max : ymax;
                
                multi->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);

                // multi->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                multi->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi->GetHistogram()->GetXaxis()->SetLimits(0.,2.4);
                multi->GetHistogram()->GetXaxis()->SetNdivisions(-3);
                leg.Draw();
                double xmin = multi->GetHistogram()->GetXaxis()->GetXmin();
                double xmax = multi->GetHistogram()->GetXaxis()->GetXmax();
                TLine line (xmin,0.0,xmax,0.0);
                line.SetLineStyle(7);
                line.Draw();
                canxiv2.Modified();
                canxiv2.SaveAs(Form("integpics/xiv2_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));

                //lets draw with just stat errors and fit
                //have to put all data into single tgraph and fit in order to find stats
                tgr_all->Draw("A");
                tgr_all->Fit("pol1","0");
                TF1* fit = (TF1*)tgr_all->GetFunction("pol1");
                fit->SetRange(GAT::limxi[0],GAT::limxi[GAT::nxi]);
                // fit->SetRange((GAT::limxi[0]+GAT::limxi[1])/2.,(GAT::limxi[GAT::nxi]+GAT::limxi[GAT::nxi])/2.);
                fit->SetLineColor(kRed-7);
                canxiv2.Modified();
                canxiv2.Update();
                TPaveStats* ps = (TPaveStats *)tgr_all->GetListOfFunctions()->FindObject("stats");
            
                multi_nosys->Draw("A");

                fit->Draw("same");

                ps->SetX1NDC(.3);
                ps->SetX2NDC(.6);
                ps->SetY1NDC(.775);
                ps->SetY2NDC(.9);
                ps->Draw();
            
                multi_nosys->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);
            
                multi_nosys->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi_nosys->GetHistogram()->GetXaxis()->SetLimits(0.,2.4);
                multi_nosys->GetHistogram()->GetXaxis()->SetNdivisions(-3);
                leg.Draw();
                line.Draw();

                /*Create a TGraphErrors to hold the confidence intervals*/
                TGraphErrors* grint = new TGraphErrors(GAT::nxi);
                for (int i=0; i<GAT::nxi; i++)
                    grint->SetPoint(i, (GAT::limxi[i+1]+GAT::limxi[i])/2., 0);
                (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,ci);
                // grint->Draw("psame");

                TGraph* grmin = new TGraphErrors(GAT::nxi);
                TGraph* grmax = new TGraphErrors(GAT::nxi);
                for (int i=0; i<GAT::nxi; i++) {
                    grmin->SetPoint(i,grint->GetX()[i],grint->GetY()[i]-grint->GetEY()[i]);
                    grmax->SetPoint(i,grint->GetX()[i],grint->GetY()[i]+grint->GetEY()[i]);
                }
            
                TGraph* grshade = new TGraph(7*2);
                for (int i=0; i<7; i++)
                {
                    grshade->SetPoint(i,i*.4,grmax->Eval(i*.4));
                    grshade->SetPoint(7+i,2.4-i*.4,grmin->Eval(2.4-i*.4));
                }

                grshade->SetFillStyle(1001);
                grshade->SetFillColor(kRed-7);
                grshade->SetFillColorAlpha(kRed-7,.2);
                grshade->Draw("fsame");

                multi_nosys->Draw();//bring points to front

                canxiv2.SaveAs(Form("integpics/xiv2_fit_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));


            }
        }
    }

    //---------- copying v2j vs xi for v2j vs ptt ----------
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        if (!isdiradj2 && iadj==2) continue;
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            // if (itype != GAT::PI0) continue;
            // if (islateadj && !(itype == GAT::DIR && iadj == 1)) continue;
            if (islateadj) continue;
            if (isdiradj2) continue;
            
            TCanvas canpttv2 ("canpttv2","",1500,800);

            for (int icent = 0; icent < nsumcent+1; icent++)
            {

                TMultiGraph* multi = new TMultiGraph();
                TMultiGraph* multi_nosys = new TMultiGraph();

                TLegend leg (0.6,0.775,0.9,0.9);

                leg.SetTextSize(.035);
                leg.SetTextAlign(22);
                leg.SetNColumns(2);

                TGraphErrors* tgr_all {};

                for (int ixi = GAT::nxi-1; ixi >= 0; ixi--)
                {
                    TGraphErrors* tgr = tgripttv2j[iadj][itype][icent][ixi];

                    tgr->SetMarkerColor(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)));
                    tgr->SetLineColor(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)));
                
                    tgr->SetMarkerSize(2.);
                    tgr->SetLineWidth(2);
                    tgr->SetMarkerStyle(markstyle[ixi]);

                    TGraphErrors* tgr_syserr = tgripttv2j_syserr[iadj][itype][icent][ixi];
                    tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)),.50);

                    multi->Add(tgr_syserr,"2");
                    multi->Add(tgr,"P");

                    multi->SetTitle(tgr->GetTitle());

                    TGraphErrors* tgr_actualx = new TGraphErrors(*tgr);
                    if (!tgr_all) tgr_all = new TGraphErrors(*tgr);
                    for (int i = 0; i < GAT::nptt; i++)
                    {
                        tgr_actualx->SetPoint(i,(GAT::limptt[i]+GAT::limptt[i+1])/2.,tgr->GetY()[i]);
                        tgr_actualx->SetPointError(i,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));

                        tgr_all->SetPoint(i+ixi*GAT::nptt,(GAT::limptt[i]+GAT::limptt[i+1])/2.,tgr->GetY()[i]);
                        tgr_all->SetPointError(i+ixi*GAT::nptt,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));
                    }

                    multi_nosys->Add(tgr_actualx,"P");
                    multi_nosys->SetTitle(tgr_actualx->GetTitle());

                }

                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    TGraphErrors* tgr = tgripttv2j[iadj][itype][icent][ixi];
                    ostringstream oss;
                    if (ixi == GAT::nxi-1)
                        leg.AddEntry((TObject*)0, "", "");
                    leg.AddEntry(tgr,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]),"P");
                }

                multi->Draw("A");
                
                double ymin = multi->GetHistogram()->GetYaxis()->GetXmin();
                double ymax = multi->GetHistogram()->GetYaxis()->GetXmax();

                ymin = (ymin < v2min) ? v2min : ymin;
                ymax = (ymax > v2max) ? v2max : ymax;
                
                multi->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);

                multi->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi->GetHistogram()->GetXaxis()->SetLimits(GAT::limptt[0],GAT::limptt[GAT::nptt]);
                multi->GetHistogram()->GetXaxis()->SetNdivisions(-GAT::limptt[0]+GAT::limptt[GAT::nptt]);
                int xlab[] = {2,4,6,7,9,10};
                for (auto ilab : xlab)
                {
                    multi->GetHistogram()->GetXaxis()->ChangeLabel(ilab,-1,-1,-1,-1,-1," ");
                }

                leg.Draw();
                double xmin = multi->GetHistogram()->GetXaxis()->GetXmin();
                double xmax = multi->GetHistogram()->GetXaxis()->GetXmax();
                TLine line (xmin,0.0,xmax,0.0);
                line.SetLineStyle(7);
                line.Draw();
                canpttv2.Modified();
                canpttv2.SaveAs(Form("integpics/pttv2_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));

                //lets draw with just stat errors and fit
                //have to put all data into single tgraph and fit in order to find stats
                tgr_all->Draw("A");
                tgr_all->Fit("pol1","0");
                TF1* fit = (TF1*)tgr_all->GetFunction("pol1");
                fit->SetRange(GAT::limptt[0],GAT::limptt[GAT::nptt]);
                // fit->SetRange((GAT::limptt[0]+GAT::limptt[1])/2.,(GAT::limptt[GAT::nptt-1]+GAT::limptt[GAT::nptt])/2.);
                fit->SetLineColor(kRed-7);
                canpttv2.Modified();
                canpttv2.Update();
                TPaveStats* ps = (TPaveStats *)tgr_all->GetListOfFunctions()->FindObject("stats");
            
                multi_nosys->Draw("A");

                fit->Draw("same");

                ps->SetX1NDC(.3);
                ps->SetX2NDC(.6);
                ps->SetY1NDC(.775);
                ps->SetY2NDC(.9);
                ps->Draw();
            
                multi_nosys->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);
            
                multi_nosys->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi_nosys->GetHistogram()->GetXaxis()->SetLimits(GAT::limptt[0],GAT::limptt[GAT::nptt]);
                multi_nosys->GetHistogram()->GetXaxis()->SetNdivisions(-GAT::limptt[0]+GAT::limptt[GAT::nptt]);
                for (auto ilab : xlab)
                {
                    multi_nosys->GetHistogram()->GetXaxis()->ChangeLabel(ilab,-1,-1,-1,-1,-1," ");
                }
                leg.Draw();
                line.Draw();

                TGraphErrors* grint = new TGraphErrors(GAT::nptt);
                for (int i=0; i<GAT::nptt; i++)
                    grint->SetPoint(i, (GAT::limptt[i+1]+GAT::limptt[i])/2., 0);
                (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,ci);
                // grint->Draw("psame");

                TGraph* grmin = new TGraphErrors(GAT::nptt);
                TGraph* grmax = new TGraphErrors(GAT::nptt);
                for (int i=0; i<GAT::nptt; i++) {
                    grmin->SetPoint(i,grint->GetX()[i],grint->GetY()[i]-grint->GetEY()[i]);
                    grmax->SetPoint(i,grint->GetX()[i],grint->GetY()[i]+grint->GetEY()[i]);
                }
            
                TGraph* grshade = new TGraph(2*(2*GAT::nptt+1));
                for (int i=0; i<2*GAT::nptt+1; i++)
                {
                    double x = GAT::limptt[i/2]+(i%2)*(GAT::limptt[i/2+1]-GAT::limptt[i/2])/2.;
                    cout<<"Setting point "<<i<<" to ("<<x<<","<<grmax->Eval(x)<<")"<<endl;
                    grshade->SetPoint(i,x,grmax->Eval(x));
                    cout<<"Setting point "<<2*GAT::nptt+1+i<<" to ("<<GAT::limptt[GAT::nptt]-(x-GAT::limptt[0])<<","<<grmin->Eval(GAT::limptt[GAT::nptt]-(x-GAT::limptt[0]))<<")"<<endl;
                    grshade->SetPoint(2*GAT::nptt+1+i,GAT::limptt[GAT::nptt]-(x-GAT::limptt[0]),grmin->Eval(GAT::limptt[GAT::nptt]-(x-GAT::limptt[0])));
                
                }

                grshade->SetFillStyle(1001);
                grshade->SetFillColor(kRed-7);
                grshade->SetFillColorAlpha(kRed-7,.2);
                grshade->Draw("fsame");

                multi_nosys->Draw();//bring points to front

                canpttv2.SaveAs(Form("integpics/pttv2_fit_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));
            }
        }
    }
    //---------- end copied v2j vs xi for v2j vs ptt ----------

    //==================== LAST BIN =========================

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        if (!isdiradj2 && iadj==2) continue;
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            // if (itype != GAT::PI0) continue;
            // if (islateadj && !(itype == GAT::DIR && iadj == 1)) continue;
            if (islateadj) continue;
            if (isdiradj2) continue;
            
            // TCanvas canxiv2 ("canxiv2","",1000,1000);
            TCanvas canxiv2 ("canxiv2","",1500,800);

            for (int icent = 0; icent < nsumcent+1; icent++)
            {

                TMultiGraph* multi = new TMultiGraph();
                TMultiGraph* multi_nosys = new TMultiGraph();
                // TLegend leg (0.3,0.775,0.9,0.9);
                TLegend leg (0.6,0.775,0.9,0.9);
                // TLegend leg (0.3,0.8375,0.9,0.9);
                leg.SetTextSize(.035);
                leg.SetTextAlign(22);
                leg.SetNColumns(2);

                TGraphErrors* tgr_all {};

                for (int iptt = GAT::nptt-1; iptt >= 0; iptt--)
                    // for (int iptt = 1; iptt >= 0; iptt--)
                {
                    TGraphErrors* tgr = tgrixiv2jlastbin[iadj][itype][iptt][icent];
                
                    // tgr->SetMarkerColor(TColor::GetColorPalette(254*iptt/(GAT::nptt-1)));
                    // tgr->SetLineColor(TColor::GetColorPalette(254*iptt/(GAT::nptt-1)));

                    tgr->SetMarkerColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                    tgr->SetLineColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));

                    // tgr->SetMarkerColor(markcol[iptt]);
                    // tgr->SetLineColor(markcol[iptt]);
                
                    tgr->SetMarkerSize(2.);
                    tgr->SetLineWidth(2);
                    tgr->SetMarkerStyle(markstyle[iptt]);

                    TGraphErrors* tgr_syserr = tgrixiv2jlastbin_syserr[iadj][itype][iptt][icent];
                    tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)),.50);

                    multi->Add(tgr_syserr,"2");
                    multi->Add(tgr,"P");

                    multi->SetTitle(tgr->GetTitle());

                    TGraphErrors* tgr_actualx = new TGraphErrors(*tgr);
                    if (!tgr_all) tgr_all = new TGraphErrors(*tgr);
                    for (int i = 0; i < GAT::nxi; i++)
                    {
                        tgr_actualx->SetPoint(i,(GAT::limxi[i]+GAT::limxi[i+1])/2.,tgr->GetY()[i]);
                        tgr_actualx->SetPointError(i,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));

                        tgr_all->SetPoint(i+iptt*GAT::nxi,(GAT::limxi[i]+GAT::limxi[i+1])/2.,tgr->GetY()[i]);
                        tgr_all->SetPointError(i+iptt*GAT::nxi,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));
                    }

                    multi_nosys->Add(tgr_actualx,"P");
                    multi_nosys->SetTitle(tgr_actualx->GetTitle());

                }

                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                    // for (int iptt = 0; iptt < 2; iptt++)
                {
                    TGraphErrors* tgr = tgrixiv2jlastbin[iadj][itype][iptt][icent];
                    ostringstream oss;
                    leg.AddEntry(tgr,Form("%.0f #leq  p_{T} < %.0f",GAT::limptt[iptt],GAT::limptt[iptt+1]),"P");
                }

            
                multi->Draw("A");
                double ymin = multi->GetHistogram()->GetYaxis()->GetXmin();
                double ymax = multi->GetHistogram()->GetYaxis()->GetXmax();

                ymin = (ymin < v2min) ? v2min : ymin;
                ymax = (ymax > v2max) ? v2max : ymax;
                
                multi->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);

                // multi->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                multi->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi->GetHistogram()->GetXaxis()->SetLimits(0.,2.4);
                multi->GetHistogram()->GetXaxis()->SetNdivisions(-3);
                leg.Draw();
                double xmin = multi->GetHistogram()->GetXaxis()->GetXmin();
                double xmax = multi->GetHistogram()->GetXaxis()->GetXmax();
                TLine line (xmin,0.0,xmax,0.0);
                line.SetLineStyle(7);
                line.Draw();
                canxiv2.Modified();
                canxiv2.SaveAs(Form("integpics/xiv2lastbin_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));

                //lets draw with just stat errors and fit
                //have to put all data into single tgraph and fit in order to find stats
                tgr_all->Draw("A");
                tgr_all->Fit("pol1","0");
                TF1* fit = (TF1*)tgr_all->GetFunction("pol1");
                fit->SetRange(GAT::limxi[0],GAT::limxi[GAT::nxi]);
                // fit->SetRange((GAT::limxi[0]+GAT::limxi[1])/2.,(GAT::limxi[GAT::nxi]+GAT::limxi[GAT::nxi])/2.);
                fit->SetLineColor(kRed-7);
                canxiv2.Modified();
                canxiv2.Update();
                TPaveStats* ps = (TPaveStats *)tgr_all->GetListOfFunctions()->FindObject("stats");
            
                multi_nosys->Draw("A");

                fit->Draw("same");

                ps->SetX1NDC(.3);
                ps->SetX2NDC(.6);
                ps->SetY1NDC(.775);
                ps->SetY2NDC(.9);
                ps->Draw();
            
                multi_nosys->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);
            
                multi_nosys->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi_nosys->GetHistogram()->GetXaxis()->SetLimits(0.,2.4);
                multi_nosys->GetHistogram()->GetXaxis()->SetNdivisions(-3);
                leg.Draw();
                line.Draw();

                /*Create a TGraphErrors to hold the confidence intervals*/
                TGraphErrors* grint = new TGraphErrors(GAT::nxi);
                for (int i=0; i<GAT::nxi; i++)
                    grint->SetPoint(i, (GAT::limxi[i+1]+GAT::limxi[i])/2., 0);
                (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,ci);
                // grint->Draw("psame");

                TGraph* grmin = new TGraphErrors(GAT::nxi);
                TGraph* grmax = new TGraphErrors(GAT::nxi);
                for (int i=0; i<GAT::nxi; i++) {
                    grmin->SetPoint(i,grint->GetX()[i],grint->GetY()[i]-grint->GetEY()[i]);
                    grmax->SetPoint(i,grint->GetX()[i],grint->GetY()[i]+grint->GetEY()[i]);
                }
            
                TGraph* grshade = new TGraph(7*2);
                for (int i=0; i<7; i++)
                {
                    grshade->SetPoint(i,i*.4,grmax->Eval(i*.4));
                    grshade->SetPoint(7+i,2.4-i*.4,grmin->Eval(2.4-i*.4));
                }

                grshade->SetFillStyle(1001);
                grshade->SetFillColor(kRed-7);
                grshade->SetFillColorAlpha(kRed-7,.2);
                grshade->Draw("fsame");

                multi_nosys->Draw();//bring points to front

                canxiv2.SaveAs(Form("integpics/xiv2lastbin_fit_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));


            }
        }
    }

    //---------- copying v2j vs xi for v2j vs ptt ----------
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        if (!isdiradj2 && iadj==2) continue;
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            // if (itype != GAT::PI0) continue;
            // if (islateadj && !(itype == GAT::DIR && iadj == 1)) continue;
            if (islateadj) continue;
            if (isdiradj2) continue;
            
            TCanvas canpttv2 ("canpttv2","",1500,800);

            for (int icent = 0; icent < nsumcent+1; icent++)
            {

                TMultiGraph* multi = new TMultiGraph();
                TMultiGraph* multi_nosys = new TMultiGraph();

                TLegend leg (0.6,0.775,0.9,0.9);

                leg.SetTextSize(.035);
                leg.SetTextAlign(22);
                leg.SetNColumns(2);

                TGraphErrors* tgr_all {};

                for (int ixi = GAT::nxi-1; ixi >= 0; ixi--)
                {
                    TGraphErrors* tgr = tgripttv2jlastbin[iadj][itype][icent][ixi];

                    tgr->SetMarkerColor(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)));
                    tgr->SetLineColor(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)));
                
                    tgr->SetMarkerSize(2.);
                    tgr->SetLineWidth(2);
                    tgr->SetMarkerStyle(markstyle[ixi]);

                    TGraphErrors* tgr_syserr = tgripttv2jlastbin_syserr[iadj][itype][icent][ixi];
                    tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)),.50);

                    multi->Add(tgr_syserr,"2");
                    multi->Add(tgr,"P");

                    multi->SetTitle(tgr->GetTitle());

                    TGraphErrors* tgr_actualx = new TGraphErrors(*tgr);
                    if (!tgr_all) tgr_all = new TGraphErrors(*tgr);
                    for (int i = 0; i < GAT::nptt; i++)
                    {
                        tgr_actualx->SetPoint(i,(GAT::limptt[i]+GAT::limptt[i+1])/2.,tgr->GetY()[i]);
                        tgr_actualx->SetPointError(i,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));

                        tgr_all->SetPoint(i+ixi*GAT::nptt,(GAT::limptt[i]+GAT::limptt[i+1])/2.,tgr->GetY()[i]);
                        tgr_all->SetPointError(i+ixi*GAT::nptt,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));
                    }

                    multi_nosys->Add(tgr_actualx,"P");
                    multi_nosys->SetTitle(tgr_actualx->GetTitle());

                }

                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    TGraphErrors* tgr = tgripttv2jlastbin[iadj][itype][icent][ixi];
                    ostringstream oss;
                    if (ixi == GAT::nxi-1)
                        leg.AddEntry((TObject*)0, "", "");
                    leg.AddEntry(tgr,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]),"P");
                }

                multi->Draw("A");
                
                double ymin = multi->GetHistogram()->GetYaxis()->GetXmin();
                double ymax = multi->GetHistogram()->GetYaxis()->GetXmax();

                ymin = (ymin < v2min) ? v2min : ymin;
                ymax = (ymax > v2max) ? v2max : ymax;
                
                multi->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);

                multi->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi->GetHistogram()->GetXaxis()->SetLimits(GAT::limptt[0],GAT::limptt[GAT::nptt]);
                multi->GetHistogram()->GetXaxis()->SetNdivisions(-GAT::limptt[0]+GAT::limptt[GAT::nptt]);
                int xlab[] = {2,4,6,7,9,10};
                for (auto ilab : xlab)
                {
                    multi->GetHistogram()->GetXaxis()->ChangeLabel(ilab,-1,-1,-1,-1,-1," ");
                }

                leg.Draw();
                double xmin = multi->GetHistogram()->GetXaxis()->GetXmin();
                double xmax = multi->GetHistogram()->GetXaxis()->GetXmax();
                TLine line (xmin,0.0,xmax,0.0);
                line.SetLineStyle(7);
                line.Draw();
                canpttv2.Modified();
                canpttv2.SaveAs(Form("integpics/pttv2lastbin_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));

                //lets draw with just stat errors and fit
                //have to put all data into single tgraph and fit in order to find stats
                tgr_all->Draw("A");
                tgr_all->Fit("pol1","0");
                TF1* fit = (TF1*)tgr_all->GetFunction("pol1");
                fit->SetRange(GAT::limptt[0],GAT::limptt[GAT::nptt]);
                // fit->SetRange((GAT::limptt[0]+GAT::limptt[1])/2.,(GAT::limptt[GAT::nptt-1]+GAT::limptt[GAT::nptt])/2.);
                fit->SetLineColor(kRed-7);
                canpttv2.Modified();
                canpttv2.Update();
                TPaveStats* ps = (TPaveStats *)tgr_all->GetListOfFunctions()->FindObject("stats");
            
                multi_nosys->Draw("A");

                fit->Draw("same");

                ps->SetX1NDC(.3);
                ps->SetX2NDC(.6);
                ps->SetY1NDC(.775);
                ps->SetY2NDC(.9);
                ps->Draw();
            
                multi_nosys->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);
            
                multi_nosys->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi_nosys->GetHistogram()->GetXaxis()->SetLimits(GAT::limptt[0],GAT::limptt[GAT::nptt]);
                multi_nosys->GetHistogram()->GetXaxis()->SetNdivisions(-GAT::limptt[0]+GAT::limptt[GAT::nptt]);
                for (auto ilab : xlab)
                {
                    multi_nosys->GetHistogram()->GetXaxis()->ChangeLabel(ilab,-1,-1,-1,-1,-1," ");
                }
                leg.Draw();
                line.Draw();

                TGraphErrors* grint = new TGraphErrors(GAT::nptt);
                for (int i=0; i<GAT::nptt; i++)
                    grint->SetPoint(i, (GAT::limptt[i+1]+GAT::limptt[i])/2., 0);
                (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,ci);
                // grint->Draw("psame");

                TGraph* grmin = new TGraphErrors(GAT::nptt);
                TGraph* grmax = new TGraphErrors(GAT::nptt);
                for (int i=0; i<GAT::nptt; i++) {
                    grmin->SetPoint(i,grint->GetX()[i],grint->GetY()[i]-grint->GetEY()[i]);
                    grmax->SetPoint(i,grint->GetX()[i],grint->GetY()[i]+grint->GetEY()[i]);
                }
            
                TGraph* grshade = new TGraph(2*(2*GAT::nptt+1));
                for (int i=0; i<2*GAT::nptt+1; i++)
                {
                    double x = GAT::limptt[i/2]+(i%2)*(GAT::limptt[i/2+1]-GAT::limptt[i/2])/2.;
                    cout<<"Setting point "<<i<<" to ("<<x<<","<<grmax->Eval(x)<<")"<<endl;
                    grshade->SetPoint(i,x,grmax->Eval(x));
                    cout<<"Setting point "<<2*GAT::nptt+1+i<<" to ("<<GAT::limptt[GAT::nptt]-(x-GAT::limptt[0])<<","<<grmin->Eval(GAT::limptt[GAT::nptt]-(x-GAT::limptt[0]))<<")"<<endl;
                    grshade->SetPoint(2*GAT::nptt+1+i,GAT::limptt[GAT::nptt]-(x-GAT::limptt[0]),grmin->Eval(GAT::limptt[GAT::nptt]-(x-GAT::limptt[0])));
                
                }

                grshade->SetFillStyle(1001);
                grshade->SetFillColor(kRed-7);
                grshade->SetFillColorAlpha(kRed-7,.2);
                grshade->Draw("fsame");

                multi_nosys->Draw();//bring points to front

                canpttv2.SaveAs(Form("integpics/pttv2lastbin_fit_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));
            }
        }
    }
    //---------- end copied v2j vs xi for v2j vs ptt ----------
    //========================= END LAST BIN v2j =========================

    //==================== FROM FIT =========================

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        if (!isdiradj2 && iadj==2) continue;
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            // if (itype != GAT::PI0) continue;
            if (islateadj && !(itype == GAT::DIR && iadj == 1)) continue;
            if (isdiradj2)
            {
                if (itype != GAT::DIR || iadj != 2)
                    continue;
            }
            
            // TCanvas canxiv2 ("canxiv2","",1000,1000);
            TCanvas canxiv2 ("canxiv2fromfit","",1500,800);

            for (int icent = 0; icent < nsumcent+1; icent++)
            {

                TMultiGraph* multi = new TMultiGraph();
                TMultiGraph* multi_nosys = new TMultiGraph();
                // TLegend leg (0.3,0.775,0.9,0.9);
                TLegend leg (0.6,0.775,0.9,0.9);
                // TLegend leg (0.3,0.8375,0.9,0.9);
                leg.SetTextSize(.035);
                leg.SetTextAlign(22);
                leg.SetNColumns(2);

                TGraphErrors* tgr_all {};

                for (int iptt = GAT::nptt-1; iptt >= 0; iptt--)
                    // for (int iptt = 1; iptt >= 0; iptt--)
                {
                    TGraphErrors* tgr = tgrixiv2jfromfit[iadj][itype][iptt][icent];
                
                    // tgr->SetMarkerColor(TColor::GetColorPalette(254*iptt/(GAT::nptt-1)));
                    // tgr->SetLineColor(TColor::GetColorPalette(254*iptt/(GAT::nptt-1)));

                    tgr->SetMarkerColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                    tgr->SetLineColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));

                    // tgr->SetMarkerColor(markcol[iptt]);
                    // tgr->SetLineColor(markcol[iptt]);
                
                    tgr->SetMarkerSize(2.);
                    tgr->SetLineWidth(2);
                    tgr->SetMarkerStyle(markstyle[iptt]);

                    TGraphErrors* tgr_syserr = tgrixiv2jfromfit_syserr[iadj][itype][iptt][icent];
                    tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)),.50);
                    // tgr_syserr->SetFillColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));

                    multi->Add(tgr_syserr,"2");
                    multi->Add(tgr,"P");

                    multi->SetTitle(tgr->GetTitle());

                    TGraphErrors* tgr_actualx = new TGraphErrors(*tgr);
                    if (!tgr_all) tgr_all = new TGraphErrors(*tgr);
                    for (int i = 0; i < GAT::nxi; i++)
                    {
                        tgr_actualx->SetPoint(i,(GAT::limxi[i]+GAT::limxi[i+1])/2.,tgr->GetY()[i]);
                        tgr_actualx->SetPointError(i,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));

                        tgr_all->SetPoint(i+iptt*GAT::nxi,(GAT::limxi[i]+GAT::limxi[i+1])/2.,tgr->GetY()[i]);
                        tgr_all->SetPointError(i+iptt*GAT::nxi,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));
                    }

                    multi_nosys->Add(tgr_actualx,"P");
                    multi_nosys->SetTitle(tgr_actualx->GetTitle());

                }

                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                    // for (int iptt = 0; iptt < 2; iptt++)
                {
                    TGraphErrors* tgr = tgrixiv2jfromfit[iadj][itype][iptt][icent];
                    ostringstream oss;
                    leg.AddEntry(tgr,Form("%.0f #leq  p_{T} < %.0f",GAT::limptt[iptt],GAT::limptt[iptt+1]),"P");
                }

            
                multi->Draw("A");
                double ymin = multi->GetHistogram()->GetYaxis()->GetXmin();
                double ymax = multi->GetHistogram()->GetYaxis()->GetXmax();

                ymin = (ymin < v2min) ? v2min : ymin;
                ymax = (ymax > v2max) ? v2max : ymax;
                
                multi->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);

                // multi->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                multi->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi->GetHistogram()->GetXaxis()->SetLimits(0.,2.4);
                multi->GetHistogram()->GetXaxis()->SetNdivisions(-3);
                leg.Draw();
                double xmin = multi->GetHistogram()->GetXaxis()->GetXmin();
                double xmax = multi->GetHistogram()->GetXaxis()->GetXmax();
                TLine line (xmin,0.0,xmax,0.0);
                line.SetLineStyle(7);
                line.Draw();
                canxiv2.Modified();
                canxiv2.SaveAs(Form("integpics/xiv2fromfit_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));

                //lets draw with just stat errors and fit
                //have to put all data into single tgraph and fit in order to find stats
                tgr_all->Draw("A");
                tgr_all->Fit("pol1","0");
                TF1* fit = (TF1*)tgr_all->GetFunction("pol1");
                fit->SetRange(GAT::limxi[0],GAT::limxi[GAT::nxi]);
                // fit->SetRange((GAT::limxi[0]+GAT::limxi[1])/2.,(GAT::limxi[GAT::nxi-1]+GAT::limxi[GAT::nxi])/2.);
                fit->SetLineColor(kRed-7);
                canxiv2.Modified();
                canxiv2.Update();
                TPaveStats* ps = (TPaveStats *)tgr_all->GetListOfFunctions()->FindObject("stats");
            
                multi_nosys->Draw("A");

                fit->Draw("same");

                ps->SetX1NDC(.3);
                ps->SetX2NDC(.6);
                ps->SetY1NDC(.775);
                ps->SetY2NDC(.9);
                ps->Draw();
            
                multi_nosys->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);
            
                multi_nosys->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi_nosys->GetHistogram()->GetXaxis()->SetLimits(0.,2.4);
                multi_nosys->GetHistogram()->GetXaxis()->SetNdivisions(-3);
                leg.Draw();
                line.Draw();

                /*Create a TGraphErrors to hold the confidence intervals*/
                TGraphErrors* grint = new TGraphErrors(GAT::nxi);
                for (int i=0; i<GAT::nxi; i++)
                    grint->SetPoint(i, (GAT::limxi[i+1]+GAT::limxi[i])/2., 0);
                (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,ci);
                // grint->Draw("psame");

                TGraph* grmin = new TGraphErrors(GAT::nxi);
                TGraph* grmax = new TGraphErrors(GAT::nxi);
                for (int i=0; i<GAT::nxi; i++) {
                    grmin->SetPoint(i,grint->GetX()[i],grint->GetY()[i]-grint->GetEY()[i]);
                    grmax->SetPoint(i,grint->GetX()[i],grint->GetY()[i]+grint->GetEY()[i]);
                }
            
                TGraph* grshade = new TGraph(7*2);
                for (int i=0; i<7; i++)
                {
                    grshade->SetPoint(i,i*.4,grmax->Eval(i*.4));
                    grshade->SetPoint(7+i,2.4-i*.4,grmin->Eval(2.4-i*.4));
                }

                grshade->SetFillStyle(1001);
                grshade->SetFillColor(kRed-7);
                grshade->SetFillColorAlpha(kRed-7,.2);
                grshade->Draw("fsame");

                multi_nosys->Draw();//bring points to front

                canxiv2.SaveAs(Form("integpics/xiv2fromfit_fit_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));


            }
        }
    }

    //---------- copying v2j vs xi for v2j vs ptt ----------
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        if (!isdiradj2 && iadj==2) continue;
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            // if (itype != GAT::PI0) continue;
            if (islateadj && !(itype == GAT::DIR && iadj == 1)) continue;
            if (isdiradj2)
            {
                if (itype != GAT::DIR || iadj != 2)
                    continue;
            }
            
            TCanvas canpttv2 ("canpttv2","",1500,800);

            for (int icent = 0; icent < nsumcent+1; icent++)
            {

                TMultiGraph* multi = new TMultiGraph();
                TMultiGraph* multi_nosys = new TMultiGraph();

                TLegend leg (0.6,0.775,0.9,0.9);

                leg.SetTextSize(.035);
                leg.SetTextAlign(22);
                leg.SetNColumns(2);

                TGraphErrors* tgr_all {};

                for (int ixi = GAT::nxi-1; ixi >= 0; ixi--)
                {
                    TGraphErrors* tgr = tgripttv2jfromfit[iadj][itype][icent][ixi];

                    tgr->SetMarkerColor(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)));
                    tgr->SetLineColor(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)));
                
                    tgr->SetMarkerSize(2.);
                    tgr->SetLineWidth(2);
                    tgr->SetMarkerStyle(markstyle[ixi]);

                    TGraphErrors* tgr_syserr = tgripttv2jfromfit_syserr[iadj][itype][icent][ixi];
                    tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)),.50);

                    multi->Add(tgr_syserr,"2");
                    multi->Add(tgr,"P");

                    multi->SetTitle(tgr->GetTitle());

                    TGraphErrors* tgr_actualx = new TGraphErrors(*tgr);
                    if (!tgr_all) tgr_all = new TGraphErrors(*tgr);
                    for (int i = 0; i < GAT::nptt; i++)
                    {
                        tgr_actualx->SetPoint(i,(GAT::limptt[i]+GAT::limptt[i+1])/2.,tgr->GetY()[i]);
                        tgr_actualx->SetPointError(i,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));

                        tgr_all->SetPoint(i+ixi*GAT::nptt,(GAT::limptt[i]+GAT::limptt[i+1])/2.,tgr->GetY()[i]);
                        tgr_all->SetPointError(i+ixi*GAT::nptt,0,sqrt(pow(tgr->GetErrorY(i),2.)+pow(tgr_syserr->GetErrorY(i),2.)));
                    }

                    multi_nosys->Add(tgr_actualx,"P");
                    multi_nosys->SetTitle(tgr_actualx->GetTitle());

                }

                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    TGraphErrors* tgr = tgripttv2jfromfit[iadj][itype][icent][ixi];
                    ostringstream oss;
                    if (ixi == GAT::nxi-1)
                        leg.AddEntry((TObject*)0, "", "");
                    leg.AddEntry(tgr,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]),"P");
                }

                multi->Draw("A");
                
                double ymin = multi->GetHistogram()->GetYaxis()->GetXmin();
                double ymax = multi->GetHistogram()->GetYaxis()->GetXmax();

                ymin = (ymin < v2min) ? v2min : ymin;
                ymax = (ymax > v2max) ? v2max : ymax;
                
                multi->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);

                multi->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi->GetHistogram()->GetXaxis()->SetLimits(GAT::limptt[0],GAT::limptt[GAT::nptt]);
                multi->GetHistogram()->GetXaxis()->SetNdivisions(-GAT::limptt[0]+GAT::limptt[GAT::nptt]);
                int xlab[] = {2,4,6,7,9,10};
                for (auto ilab : xlab)
                {
                    multi->GetHistogram()->GetXaxis()->ChangeLabel(ilab,-1,-1,-1,-1,-1," ");
                }

                leg.Draw();
                double xmin = multi->GetHistogram()->GetXaxis()->GetXmin();
                double xmax = multi->GetHistogram()->GetXaxis()->GetXmax();
                TLine line (xmin,0.0,xmax,0.0);
                line.SetLineStyle(7);
                line.Draw();
                canpttv2.Modified();
                canpttv2.SaveAs(Form("integpics/pttv2fromfit_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));

                //lets draw with just stat errors and fit
                //have to put all data into single tgraph and fit in order to find stats
                tgr_all->Draw("A");
                tgr_all->Fit("pol1","0");
                TF1* fit = (TF1*)tgr_all->GetFunction("pol1");
                fit->SetRange(GAT::limptt[0],GAT::limptt[GAT::nptt]);
                // fit->SetRange((GAT::limptt[0]+GAT::limptt[1])/2.,(GAT::limptt[GAT::nptt-1]+GAT::limptt[GAT::nptt])/2.);
                fit->SetLineColor(kRed-7);
                canpttv2.Modified();
                canpttv2.Update();
                TPaveStats* ps = (TPaveStats *)tgr_all->GetListOfFunctions()->FindObject("stats");
            
                multi_nosys->Draw("A");

                fit->Draw("same");

                ps->SetX1NDC(.3);
                ps->SetX2NDC(.6);
                ps->SetY1NDC(.775);
                ps->SetY2NDC(.9);
                ps->Draw();
            
                multi_nosys->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.2);
            
                multi_nosys->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                multi_nosys->GetHistogram()->GetXaxis()->SetLimits(GAT::limptt[0],GAT::limptt[GAT::nptt]);
                multi_nosys->GetHistogram()->GetXaxis()->SetNdivisions(-GAT::limptt[0]+GAT::limptt[GAT::nptt]);
                for (auto ilab : xlab)
                {
                    multi_nosys->GetHistogram()->GetXaxis()->ChangeLabel(ilab,-1,-1,-1,-1,-1," ");
                }
                leg.Draw();
                line.Draw();

                TGraphErrors* grint = new TGraphErrors(GAT::nptt);
                for (int i=0; i<GAT::nptt; i++)
                    grint->SetPoint(i, (GAT::limptt[i+1]+GAT::limptt[i])/2., 0);
                (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,ci);
                // grint->Draw("psame");

                TGraph* grmin = new TGraphErrors(GAT::nptt);
                TGraph* grmax = new TGraphErrors(GAT::nptt);
                for (int i=0; i<GAT::nptt; i++) {
                    grmin->SetPoint(i,grint->GetX()[i],grint->GetY()[i]-grint->GetEY()[i]);
                    grmax->SetPoint(i,grint->GetX()[i],grint->GetY()[i]+grint->GetEY()[i]);
                }
            
                TGraph* grshade = new TGraph(2*(2*GAT::nptt+1));
                for (int i=0; i<2*GAT::nptt+1; i++)
                {
                    double x = GAT::limptt[i/2]+(i%2)*(GAT::limptt[i/2+1]-GAT::limptt[i/2])/2.;
                    cout<<"Setting point "<<i<<" to ("<<x<<","<<grmax->Eval(x)<<")"<<endl;
                    grshade->SetPoint(i,x,grmax->Eval(x));
                    cout<<"Setting point "<<2*GAT::nptt+1+i<<" to ("<<GAT::limptt[GAT::nptt]-(x-GAT::limptt[0])<<","<<grmin->Eval(GAT::limptt[GAT::nptt]-(x-GAT::limptt[0]))<<")"<<endl;
                    grshade->SetPoint(2*GAT::nptt+1+i,GAT::limptt[GAT::nptt]-(x-GAT::limptt[0]),grmin->Eval(GAT::limptt[GAT::nptt]-(x-GAT::limptt[0])));
                
                }

                grshade->SetFillStyle(1001);
                grshade->SetFillColor(kRed-7);
                grshade->SetFillColorAlpha(kRed-7,.2);
                grshade->Draw("fsame");

                multi_nosys->Draw();//bring points to front

                canpttv2.SaveAs(Form("integpics/pttv2fromfit_fit_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));
            }
        }
    }
    //---------- end copied v2j vs xi for v2j vs ptt ----------
    //========================= END FROM FIT v2j =========================



    //Draw integrated away side yield
    //function of xi
    //A) all-plane. new plot for each pt. same plot for each cent.
    //B) all-plane. new plot for each cent. same plot for each pt
    //C) in-out-all plane cases same plot. new plot for each pt, cent
    gStyle->SetEndErrorSize(5);
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        if (!isdiradj2 && iadj==2) continue;
        
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            // if (itype != GAT::PI0) continue;
            if (islateadj && !(itype == GAT::DIR && iadj == 1)) continue;

            if (isdiradj2)
            {
                if (itype != GAT::DIR || iadj != 2)
                    continue;
            }

            
            // vector<TMultiGraph> mgra (GAT::nptt,TMultiGraph());
            // vector<TMultiGraph> mgrb (GAT::ncent,TMultiGraph());
            // vector<TMultiGraph> mgrc (GAT::nptt*GAT::ncent,TMultiGraph());

            // vector<TLegend> lega (GAT::nptt,TLegend(0.15,0.7,.45,.875));
            // vector<TLegend> legb (GAT::ncent,TLegend(0.15,0.7,.45,.875));
            // vector<TLegend> legc (GAT::nptt*GAT::ncent,TLegend(0.15,0.7,.45,.875));

            //B) all-plane. new plot for each cent. same plot for each pt
            for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
            {
                if (islateadj && idpsi != GAT::ndpsi) continue;
                // if (isdiradj2 && idpsi != GAT::ndpsi) continue;
                // if (isdiradj2) break;//y vs xi not available
                
                TCanvas* can[nsumcent+1];
                TMultiGraph* mgr[nsumcent+1];
                double ymaxer = numeric_limits<double>::min();
                double yminer = numeric_limits<double>::max();
                for (int icent = 0; icent < nsumcent+1; icent++)
                {
                    // TCanvas can ("can","",1000,1000);
                    // TCanvas can ("can","",800,900);
                    // TCanvas can ("can","",1100,1450);
                    can[icent] = new TCanvas(Form("can%i",icent),"",800,900);
                    // can[icent]->SetGrid(0,1);
                    can[icent]->SetLogy();
                    // TMultiGraph* mgr = new TMultiGraph();
                    mgr[icent] = new TMultiGraph();
                    // TLegend leg (.15,.7,.45,.875);
                    // TLegend leg (.55,.15,.85,.325);
                    TLegend* leg = new TLegend(0.3,0.775,0.9,0.9);
                    // TLegend leg (0.3,0.8,0.9,0.9);
                    leg->SetTextSize(.035);
                    leg->SetTextAlign(22);
                    leg->SetNColumns(2);
                    for (int iptt = GAT::nptt-1; iptt >= 0; iptt--)             
                    {
                        TGraphErrors* tgr = tgrixiyield[iadj][itype][iptt][icent][idpsi];

                        tgr->SetMarkerColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                        tgr->SetLineColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                
                        tgr->SetMarkerSize(2.);
                        tgr->SetLineWidth(3);
                        tgr->SetMarkerStyle(markstyle[iptt]);
                        // tgr->GetYaxis()->SetTitleSize(0.5);
                        // tgr->GetXaxis()->SetTitleSize(0.5);

                        // TGraphErrors* tgr_syserr = tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi];
                        TGraphErrors* tgr_syserr = tgrixiyield_fullsyserr[iadj][itype][iptt][icent][idpsi];
                        // cout<<"Using sys errors from tgraph tgrixiyield_fullsyserr["<<iadj<<"]["<<itype<<"]["<<iptt<<"]["<<icent<<"]["<<idpsi<<"]"<<endl;
                        tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)),.50);

                        mgr[icent]->Add(tgr_syserr,"2");
                        mgr[icent]->Add(tgr,"P");

                        mgr[icent]->SetTitle(tgr->GetTitle());
                    }

                    for (int iptt = 0; iptt < GAT::nptt; iptt++)
                    {
                        TGraphErrors* tgr = tgrixiyield[iadj][itype][iptt][icent][idpsi];
                        ostringstream oss;
                        leg->AddEntry(tgr,Form("%.0f #leq  p_{T} < %.0f",GAT::limptt[iptt],GAT::limptt[iptt+1]),"P");
                    }
            
                    // mgr[icent]->GetHistogram()->GetYaxis()->SetRangeUser(1.e-2,3.);
                    mgr[icent]->Draw("A");
                
                    // if (icent == 0)
                    //     mgr[icent]->GetHistogram()->GetYaxis()->SetLimits(-1.4,2.2);
                    // else if (icent == 1)
                    //     mgr[icent]->GetHistogram()->GetYaxis()->SetLimits(-1.0,1.0);

                    double ymin = mgr[icent]->GetHistogram()->GetYaxis()->GetXmin();
                    double ymax = mgr[icent]->GetHistogram()->GetYaxis()->GetXmax();
                    mgr[icent]->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*5.);

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

                    TLatex* txtcent;
                    if (icent != nsumcent)
                    {
                        if (idpsi == IN)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV In-Plane",limsumcent[icent],limsumcent[icent+1]));
                        else if (idpsi == OUT)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV Out-of-Plane",limsumcent[icent],limsumcent[icent+1]));
                        else
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV",limsumcent[icent],limsumcent[icent+1]));
                    }
                    else
                    {
                        if (idpsi == IN)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV In-Plane",limsumcent[0],limsumcent[nsumcent]));
                        else if (idpsi == OUT)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV Out-of-Plane",limsumcent[0],limsumcent[nsumcent]));
                        else
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV",limsumcent[0],limsumcent[nsumcent]));
                    }
                    txtcent->SetTextAlign(22);
                    txtcent->SetNDC();
                    txtcent->SetTextSize(0.035);
                    txtcent->Draw();
                }

                for (int icent = 0; icent < nsumcent+1; icent++)
                {
                    mgr[icent]->GetHistogram()->GetYaxis()->SetRangeUser(yminer,ymaxer+(ymaxer-yminer)*5.);
                    can[icent]->Modified();
                    can[icent]->SaveAs(Form("integpics/y_adj%i_%s_cent%i_dpsi%i.pdf",iadj,GAT::spec[itype],icent,idpsi));

                    //make wider version too. maybe this will work.
                    can[icent]->SetWindowSize(1600,900);
                    can[icent]->Modified();
                    can[icent]->SaveAs(Form("integpics/y_wide_adj%i_%s_cent%i_dpsi%i.pdf",iadj,GAT::spec[itype],icent,idpsi));
                }
            }


            //B) all-plane. new plot for each cent. same plot for each pt
            for (int icent = 0; icent < nsumcent+1; icent++)
            {
                // if (isdiradj2) break; // dont have these plots
                
                TCanvas* can[GAT::ndpsi+1];
                TMultiGraph* mgr[GAT::ndpsi+1];
                double ymaxer = numeric_limits<double>::min();
                double yminer = numeric_limits<double>::max();

                for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                {
                    if (islateadj && idpsi != GAT::ndpsi) continue;
                    
                    // TCanvas can ("can","",1000,1000);
                    // TCanvas can ("can","",800,900);
                    // TCanvas can ("can","",1100,1450);
                    can[idpsi] = new TCanvas(Form("can_dpsi%i",idpsi),"",800,900);
                    // can[idpsi]->SetGrid(0,1);
                    can[idpsi]->SetLogy();
                    // TMultiGraph* mgr = new TMultiGraph();
                    mgr[idpsi] = new TMultiGraph();
                    // TLegend leg (.15,.7,.45,.875);
                    // TLegend leg (.55,.15,.85,.325);
                    TLegend* leg = new TLegend(0.3,0.775,0.9,0.9);
                    // TLegend leg (0.3,0.8,0.9,0.9);
                    leg->SetTextSize(.035);
                    leg->SetTextAlign(22);
                    leg->SetNColumns(2);
                    for (int iptt = GAT::nptt-1; iptt >= 0; iptt--)             
                    {
                        TGraphErrors* tgr = tgrixiyield[iadj][itype][iptt][icent][idpsi];

                        tgr->SetMarkerColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                        tgr->SetLineColor(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)));
                
                        tgr->SetMarkerSize(2.);
                        tgr->SetLineWidth(3);
                        tgr->SetMarkerStyle(markstyle[iptt]);
                        // tgr->GetYaxis()->SetTitleSize(0.5);
                        // tgr->GetXaxis()->SetTitleSize(0.5);

                        // TGraphErrors* tgr_syserr = tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi];
                        TGraphErrors* tgr_syserr = tgrixiyield_fullsyserr[iadj][itype][iptt][icent][idpsi];
                        tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*iptt/(GAT::nptt-1)),.50);

                        mgr[idpsi]->Add(tgr_syserr,"2");
                        mgr[idpsi]->Add(tgr,"P");

                        mgr[idpsi]->SetTitle(tgr->GetTitle());
                    }

                    for (int iptt = 0; iptt < GAT::nptt; iptt++)
                    {
                        TGraphErrors* tgr = tgrixiyield[iadj][itype][iptt][icent][idpsi];
                        ostringstream oss;
                        leg->AddEntry(tgr,Form("%.0f #leq  p_{T} < %.0f",GAT::limptt[iptt],GAT::limptt[iptt+1]),"P");
                    }
            
                    mgr[idpsi]->Draw("A");          

                    double ymin = mgr[idpsi]->GetHistogram()->GetYaxis()->GetXmin();
                    double ymax = mgr[idpsi]->GetHistogram()->GetYaxis()->GetXmax();
                    mgr[idpsi]->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*5.);

                    ymaxer = max(ymax,ymaxer);
                    yminer = min(ymin,yminer);
            
                    mgr[idpsi]->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                    mgr[idpsi]->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                    mgr[idpsi]->GetHistogram()->GetXaxis()->SetLimits(0.,2.4);
                    mgr[idpsi]->GetHistogram()->GetXaxis()->SetNdivisions(-3);
                
                    leg->Draw();

                    double xmin = mgr[idpsi]->GetHistogram()->GetXaxis()->GetXmin();
                    double xmax = mgr[idpsi]->GetHistogram()->GetXaxis()->GetXmax();
                    TLine* line = new TLine(xmin,0.0,xmax,0.0);
                    line->SetLineStyle(7);
                    line->Draw();

                    TLatex* txtcent;
                    if (icent != nsumcent)
                    {
                        if (idpsi == IN)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV In-Plane",limsumcent[icent],limsumcent[icent+1]));
                        else if (idpsi == OUT)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV Out-of-Plane",limsumcent[icent],limsumcent[icent+1]));
                        else
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV",limsumcent[icent],limsumcent[icent+1]));
                    }
                    else
                    {
                        if (idpsi == IN)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV In-Plane",limsumcent[0],limsumcent[nsumcent]));
                        else if (idpsi == OUT)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV Out-of-Plane",limsumcent[0],limsumcent[nsumcent]));
                        else
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV",limsumcent[0],limsumcent[nsumcent]));
                    }
                    txtcent->SetTextAlign(22);
                    txtcent->SetNDC();
                    txtcent->SetTextSize(0.035);
                    txtcent->Draw();
                }

                for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                {
                    if (islateadj && idpsi != GAT::ndpsi) continue;
                    
                    mgr[idpsi]->GetHistogram()->GetYaxis()->SetRangeUser(yminer,ymaxer+(ymaxer-yminer)*5.);
                    can[idpsi]->Modified();
                    can[idpsi]->SaveAs(Form("integpics/y_adj%i_%s_dpsi%i_cent%i.pdf",iadj,GAT::spec[itype],idpsi,icent));
                }
            }
            //finish copied section to match limits across dpsi not cent

            //---------- even worse coding. copying both above and adapting to ptt ----------
            for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
            {
                if (islateadj && idpsi != GAT::ndpsi) continue;
                // if (isdiradj2 && idpsi != GAT::ndpsi) continue;
                
                TCanvas* can[nsumcent+1];
                TMultiGraph* mgr[nsumcent+1];
                double ymaxer = numeric_limits<double>::min();
                double yminer = numeric_limits<double>::max();
                for (int icent = 0; icent < nsumcent+1; icent++)
                {
                    // TCanvas can ("can","",1000,1000);
                    // TCanvas can ("can","",800,900);
                    // TCanvas can ("can","",1100,1450);
                    can[icent] = new TCanvas(Form("can%i",icent),"",800,900);
                    // can[icent]->SetGrid(0,1);
                    can[icent]->SetLogy();
                    // TMultiGraph* mgr = new TMultiGraph();
                    mgr[icent] = new TMultiGraph();
                    // TLegend leg (.15,.7,.45,.875);
                    // TLegend leg (.55,.15,.85,.325);
                    TLegend* leg = new TLegend(0.3,0.775,0.9,0.9);
                    // TLegend leg (0.3,0.8,0.9,0.9);
                    leg->SetTextSize(.035);
                    leg->SetTextAlign(22);
                    leg->SetNColumns(2);
                    for (int ixi = GAT::nxi-1; ixi >= 0; ixi--)         
                    {
                        TGraphErrors* tgr = tgripttyield[iadj][itype][icent][idpsi][ixi];

                        tgr->SetMarkerColor(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)));
                        tgr->SetLineColor(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)));
                
                        tgr->SetMarkerSize(2.);
                        tgr->SetLineWidth(3);
                        tgr->SetMarkerStyle(markstyle[ixi]);
                        // tgr->GetYaxis()->SetTitleSize(0.5);
                        // tgr->GetXaxis()->SetTitleSize(0.5);

                        // TGraphErrors* tgr_syserr = tgripttyield_syserr[iadj][itype][icent][idpsi][ixi];
                        TGraphErrors* tgr_syserr = tgripttyield_fullsyserr[iadj][itype][icent][idpsi][ixi];

                        tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)),.50);
                        mgr[icent]->Add(tgr_syserr,"2");
                        
                        mgr[icent]->Add(tgr,"P");

                        mgr[icent]->SetTitle(tgr->GetTitle());
                    }

                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        TGraphErrors* tgr = tgripttyield[iadj][itype][icent][idpsi][ixi];
                        ostringstream oss;
                        if (ixi == GAT::nxi-1)
                            leg->AddEntry((TObject*)0, "", "");
                        leg->AddEntry(tgr,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]),"P");
                    }

            
                    // mgr[icent]->GetHistogram()->GetYaxis()->SetRangeUser(1.e-2,3.);
                    mgr[icent]->Draw("A");
                
                    // if (icent == 0)
                    //     mgr[icent]->GetHistogram()->GetYaxis()->SetLimits(-1.4,2.2);
                    // else if (icent == 1)
                    //     mgr[icent]->GetHistogram()->GetYaxis()->SetLimits(-1.0,1.0);

                    double ymin = mgr[icent]->GetHistogram()->GetYaxis()->GetXmin();
                    double ymax = mgr[icent]->GetHistogram()->GetYaxis()->GetXmax();
                    mgr[icent]->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*5.);

                    ymaxer = max(ymax,ymaxer);
                    yminer = min(ymin,yminer);
            
                    mgr[icent]->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                    mgr[icent]->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                    mgr[icent]->GetHistogram()->GetXaxis()->SetLimits(GAT::limptt[0],GAT::limptt[GAT::nptt]);
                    mgr[icent]->GetHistogram()->GetXaxis()->SetNdivisions(-GAT::limptt[0]+GAT::limptt[GAT::nptt]);
                    int xlab[] = {2,4,6,7,9,10};
                    for (auto ilab : xlab)
                    {
                        mgr[icent]->GetHistogram()->GetXaxis()->ChangeLabel(ilab,-1,-1,-1,-1,-1," ");
                    }
                
                    leg->Draw();

                    double xmin = mgr[icent]->GetHistogram()->GetXaxis()->GetXmin();
                    double xmax = mgr[icent]->GetHistogram()->GetXaxis()->GetXmax();
                    TLine* line = new TLine(xmin,0.0,xmax,0.0);
                    line->SetLineStyle(7);
                    line->Draw();

                    TLatex* txtcent;
                    if (icent != nsumcent)
                    {
                        if (idpsi == IN)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV In-Plane",limsumcent[icent],limsumcent[icent+1]));
                        else if (idpsi == OUT)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV Out-of-Plane",limsumcent[icent],limsumcent[icent+1]));
                        else
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV",limsumcent[icent],limsumcent[icent+1]));
                    }
                    else
                    {
                        if (idpsi == IN)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV In-Plane",limsumcent[0],limsumcent[nsumcent]));
                        else if (idpsi == OUT)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV Out-of-Plane",limsumcent[0],limsumcent[nsumcent]));
                        else
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV",limsumcent[0],limsumcent[nsumcent]));
                    }
                    txtcent->SetTextAlign(22);
                    txtcent->SetNDC();
                    txtcent->SetTextSize(0.035);
                    txtcent->Draw();
                }

                for (int icent = 0; icent < nsumcent+1; icent++)
                {
                    mgr[icent]->GetHistogram()->GetYaxis()->SetRangeUser(yminer,ymaxer+(ymaxer-yminer)*5.);
                    can[icent]->Modified();
                    can[icent]->SaveAs(Form("integpics/yptt_adj%i_%s_cent%i_dpsi%i.pdf",iadj,GAT::spec[itype],icent,idpsi));
                }
            }

            for (int icent = 0; icent < nsumcent+1; icent++)
            {
                for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                {
                    
                    // int icent = nsumcent;
                    // int idpsi = ALL;
                    // int itype = GAT::PI0;
                    // int iadj = 0;

                    TCanvas can ("can","",1500,800);
                    //can.SetLogy();
                    gStyle->SetOptFit(0);
                        
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        TMultiGraph* mgr = new TMultiGraph();
                        TGraphErrors* tgr = tgripttyield[iadj][itype][icent][idpsi][ixi];
                        TGraphErrors* tgr_syserr = tgripttyield_syserr[iadj][itype][icent][idpsi][ixi];

                        // TGraphErrors* tgrtmp = new TGraphErrors();
                        TGraphErrors* tgrtmp = (TGraphErrors*) tgr->Clone(); 
                        for (int iptt = 0; iptt < GAT::nptt; iptt++)
                        {
                            tgrtmp->SetPoint(iptt,(GAT::limptt[iptt]+GAT::limptt[iptt+1])/2.,tgr->GetY()[iptt]);
                            // tgrtmp->SetPoint(iptt,GAT::limptt[iptt],tgr->GetY()[iptt]);

                            double err = pow(tgr->GetEY()[iptt],2.) + pow(tgr_syserr->GetEY()[iptt],2.);
                            err = sqrt(err);
                            tgrtmp->SetPointError(iptt,0,err);
                        }

                        // mgr->Add(tgr_syserr,"2");
                        mgr->Add(tgrtmp,"P");

                        mgr->SetTitle(tgrtmp->GetTitle());

                        mgr->Draw("A");

                        mgr->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                        mgr->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                        // mgr->GetHistogram()->GetXaxis()->SetLimits(GAT::limptt[0],GAT::limptt[GAT::nptt]);
                        // mgr->GetHistogram()->GetXaxis()->SetNdivisions(-GAT::limptt[0]+GAT::limptt[GAT::nptt]);
                        // int xlab[] = {2,4,6,7,9,10};
                        // for (auto ilab : xlab)
                        // {
                        //      mgr->GetHistogram()->GetXaxis()->ChangeLabel(ilab,-1,-1,-1,-1,-1," ");
                        // }
                        mgr->GetHistogram()->GetXaxis()->SetLimits(GAT::limptt[0],GAT::limptt[GAT::nptt]+1.);

                        double ymin = mgr->GetHistogram()->GetYaxis()->GetXmin();
                        double ymax = mgr->GetHistogram()->GetYaxis()->GetXmax();
                        // ymaxer = max(ymax,ymaxer);
                        // yminer = min(ymin,yminer);
                        mgr->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*.1);
                        if (icent == nsumcent && idpsi == ALL && itype == GAT::PI0 && iadj == 0)
                        {
                            // if (ixi == 0)
                            //  mgr->GetHistogram()->GetYaxis()->SetRangeUser(.0055,.027);
                            // if (ixi == 1)
                            //  mgr->GetHistogram()->GetYaxis()->SetRangeUser(.056,.083);
                            // if (ixi == 2)
                            //  mgr->GetHistogram()->GetYaxis()->SetRangeUser(.06,.41);

                            if (ixi == 0)
                                mgr->GetHistogram()->GetYaxis()->SetRangeUser(.014,.068);
                            // if (ixi == 1)
                            //  mgr->GetHistogram()->GetYaxis()->SetRangeUser(.056,.083);
                            if (ixi == 2)
                                mgr->GetHistogram()->GetYaxis()->SetRangeUser(.22,.94);
                        }
                    
                        can.Update();

                        mgr->Fit("pol1");

                        TGraphErrors* grint = new TGraphErrors(GAT::nptt);
                        for (int i=0; i<GAT::nptt; i++)
                            grint->SetPoint(i, tgr->GetX()[i], 0);
                        // (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,ci);
                        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,.68);//1 sigma
                        // grint->Draw("psame");

                        TGraph* grmin = new TGraphErrors(GAT::nptt);
                        TGraph* grmax = new TGraphErrors(GAT::nptt);
                        for (int i=0; i<GAT::nptt; i++) {
                            grmin->SetPoint(i,grint->GetX()[i],grint->GetY()[i]-grint->GetEY()[i]);
                            grmax->SetPoint(i,grint->GetX()[i],grint->GetY()[i]+grint->GetEY()[i]);
                        }

                        double limptt[GAT::nptt+1] = {5.,7.,9.,12.,16.};//easiest way to extend drawing of ci to 16
                        TGraph* grshade = new TGraph(2*(2*GAT::nptt+1));
                        for (int i=0; i<2*GAT::nptt+1; i++)
                        {
                            double x = limptt[i/2]+(i%2)*(limptt[i/2+1]-limptt[i/2])/2.;
                            // cout<<"Setting point "<<i<<" to ("<<x<<","<<grmax->Eval(x)<<")"<<endl;
                            grshade->SetPoint(i,x,grmax->Eval(x));
                            // cout<<"Setting point "<<2*GAT::nptt+1+i<<" to ("<<limptt[GAT::nptt]-(x-limptt[0])<<","<<grmin->Eval(limptt[GAT::nptt]-(x-limptt[0]))<<")"<<endl;
                            grshade->SetPoint(2*GAT::nptt+1+i,limptt[GAT::nptt]-(x-limptt[0]),grmin->Eval(limptt[GAT::nptt]-(x-limptt[0])));
                
                        }

                        grshade->SetFillStyle(1001);
                        grshade->SetFillColor(kRed-7);
                        grshade->SetFillColorAlpha(kRed-7,.2);
                        grshade->Draw("fsame");

                        //draw coordinates
                        TF1* fit = mgr->GetFunction("pol1");
                        // fit->Print();

                        if (icent == nsumcent && idpsi == ALL && itype == GAT::PI0 && iadj == 0)
                        {
                            for (int iptt = 0; iptt < GAT::nptt; iptt++)
                            {
                                double misx = pt_missed[iptt];
                                double misy = fit->Eval(misx);
                                TLatex* mistxt;
                                if (ixi == 0)
                                    mistxt = new TLatex(misx,misy,Form("#rightarrow  (%.2f,%.4f)",misx,misy));
                                if (ixi == 1)
                                    mistxt = new TLatex(misx,misy,Form("#rightarrow  (%.2f,%.4f)",misx,misy));
                                if (ixi == 2)
                                    mistxt = new TLatex(misx,misy,Form("#rightarrow  (%.2f,%.3f)",misx,misy));
                                mistxt->SetTextAngle(90);
                                mistxt->SetTextAlign(12);
                                mistxt->SetTextSize(.03);
                                mistxt->SetTextColor(TColor::GetColorPalette(220*.2/(GAT::nxi-1)));
                                mistxt->Draw();

                                double tagx = pt_tagged[iptt];
                                double tagy = fit->Eval(tagx);
                                TLatex* tagtxt;
                                if (ixi == 0)
                                    tagtxt = new TLatex(tagx,tagy,Form("(%.2f,%.4f)  #leftarrow",tagx,tagy));
                                if (ixi == 1)
                                    tagtxt = new TLatex(tagx,tagy,Form("(%.2f,%.4f)  #leftarrow",tagx,tagy));
                                if (ixi == 2)
                                    tagtxt = new TLatex(tagx,tagy,Form("(%.2f,%.3f)  #leftarrow",tagx,tagy));
                                tagtxt->SetTextAngle(90);
                                tagtxt->SetTextAlign(32);
                                tagtxt->SetTextSize(.03);
                                tagtxt->SetTextColor(TColor::GetColorPalette(220*.8/(GAT::nxi-1)));
                                tagtxt->Draw();
                            }
                        }

                        TLegend* leg = new TLegend(0.7,0.8375,0.9,0.9);
                        leg->SetTextSize(.035);
                        leg->SetTextAlign(22);
                        leg->AddEntry(tgrtmp,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]),"P");
                        leg->Draw();
                            
                        can.SaveAs(Form("integpics/yptt_adj%i_%s_cent%i_dpsi%i_xi%i.pdf",iadj,GAT::spec[itype],icent,idpsi,ixi));
                    }
                }
                gStyle->SetOptFit();
            }       


            //B) all-plane. new plot for each cent. same plot for each pt
            for (int icent = 0; icent < nsumcent+1; icent++)
            {
                TCanvas* can[GAT::ndpsi+1];
                TMultiGraph* mgr[GAT::ndpsi+1];
                double ymaxer = numeric_limits<double>::min();
                double yminer = numeric_limits<double>::max();

                // ---------- introducing my big canvas ----------
                TCanvas* bigcan = new TCanvas("bigcan","",3000,3000);
                TPad* bigpads[3][3] = {};
                TMultiGraph* lilmgr[3][3] = {};
                TGraphErrors* tgr_copy[3][3] = {};
                TGraphErrors* tgr_syserr_copy[3][3] = {};
                TGraphErrors* tgr_efferr[3][3] = {};
                TGraphErrors* tgr_smear[3][3] = {};
                TGraphErrors* tgr_total[3][3] = {};
                double ypos = 0.;
                for (int iy = 0; iy < 3; iy++)
                {

                    double top = 0.;
                    if (iy == 2) top = .3;
                    double bot = 0.;
                    if (iy == 0) bot = .3;

                    // 1 = h0 + h1 + h2
                    // h1 = .7h0 = .7h2 = (1.-top-bot)hi
                    // 1 = 1./.7 h1 + h1 + 1./.7 h1 = sum 1./(1.-top-bot) hi
                    // h1 = 1./(2.7/.7) = .7/2.7;
                    double ht = .7/2.7;
                    if (iy != 1) ht = 1./2.7;

                    double xpos = 0.;
                    for (int ix = 0; ix < 3; ix++)
                    {
                        double left = 0.;
                        if (ix == 0) left = .25;
                        double right = 0.;
                        if (ix == 2) right = .2;
                        // 1 = w0 + w1 + w2
                        // w1 = .75w0 = .8w2 = (1.-top-bot)wi
                        // 1 = 1./.75 w1 + w1 + 1./.8 w1 = sum 1./(1.-top-bot) wi
                        // w1 = 1./((15.5/6.) + 1.) = 1./21.5/6. = 6./21.5
                        double wd = 6./21.5;
                        if (ix == 0) wd/=.75;
                        if (ix == 2) wd/=.8;

                        bigpads[ix][iy] = new TPad(Form("pad%i%i",ix,iy),"",xpos,ypos,xpos+wd,ypos+ht);
                        bigpads[ix][iy]->SetTopMargin(top);
                        bigpads[ix][iy]->SetBottomMargin(bot);
                        bigpads[ix][iy]->SetRightMargin(right);
                        bigpads[ix][iy]->SetLeftMargin(left);
                        bigpads[ix][iy]->Draw();

                        if (ix != 2) bigpads[ix][iy]->SetTicky(2);
                        if (iy != 2) bigpads[ix][iy]->SetTickx(2);
                        
                        xpos+=wd;
                    }
                    ypos+=ht;
                }

                //---------- big canvas continuing----------
                gStyle->SetEndErrorSize(12);
                bool gobig = true;
                for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                {
                    if (islateadj && idpsi != GAT::ndpsi) continue;
                    
                    for (int ixi = GAT::nxi-1; ixi >= 0; ixi--)         
                    {
                        TGraphErrors* tgr = tgripttyield[iadj][itype][icent][idpsi][ixi];
                        TGraphErrors* tgr_syserr = tgripttyield_syserr[iadj][itype][icent][idpsi][ixi];

                        bigpads[idpsi][ixi]->cd();
                        // make copy of tgraphs without shift in pt
                        vector<double> vptt, vptt_left, vptt_right;
                        vector<double> vex, vexbig;
                        vector<double> vefferr;
                        vector<double> vsmr;
                        vector<double> vtot;
                        for (int iptt = 0; iptt < GAT::nptt; iptt++)
                        {
                            double ew = .25;
                            if (gobig) ew = .2;//maybe use different error width
                            
                            vptt.push_back((GAT::limptt[iptt]+GAT::limptt[iptt+1])/2.);
                            vptt_left.push_back((GAT::limptt[iptt]+GAT::limptt[iptt+1])/2. - 2.*ew);
                            vptt_right.push_back((GAT::limptt[iptt]+GAT::limptt[iptt+1])/2. + 2.*ew);

                            vex.push_back(ew);
                            vexbig.push_back(3.*ew);                        

                            double efer = efferr*tgr->GetY()[iptt];
                            vefferr.push_back(efer);

                            //take care: smearing ratio binned in xi
                            // double smr = fabs(tgrsmear_ratio[iadj][itype][iptt]->GetY()[ixi]-1.)*tgr->GetY()[iptt];
                            double smr = -9999.;
                            if (isdiradj2)
                                smr = fabs(tgrsmear_ratio[0][GAT::INC][iptt]->GetY()[ixi]-1.)*tgr->GetY()[iptt];
                            else
                                smr = fabs(tgrsmear_ratio[iadj][itype][iptt]->GetY()[ixi]-1.)*tgr->GetY()[iptt];
                            vsmr.push_back(smr);

                            double se = tgr_syserr->GetEY()[iptt];
                            
                            double tot = sqrt(se*se + smr*smr + efer*efer);
                            vtot.push_back(tot);
                        }
                        // TGraphErrors* tgr_copy;
                        // TGraphErrors* tgr_syserr_copy;
                        // TGraphErrors* tgr_efferr;
                        // TGraphErrors* tgr_smear;
                        // TGraphErrors* tgr_total;

                        if (!gobig)
                        {
                            tgr_copy[idpsi][ixi] = new TGraphErrors(GAT::nptt,vptt.data(),tgr->GetY(),tgr->GetEX(),tgr->GetEY());
                            tgr_syserr_copy[idpsi][ixi] = new TGraphErrors(GAT::nptt,vptt.data(),tgr_syserr->GetY(),vex.data(),tgr_syserr->GetEY());
                            tgr_efferr[idpsi][ixi] = new TGraphErrors(GAT::nptt,vptt.data(),tgr->GetY(),vex.data(),vefferr.data());
                            tgr_smear[idpsi][ixi] = new TGraphErrors(GAT::nptt,vptt.data(),tgr->GetY(),vex.data(),vsmr.data());
                            tgr_total[idpsi][ixi] = new TGraphErrors(GAT::nptt,vptt.data(),tgr->GetY(),vex.data(),vtot.data());
                        }
                        else
                        {
                            tgr_copy[idpsi][ixi] = new TGraphErrors(GAT::nptt,vptt.data(),tgr->GetY(),tgr->GetEX(),tgr->GetEY());
                            tgr_syserr_copy[idpsi][ixi] = new TGraphErrors(GAT::nptt,vptt.data(),tgr_syserr->GetY(),vex.data(),tgr_syserr->GetEY());
                            tgr_efferr[idpsi][ixi] = new TGraphErrors(GAT::nptt,vptt_right.data(),tgr->GetY(),vex.data(),vefferr.data());
                            tgr_smear[idpsi][ixi] = new TGraphErrors(GAT::nptt,vptt_left.data(),tgr->GetY(),vex.data(),vsmr.data());
                            tgr_total[idpsi][ixi] = new TGraphErrors(GAT::nptt,vptt.data(),tgr->GetY(),vexbig.data(),vtot.data());                          
                        }
                        
                        tgr_copy[idpsi][ixi]->SetMarkerColor(kBlack);
                        tgr_copy[idpsi][ixi]->SetLineColor(kBlack);
                        tgr_copy[idpsi][ixi]->SetMarkerSize(2);
                        tgr_copy[idpsi][ixi]->SetLineWidth(1);
                        tgr_copy[idpsi][ixi]->SetMarkerStyle(kFullCircle);

                        tgr_total[idpsi][ixi]->SetFillColorAlpha(TColor::GetColorPalette(220*3/3),.5);
                        tgr_efferr[idpsi][ixi]->SetFillColorAlpha(TColor::GetColorPalette(220*2/3),.5);
                        tgr_syserr_copy[idpsi][ixi]->SetFillColorAlpha(TColor::GetColorPalette(220*1/3),.50);
                        tgr_smear[idpsi][ixi]->SetFillColorAlpha(TColor::GetColorPalette(220*0/3),.50);
                        
                        // TMultiGraph* lilmgr = new TMultiGraph();
                        lilmgr[idpsi][ixi] = new TMultiGraph();
                        lilmgr[idpsi][ixi]->Add(tgr_total[idpsi][ixi],"2");
                        lilmgr[idpsi][ixi]->Add(tgr_efferr[idpsi][ixi],"2");
                        lilmgr[idpsi][ixi]->Add(tgr_syserr_copy[idpsi][ixi],"2");
                        lilmgr[idpsi][ixi]->Add(tgr_smear[idpsi][ixi],"2");
                        lilmgr[idpsi][ixi]->Add(tgr_copy[idpsi][ixi],"P");

                        lilmgr[idpsi][ixi]->Draw("A");

                        if (ixi == 2)
                        {
                            TText* txt;
                            double rt = bigpads[idpsi][ixi]->GetRightMargin();
                            double lt = bigpads[idpsi][ixi]->GetLeftMargin();
                            double xpos = ((1.-rt)+lt)/2.;
                            if (idpsi == IN)
                            {
                                txt = new TText(xpos,.80,"in-plane");
                                // lilmgr[idpsi][ixi]->SetTitle("in-plane");
                            }
                            if (idpsi == OUT)
                            {
                                txt = new TText(xpos,.80,"out-of-plane");
                                // lilmgr[idpsi][ixi]->SetTitle("out-of-plane");
                            }                           
                            if (idpsi == ALL)
                            {
                                txt = new TText(xpos,.80,"all-plane");
                                // lilmgr[idpsi][ixi]->SetTitle("all-plane");
                            }
                            txt->SetNDC();
                            txt->SetTextAlign(22);
                            txt->SetTextFont(43);
                            txt->SetTextSize(60);
                            txt->Draw();
                        }

                        if (idpsi == 2)
                        {
                            TLatex* txt;
                            double tp = bigpads[idpsi][ixi]->GetTopMargin();
                            double bt = bigpads[idpsi][ixi]->GetBottomMargin();
                            double ypos = ((1.-tp)+bt)/2.;

                            txt = new TLatex(.9,ypos,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]));
                            
                            txt->SetNDC();
                            txt->SetTextAlign(22);
                            txt->SetTextAngle(90);
                            // txt->SetTextAngle(270);
                            txt->SetTextFont(43);
                            txt->SetTextSize(60);
                            txt->Draw();
                        }

                        if (ixi == 0)
                        {
                            double rt = bigpads[idpsi][ixi]->GetRightMargin();
                            double lt = bigpads[idpsi][ixi]->GetLeftMargin();
                            // double xpos = ((1.-rt)+lt)/2.;
                            double xpos = .5;
                            double legw = 1.;

                            // TLegend* leg = new TLegend(xpos-legw/2.,.05,xpos+legw/2.,.15);
                            TLegend* leg = new TLegend(lt,.05,1.-rt,.15);
                            // leg->SetMargin(.1/(1.-rt-lt));
                            leg->SetMargin(.15);
                            if (idpsi == IN)
                            {
                                leg->AddEntry(tgr_smear[idpsi][ixi]," momentum resolution","f");
                            }
                            if (idpsi == OUT)
                            {
                                leg->AddEntry(tgr_syserr_copy[idpsi][ixi]," background subtraction","f");
                            }
                            if (idpsi == ALL)
                            {
                                leg->SetNColumns(2);
                                leg->SetMargin(.3);
                                leg->AddEntry(tgr_efferr[idpsi][ixi]," hadron efficiency","f");
                                leg->AddEntry(tgr_total[idpsi][ixi]," total","f");
                            }

                            // leg->SetTextAlign(22);
                            leg->SetBorderSize(0);
                            leg->SetTextFont(43);
                            leg->SetTextSize(55);

                            leg->Draw();
                        }

                        //tgr has x shifted. should copy to new tgraph without x shift
                        // lilmgr[idpsi][ixi]->GetHistogram()->GetXaxis()->SetLimits(4.8,15.2);
                        // lilmgr[idpsi][ixi]->GetHistogram()->GetXaxis()->SetLimits(4.75,14.75);
                        lilmgr[idpsi][ixi]->GetHistogram()->GetXaxis()->SetLimits(4.5,15.);
                        // lilmgr[idpsi][ixi]->GetHistogram()->GetXaxis()->SetLimits(5.,14.5);
                        if (biglims[iadj][itype][icent][ixi][0] != biglims[iadj][itype][icent][ixi][1])
                            lilmgr[idpsi][ixi]->GetHistogram()->GetYaxis()->SetRangeUser(biglims[iadj][itype][icent][ixi][0],biglims[iadj][itype][icent][ixi][1]);

                        cout<<"WARNING!! Must set big lims for other trigger types!"<<endl;
                        
                        if (ixi == 0)
                        {
                            lilmgr[idpsi][ixi]->GetHistogram()->GetXaxis()->SetLabelFont(43);
                            lilmgr[idpsi][ixi]->GetHistogram()->GetXaxis()->SetLabelSize(50);
                            
                            lilmgr[idpsi][ixi]->GetHistogram()->GetXaxis()->SetTitleFont(43);
                            lilmgr[idpsi][ixi]->GetHistogram()->GetXaxis()->SetTitleSize(60);                       
                            lilmgr[idpsi][ixi]->GetHistogram()->GetXaxis()->SetTitleOffset(3.0);
                            lilmgr[idpsi][ixi]->GetHistogram()->GetXaxis()->SetTitle("p_{T t} ");
                        }
                        lilmgr[idpsi][ixi]->GetHistogram()->GetXaxis()->SetNdivisions(506);

                        if (idpsi == 0)
                        {
                            lilmgr[idpsi][ixi]->GetHistogram()->GetYaxis()->SetLabelFont(43);
                            lilmgr[idpsi][ixi]->GetHistogram()->GetYaxis()->SetLabelSize(50);

                            lilmgr[idpsi][ixi]->GetHistogram()->GetYaxis()->SetTitleFont(43);
                            lilmgr[idpsi][ixi]->GetHistogram()->GetYaxis()->SetTitleSize(60);                       
                            lilmgr[idpsi][ixi]->GetHistogram()->GetYaxis()->SetTitleOffset(5.0);                            
                            lilmgr[idpsi][ixi]->GetHistogram()->GetYaxis()->SetTitle("1/N_{t} dN_{h} /d#xi ");
                        }
                        lilmgr[idpsi][ixi]->GetHistogram()->GetYaxis()->SetNdivisions(508);

                    }
                }

                bigcan->SaveAs(Form("integpics/bigcan_yptt_adj%i_%s_cent%i.pdf",iadj,GAT::spec[itype],icent));
                        
                gStyle->SetEndErrorSize(5);
                for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                {
                    if (islateadj && idpsi != GAT::ndpsi) continue;
                    
                    // TCanvas can ("can","",1000,1000);
                    // TCanvas can ("can","",800,900);
                    // TCanvas can ("can","",1100,1450);
                    can[idpsi] = new TCanvas(Form("can_dpsi%i",idpsi),"",800,900);
                    // can[idpsi]->SetGrid(0,1);
                    can[idpsi]->SetLogy();
                    // TMultiGraph* mgr = new TMultiGraph();
                    mgr[idpsi] = new TMultiGraph();
                    // TLegend leg (.15,.7,.45,.875);
                    // TLegend leg (.55,.15,.85,.325);
                    TLegend* leg = new TLegend(0.3,0.775,0.9,0.9);
                    // TLegend leg (0.3,0.8,0.9,0.9);
                    leg->SetTextSize(.035);
                    leg->SetTextAlign(22);
                    leg->SetNColumns(2);
                    for (int ixi = GAT::nxi-1; ixi >= 0; ixi--)         
                    {
                        TGraphErrors* tgr = tgripttyield[iadj][itype][icent][idpsi][ixi];

                        tgr->SetMarkerColor(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)));
                        tgr->SetLineColor(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)));
                
                        tgr->SetMarkerSize(2.);
                        tgr->SetLineWidth(3);
                        tgr->SetMarkerStyle(markstyle[ixi]);
                        // tgr->GetYaxis()->SetTitleSize(0.5);
                        // tgr->GetXaxis()->SetTitleSize(0.5);

                        // TGraphErrors* tgr_syserr = tgripttyield_syserr[iadj][itype][icent][idpsi][ixi];
                        TGraphErrors* tgr_syserr = tgripttyield_fullsyserr[iadj][itype][icent][idpsi][ixi];
                        tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*ixi/(GAT::nxi-1)),.50);

                        mgr[idpsi]->Add(tgr_syserr,"2");
                        mgr[idpsi]->Add(tgr,"P");

                        mgr[idpsi]->SetTitle(tgr->GetTitle());

                    }

                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        TGraphErrors* tgr = tgripttyield[iadj][itype][icent][idpsi][ixi];
                        ostringstream oss;
                        if (ixi == GAT::nxi-1)
                            leg->AddEntry((TObject*)0, "", "");
                        leg->AddEntry(tgr,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]),"P");
                    }

                    can[idpsi]->cd();
                    mgr[idpsi]->Draw("A");          

                    double ymin = mgr[idpsi]->GetHistogram()->GetYaxis()->GetXmin();
                    double ymax = mgr[idpsi]->GetHistogram()->GetYaxis()->GetXmax();
                    mgr[idpsi]->GetHistogram()->GetYaxis()->SetRangeUser(ymin,ymax+(ymax-ymin)*5.);

                    ymaxer = max(ymax,ymaxer);
                    yminer = min(ymin,yminer);
            
                    mgr[idpsi]->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                    mgr[idpsi]->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                    mgr[idpsi]->GetHistogram()->GetXaxis()->SetLimits(GAT::limptt[0],GAT::limptt[GAT::nptt]);
                    mgr[idpsi]->GetHistogram()->GetXaxis()->SetNdivisions(-GAT::limptt[0]+GAT::limptt[GAT::nptt]);
                    int xlab[] = {2,4,6,7,9,10};
                    for (auto ilab : xlab)
                    {
                        mgr[idpsi]->GetHistogram()->GetXaxis()->ChangeLabel(ilab,-1,-1,-1,-1,-1," ");
                    }
                
                    leg->Draw();

                    double xmin = mgr[idpsi]->GetHistogram()->GetXaxis()->GetXmin();
                    double xmax = mgr[idpsi]->GetHistogram()->GetXaxis()->GetXmax();
                    TLine* line = new TLine(xmin,0.0,xmax,0.0);
                    line->SetLineStyle(7);
                    line->Draw();

                    TLatex* txtcent;
                    if (icent != nsumcent)
                    {
                        if (idpsi == IN)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV In-Plane",limsumcent[icent],limsumcent[icent+1]));
                        else if (idpsi == OUT)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV Out-of-Plane",limsumcent[icent],limsumcent[icent+1]));
                        else
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV",limsumcent[icent],limsumcent[icent+1]));
                    }
                    else
                    {
                        if (idpsi == IN)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV In-Plane",limsumcent[0],limsumcent[nsumcent]));
                        else if (idpsi == OUT)
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV Out-of-Plane",limsumcent[0],limsumcent[nsumcent]));
                        else
                            txtcent = new TLatex(0.5,.72,Form("Au+Au %i-%i%% 200GeV",limsumcent[0],limsumcent[nsumcent]));
                    }
                    txtcent->SetTextAlign(22);
                    txtcent->SetNDC();
                    txtcent->SetTextSize(0.035);
                    txtcent->Draw();
                }

                for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                {
                    if (islateadj && idpsi != GAT::ndpsi) continue;
                    // if (isdiradj2 && idpsi != GAT::ndpsi) continue;
                    mgr[idpsi]->GetHistogram()->GetYaxis()->SetRangeUser(yminer,ymaxer+(ymaxer-yminer)*5.);
                    can[idpsi]->Modified();
                    can[idpsi]->SaveAs(Form("integpics/yptt_adj%i_%s_dpsi%i_cent%i.pdf",iadj,GAT::spec[itype],idpsi,icent));
                }
            }
            //finish copied section to match limits across dpsi not cent
            //---------- end copied ixi to iptt plotting ----------
            
            //A) all-plane. new plot for each pt. same plot for each cent.

            for (int iptt = 0; iptt < GAT::nptt; iptt++)                    
            {
                // if (isdiradj2) break; // dont have these plots
                
                TCanvas can ("can","",1000,1000);
                can.SetLogy();
                TMultiGraph* mgr = new TMultiGraph();
                TLegend leg (.15,.78,.35,.875);

                for (int icent = 0; icent < nsumcent; icent++)
                {
                    TGraphErrors* tgr = tgrixiyield[iadj][itype][iptt][icent][ALL];

                    tgr->SetMarkerColor(TColor::GetColorPalette(220*icent/(nsumcent-1)));
                    tgr->SetLineColor(TColor::GetColorPalette(220*icent/(nsumcent-1)));
                
                    tgr->SetMarkerSize(2.);
                    tgr->SetLineWidth(3);
                    tgr->SetMarkerStyle(markstyle[icent]);

                    // TGraphErrors* tgr_syserr = tgrixiyield_syserr[iadj][itype][iptt][icent][ALL];
                    TGraphErrors* tgr_syserr = tgrixiyield_fullsyserr[iadj][itype][iptt][icent][ALL];
                    tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*icent/(nsumcent-1)),.50);
                    
                    mgr->Add(tgr_syserr,"2");
                    mgr->Add(tgr,"P");

                    ostringstream oss;
                    if (icent == nsumcent)
                        leg.AddEntry(tgr,Form("%i-%i%%",limsumcent[0],limsumcent[nsumcent]),"P");
                    else
                        leg.AddEntry(tgr,Form("%i-%i%%",limsumcent[icent],limsumcent[icent+1]),"P");

                    mgr->SetTitle(tgr->GetTitle());
                }

            
                mgr->Draw("A");
                // if (icent == 0)
                //     mgr->GetHistogram()->GetYaxis()->SetLimits(-1.4,2.2);
                // else if (icent == 1)
                //     mgr->GetHistogram()->GetYaxis()->SetLimits(-1.0,1.0);
                mgr->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                mgr->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                mgr->GetHistogram()->GetXaxis()->SetLimits(0.,2.4);
                mgr->GetHistogram()->GetXaxis()->SetNdivisions(-3);
                
                leg.Draw();

                double xmin = mgr->GetHistogram()->GetXaxis()->GetXmin();
                double xmax = mgr->GetHistogram()->GetXaxis()->GetXmax();
                TLine* line = new TLine(xmin,0.0,xmax,0.0);
                line->SetLineStyle(7);
                line->Draw();

                can.Modified();
                can.SaveAs(Form("integpics/y_adj%i_%s_ptt%i.pdf",iadj,GAT::spec[itype],iptt));
            }

            //C) in-out-all plane cases same plot. new plot for each pt, cent
            for (int iptt = 0; iptt < GAT::nptt; iptt++)                    
            {
                // if (isdiradj2) break; // dont have these plots
                
                for (int icent = 0; icent < nsumcent+1; icent++)
                {
                    TCanvas can ("can","",1000,1000);
                    can.SetLogy();
                    TMultiGraph* mgr = new TMultiGraph();
                    TLegend leg (.15,.7875,.40,.875);

                    for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                    {
                        if (islateadj && idpsi != GAT::ndpsi) continue;
                        // if (isdiradj2 && idpsi != GAT::ndpsi) continue;
                        
                        if (idpsi == GAT::ndpsi) continue;//skip all plane case
                        TGraphErrors* tgr = tgrixiyield[iadj][itype][iptt][icent][idpsi];
                

                        tgr->SetMarkerColor(TColor::GetColorPalette(220*idpsi/(GAT::ndpsi-1)));
                        tgr->SetLineColor(TColor::GetColorPalette(220*idpsi/(GAT::ndpsi-1)));
                
                        tgr->SetMarkerSize(2.);
                        tgr->SetLineWidth(3);
                        tgr->SetMarkerStyle(markstyle[idpsi]);

                        // TGraphErrors* tgr_syserr = tgrixiyield_syserr[iadj][itype][iptt][icent][idpsi];
                        TGraphErrors* tgr_syserr = tgrixiyield_fullsyserr[iadj][itype][iptt][icent][idpsi];
                        tgr_syserr->SetFillColorAlpha(TColor::GetColorPalette(220*idpsi/(GAT::ndpsi-1)),.50);

                        mgr->Add(tgr_syserr,"2");//add sys error box first
                        mgr->Add(tgr,"P");

                        ostringstream oss;
                        leg.AddEntry(tgr,Form("%s",strplane[idpsi].c_str()),"P");

                        mgr->SetTitle(tgr->GetTitle());
                    }

            
                    mgr->Draw("A");
                    // if (icent == 0)
                    //     mgr->GetHistogram()->GetYaxis()->SetLimits(-1.4,2.2);
                    // else if (icent == 1)
                    //     mgr->GetHistogram()->GetYaxis()->SetLimits(-1.0,1.0);
                    mgr->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                    mgr->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
                    mgr->GetHistogram()->GetXaxis()->SetLimits(0.,2.4);
                    mgr->GetHistogram()->GetXaxis()->SetNdivisions(-3);
                
                    leg.Draw();

                    double xmin = mgr->GetHistogram()->GetXaxis()->GetXmin();
                    double xmax = mgr->GetHistogram()->GetXaxis()->GetXmax();
                    TLine* line = new TLine(xmin,0.0,xmax,0.0);
                    line->SetLineStyle(7);
                    line->Draw();

                    can.Modified();
                    can.SaveAs(Form("integpics/y_adj%i_%s_ptt%i_cent%i.pdf",iadj,GAT::spec[itype],iptt,icent));

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



}

bool FindObject(TFile* input, TString name, TGraphErrors*& obj)
{
    // cout<<"Looking for "<<name<<" in "<<input->GetName()<<endl;
    
    input->GetObject(name,obj);

    if (!obj)
    {
        cerr<<"Could not find "<<name<<". Exiting..."<<endl;
        return false;
    }

    // cout<<"Successfully found "<<name<<endl<<endl;
    return true;
}

bool FindObject(TFile* input, TString name, TGraph*& obj)
{
    // cout<<"Looking for "<<name<<" in "<<input->GetName()<<endl;
    
    input->GetObject(name,obj);

    if (!obj)
    {
        cerr<<"Could not find "<<name<<". Exiting..."<<endl;
        return false;
    }

    // cout<<"Successfully found "<<name<<endl<<endl;
    return true;
}

//this method has not been configured for islateadj
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
                        tgripttyield_fullsyserr[iadj][itype][icent][idpsi][ixi]->SetTitle(";p_{Tt};1/N_{t} dN/d#xi");
                        tgripttyield_fullsyserr[iadj][itype][icent][idpsi][ixi]->SetName(Form("ipttyield_fullsyserr_adj%i_%s_cent%i_dpsi%i_xi%i",iadj,GAT::spec[itype],icent,idpsi,ixi));
                    }

                    for (int iptt = 0; iptt < GAT::nptt; iptt++)
                    {
                        tgrixiyield_fullsyserr[iadj][itype][iptt][icent][idpsi] = new TGraphErrors(GAT::nxi);
                        tgrixiyield_fullsyserr[iadj][itype][iptt][icent][idpsi]->SetTitle(";#xi;1/N_{t} dN/d#xi");
                        tgrixiyield_fullsyserr[iadj][itype][iptt][icent][idpsi]->SetName(Form("ixiyield_fullsyserr_adj%i_%s_ptt%i_cent%i_dpsi%i",iadj,GAT::spec[itype],iptt,icent,idpsi));
                    }
                    
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        for (int iptt = 0; iptt < GAT::nptt; iptt++)
                        {
                            {
                                TGraphErrors* tgr = tgripttyield_syserr[iadj][itype][icent][idpsi][ixi];

                                double efer = efferr*tgr->GetY()[iptt];

                                //for adj2_dir smearing we will be using adj0_inc
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

                                //for adj2_dir smearing we will be using adj0_inc
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
