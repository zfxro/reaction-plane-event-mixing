//make some spiffy plots with kcomponent output
//much is copied directly from calcKComponents.C....

//drawKComps is being overwritten by something simpler for thesis
//see backuped version 1 or commented regions for previous method

const int nadj = 2;

// const int npartmin = 2;
// const int npartmax = 390;
int npartmin = -999;
int npartmax = -999;

const float PI = (float) TMath::Pi();

namespace GAT {
    #include "../src/gamAnaTree/gamAnaTree.h"
}

TLatex* txtmom[GAT::nptt];
TLatex* txtxi[GAT::nxi];
TLatex* txtmomlow[GAT::nptt];
TLatex* txtxilow[GAT::nxi];
// const string kname[4] {"#kappa_{a}","#kappa_{b}","#kappa_{#gammah}","#kappa_{h#gamma}"};
const string kname[6] {"#kappa_{a}","#kappa_{b}","#kappa_{#gammah}","#kappa_{h#gamma}","#kappa_{b#gamma#sigma}","#kappa_{bh#sigma}"};

int GetNpartDists(TH1F* hnpart[GAT::ncent], vector<double>& vnpart, vector<double>& vnparte);

TGraphErrors* MakeProd(TGraphErrors* tgra, TGraphErrors* tgrae, TGraphErrors* tgrb, TGraphErrors* tgrbe, string name);

// void DrawKComps(TGraphErrors* gra, TGraphErrors* grae, string laba, TGraphErrors* grb, TGraphErrors* grbe, string labb, TGraphErrors* grp, string labp, TH1F* hnpart[GAT::ncent], string kname = "", TGraphErrors* grk = nullptr, TGraphErrors* gra0 = nullptr, TGraphErrors* grae0 = nullptr, string laba0 = "", TGraphErrors* grb0 = nullptr, TGraphErrors* grbe0 = nullptr, string labb0 = "");

void DrawKComps(TGraphErrors* gra, TGraphErrors* grae, string laba, TGraphErrors* grb, TGraphErrors* grbe, string labb, TGraphErrors* grp, string labp, TH1F* hnpart[GAT::ncent], string kname = "", TGraphErrors* grk = nullptr);

// void DrawKCompsSmall(TGraphErrors* gra, TGraphErrors* grae, string laba, TGraphErrors* grb, TGraphErrors* grbe, string labb, TGraphErrors* grp, string labp, TH1F* hnpart[GAT::ncent], TGraphErrors* grk, int ik, int iptt, int ixi);

void DrawKCompsSmall(TGraphErrors* gra, TGraphErrors* grae, string laba, TGraphErrors* grb, TGraphErrors* grbe, string labb, TGraphErrors* grp, string labp, TH1F* hnpart[GAT::ncent], TGraphErrors* grk, int ik, int iptt, int ixi,TGraphErrors* grad = nullptr, TGraphErrors* grade = nullptr, TGraphErrors* grbd = nullptr, TGraphErrors* grbde = nullptr, TGraphErrors* sig = nullptr);

void CalcExp(TF1* fa, TGraphErrors* grae, TH1F* hnpart, double& expa, double& expae);

void drawKComponents(
    const char* inname = "kComponents.root"
    )
{

    //Get Npart distributions
    TH1F* hnpart[GAT::ncent];
    vector<double> vnpart;
    vector<double> vnparte;
    if (GetNpartDists(hnpart,vnpart,vnparte)<0) return;

    npartmin = hnpart[GAT::ncent-1]->GetXaxis()->GetBinLowEdge(hnpart[GAT::ncent-1]->FindFirstBinAbove(0.0));
    npartmax = hnpart[GAT::ncent-1]->GetXaxis()->GetBinUpEdge(hnpart[0]->FindLastBinAbove(0.0));

    //Fill vector of all npart values. Used for TGraphing products as function of Npart
    // vector<double> vnpartall;
    // for (int npart = npartmin; npart <= npartmax; npart++)
    // {
    //  vnpartall.push_back(npart);
    // }

    //tlatexs
    for (int iptt = 0; iptt < GAT::nptt; iptt++)
    {   
        txtmom[iptt] = new TLatex (0.5,.85,Form("%.1f #leq p_{T} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
        txtmom[iptt]->SetTextAlign(22);
        txtmom[iptt]->SetNDC();
        txtmom[iptt]->SetTextSize(0.05);

        txtmomlow[iptt] = new TLatex (0.5,.3,Form("%.1f #leq p_{T} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
        txtmomlow[iptt]->SetTextAlign(22);
        txtmomlow[iptt]->SetNDC();
        txtmomlow[iptt]->SetTextSize(0.05);
    }

    for (int ixi = 0; ixi < GAT::nxi; ixi++)
    {
        txtxi[ixi] = new TLatex (0.5,.775,Form("%.1f #leq #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]));
        txtxi[ixi]->SetTextAlign(22);
        txtxi[ixi]->SetNDC();
        txtxi[ixi]->SetTextSize(0.05);

        txtxilow[ixi] = new TLatex (0.5,.225,Form("%.1f #leq #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]));
        txtxilow[ixi]->SetTextAlign(22);
        txtxilow[ixi]->SetNDC();
        txtxilow[ixi]->SetTextSize(0.05);
    }


    //open input
    TFile* input = new TFile(inname,"READ");
    if (!input || input->IsZombie())
    {
        cerr<<"Could not open "<<inname<<endl;
        return;
    }

    //get sigma (rp resolution)
    TGraphErrors* sig = (TGraphErrors*) gROOT->FindObject("npartsigma");

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            // if (itype != GAT::PI0) continue;
            
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                ostringstream sufg;
                sufg<<"adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt;
                
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    //manually inserting exponent for this graph. skipping all others. search for this warning message to see other continues

                    ostringstream sufh;
                    sufh<<sufg.str()<<"_xi"<<ixi;
                    
                    //---------- k_a ----------

                    TGraphErrors* ag = (TGraphErrors*) gROOT->FindObject(Form("ag_Np_%s",sufg.str().c_str()));
                    TGraphErrors* age = (TGraphErrors*) gROOT->FindObject(Form("age_Np_%s",sufg.str().c_str()));
                    TGraphErrors* ah = (TGraphErrors*) gROOT->FindObject(Form("ah_Np_%s",sufh.str().c_str()));
                    TGraphErrors* ahe = (TGraphErrors*) gROOT->FindObject(Form("ahe_Np_%s",sufh.str().c_str()));
                    
                    ostringstream ossagah;
                    ossagah<<"agah_"<<sufh.str();
                    TGraphErrors* agah = MakeProd(ag,age,ah,ahe,ossagah.str());

                    TGraphErrors* ka = (TGraphErrors*) gROOT->FindObject(Form("ka_Np_%s",sufh.str().c_str()));
                    
                    // DrawKComps(ag,age,"a_{#gamma}",ah,ahe,"a_{h}",agah,"a_{#gamma} a_{h}",hnpart,"k_{a}",ka);
                    DrawKCompsSmall(ag,age,"a_{#gamma}",ah,ahe,"a_{h}",agah,"a_{#gamma} a_{h}",hnpart,ka,0,iptt,ixi);

                    //---------- k_b ----------            
                    
                    TGraphErrors* bg = (TGraphErrors*) gROOT->FindObject(Form("bg_Np_%s",sufh.str().c_str()));
                    TGraphErrors* bge = (TGraphErrors*) gROOT->FindObject(Form("bge_Np_%s",sufh.str().c_str()));
                    TGraphErrors* bh = (TGraphErrors*) gROOT->FindObject(Form("bh_Np_%s",sufh.str().c_str()));
                    TGraphErrors* bhe = (TGraphErrors*) gROOT->FindObject(Form("bhe_Np_%s",sufh.str().c_str()));

                    TGraphErrors* bgs = (TGraphErrors*) gROOT->FindObject(Form("bgs_Np_%s",sufh.str().c_str()));
                    TGraphErrors* bgse = (TGraphErrors*) gROOT->FindObject(Form("bgse_Np_%s",sufh.str().c_str()));
                    TGraphErrors* bhs = (TGraphErrors*) gROOT->FindObject(Form("bhs_Np_%s",sufh.str().c_str()));
                    TGraphErrors* bhse = (TGraphErrors*) gROOT->FindObject(Form("bhse_Np_%s",sufh.str().c_str()));

                    ostringstream ossbgbh;
                    ossbgbh<<"bgbh_"<<sufh.str();
                    TGraphErrors* bgbh = MakeProd(bg,bge,bh,bhe,ossbgbh.str());

                    TGraphErrors* kb = (TGraphErrors*) gROOT->FindObject(Form("kb_Np_%s",sufh.str().c_str()));

                    TGraphErrors* kb2 = (TGraphErrors*) gROOT->FindObject(Form("kb2_Np_%s",sufh.str().c_str()));

                    //this drawing wasnt being used anyway...
                    // DrawKComps(bg,bge,"b_{#gamma}",bh,bhe,"b_{h}",bgbh,"b_{#gamma} b_{h}",hnpart,"k_{b}",kb,bgs,bgse,"b_{#gamma}#sigma",bhs,bhse,"b_{h}#sigma");

                    // DrawKCompsSmall(bg,bge,"b_{#gamma}",bh,bhe,"b_{h}",bgbh,"b_{#gamma} b_{h}",hnpart,kb2,1,iptt,ixi);

                    // void DrawKCompsSmall(TGraphErrors* gra, TGraphErrors* grae, string laba, TGraphErrors* grb, TGraphErrors* grbe, string labb, TGraphErrors* grp, string labp, TH1F* hnpart[GAT::ncent], TGraphErrors* grk, int ik, int iptt, int ixi,TGraphErrors* grad = nullptr, TGraphErrors* grade = nullptr, TGraphErrors* grbd = nullptr, TGraphErrors* grbe = nullptr);
                    
                    DrawKCompsSmall(bg,bge,"b_{#gamma}",bh,bhe,"b_{h}",bgbh,"b_{#gamma} b_{h}",hnpart,kb2,1,iptt,ixi,bgs,bgse,bhs,bhse,sig);

                    //---------- k_hg ----------

                    ostringstream ossahbgs;
                    ossahbgs<<"ahbgs_"<<sufh.str();
                    TGraphErrors* ahbgs = MakeProd(ah,ahe,bgs,bgse,ossahbgs.str());

                    TGraphErrors* khg = (TGraphErrors*) gROOT->FindObject(Form("khg_Np_%s",sufh.str().c_str()));
                    
                    if (iadj == 0 && itype == GAT::PI0 && iptt == 0 && ixi == 1)
                        DrawKComps(ah,ahe,"a_{h}",bgs,bgse,"b_{#gamma} #sigma",ahbgs,"a_{h} b_{#gamma} #sigma",hnpart,"k_{h #gamma}",khg);
                    
                    DrawKCompsSmall(ah,ahe,"a_{h}",bgs,bgse,"b_{#gamma} #sigma",ahbgs,"a_{h} b_{#gamma} #sigma",hnpart,khg,3,iptt,ixi);

                    //---------- k_gh ----------

                    ostringstream ossagbhs;
                    ossagbhs<<"agbhs_"<<sufh.str();
                    TGraphErrors* agbhs = MakeProd(ag,age,bhs,bhse,ossagbhs.str());

                    TGraphErrors* kgh = (TGraphErrors*) gROOT->FindObject(Form("kgh_Np_%s",sufh.str().c_str()));

                    DrawKCompsSmall(ag,age,"a_{h}",bhs,bhse,"b_{h} #sigma",agbhs,"a_{#gamma} b_{h} #sigma",hnpart,kgh,2,iptt,ixi);
                }
            }
        }
    }
    cout<<"All Done!"<<endl;
}

TGraphErrors* MakeProd(TGraphErrors* tgra, TGraphErrors* tgrae, TGraphErrors* tgrb, TGraphErrors* tgrbe, string name)
{
    vector<double> vprod;
    vector<double> vprode;
    vector<double> vnpart;
    for (int npart = npartmin; npart <= npartmax; npart++)
    {
        double a, ae, b, be;
        
        a = tgra->Eval(npart);
        if (a < 0.)
        {
            ostringstream oss;
            oss<<tgra->GetName();
            string str = oss.str();
            if (str[0] == 'a')
            {
                cout<<tgra->GetName()<<"->Eval("<<npart<<") evaluated to "<<a<<". Setting to 0."<<endl;
                a = 0;
            }
        }

        if (tgrae)
        {
            ae = tgrae->Eval(npart);
            if (ae < 0.) ae = 0.;
        }
        else ae = 0.;
                
        b = tgrb->Eval(npart);
        if (b < 0.)
        {
            ostringstream oss;
            oss<<tgrb->GetName();
            string str = oss.str();
            if (str[0] == 'a')
            {
                cout<<tgrb->GetName()<<"->Eval("<<npart<<") evaluated to "<<b<<". Setting to 0."<<endl;
                b = 0;
            }
        }

        if (tgrbe)
        {
            be = tgrbe->Eval(npart);
            if (be < 0.) be = 0.;
        }
        else be = 0.;
        
        double prod = a*b;
        double prode = sqrt(pow(a*be,2.) + pow(ae*b,2.));

        vprod.push_back(prod);
        vprode.push_back(prode);
        vnpart.push_back(npart);
    }

    TGraphErrors* tgr = new TGraphErrors(vprod.size(),&vnpart[0],&vprod[0],0,&vprode[0]);
    tgr->SetName(name.c_str());

    return tgr;
}

int GetNpartDists(TH1F* hnpart[GAT::ncent], vector<double>& vnpart, vector<double>& vnparte)
{
    const char* centdistsname = "CentDists_9bins.root";
    TFile* fcentdists = new TFile(centdistsname,"READ");
    if (!fcentdists || fcentdists->IsZombie())
    {
        cerr<<"Could not open "<<centdistsname<<endl;
        return -1;
    }
    
    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        hnpart[icent] = (TH1F*) gROOT->FindObject(Form("AuAu_Npart_%i_%i",GAT::limcent[icent],GAT::limcent[icent+1]));
        if (!hnpart[icent])
        {
            return -1;
        }
        hnpart[icent]->SetDirectory(0);
        vnpart.push_back(hnpart[icent]->GetMean());
        vnparte.push_back(hnpart[icent]->GetRMS());
    }
    fcentdists->Close();
    return 0;
}

void DrawKComps(TGraphErrors* gra, TGraphErrors* grae, string laba, TGraphErrors* grb, TGraphErrors* grbe, string labb, TGraphErrors* grp, string labp, TH1F* hnpart[GAT::ncent], string kname = "", TGraphErrors* grk = nullptr)
{

    double xmin = npartmin;
    double xmax = npartmax;

    // gROOT->ForceStyle(0);
    gStyle->SetEndErrorSize(5);
    gStyle->SetPalette(kBird);
    // gStyle->SetLabelSize(.05);
    gStyle->SetLegendBorderSize(0);

    // TGaxis::SetExponentOffset(-.06,-.05,"y");
    TGaxis::SetMaxDigits(3);
    // TGaxis::SetExponentOffset(-.06,-.04,"y");
    // TGaxis::SetExponentOffset(-.04,0,"y");

    TF1* fa = new TF1(Form("%s_eval",gra->GetName()),[=](double* x, double* par){return gra->Eval(x[0]);},xmin,xmax,0);
    fa->SetLineColor(TColor::GetColorPalette(254*3/5));
    fa->SetLineWidth(2);

    TF1* fb = new TF1("fb",[=](double* x, double* par){return grb->Eval(x[0]);},xmin,xmax,0);
    // fb->SetLineColor(kRed);
    fb->SetLineColor(TColor::GetColorPalette(254*1/5));
    fb->SetLineWidth(2);

    vector<double> vmean, vrms, vexpa, vexpae, vexpb, vexpbe, vprod, vprode;    

    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        double mean = hnpart[icent]->GetMean();
        vmean.push_back(mean);
        
        double rms = hnpart[icent]->GetRMS();
        vrms.push_back(rms);

        //Recalculate <a> and plot as single point tgrapherror. copying from calcKComponents
        double expa = -9999.;
        double expae = -9999.;

        CalcExp(fa,grae,hnpart[icent],expa,expae);

        vexpa.push_back(expa);
        vexpae.push_back(expae);

        double expb = -9999.;
        double expbe = -9999.;       

        CalcExp(fb,grbe,hnpart[icent],expb,expbe);

        vexpb.push_back(expb);
        vexpbe.push_back(expbe);
                
        double prod = expa*expb;
        vprod.push_back(prod);
        
        double prode = sqrt(pow(expa*expbe,2.) + pow(expae*expb,2.));
        vprode.push_back(prode);
    }
    
    // TCanvas* can = new TCanvas("can","",1200,1350);
    TCanvas* can = new TCanvas("can","",1700,1900);
    TPad *pad1 = new TPad("pad1","",0,0.67,1,1);
    pad1->SetBottomMargin(0.0);
    pad1->SetTopMargin(0.005);
    pad1->Draw();
    pad1->cd();

    double ymin = gra->Eval(xmin);
    double ymax = gra->Eval(xmin);
    for (int ix = xmin; ix < xmax; ++ix)
    {
        double eval = gra->Eval(ix);
        if (eval < ymin) ymin = eval;
        if (eval > ymax) ymax = eval;
    }

    ymax = ymax + fabs(ymax)*.1;
    ymin = ymin - fabs(ymin)*.1;
    ymin = ymin - (ymax - ymin)*.1;

    TH1F *frame1 = can->DrawFrame(xmin,ymin,xmax,ymax);
    frame1->SetTitle("");
    frame1->GetYaxis()->SetLabelSize(0.04);
    frame1->GetYaxis()->SetAxisColor(TColor::GetColorPalette(254*3/5));
    frame1->GetYaxis()->SetTitle(Form("%s      ",laba.c_str()));
    // frame1->GetYaxis()->SetTitle(laba.c_str());
    frame1->GetYaxis()->SetLabelSize(0.04);
    frame1->GetYaxis()->SetTitleOffset(0.65);
    frame1->GetYaxis()->SetTitleSize(0.0475);
    int font = frame1->GetYaxis()->GetLabelFont();
    // frame1->SetTitle(grp->GetName());

    fa->SetTitle("");
    // fa->GetXaxis()->SetRangeUser(xmin,xmax);
    fa->GetXaxis()->SetLimits(xmin,xmax);
    fa->Draw("same");


    //create a transparent pad drawn on top of the main pad
    can->cd();
    TPad *overlay1 = new TPad("overlay1","",0,0.67,1,1);
    overlay1->SetBottomMargin(0.0);
    overlay1->SetTopMargin(0.005);
    overlay1->SetFillStyle(0);
    overlay1->SetFillColor(0);
    overlay1->SetFrameFillStyle(0);
    // overlay1->Draw("FA");
    overlay1->Draw();
    overlay1->cd();

    
    ymin = grb->Eval(xmin);
    ymax = grb->Eval(xmin);
    for (int ix = xmin; ix < xmax; ++ix)
    {
        double eval = grb->Eval(ix);
        if (eval < ymin) ymin = eval;
        if (eval > ymax) ymax = eval;
    }

    ymax = ymax + fabs(ymax)*.1;
    ymin = ymin - fabs(ymin)*.1;
    ymin = ymin - (ymax - ymin)*.1;

    TH1F *frame2 = overlay1->DrawFrame(xmin,ymin,xmax,ymax);
    frame2->GetXaxis()->SetLabelOffset(99);
    frame2->GetYaxis()->SetTickLength(0);
    frame2->GetYaxis()->SetLabelOffset(99);
    // grb->Draw("LP");

    fb->SetTitle("");
    // fb->GetXaxis()->SetRangeUser(xmin,xmax);
    fb->GetXaxis()->SetLimits(xmin,xmax);
    fb->Draw("same");

    overlay1->Update();

    //Draw an axis on the right side
    // TGaxis::SetExponentOffset(.01,-.01,"y");
    TGaxis *axis1 = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
    // TAxis *axis1 = new TAxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
    axis1->SetLineColor(TColor::GetColorPalette(254*1/5));
    // axis1->SetAxisColor(TColor::GetColorPalette(254*1/5));
    axis1->SetTitle(Form("%s      ",labb.c_str()));
    // axis1->SetTitle(labb.c_str());
    // axis1->fAxis->RotateTitle();//fAxis is protected. TGaxis has no method rotate
    axis1->SetTitleFont(font);
    axis1->SetTitleSize(.0475);
    axis1->SetTitleOffset(.65);
    axis1->SetLabelSize(.04);
    axis1->SetLabelFont(font);
    // axis1->SetExponentOffset(.06,.04);
    // axis1->SetNoExponent();
    axis1->Draw();
    // TGaxis::SetExponentOffset(-.01,-.01,"y");

    cerr<<"WARNING! Must skip if not ahbgs_adj0_pi0_ptt0_xi1.pdf"<<endl;                   
    // ah : no exponent
    // bgs: 10^-6
    // ahbgs: 10^-6
    TLatex latex;
    latex.SetTextAlign(03);
    latex.SetTextFont(font);
    latex.SetTextSize(.04);
    latex.DrawLatexNDC(.925,.99,"#times10^{-6}");
    
    TLegend* leg1 = new TLegend(0.75,0.11,0.81,0.285);
    leg1->SetMargin(.35);
    leg1->SetTextSize(.0475);
    leg1->AddEntry(fa,laba.c_str(),"l");
    leg1->AddEntry(fb,labb.c_str(),"l");
    leg1->Draw();

    //---------- center pad ----------

    can->cd();

    TPad *pad2 = new TPad("pad2","",0,.34,1,0.67);
    pad2->SetTickx(2);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.);
    // pad2->SetFillColor(42);
    // pad2->SetGrid();
    pad2->Draw();
    pad2->cd();

    gStyle->SetOptStat(0);
    hnpart[GAT::ncent-1]->SetTitle("");
    hnpart[GAT::ncent-1]->SetLineColor(TColor::GetColorPalette(254*4/5));
    hnpart[GAT::ncent-1]->Draw();
    hnpart[GAT::ncent-1]->GetYaxis()->SetTickLength(0);
    hnpart[GAT::ncent-1]->GetYaxis()->SetLabelOffset(99);
    ymax = hnpart[GAT::ncent-1]->GetBinContent(hnpart[GAT::ncent-1]->GetMaximumBin()) * 2;
    hnpart[GAT::ncent-1]->SetMaximum(ymax);
    hnpart[GAT::ncent-1]->SetMinimum(0.);
    // hnpart[GAT::ncent-1]->GetXaxis()->SetRangeUser(xmin,xmax);
    hnpart[GAT::ncent-1]->GetXaxis()->SetLimits(xmin,xmax);
    pad2->Update();

    for (int icent = 0; icent < GAT::ncent-1; icent++)
    {
        hnpart[icent]->SetLineColor(TColor::GetColorPalette(254*4/5));
        hnpart[icent]->Draw("same");
    }
    pad2->Update();
    
    
    //create a transparent pad drawn on top of the main pad
    can->cd();
    TPad *overlay2 = new TPad("overlay2","",0,.34,1,0.67);
    overlay2->SetTopMargin(0.0);
    overlay2->SetBottomMargin(0.);
    overlay2->SetFillStyle(0);
    overlay2->SetFillColor(0);
    overlay2->SetFrameFillStyle(0);
    overlay2->Draw("FA");
    overlay2->cd();

    ymin = grp->Eval(xmin);
    ymax = grp->Eval(xmin);
    for (int ix = xmin; ix < xmax; ++ix)
    {
        double eval = grp->Eval(ix);
        if (eval < ymin) ymin = eval;
        if (eval > ymax) ymax = eval;
    }
    
    ymin = ymin - fabs(ymin)*.1;
    ymax = ymax + fabs(ymax)*.1;
    ymin = ymin - (ymax - ymin)*.1;

    TH1F *frame3 = overlay2->DrawFrame(xmin,ymin,xmax,ymax);
    frame3->GetXaxis()->SetLabelOffset(99);
    frame3->GetYaxis()->SetTickLength(0);
    frame3->GetYaxis()->SetLabelOffset(99);
    grp->SetMarkerSize(0);
    grp->SetLineColor(TColor::GetColorPalette(254*2/5));
    grp->SetLineWidth(2);
    // grp->GetXaxis()->SetRangeUser(xmin,xmax);
    grp->GetXaxis()->SetLimits(xmin,xmax);
    grp->Draw("CX");

    TGraphErrors* grprodall = new TGraphErrors(GAT::ncent,&vmean[0],&vprod[0],&vrms[0],&vprode[0]);
    grprodall->SetMarkerStyle(kFullCircle);
    // grprodall->SetMarkerStyle(21);
    grprodall->SetMarkerSize(2);
    // grprodall->SetMarkerColor(kMagenta);
    // grprodall->SetLineColor(kMagenta);
    grprodall->SetMarkerColor(TColor::GetColorPalette(254*0/5));
    grprodall->SetLineColor(TColor::GetColorPalette(254*0/5));
    grprodall->SetLineWidth(1);
    // grprodall->GetXaxis()->SetRangeUser(xmin,xmax);
    grprodall->GetXaxis()->SetLimits(xmin,xmax);
    grprodall->Draw("P same");

    TGaxis *axis2 = new TGaxis(xmin,ymin,xmin,ymax, ymin,ymax,510);
    axis2->SetTitle(Form("%s      ",labp.c_str()));
    // axis2->SetTitle(labp.c_str());
    axis2->SetTitleFont(font);
    axis2->SetTitleSize(.0475);
    axis2->SetTitleOffset(.65);
    axis2->SetLabelSize(.04);
    axis2->SetLabelFont(font);
    axis2->Draw();

    latex.DrawLatexNDC(.055,.99,"#times10^{-6}");

    TLegend* leg2 = new TLegend(0.20,0.65,0.3,0.85);
    leg2->SetMargin(.225);
    leg2->SetTextSize(0.0475);
    // leg2->SetFillStyle(0);
    // leg2->SetBorderSize(0);
    leg2->AddEntry(grp,labp.c_str(),"l");

    leg2->AddEntry(grprodall,Form("#LT%s#GT#LT%s#GT",laba.c_str(),labb.c_str()),"p");

    leg2->Draw();

    //-------------------- bottom pad --------------------

    can->cd();
    TPad* pad3 = new TPad("pad3","",0,0,1,.34);
    pad3->SetTopMargin(0.);
    pad3->SetTickx(2);
    pad3->Draw();
    pad3->cd();


    // grk->SetTitle(Form(";N_{part};#frac{<%s%s>}{<%s><%s>}",laba.c_str(),labb.c_str(),laba.c_str(),labb.c_str()));
    // grk->SetTitle(Form(";N_{part};#frac{#LT%s%s#GT}{#LT%s#GT#LT%s#GT}",laba.c_str(),labb.c_str(),laba.c_str(),labb.c_str()));
    // grk->SetTitle(Form(";N_{part};#LT%s%s#GT/#LT%s#GT#LT%s#GT",laba.c_str(),labb.c_str(),laba.c_str(),labb.c_str()));
    grk->SetTitle(";N_{part};#kappa_{h#gamma}      ");
    
    grk->SetMarkerStyle(kFullCircle);
    grk->SetMarkerColor(TColor::GetColorPalette(254*1/5));
    grk->SetMarkerSize(2);
    grk->SetLineColor(TColor::GetColorPalette(254*1/5));
    grk->SetLineWidth(1);
    
    grk->GetXaxis()->SetLimits(xmin,xmax);
    // grk->GetXaxis()->SetRangeUser(xmin,xmax);
    grk->Draw("AP");
    pad3->Update();
    
    // ymin = grk->GetHistogram()->GetMinimum();
    // ymax = grk->GetHistogram()->GetMaximum();
    // ymin = ymin - (ymax-ymin)*.1;
    // if (ymin > 1) ymin = 1;

    // grk->GetYaxis()->SetRangeUser(ymin,ymax);
    grk->GetYaxis()->SetLabelSize(.04);
    grk->GetYaxis()->SetTitleSize(.0475);
    grk->GetYaxis()->SetTitleOffset(.65);
    
    grk->GetXaxis()->SetLabelSize(.04);
    grk->GetXaxis()->SetTitleSize(.0475);
    grk->GetXaxis()->SetTitleOffset(.8);
    
    // grk->GetYaxis()->LabelsOption("h");//neither v nor h have an effect
    TLine* tl = new TLine(xmin,1.0,xmax,1.0);
    tl->SetLineStyle(9);
    tl->Draw();
    pad3->Update();

    can->SaveAs(Form("kdraw/%s.pdf",grp->GetName()));

    
    delete overlay1;
    
    delete overlay2;
    
    delete pad1;
    
    delete pad2;
    
    delete pad3;
    
    delete can;
    
}

//this method is an attempt to merge DrawKComp with the style from drawKComponents_terms.C DrawKPlot
// void DrawKCompsSmall(TGraphErrors* gra, TGraphErrors* grae, string laba, TGraphErrors* grb, TGraphErrors* grbe, string labb, TGraphErrors* grp, string labp, TH1F* hnpart[GAT::ncent], TGraphErrors* grk, int ik, int iptt, int ixi)
void DrawKCompsSmall(TGraphErrors* gra, TGraphErrors* grae, string laba, TGraphErrors* grb, TGraphErrors* grbe, string labb, TGraphErrors* grp, string labp, TH1F* hnpart[GAT::ncent], TGraphErrors* grk, int ik, int iptt, int ixi,TGraphErrors* grad = nullptr, TGraphErrors* grade = nullptr, TGraphErrors* grbd = nullptr, TGraphErrors* grbde = nullptr, TGraphErrors* sig)
{

    double xmin = npartmin;
    double xmax = npartmax;

    // gROOT->ForceStyle(0);
    gStyle->SetEndErrorSize(15);
    gStyle->SetPalette(kBird);
    // gStyle->SetLabelSize(.05);
    gStyle->SetLegendBorderSize(0);

    // TGaxis::SetExponentOffset(-.06,-.05,"y");
    TGaxis::SetMaxDigits(3);
    // TGaxis::SetExponentOffset(-.06,-.04,"y");
    // TGaxis::SetExponentOffset(-.04,0,"y");

    TF1* fa = new TF1(Form("%s_eval",gra->GetName()),[=](double* x, double* par){return gra->Eval(x[0]);},xmin,xmax,0);
    fa->SetLineColor(TColor::GetColorPalette(254*3/5));
    fa->SetLineWidth(2);

    TF1* fb = new TF1("fb",[=](double* x, double* par){return grb->Eval(x[0]);},xmin,xmax,0);
    // fb->SetLineColor(kRed);
    fb->SetLineColor(TColor::GetColorPalette(254*1/5));
    fb->SetLineWidth(2);

    //functions for kb only
    TF1* fad = nullptr;
    TF1* fbd = nullptr;
    if (grad && grbd)
    {
        //these functions are not being drawn. no need to set color etc.
        fad = new TF1("fad",[=](double* x, double* par){return grad->Eval(x[0]);},xmin,xmax,0);
        // fad->SetLineColor(TColor::GetColorPalette(254*3/5));
        // fad->SetLineWidth(2);

        fbd = new TF1("fbd",[=](double* x, double* par){return grbd->Eval(x[0]);},xmin,xmax,0);
        // fbd->SetLineColor(TColor::GetColorPalette(254*1/5));
        // fbd->SetLineWidth(2);
    }   

    //so vexp* aren't being used
    // vector<double> vmean, vrms, vexpa, vexpae, vexpb, vexpbe, vprod, vprode;
    vector<double> vmean, vrms, vprod, vprode;    

    //this loop is what must be altered for kb.
    //we want to recalculate <bgs> and <bhs> as is done in calculation of kb.
    //although sigma is not being recalculated...
    //honestly the recalculation prob not nec for this drawing..
    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        double mean = hnpart[icent]->GetMean();
        vmean.push_back(mean);
        
        double rms = hnpart[icent]->GetRMS();
        vrms.push_back(rms);

        //Recalculate <a> and plot as single point tgrapherror. copying from calcKComponents
        double expa = -9999.;
        double expae = -9999.;

        if (fad)
            CalcExp(fad,grade,hnpart[icent],expa,expae);
        else
            CalcExp(fa,grae,hnpart[icent],expa,expae);
        
        // vexpa.push_back(expa);
        // vexpae.push_back(expae);

        double expb = -9999.;
        double expbe = -9999.;       

        if (fbd)
            CalcExp(fbd,grbde,hnpart[icent],expb,expbe);
        else
            CalcExp(fb,grbe,hnpart[icent],expb,expbe);

        // vexpb.push_back(expb);
        // vexpbe.push_back(expbe);
                
        double prod = expa*expb;
        double prode = sqrt(pow(expa*expbe,2.) + pow(expae*expb,2.));
        
        if (ik == 1)//if kb
        {
            //just check
            if (!fbd || !fad)
            {
                cerr<<"Something wen't wrong"<<endl;
                return;
            }
            prod /= sig->GetY()[icent]*sig->GetY()[icent];
            prode /= sig->GetY()[icent]*sig->GetY()[icent];
        }
        
        vprod.push_back(prod);  
        vprode.push_back(prode);
    }
    
    // TCanvas* can = new TCanvas("can","",1700,1900);
    // TCanvas* can = new TCanvas("can","",1200,4500);
    TCanvas* can = new TCanvas("can","",1200,4000);
    TPad *pad1 = new TPad("pad1","",0,0.67,1,1);
    pad1->SetBottomMargin(0.005);
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.15);
    pad1->SetTopMargin(0.05);
    pad1->Draw();
    pad1->cd();

    double ymin = gra->Eval(xmin);
    double ymax = gra->Eval(xmin);
    for (int ix = xmin; ix < xmax; ++ix)
    {
        double eval = gra->Eval(ix);
        if (eval < ymin) ymin = eval;
        if (eval > ymax) ymax = eval;
    }

    double ydiff = ymax - ymin;
    ymax = ymax + .15*ydiff;
    ymin = ymin - .05*ydiff;

    TH1F *frame1 = pad1->DrawFrame(xmin,ymin,xmax,ymax);
    frame1->SetTitle("");
    // frame1->GetYaxis()->SetLabelSize(.045);
    // frame1->SetTitle(grp->GetName());
    frame1->GetYaxis()->SetAxisColor(TColor::GetColorPalette(254*3/5));
    frame1->GetYaxis()->SetTitle(Form("%s",laba.c_str()));
    // frame1->GetYaxis()->SetTitle(laba.c_str());
    frame1->GetYaxis()->SetLabelSize(.045);
    frame1->GetYaxis()->SetTitleOffset(1.4);
    frame1->GetYaxis()->SetTitleSize(.05);
    frame1->GetYaxis()->SetMaxDigits(2);
    frame1->GetXaxis()->SetLabelSize(0.0);

    int font = fa->GetYaxis()->GetLabelFont();

    pad1->Update();
        
    // fa->SetTitle("");
    fa->GetXaxis()->SetLimits(xmin,xmax);

    // fa->GetYaxis()->SetAxisColor(TColor::GetColorPalette(254*3/5));
    // fa->GetYaxis()->SetTitle(Form("%s",laba.c_str()));
    // // fa->GetYaxis()->SetTitle(laba.c_str());
    // fa->GetYaxis()->SetLabelSize(.045);
    // fa->GetYaxis()->SetTitleOffset(1.4);
    // fa->GetYaxis()->SetTitleSize(.05);
    // fa->GetYaxis()->SetMaxDigits(2);
    // fa->GetXaxis()->SetLabelSize(0.0);
    
    fa->Draw("same");

    pad1->Update();
    

    //create a transparent pad drawn on top of the main pad
    can->cd();
    TPad *overlay1 = new TPad("overlay1","",0,0.67,1,1);

    overlay1->SetBottomMargin(0.005);
    overlay1->SetLeftMargin(0.15);
    overlay1->SetRightMargin(0.15);
    overlay1->SetTopMargin(0.05);

    // overlay1->SetBottomMargin(0.0);
    // overlay1->SetTopMargin(0.005);
    overlay1->SetFillStyle(0);
    
    overlay1->SetFillColor(0);
    overlay1->SetFrameFillStyle(0);
    
    // overlay1->Draw("FA");
    overlay1->Draw();
    overlay1->cd();
    
    ymin = grb->Eval(xmin);
    ymax = grb->Eval(xmin);
    for (int ix = xmin; ix < xmax; ++ix)
    {
        double eval = grb->Eval(ix);
        if (eval < ymin) ymin = eval;
        if (eval > ymax) ymax = eval;
    }

    ydiff = ymax - ymin;
    ymax = ymax + .15*ydiff;
    ymin = ymin - .05*ydiff;

    // ymax = ymax + fabs(ymax)*.1;
    // ymin = ymin - fabs(ymin)*.1;
    // ymin = ymin - (ymax - ymin)*.1;

    TH1F *frame2 = overlay1->DrawFrame(xmin,ymin,xmax,ymax);
    frame2->GetXaxis()->SetLabelOffset(99);
    frame2->GetYaxis()->SetTickLength(0);
    frame2->GetYaxis()->SetLabelOffset(99);
    // grb->Draw("LP");

    fb->SetTitle("");
    fb->GetXaxis()->SetLimits(xmin,xmax);
    // fb->GetXaxis()->SetRangeUser(xmin,xmax);
    fb->Draw("same");
    overlay1->Update();

    //Draw an axis on the right side
    // TGaxis::SetExponentOffset(.01,-.01,"y");
    TGaxis *axis1 = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
    // TAxis *axis1 = new TAxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
    axis1->SetLineColor(TColor::GetColorPalette(254*1/5));
    // axis1->SetAxisColor(TColor::GetColorPalette(254*1/5));
    axis1->SetTitle(Form("%s",labb.c_str()));
    // axis1->SetTitle(labb.c_str());
    // axis1->fAxis->RotateTitle();//fAxis is protected. TGaxis has no method rotate
    axis1->SetTitleFont(font);
    axis1->SetTitleSize(.05);
    axis1->SetTitleOffset(1.4);
    axis1->SetLabelSize(.045);
    axis1->SetLabelFont(font);
    // axis1->SetExponentOffset(.06,.04);
    // axis1->SetNoExponent();
    axis1->Draw();
    // TGaxis::SetExponentOffset(-.01,-.01,"y");
    
    // TLegend* leg1 = new TLegend(0.75,0.11,0.81,0.285);
    // leg1->SetMargin(.35);
    TLegend* leg1;
    if (ik == 2 || ik == 3)
        leg1 = new TLegend(0.20,0.825,0.35,0.95);
    else
        leg1 = new TLegend(0.20,0.825,0.325,0.95);

    leg1->SetFillStyle(0);
    // leg1->SetLegendBorderSize(0);
    leg1->SetMargin(.35);
    leg1->SetTextSize(.05);
    leg1->AddEntry(fa,laba.c_str(),"l");
    leg1->AddEntry(fb,labb.c_str(),"l");
    leg1->Draw();

    //---------- center pad ----------

    can->cd();

    TPad *pad2 = new TPad("pad2","",0,.34,1,0.67);
    pad2->SetBottomMargin(0.02);
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.15);
    pad2->SetTopMargin(0.05);
    // pad2->SetFillColor(42);
    // pad2->SetGrid();
    pad2->Draw();
    pad2->cd();
    
    ymin = grp->Eval(xmin);
    ymax = grp->Eval(xmin);
    for (int ix = xmin; ix < xmax; ++ix)
    {
        double eval = grp->Eval(ix);
        if (eval < ymin) ymin = eval;
        if (eval > ymax) ymax = eval;
    }
    
    ydiff = ymax - ymin;
    ymax = ymax + .15*ydiff;
    ymin = ymin - .05*ydiff;

    // ymin = ymin - fabs(ymin)*.1;
    // ymax = ymax + fabs(ymax)*.1;
    // ymin = ymin - (ymax - ymin)*.1;

    grp->SetMarkerSize(0);
    grp->SetLineColor(TColor::GetColorPalette(254*2/5));
    grp->SetLineWidth(2);
    // grp->GetXaxis()->SetRangeUser(xmin,xmax);
    grp->GetXaxis()->SetLimits(xmin,xmax);

    grp->GetYaxis()->SetTitleSize(.05);
    grp->GetYaxis()->SetTitleOffset(1.4);
    grp->GetYaxis()->SetLabelSize(.045);

    grp->GetXaxis()->SetLabelSize(0.);
    
    grp->SetTitle(Form(";;%s",labp.c_str()));
    grp->Draw("ACX");

    grp->GetHistogram()->SetMinimum(ymin);
    grp->GetHistogram()->SetMaximum(ymax);
    pad2->Update();

    TGraphErrors* grprodall = new TGraphErrors(GAT::ncent,&vmean[0],&vprod[0],&vrms[0],&vprode[0]);
    grprodall->SetMarkerStyle(kFullCircle);
    // grprodall->SetMarkerStyle(21);
    grprodall->SetMarkerSize(4);
    // grprodall->SetMarkerColor(kMagenta);
    // grprodall->SetLineColor(kMagenta);
    grprodall->SetMarkerColor(TColor::GetColorPalette(254*0/5));
    grprodall->SetLineColor(TColor::GetColorPalette(254*0/5));
    grprodall->SetLineWidth(1);
    // grprodall->GetXaxis()->SetRangeUser(xmin,xmax);
    grprodall->GetXaxis()->SetLimits(xmin,xmax);
    grprodall->Draw("P same");

    TLegend* leg2;
    if (ik == 2 || ik == 3)
        leg2 = new TLegend(0.20,0.8,0.45,0.95);
    else
        leg2 = new TLegend(0.20,0.8,0.40,0.95);
    // leg2->SetTextAlign(12);//default anyway
    // leg2->SetLegendBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetMargin(.225);
    leg2->SetTextSize(0.05);
    // leg2->SetFillStyle(0);
    // leg2->SetBorderSize(0);

    //needs to be changed if ik == 1. legend may need resizing as well
    if (ik == 1)
        leg2->AddEntry(grprodall,"#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT/#LT#sigma#GT^{2}");
    // leg2->AddEntry(grprodall,"#LT#sigma#GT^{2}/(#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT)");
    // leg2->AddEntry(grprodall,"#frac{#LT#sigma#GT^{2}}{#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT}");
    else
        leg2->AddEntry(grprodall,Form("#LT%s#GT#LT%s#GT",laba.c_str(),labb.c_str()),"p");
    
    leg2->AddEntry(grp,labp.c_str(),"l");

    leg2->Draw();

    //-------------------- bottom pad --------------------

    can->cd();
    TPad* pad3 = new TPad("pad3","",0,0,1,.34);
    pad3->SetBottomMargin(0.11);
    pad3->SetLeftMargin(0.15);
    pad3->SetRightMargin(0.15);
    pad3->SetTopMargin(0.01);
    pad3->Draw();
    pad3->cd();


    // grk->SetTitle(Form(";N_{part};#frac{<%s%s>}{<%s><%s>}",laba.c_str(),labb.c_str(),laba.c_str(),labb.c_str()));
    // grk->SetTitle(Form(";N_{part};#frac{#LT%s%s#GT}{#LT%s#GT#LT%s#GT}",laba.c_str(),labb.c_str(),laba.c_str(),labb.c_str()));
    // grk->SetTitle(Form(";N_{part};#LT%s%s#GT/#LT%s#GT#LT%s#GT",laba.c_str(),labb.c_str(),laba.c_str(),labb.c_str()));
    grk->SetTitle(Form(";N_{part};%s",kname[ik].c_str()));
    
    grk->SetMarkerStyle(kFullCircle);
    grk->SetMarkerColor(TColor::GetColorPalette(254*1/5));
    grk->SetMarkerSize(4);
    grk->SetLineColor(TColor::GetColorPalette(254*1/5));
    grk->SetLineWidth(1);
    
    // grk->GetXaxis()->SetRangeUser(xmin,xmax);
    grk->GetXaxis()->SetLimits(xmin,xmax);
    grk->Draw("AP");
    pad3->Update();
    
    ymin = grk->GetHistogram()->GetMinimum();
    ymax = grk->GetHistogram()->GetMaximum();
    grk->GetYaxis()->SetNdivisions(505);
    
    grk->GetYaxis()->SetLabelSize(.045);
    grk->GetYaxis()->SetTitleSize(.05);
    grk->GetYaxis()->SetTitleOffset(1.4);
    
    grk->GetXaxis()->SetLabelSize(.045);
    grk->GetXaxis()->SetTitleSize(.05);
    grk->GetXaxis()->SetTitleOffset(1.0);

    if (ik == 4 || ik == 5 || ik == 1)
    {
        txtmomlow[iptt]->Draw();
        txtxilow[ixi]->Draw();
    }
    else        
    {
        txtmom[iptt]->Draw();
        txtxi[ixi]->Draw();
    }
    
    // grk->GetYaxis()->LabelsOption("h");//neither v nor h have an effect
    TLine* tl = new TLine(xmin,1.0,xmax,1.0);
    tl->SetLineStyle(7);
    tl->SetLineWidth(1);
    tl->Draw();
    pad3->Update();

    can->SaveAs(Form("kdraw/%s_small.pdf",grp->GetName()));
    
    delete overlay1;    
    delete pad1;    
    delete pad2;    
    delete pad3;    
    delete can;    
}

void CalcExp(TF1* fa, TGraphErrors* grae, TH1F* hnpart, double& expa, double& expae)
{
    int nbinx = hnpart->GetNbinsX();
        
    // double expa = 0.;
    // double expae = 0.;

    expa = 0.;
    expae = 0.;
        
    for (int ibin = 1; ibin <= nbinx; ibin++)
    {
        double w = hnpart->GetBinContent(ibin);
        // wtotal += w;
        if (w <= 0.0) continue;
            
        double we = hnpart->GetBinError(ibin);
            
        double npart = hnpart->GetBinCenter(ibin);          
            
        double a = fa->Eval(npart);
        if (a < 0.)
        {
            ostringstream oss;
            oss<<fa->GetName();
            string str = oss.str();
            if (str[0] == 'a')
            {
                cout<<fa->GetName()<<"->Eval("<<npart<<") evaluated to "<<a<<". Setting to 0."<<endl;
                continue;
                // a = 0;
            }
        }

        double ae = 0.;
        if (grae)
            ae = grae->Eval(npart);
        if (ae < 0.) ae = 0.;                       
            
        // cout<<"<a> = "<<expa<<" + ("<<w<<"*"<<a<<" = "<<w*a<<") = ";
        expa += w*a;
        // cout<<expa<<endl;

        expae += w*ae;

    }

    return;
}
