//make some spiffy plots with kcomponent output
//much is copied directly from calcKComponents.C....

enum RPO {IN,OUT,ALL};

const int nadj = 2;

const int npartmin = 2;
const int npartmax = 390;

const float PI = (float) TMath::Pi();

// const int marksty[4] = {20,21,33,47};
// const int absmarksty[4] = {24,25,27,46};
const int marksty[4] = {20,21,22,23};
const int absmarksty[4] = {24,25,26,32};
const int markcol[4] = {1,4,2,3};

const string klab[4][3] = {
    {"a_{t} a_{h}","#LTa_{t} a_{h}#GT","#LTa_{t}#GT#LTa_{h}#GT"},
    {"b_{t} b_{h}","#LTb_{t} b_{h}#GT","#LTb_{t}#sigma#GT#LTb_{h}#sigma#GT"},
    {"a_{t} b_{h} #sigma","#LTa_{t} b_{h} #sigma#GT","#LTa_{t}#GT#LTb_{h}#sigma#GT"},
    {"a_{h} b_{t} #sigma","#LTa_{h} b_{t} #sigma#GT","#LTa_{h}#GT#LTb_{t}#sigma#GT"}};

const string kname[4] {"#kappa_{a}","#kappa_{b}","#kappa_{th}","#kappa_{ht}"};

const double lt = .15;
const double rt = .05;
const double legw = .275;
const double legw2 = .625;//legend width for 2 columns

const int errsize = 20;

namespace GAT {
    #include "../src/gamAnaTree/gamAnaTree.h"
}

int GetNpartDists(TH1F* hnpart[GAT::ncent], vector<double>& vnpart, vector<double>& vnparte);

void AddToMultiGraph(TGraphErrors* tgr, TMultiGraph& mgr, TLegend& leg, int igr, bool isxi = true, bool drawl = false, bool doabs = true);

void DrawKPlot(TGraphErrors* tgrnum, TGraphErrors* tgrden, TGraphErrors* tgrk, const char* name, const int ik, const int iptt, const int ixi);

void PrintData(TGraphErrors* tgr[nadj][GAT::ntype+1][GAT::nptt][GAT::nxi],string pref);

void DrawMultiGraph(TCanvas& can, TMultiGraph& mgr, bool setlogy = true, bool noex = false);

TLatex* txtmom[GAT::nptt];
TLatex* txtmom_vh[GAT::nptt];
TLatex* txtxi[GAT::nxi];
// TLatex* txtxi3[GAT::nxi];
TLatex* txtmom2[GAT::nptt];
// TLatex* txtmom3[GAT::nptt];

void drawKComponents_terms(
    const char* inname = "kComponents.root"
    )
{
    gStyle->SetEndErrorSize(errsize);
    // gStyle->SetMarkerSize(4.5);
    // gStyle->SetLineWidth(5);
    gStyle->SetTitleOffset(1.7,"Y");
    gStyle->SetTitleOffset(1.05,"X");
    gStyle->SetImageScaling(3.);

    for (int iptt = 0; iptt < GAT::nptt; iptt++)
    {
        // txtmom[iptt] = new TLatex (0.475,.85,Form("%.1f #leq p_{T} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
        // txtmom[iptt] = new TLatex (0.4,.825,Form("%.1f #leq p_{T} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
        txtmom[iptt] = new TLatex ((lt+1.-rt-legw)/2.,.825,Form("%.1f #leq p_{T} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
        txtmom[iptt]->SetTextAlign(22);
        txtmom[iptt]->SetNDC();
        txtmom[iptt]->SetTextSize(0.0375);

        txtmom_vh[iptt] = new TLatex (0.45,.85,Form("%.1f #leq p_{T} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
        txtmom_vh[iptt]->SetTextAlign(22);
        txtmom_vh[iptt]->SetNDC();
        txtmom_vh[iptt]->SetTextSize(0.0375);
        
        txtmom2[iptt] = new TLatex (0.55,.85,Form("%.1f #leq p_{T} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
        txtmom2[iptt]->SetTextAlign(22);
        txtmom2[iptt]->SetNDC();
        txtmom2[iptt]->SetTextSize(0.05);

        // txtmom3[iptt] = new TLatex (0.55,.35,Form("%.1f #leq p_{T} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
        // txtmom3[iptt]->SetTextAlign(22);
        // txtmom3[iptt]->SetNDC();
        // txtmom3[iptt]->SetTextSize(0.05);
    }

    for (int ixi = 0; ixi < GAT::nxi; ixi++)
    {
        txtxi[ixi] = new TLatex (0.55,.775,Form("%.1f #leq #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]));
        txtxi[ixi]->SetTextAlign(22);
        txtxi[ixi]->SetNDC();
        txtxi[ixi]->SetTextSize(0.05);

        // txtxi3[ixi] = new TLatex (0.55,.175,Form("%.1f #leq #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]));
        // txtxi3[ixi]->SetTextAlign(22);
        // txtxi3[ixi]->SetNDC();
        // txtxi3[ixi]->SetTextSize(0.05);
    }
        
    //Get Npart distributions
    TH1F* hnpart[GAT::ncent];
    vector<double> vnpart;
    vector<double> vnparte;
    if (GetNpartDists(hnpart,vnpart,vnparte)<0) return;

    //open input
    TFile* input = new TFile(inname,"READ");
    if (!input || input->IsZombie())
    {
        cerr<<"Could not open "<<inname<<endl;
        return;
    }

    TGraphErrors* agah4pi[nadj][GAT::ntype+1][GAT::nptt][GAT::nxi];
    TGraphErrors* bgsbhs2pi[nadj][GAT::ntype+1][GAT::nptt][GAT::nxi];
    TGraphErrors* agbhs8[nadj][GAT::ntype+1][GAT::nptt][GAT::nxi];
    TGraphErrors* ahbgs8[nadj][GAT::ntype+1][GAT::nptt][GAT::nxi];

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            // if (itype != GAT::PI0) continue;

            TMultiGraph mgr_ag ("mgr_ag",";#LTN_{part}#GT;a_{t}");
            // TLegend leg_ag (0.55,0.15,0.85,0.4);
            TLegend leg_ag (1.-rt-legw2,0.8,1.-rt,0.9);
            leg_ag.SetTextSize(0.0375);

            TMultiGraph mgr_bg[GAT::nxi];
            TLegend leg_bg[GAT::nxi];
            for (int ixi = 0; ixi < GAT::nxi; ixi++)
            {
                mgr_bg[ixi].SetTitle(";#LTN_{part}#GT;b_{t}");
                // leg_bg[ixi].SetX1(0.275);
                leg_bg[ixi].SetX1(1.-rt-legw2);
                leg_bg[ixi].SetY1(0.8);
                leg_bg[ixi].SetX2(01.-rt);
                leg_bg[ixi].SetY2(0.9);
                leg_bg[ixi].SetTextSize(0.0375);
            }
            
            TMultiGraph mgr_bgs[GAT::nxi];
            TLegend leg_bgs[GAT::nxi];
            for (int ixi = 0; ixi < GAT::nxi; ixi++)
            {
                mgr_bgs[ixi].SetTitle(";#LTN_{part}#GT;b_{t}#sigma");
                leg_bgs[ixi].SetX1(1.-rt-legw2);
                leg_bgs[ixi].SetY1(0.8);
                leg_bgs[ixi].SetX2(01.-rt);
                leg_bgs[ixi].SetY2(0.9);
                leg_bgs[ixi].SetTextSize(0.0375);
            }
            
            TMultiGraph mgr_vg[GAT::nxi];
            TLegend leg_vg[GAT::nxi];
            for (int ixi = 0; ixi < GAT::nxi; ixi++)
            {
                if (itype == GAT::PI0)
                    mgr_vg[ixi].SetTitle(";#LTN_{part}#GT;v_{2,#pi0}");
                else
                    mgr_vg[ixi].SetTitle(";#LTN_{part}#GT;v_{2,t}");
                leg_vg[ixi].SetX1(01.-rt-legw2);
                leg_vg[ixi].SetY1(0.8);
                leg_vg[ixi].SetX2(01.-rt);
                leg_vg[ixi].SetY2(0.9);
                leg_vg[ixi].SetTextSize(0.0375);
            }
            
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                ostringstream sufg;
                sufg<<"adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt;

                TGraphErrors* tgrag = (TGraphErrors*) gROOT->FindObject(Form("ag_Np_%s",sufg.str().c_str()));
                TGraphErrors* tgrage = (TGraphErrors*) gROOT->FindObject(Form("age_Np_%s",sufg.str().c_str()));

                AddToMultiGraph(tgrag,mgr_ag,leg_ag,iptt,false);

                TMultiGraph mgr_flat0 ("mgr_flat0",";#LTN_{part}#GT;4#pi#LTa_{t}#GT#LTa_{h}#GT");
                // TLegend leg_flat0 (0.58,0.13,0.88,0.33);
                TLegend leg_flat0 (1.-rt-legw,0.75,1.-rt,0.9);
                leg_flat0.SetTextSize(0.0375);

                TMultiGraph mgr_cos0 ("mgr_cos0",";#LTN_{part}#GT;2#pi#LTb_{t}#sigma#GT#LTb_{h}#sigma#GT");;
                TLegend leg_cos0 (1.-rt-legw,0.75,1.-rt,0.9);

                leg_cos0.SetTextSize(0.0375);
                TMultiGraph mgr_flat1 ("mgr_flat1",";#LTN_{part}#GT;8#LTa_{h}#GT#LTb_{t}#sigma#GT");
                TLegend leg_flat1 (1.-rt-legw,0.75,1.-rt,0.9);
                leg_flat1.SetTextSize(0.0375);
                
                TMultiGraph mgr_cos1 ("mgr_cos1",";#LTN_{part}#GT;8#LTa_{t}#GT#LTb_{h}#sigma#GT");
                TLegend leg_cos1 (1.-rt-legw,0.75,1.-rt,0.9);
                leg_cos1.SetTextSize(0.0375);

                TMultiGraph mgr_ah ("mgr_ah",";#LTN_{part}#GT;a_{h}");
                // TLegend leg_ah (0.58,0.13,0.88,0.33);
                TLegend leg_ah (1.-rt-legw,0.75,1.-rt,0.9);
                leg_ah.SetTextSize(0.0375);
                
                TMultiGraph mgr_vh ("mgr_vh",";#LTN_{part}#GT;v_{h}");
                TLegend leg_vh (1.-rt-legw,0.75,1.-rt,0.9);
                leg_vh.SetTextSize(0.0375);
                
                TMultiGraph mgr_bh ("mgr_bh",";#LTN_{part}#GT;b_{h}");
                // TLegend leg_bh (0.48,0.45,0.78,0.65);
                // TLegend leg_bh (0.48,0.13,0.78,0.33);
                TLegend leg_bh (1.-rt-legw,0.75,1.-rt,0.9);
                leg_bh.SetTextSize(0.0375);
                
                TMultiGraph mgr_bhs ("mgr_bhs",";#LTN_{part}#GT;b_{h}#sigma");
                // TLegend leg_bhs (0.48,0.45,0.78,0.65);
                // TLegend leg_bhs (0.48,0.13,0.78,0.33);
                TLegend leg_bhs (1.-rt-legw,0.75,1.-rt,0.9);
                leg_bhs.SetTextSize(0.0375);
                
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    ostringstream sufh;
                    sufh<<sufg.str()<<"_xi"<<ixi;
                    
                    //---------- first flat term ----------

                    TGraphErrors* tgrah = (TGraphErrors*) gROOT->FindObject(Form("ah_Np_%s",sufh.str().c_str()));
                    AddToMultiGraph(tgrah,mgr_ah,leg_ah,ixi);
                    
                    TGraphErrors* tgrahe = (TGraphErrors*) gROOT->FindObject(Form("ahe_Np_%s",sufh.str().c_str()));

                    ostringstream ossagah;
                    ossagah<<"agah4pi_"<<sufh.str();
                    agah4pi[iadj][itype][iptt][ixi] = new TGraphErrors();
                    for (int i = 0; i < tgrag->GetN(); i++)
                    {
                        double np = tgrag->GetX()[i];
                        double npe = tgrag->GetEX()[i];
                        
                        double ag = tgrag->GetY()[i];
                        double age = tgrag->GetEY()[i];
                        double ah = tgrah->GetY()[i];
                        double ahe = tgrah->GetEY()[i];
                        
                        agah4pi[iadj][itype][iptt][ixi]->SetPoint(i,np,ag*ah*4*PI);
                        agah4pi[iadj][itype][iptt][ixi]->SetPointError(i,npe,4*PI*sqrt(pow(ag*ahe,2) + pow(age*ah,2)));
                    }
                    agah4pi[iadj][itype][iptt][ixi]->SetName(ossagah.str().c_str());

                    AddToMultiGraph(agah4pi[iadj][itype][iptt][ixi],mgr_flat0,leg_flat0,ixi,true,false);
                    
                    //---------- first cosine term ----------               

                    bgsbhs2pi[iadj][itype][iptt][ixi] = (TGraphErrors*) gROOT->FindObject(Form("bgsbhs_Np_%s",sufh.str().c_str()));
                    for (int i = 0; i < bgsbhs2pi[iadj][itype][iptt][ixi]->GetN(); i++)
                    {
                        bgsbhs2pi[iadj][itype][iptt][ixi]->GetY()[i] *= 2*PI;
                    }

                    AddToMultiGraph(bgsbhs2pi[iadj][itype][iptt][ixi],mgr_cos0,leg_cos0,ixi,true,false,true);

                    //---------- second flat term ----------

                    TGraphErrors* tgrvg = (TGraphErrors*) gROOT->FindObject(Form("vg_Np_%s",sufh.str().c_str()));
                    AddToMultiGraph(tgrvg,mgr_vg[ixi],leg_vg[ixi],iptt,false);

                    TGraphErrors* tgrbg = (TGraphErrors*) gROOT->FindObject(Form("bg_Np_%s",sufh.str().c_str()));
                    AddToMultiGraph(tgrbg,mgr_bg[ixi],leg_bg[ixi],iptt,false);

                    TGraphErrors* tgrbgs = (TGraphErrors*) gROOT->FindObject(Form("bgs_Np_%s",sufh.str().c_str()));
                    AddToMultiGraph(tgrbgs,mgr_bgs[ixi],leg_bgs[ixi],iptt,false);

                    TGraphErrors* tgrbgse = (TGraphErrors*) gROOT->FindObject(Form("bgse_Np_%s",sufh.str().c_str()));

                    ostringstream ossahbgs;
                    ossahbgs<<"ahbgs8_"<<sufh.str();

                    ahbgs8[iadj][itype][iptt][ixi] = new TGraphErrors();
                    for (int i = 0; i < tgrah->GetN(); i++)
                    {
                        double np = tgrah->GetX()[i];
                        double npe = tgrah->GetEX()[i];
                        
                        double ah = tgrah->GetY()[i];
                        double ahe = tgrah->GetEY()[i];
                        double bgs = tgrbgs->GetY()[i];
                        double bgse = tgrbgs->GetEY()[i];

                        ahbgs8[iadj][itype][iptt][ixi]->SetPoint(i,np,ah*bgs*8);
                        ahbgs8[iadj][itype][iptt][ixi]->SetPointError(i,npe,8*sqrt(pow(ah*bgse,2) + pow(ahe*bgs,2)));
                    }
                    ahbgs8[iadj][itype][iptt][ixi]->SetName(ossahbgs.str().c_str());

                    AddToMultiGraph(ahbgs8[iadj][itype][iptt][ixi],mgr_flat1,leg_flat1,ixi,true,false,true);
                    
                    //---------- second cosine term ----------

                    TGraphErrors* tgrvh = (TGraphErrors*) gROOT->FindObject(Form("vh_Np_%s",sufh.str().c_str()));
                    AddToMultiGraph(tgrvh,mgr_vh,leg_vh,ixi);

                    TGraphErrors* tgrbh = (TGraphErrors*) gROOT->FindObject(Form("bh_Np_%s",sufh.str().c_str()));
                    AddToMultiGraph(tgrbh,mgr_bh,leg_bh,ixi);

                    TGraphErrors* tgrbhs = (TGraphErrors*) gROOT->FindObject(Form("bhs_Np_%s",sufh.str().c_str()));
                    AddToMultiGraph(tgrbhs,mgr_bhs,leg_bhs,ixi);
                    
                    TGraphErrors* tgrbhse = (TGraphErrors*) gROOT->FindObject(Form("bhse_Np_%s",sufh.str().c_str()));

                    ostringstream ossagbhs;
                    ossagbhs<<"agbhs8_"<<sufh.str();

                    agbhs8[iadj][itype][iptt][ixi] = new TGraphErrors();
                    for (int i = 0; i < tgrag->GetN(); i++)
                    {
                        double np = tgrag->GetX()[i];
                        double npe = tgrag->GetEX()[i];
                        
                        double ag = tgrag->GetY()[i];
                        double age = tgrag->GetEY()[i];
                        double bhs = tgrbhs->GetY()[i];
                        double bhse = tgrbhs->GetEY()[i];
                        
                        agbhs8[iadj][itype][iptt][ixi]->SetPoint(i,np,ag*bhs*8);
                        agbhs8[iadj][itype][iptt][ixi]->SetPointError(i,npe,8*sqrt(pow(ag*bhse,2) + pow(age*bhs,2)));
                    }
                    agbhs8[iadj][itype][iptt][ixi]->SetName(ossagbhs.str().c_str());

                    AddToMultiGraph(agbhs8[iadj][itype][iptt][ixi],mgr_cos1,leg_cos1,ixi,true,false,true);

                    //==================================================
                    // draw k num and den and ratio that is k
                    //ka, kb, kgh, khg
                    //==================================================

                    //---------- ka ----------
                    TGraphErrors* tgrkanum = (TGraphErrors*) gROOT->FindObject(Form("kanum_Np_%s",sufh.str().c_str()));
                    TGraphErrors* tgrkaden = (TGraphErrors*) gROOT->FindObject(Form("kaden_Np_%s",sufh.str().c_str()));
                    TGraphErrors* tgrka = (TGraphErrors*) gROOT->FindObject(Form("ka_Np_%s",sufh.str().c_str()));
                    ostringstream kaoss;
                    kaoss<<"termpics/ka_"<<sufh.str()<<".pdf";
                    // cout<<"Calling DrawKPlot() to create "<<kaoss.str()<<endl;
                    DrawKPlot(tgrkanum,tgrkaden,tgrka,kaoss.str().c_str(),0,iptt,ixi);
                    // DrawKComps(tgrag,age,"a_{t}",ah,ahe,"a_{h}",agah,"a_{t} a_{h}",ka,ik);

                    //---------- kb ----------
                    TGraphErrors* tgrkbnum = (TGraphErrors*) gROOT->FindObject(Form("kbnum_Np_%s",sufh.str().c_str()));
                    TGraphErrors* tgrkbden = (TGraphErrors*) gROOT->FindObject(Form("kbden_Np_%s",sufh.str().c_str()));
                    TGraphErrors* tgrkb = (TGraphErrors*) gROOT->FindObject(Form("kb_Np_%s",sufh.str().c_str()));
                    ostringstream kboss;
                    kboss<<"termpics/kb_"<<sufh.str()<<".pdf";
                    // cout<<"Calling DrawKPlot() to create "<<kboss.str()<<endl;
                    DrawKPlot(tgrkbnum,tgrkbden,tgrkb,kboss.str().c_str(),1,iptt,ixi);

                    //---------- kgh ----------
                    TGraphErrors* tgrkghnum = (TGraphErrors*) gROOT->FindObject(Form("kghnum_Np_%s",sufh.str().c_str()));
                    TGraphErrors* tgrkghden = (TGraphErrors*) gROOT->FindObject(Form("kghden_Np_%s",sufh.str().c_str()));
                    TGraphErrors* tgrkgh = (TGraphErrors*) gROOT->FindObject(Form("kgh_Np_%s",sufh.str().c_str()));
                    ostringstream kghoss;
                    kghoss<<"termpics/kgh_"<<sufh.str()<<".pdf";
                    // cout<<"Calling DrawKPlot() to create "<<kghoss.str()<<endl;
                    DrawKPlot(tgrkghnum,tgrkghden,tgrkgh,kghoss.str().c_str(),2,iptt,ixi);

                    //---------- khg ----------
                    TGraphErrors* tgrkhgnum = (TGraphErrors*) gROOT->FindObject(Form("khgnum_Np_%s",sufh.str().c_str()));
                    TGraphErrors* tgrkhgden = (TGraphErrors*) gROOT->FindObject(Form("khgden_Np_%s",sufh.str().c_str()));
                    TGraphErrors* tgrkhg = (TGraphErrors*) gROOT->FindObject(Form("khg_Np_%s",sufh.str().c_str()));
                    ostringstream khgoss;
                    khgoss<<"termpics/khg_"<<sufh.str()<<".pdf";
                    // cout<<"Calling DrawKPlot() to create "<<khgoss.str()<<endl;
                    DrawKPlot(tgrkhgnum,tgrkhgden,tgrkhg,khgoss.str().c_str(),3,iptt,ixi);
                    
                    
                }

                TCanvas can_flat0 ("can_flat0","",2000,2500);
                DrawMultiGraph(can_flat0,mgr_flat0);
                // double max_flat0 = mgr_flat0.GetHistogram()->GetMaximum();
                // mgr_flat0.GetHistogram()->SetMaximum(max_flat0*50);
                leg_flat0.Draw();
                txtmom[iptt]->Draw();
                can_flat0.SaveAs(Form("termpics/agah4pi_%s.pdf",sufg.str().c_str()));

                TCanvas can_cos0 ("can_cos0","",2000,2500);
                DrawMultiGraph(can_cos0,mgr_cos0);
                leg_cos0.Draw();
                txtmom[iptt]->Draw();
                can_cos0.SaveAs(Form("termpics/bgsbhs2pi_%s.pdf",sufg.str().c_str()));                      

                TCanvas can_flat1 ("can_flat1","",2000,2500);
                DrawMultiGraph(can_flat1,mgr_flat1);
                // double max_flat1 = mgr_flat1.GetHistogram()->GetMaximum();
                // mgr_flat1.GetHistogram()->SetMaximum(max_flat1*50);
                leg_flat1.Draw();
                txtmom[iptt]->Draw();
                can_flat1.SaveAs(Form("termpics/ahbgs8_%s.pdf",sufg.str().c_str()));

                TCanvas can_cos1 ("can_cos1","",2000,2500);
                DrawMultiGraph(can_cos1,mgr_cos1);
                // double max_cos1 = mgr_cos1.GetHistogram()->GetMaximum();
                // mgr_cos1.GetHistogram()->SetMaximum(max_cos1*50);            
                leg_cos1.Draw();
                txtmom[iptt]->Draw();
                can_cos1.SaveAs(Form("termpics/agbhs8_%s.pdf",sufg.str().c_str()));

                TCanvas can_ah ("can_ah","",2000,2500);
                DrawMultiGraph(can_ah,mgr_ah);          
                // double max_ah = mgr_ah.GetHistogram()->GetMaximum();
                // mgr_ah.GetHistogram()->SetMaximum(max_ah*50);
                leg_ah.Draw();
                txtmom[iptt]->Draw();
                can_ah.SaveAs(Form("termpics/ah_%s.pdf",sufg.str().c_str()));

                TCanvas can_vh ("can_vh","",2000,2500);
                DrawMultiGraph(can_vh,mgr_vh,false);
                // mgr_vh.GetYaxis()->SetRangeUser(.040,.235);
                leg_vh.Draw();
                txtmom_vh[iptt]->Draw();
                can_vh.SaveAs(Form("termpics/vh_%s.pdf",sufg.str().c_str()));

                TCanvas can_bh ("can_bh","",2000,2500);
                DrawMultiGraph(can_bh,mgr_bh);
                // double max_bh = mgr_bh.GetHistogram()->GetMaximum();
                // mgr_bh.GetHistogram()->SetMaximum(max_bh*50);                
                leg_bh.Draw();
                txtmom[iptt]->Draw();                           
                can_bh.SaveAs(Form("termpics/bh_%s.pdf",sufg.str().c_str()));

                TCanvas can_bhs ("can_bhs","",2000,2500);
                DrawMultiGraph(can_bhs,mgr_bhs);
                // double max_bhs = mgr_bhs.GetHistogram()->GetMaximum();
                // mgr_bhs.GetHistogram()->SetMaximum(max_bhs*50);              
                leg_bhs.Draw();
                txtmom[iptt]->Draw();
                can_bhs.SaveAs(Form("termpics/bhs_%s.pdf",sufg.str().c_str()));
            }

            TCanvas can_ag ("can_ag","",2000,2500);
            DrawMultiGraph(can_ag,mgr_ag);
            // mgr_ag.GetYaxis()->SetRangeUser(2e-8,6e-4);
            // mgr_ag.GetYaxis()->SetTitleOffset(1.4);
            leg_ag.SetNColumns(2);
            leg_ag.Draw();
            can_ag.SaveAs(Form("termpics/ag_adj%i_%s.pdf",iadj,GAT::spec[itype]));

            for (int ixi = 0; ixi < GAT::nxi; ixi++)
            {
                TCanvas can_vg (Form("can_vg_xi%i",ixi),"",2000,2500);
                DrawMultiGraph(can_vg,mgr_vg[ixi],false);
                // mgr_vg[ixi].GetYaxis()->SetRangeUser(.01,.19);
                leg_vg[ixi].SetNColumns(2);
                leg_vg[ixi].Draw();
                can_vg.SaveAs(Form("termpics/vg_adj%i_%s_xi%i.pdf",iadj,GAT::spec[itype],ixi));

                TCanvas can_bg (Form("can_bg_xi%i",ixi),"",2000,2500);
                DrawMultiGraph(can_bg,mgr_bg[ixi]);
                // mgr_bg[ixi].GetYaxis()->SetRangeUser(3e-9,8e-5);
                leg_bg[ixi].SetNColumns(2);
                leg_bg[ixi].Draw();
                can_bg.SaveAs(Form("termpics/bg_adj%i_%s_xi%i.pdf",iadj,GAT::spec[itype],ixi));

                TCanvas can_bgs (Form("can_bgs_xi%i",ixi),"",2000,2500);
                DrawMultiGraph(can_bgs,mgr_bgs[ixi]);
                // mgr_bgs[ixi].GetYaxis()->SetRangeUser(3e-9,8e-5);
                leg_bgs[ixi].SetNColumns(2);
                leg_bgs[ixi].Draw();
                can_bgs.SaveAs(Form("termpics/bgs_adj%i_%s_xi%i.pdf",iadj,GAT::spec[itype],ixi));
            }
                
        }
    }

    PrintData(bgsbhs2pi,"bgsbhs2pi");

    cout<<"Completed drawKComponents_terms!"<<endl;
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

void AddToMultiGraph(TGraphErrors* tgr, TMultiGraph& mgr, TLegend& leg, int igr, bool isxi = true, bool drawl = false, bool doabs = true)
{
    double marksize = 6.0;

    //this way multigraph doesn't delete my tgraph dammit
    TGraphErrors* tgrclone0 = (TGraphErrors*) tgr->Clone();

    // used to show end errors when error less than marker size
    // TGraphErrors* tgrclone1 = (TGraphErrors*) tgr->Clone();
    // tgrclone1->SetMarkerSize(0.);
    // mgr.Add(tgrclone1);    

    if (doabs)
    {
        TGraphErrors* tgrabs = new TGraphErrors();
        TGraphErrors* tgrpos = new TGraphErrors();//positive values only
        int iptabs = 0;
        int iptpos = 0;
        for (int ipt = 0; ipt < tgr->GetN(); ipt++)
        {
            double x = tgr->GetX()[ipt];
            double y = tgr->GetY()[ipt];
            double ex = tgr->GetEX()[ipt];
            double ey = tgr->GetEY()[ipt];
            if (y >= 0)
            {
                tgrpos->SetPoint(iptpos,x,y);
                tgrpos->SetPointError(iptpos,ex,ey);

                iptpos++;
            }
            else
            {
                tgrabs->SetPoint(iptabs,x,-y);
                tgrabs->SetPointError(iptabs,ex,ey);

                iptabs++;
            }       
        }

        if (iptabs > 0)
        {
            tgrabs->SetMarkerStyle(absmarksty[igr]);
            tgrabs->SetMarkerColor(TColor::GetColorPalette(220*igr/3));
            tgrabs->SetLineColor(TColor::GetColorPalette(220*igr/3));
            tgrabs->SetLineWidth(2);
            tgrabs->SetMarkerSize(marksize);
    
            mgr.Add(tgrabs,"p");

            tgrclone0 = tgrpos;
        }
    }

    tgrclone0->SetMarkerStyle(marksty[igr]);
    // tgrclone0->SetMarkerColor(markcol[igr]);
    // tgrclone0->SetLineColor(markcol[igr]);
    tgrclone0->SetMarkerColor(TColor::GetColorPalette(220*igr/3));
    tgrclone0->SetLineColor(TColor::GetColorPalette(220*igr/3));
    tgrclone0->SetLineWidth(2);
    tgrclone0->SetMarkerSize(marksize);

    if (drawl)
        mgr.Add(tgrclone0,"lp");
    else
        mgr.Add(tgrclone0,"p"); 
    if (isxi)
        leg.AddEntry(tgrclone0,Form("#xi #in [%1.1f,%1.1f)",GAT::limxi[igr],GAT::limxi[igr+1]),"p");
    else
        leg.AddEntry(tgrclone0,Form("p_{T} #in [%1.1f,%1.1f)",GAT::limptt[igr],GAT::limptt[igr+1]),"p");
}

void DrawKPlot(TGraphErrors* tgrnum, TGraphErrors* tgrden, TGraphErrors* tgrk, const char* name, const int ik, const int iptt, const int ixi)
{
    // cout<<"Drawing k plot "<<name<<endl;
    gStyle->SetEndErrorSize(20);
    
    tgrnum->SetMarkerStyle(kFullCircle);
    // tgrnum->SetMarkerStyle(34);//cross
    // tgrnum->SetMarkerColor(kBlue);
    // tgrnum->SetLineColor(kBlue);
    tgrnum->SetMarkerColor(TColor::GetColorPalette(220*0/3));
    tgrnum->SetLineColor(TColor::GetColorPalette(220*0/3));
    tgrnum->SetMarkerSize(6.);
    tgrnum->SetLineWidth(2);

    tgrden->SetMarkerStyle(kFullSquare);
    // tgrden->SetMarkerStyle(47);//x
    // tgrden->SetMarkerColor(kRed);
    // tgrden->SetLineColor(kRed);
    tgrden->SetMarkerColor(TColor::GetColorPalette(220*1/3));
    tgrden->SetLineColor(TColor::GetColorPalette(220*1/3));
    tgrden->SetMarkerSize(6.);
    tgrden->SetLineWidth(2);
    
    tgrk->SetMarkerStyle(kFullCircle);
    // tgrk->SetMarkerColor(kBlack);
    // tgrk->SetLineColor(kBlack);
    tgrnum->SetMarkerColor(TColor::GetColorPalette(220*2/3));
    tgrnum->SetLineColor(TColor::GetColorPalette(220*2/3));
    tgrk->SetMarkerSize(6.);
    tgrk->SetLineWidth(2);
    
    TMultiGraph mgr ("mgr",Form(";#LTN_{part}#GT;%s",klab[ik][0].c_str()));
    mgr.Add(tgrnum,"p");
    mgr.Add(tgrden,"p");

    TLegend leg (.675,.75,.95,.9);
    // TLegend leg (.675,.75,.9,.9);
    leg.AddEntry(tgrnum,klab[ik][1].c_str(),"p");
    leg.AddEntry(tgrden,klab[ik][2].c_str(),"p");

    // TCanvas can ("can","",2000,2500);
    TCanvas can ("can","",2000,5000);
    can.cd();
    
    TPad hipad ("hipad","",0,.5,1,1);
    hipad.SetBottomMargin(0.02);
    hipad.SetLeftMargin(0.15);
    hipad.SetRightMargin(0.05);
    // hipad.SetGridx();
    hipad.Draw();
    hipad.cd();
    mgr.Draw("A");
    mgr.GetXaxis()->SetLimits(0.,390.);
    mgr.GetXaxis()->SetLabelSize(0.0);
    mgr.GetYaxis()->SetLabelSize(0.05);
    mgr.GetYaxis()->SetTitleSize(0.05);
    mgr.GetYaxis()->SetTitleOffset(1.4);
    mgr.GetYaxis()->SetMaxDigits(2);
    mgr.SetMaximum(1.25*mgr.GetHistogram()->GetMaximum() - .25*mgr.GetHistogram()->GetMinimum());
    hipad.Update();

    leg.SetTextSize(.05);
    leg.Draw();

    can.cd();
    
    TPad lopad ("lopad","",0,0,1,.5);
    lopad.SetTopMargin(0.01);
    lopad.SetBottomMargin(0.11);
    lopad.SetLeftMargin(0.15);
    lopad.SetRightMargin(0.05);
    // lopad.SetGridx();
    lopad.SetTickx(2);
    lopad.Draw();
    lopad.cd();
    tgrk->SetTitle(Form(";#LTN_{part}#GT;%s",kname[ik].c_str()));
    tgrk->GetXaxis()->SetLimits(0.,390.);
    tgrk->Draw("AP");
    tgrk->GetXaxis()->SetLabelSize(0.05);
    tgrk->GetYaxis()->SetLabelSize(0.05);
    tgrk->GetYaxis()->SetTitleSize(0.05);
    tgrk->GetYaxis()->SetTitleOffset(1.4);
    tgrk->GetXaxis()->SetTitleSize(0.05);
    tgrk->GetXaxis()->SetTitleOffset(1.0);
    lopad.Update();

    txtmom2[iptt]->Draw();
    txtxi[ixi]->Draw();
    
    TLine line (0.,1.,390.,1.);
    line.SetLineStyle(7);
    line.SetLineWidth(1);

    if (ik != 1)
    {
        line.Draw();
    }
    
    can.Update();
    can.SaveAs(name);

    // cout<<"Finished drawing k plot "<<name<<endl;

    gStyle->SetEndErrorSize(errsize);//reset end error size
}


void PrintData(TGraphErrors* tgr[nadj][GAT::ntype+1][GAT::nptt][GAT::nxi],string pref)
{
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < GAT::ntype+1; itype++)      
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;
            // if (itype != GAT::PI0) continue;
            
            ostringstream ossname0;
            ossname0<<"termtxt/"<<pref<<"_adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt01.txt";
            ostringstream ossname1;
            ossname1<<"termtxt/"<<pref<<"_adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt23.txt";
            ofstream ofs[2] {ofstream(ossname0.str()),ofstream(ossname1.str())};

            for (int ipoint = 0; ipoint < tgr[iadj][itype][0][0]->GetN(); ipoint++)
            {
                for (int iptt = 0; iptt<GAT::nptt; iptt++)
                {
                    ofstream* ofsptr = &(ofs[0]);
                    if (iptt>=2) ofsptr = &(ofs[1]);
                    
                    for (int ixi = 0; ixi<GAT::nxi; ixi++)
                    {
                        if ((iptt==0 || iptt==2) && ixi==0)
                        {
                            (*ofsptr)<<fixed<<setprecision(1);
                            (*ofsptr)<<"$"<<(tgr[iadj][itype][iptt][ixi]->GetX()[ipoint])<<"$";
                        }
                        // (*ofsptr)<<setprecision(3)<<scientific;
                        // (*ofsptr)<<" & \\makecell{$\\num{"<<(tgr[iadj][itype][iptt][ixi]->GetY()[ipoint])<<"}$ \\\\ $\\pm \\num{"<<(tgr[iadj][itype][iptt][ixi]->GetEY()[ipoint])<<"}$}";
                        (*ofsptr)<<scientific<<setprecision(2);
                        (*ofsptr)<<" & \\makecell{$\\num{"<<(tgr[iadj][itype][iptt][ixi]->GetY()[ipoint])<<"}$ \\\\ $\\pm \\num{";
                        (*ofsptr)<<setprecision(1);
                        (*ofsptr)<<(tgr[iadj][itype][iptt][ixi]->GetEY()[ipoint])<<"}$}";

                        // if ((iptt==0 || iptt==2) && ixi==0)
                        // {
                        //     (*ofsptr)<<"$"<<setprecision(4)<<(tgr[iadj][itype][iptt][ixi]->GetX()[ipoint])<<"$";
                        // }
                        // (*ofsptr)<<" & \\makecell{$\\num{"<<setprecision(3)<<(tgr[iadj][itype][iptt][ixi]->GetY()[ipoint])<<"}$ \\\\ $\\pm \\num{"<<setprecision(3)<<(tgr[iadj][itype][iptt][ixi]->GetEY()[ipoint])<<"}$}";
                    }
                    if (iptt==1 || iptt==3)
                        (*ofsptr)<<"\\\\"<<"\n"<<"\\hline"<<"\n";
                }
            }
            ofs[0].close();
            ofs[1].close();
        }
    }
}

void DrawMultiGraph(TCanvas& can, TMultiGraph& mgr, bool setlogy = true, bool noex = false)
{
    can.SetLeftMargin(lt);
    can.SetRightMargin(rt);

    if (setlogy)
        can.SetLogy();
    mgr.Draw("A");
    // mgr.GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
                    
    double ymin = numeric_limits<double>::max();
    double ymax = numeric_limits<double>::min();

    for (const auto&& obj : *(mgr.GetListOfGraphs()))
    {
        TGraphErrors* gr = (TGraphErrors*) obj;
        for (int ipt = 0; ipt < gr->GetN(); ipt++)
        {
            ymin = min(ymin,fabs(gr->GetY()[ipt]));
            ymax = max(ymax,fabs(gr->GetY()[ipt]));

            if (noex)
                gr->SetPointError(ipt,0.,gr->GetEY()[ipt]);
        }

    }
    // cout<<"ymin = "<<ymin<<"; ymax = "<<ymax<<endl;

    if (setlogy)
        mgr.GetYaxis()->SetRangeUser(ymin/5.,ymax*30.);
    else
    {
        double diff = ymax-ymin;
        mgr.GetYaxis()->SetRangeUser(ymin-.1*diff,ymax+.2*diff);
    }
    // double max = mgr.GetHistogram()->GetMaximum();
    // mgr.GetHistogram()->SetMaximum(max*50);

    mgr.GetYaxis()->SetTitleSize(.0425);
    mgr.GetYaxis()->SetLabelSize(.04);
    mgr.GetXaxis()->SetTitleSize(.0425);
    mgr.GetXaxis()->SetLabelSize(.04);

    can.Modified();
    can.Update();
}
