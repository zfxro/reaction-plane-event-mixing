
//main drawing functions are SmallDrawer and BottomDrawer
//leaving old drawing -- can cause confusion when reading code though

const int ncent = 2;
const int limcent[ncent+1] = {0,20,60};


enum RPO {IN,OUT,ALL};

const int nadj = 2;

const float PI = (float) TMath::Pi();

struct errorinfo_t {
    Int_t iadj, itype, iptt, ixi, icent, idpsi;
    Double_t ntrg;
    Double_t flate[2] = {};//flat term error
    Double_t cose[2] = {};//cosine term error
};

const EColor colors[4] = {kRed,kBlue,kGreen,kViolet};
const short styles[4] = {1,9,7,3};

const double yoff = 1.4;

namespace GAT {
    #include "gamAnaTree.h"
}


int ReadFile(const char* inname, TFile*& file);
double GetNTrg(TH1D* h1cent, int icent);
vector<double> GetNevt(TFile* file);
void OpenTheDrawer(TPad* pad, TH1* h1, TF1* f1, TF1* f2 = nullptr);
void SmallDrawer(TH1F* h1, TF1* f1, const char* fname, const char* label, const char* filename, double err);
void BottomDrawer(TH1F* h1, const char* label, const char* filename, double fe = -999., double ce = -999.);
void BigDrawer(TH1F* h1, const list<pair<TF1*,string>>& kdmap, const char* filename);
int GetCentBin(int ic);

void AddToHisto(TH1F* h1, TF1* f1, double wt = 1.);

TFile* output;

void subtractKWeightedBG(
    const char* namegt = "gamtree_13466_acceptcorrect.root",
    const char* namek = "kComponents.root",
    const char* nameout = "weighted_bg_subtracted.root",
    bool usekerr = true,
    bool usetermerr = true,
    bool addlinear = false,
    bool usekmc = true
    )
{
    ofstream ofsntrgall("ntrgall.txt");
    ofstream ofsntrgin("ntrgin.txt");
    ofstream ofsntrgout("ntrgout.txt");
    
    //style optimized for smalldrawer canvas
    //comment out most if interested in other plots

    gStyle->SetLineScalePS(1);
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0.);
    gStyle->SetEndErrorSize(20);
    gStyle->SetMarkerStyle(kFullCircle);
    gStyle->SetMarkerSize(6);
    gStyle->SetHistLineWidth(2);
    gStyle->SetPalette(kBird);
    gStyle->SetMarkerColor(TColor::GetColorPalette(254*0/5));
    gStyle->SetHistLineColor(TColor::GetColorPalette(254*0/5));
    //if change these settings please change in calcKComponents.C as well
    gStyle->SetTitleOffset(yoff,"Y");
    gStyle->SetLabelSize(.045,"Y");
    gStyle->SetTitleSize(.05,"Y");
    gStyle->SetLabelSize(.045,"X");
    gStyle->SetTitleSize(.05,"X");
    gROOT->ForceStyle();


    TFile* ingt = nullptr;
    if (ReadFile(namegt,ingt)<0) return;

    TFile* ink = nullptr;
    if (ReadFile(namek,ink)<0) return;

    TGraphErrors* grsigma = nullptr;
    ink->GetObject("npartsigma",grsigma);

    output = new TFile(nameout,"RECREATE");
    output->WriteTObject(grsigma);

    errorinfo_t info;
    TTree tree("Terrorinfo","Tree with information on error of subtracted background");
    tree.Branch("info",&info,"iadj/I:itype:iptt:ixi:icent:idpsi:ntrg/D:flate[2]:cose[2]");

    vector<double> vnevt = GetNevt(ink);

    TCanvas* canio = new TCanvas("canio","",1200,1350);
    canio->Divide(2,3);

    TCanvas* canall = new TCanvas("canall","",600,1350);
    canall->Divide(1,4);

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            //saving too many histograms. skip these

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                //suffix of trigger related objects (histos, tgraphs, etc.)
                ostringstream sufg;
                sufg<<"adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt;

                //--------------- get ka correction ----------
                TGraphErrors* grag = nullptr;
                ink->GetObject(Form("ag_Np_%s",sufg.str().c_str()),grag);

                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    //suffix of hadron related objects (histos, tgraphs, etc.)
                    ostringstream sufh;
                    sufh<<sufg.str()<<"_xi"<<ixi;

                    // ----- histograms that encompass multiple centrality bins -----

                    TH1F* h1sumall[ncent+1] = {};
                    TH1F* h1sumin[ncent+1] = {};
                    TH1F* h1sumout[ncent+1] = {};

                    double f0sum[ncent+1] = {};
                    double kf0sum[ncent+1] = {};
                    double kf0errsum[ncent+1] = {};

                    double c0sum[ncent+1] = {};
                    double kc0sum[ncent+1] = {};
                    double kc0errsum[ncent+1] = {};

                    double f1sum[ncent+1] = {};
                    double kf1sum[ncent+1] = {};
                    double kf1errsum[ncent+1] = {};

                    double c1sum[ncent+1] = {};
                    double kc1sum[ncent+1] = {};
                    double kc1errsum[ncent+1] = {};


                    double ntrgall[ncent+1] = {};
                    double ntrgin[ncent+1] = {};
                    double ntrgout[ncent+1] = {};


                    //----- get remaining contributions and k corrections -----

                    //--- terms... to be deprecated ---
                    TGraphErrors* grah = nullptr;
                    ink->GetObject(Form("ah_Np_%s",sufh.str().c_str()),grah);

                    TGraphErrors* grbgs = nullptr;
                    ink->GetObject(Form("bgs_Np_%s",sufh.str().c_str()),grbgs);

                    TGraphErrors* grbhs = nullptr;
                    ink->GetObject(Form("bhs_Np_%s",sufh.str().c_str()),grbhs);

                    //--- contributions. will be using these instead of terms ---

                    TGraphErrors* gragah = nullptr;
                    ink->GetObject(Form("agah_Np_%s",sufh.str().c_str()),gragah);

                    TGraphErrors* grahbgs = nullptr;
                    ink->GetObject(Form("ahbgs_Np_%s",sufh.str().c_str()),grahbgs);

                    TGraphErrors* gragbhs = nullptr;
                    ink->GetObject(Form("agbhs_Np_%s",sufh.str().c_str()),gragbhs);

                    TGraphErrors* grbgsbhs = nullptr;
                    ink->GetObject(Form("bgsbhs_Np_%s",sufh.str().c_str()),grbgsbhs);

                    //--- k corrections ---
                    TGraphErrors* grka = nullptr;
                    TGraphErrors* grkb = nullptr;
                    TGraphErrors* grkgh = nullptr;
                    TGraphErrors* grkhg = nullptr;

                    if (usekmc)
                    {
                        ink->GetObject(Form("kanew_Np_%s",sufh.str().c_str()),grka);
                        ink->GetObject(Form("kbnew_Np_%s",sufh.str().c_str()),grkb);
                        ink->GetObject(Form("kghnew_Np_%s",sufh.str().c_str()),grkgh);
                        ink->GetObject(Form("khgnew_Np_%s",sufh.str().c_str()),grkhg);
                    }
                    else
                    {
                        ink->GetObject(Form("ka_Np_%s",sufh.str().c_str()),grka);
                        ink->GetObject(Form("kb_Np_%s",sufh.str().c_str()),grkb);
                        ink->GetObject(Form("kgh_Np_%s",sufh.str().c_str()),grkgh);
                        ink->GetObject(Form("khg_Np_%s",sufh.str().c_str()),grkhg);
                    }
                    
                    //------------------------------------------------------------
                    int curcentbin = -999;

                    for (int icent = 0; icent < GAT::ncent; icent++)
                    {

                        int centbin = GetCentBin(icent);
                        cout<<"curcentbin "<<curcentbin<<". newcentbin "<<centbin<<endl;

                        //nevt used to rescale nevt normalized contributions
                        double nevt = vnevt[icent];

                        //reaction plane resolution
                        double sigma = grsigma->GetY()[icent];

                        //---------- get contributions for this centrality ----------
                        double ag = grag->GetY()[icent];

                        double bgs = grbgs->GetY()[icent];

                        double ah = grah->GetY()[icent];

                        double bhs = grbhs->GetY()[icent];

                        double bgsbhs = grbgsbhs->GetY()[icent];

                        //---------- get background terms for this centrality ----------

                        //switching to using terms not components                       
                        double agah = gragah->GetY()[icent];
                        double ahbgs = grahbgs->GetY()[icent];
                        double agbhs = gragbhs->GetY()[icent];

                        //errors
                        double agahe = gragah->GetEY()[icent];
                        double ahbgse = grahbgs->GetEY()[icent];
                        double agbhse = gragbhs->GetEY()[icent];
                        double bgsbhse = grbgsbhs->GetEY()[icent];
                        
                        //---------- get corrections for this centrality ----------
                        double ka = grka->GetY()[icent];
                        double kae = grka->GetEY()[icent];

                        double kb = grkb->GetY()[icent];
                        double kbe = grkb->GetEY()[icent];

                        double kgh = grkgh->GetY()[icent];
                        double kghe = grkgh->GetEY()[icent];

                        double khg = grkhg->GetY()[icent];
                        double khge = grkhg->GetEY()[icent];


                        //-------------------- all-plane case --------------------
                        TH1D* h1trgcent = nullptr;
                        ingt->GetObject(Form("cent_adj%i_%s_ptt%i_dpsi%i",iadj,GAT::spec[itype],iptt,ALL),h1trgcent);
                        double ntrg = GetNTrg(h1trgcent,icent);
                        if (centbin != curcentbin)
                        {
                            if (ntrgall[centbin] != 0.)
                            {
                                cerr<<"Error: unexpected state"<<endl;
                                return;
                            }
                        }
                        ntrgall[centbin] += ntrg;
                        ntrgall[ncent] += ntrg;

                        //get signal
                        ostringstream oss;
                        oss<<"dphi";
                        if (iadj)
                            oss<<"_adj";
                        oss<<"_"<<GAT::spec[itype]<<"_mix0_ptt"<<iptt<<"_cent"<<icent<<"_dpsi"<<ALL<<"_xi"<<ixi;

                        TH1F* h1dphiall = nullptr;
                        ingt->GetObject(oss.str().c_str(),h1dphiall);
                        h1dphiall->SetTitle(";#Delta#phi/#pi;1/N_{t}  dN_{h}/d#Delta#phi");
                        h1dphiall->Scale(1./(h1dphiall->GetBinWidth(1)*PI));

                        TH1F* h1dphiall_nok = (TH1F*) h1dphiall->Clone(Form("%s_nok",oss.str().c_str()));

                        if (centbin != curcentbin)
                        {
                            ostringstream osssum;
                            osssum<<"dphi";
                            if (iadj)
                                osssum<<"_adj";
                            osssum<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_sumcent"<<centbin<<"_dpsi"<<ALL<<"_xi"<<ixi;

                            cout<<"Cloning "<<h1dphiall->GetName()<<" to "<<osssum.str()<<endl;
                            h1sumall[centbin] = (TH1F*) h1dphiall->Clone(osssum.str().c_str());

                            if (centbin == 0)
                            {
                                ostringstream osssum;
                                osssum<<"dphi";
                                if (iadj)
                                    osssum<<"_adj";
                                osssum<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_sumcent"<<ncent<<"_dpsi"<<ALL<<"_xi"<<ixi;

                                cout<<"Cloning "<<h1dphiall->GetName()<<" to "<<osssum.str()<<endl;
                                h1sumall[ncent] = (TH1F*) h1dphiall->Clone(osssum.str().c_str());
                            }
                            else
                            {
                                h1sumall[ncent]->Add(h1dphiall);
                            }
                        }
                        else
                        {
                            h1sumall[centbin]->Add(h1dphiall);
                            h1sumall[ncent]->Add(h1dphiall);
                        }

                        //normalize
                        h1dphiall->Scale(1./ntrg);
                        h1dphiall_nok->Scale(1./ntrg);

                        //---------- initial flat contribution ----------

                        double f0scale = 4.*PI*agah;
                        f0scale*=nevt;
                        f0sum[centbin] += f0scale;
                        f0sum[ncent] += f0scale;
                        TF1* f0 = new TF1("f0","[0]",0.,1.);
                        f0->SetParameter(0,f0scale/ntrg);

                        double kf0scale = ka*f0scale;
                        kf0sum[centbin] += kf0scale;
                        kf0sum[ncent] += kf0scale;
                        TF1* kf0 = new TF1("kf0","[0]",0.,1.);
                        kf0->SetParameter(0,kf0scale/ntrg);

                        //squared (if not addlinear) error for k corrected bg term
                        double kf0err = 0.;                     
                        if (usekerr)
                        {
                            if (addlinear)
                                kf0err += fabs(kae*kf0scale/ka);
                            else
                                kf0err += pow(kae*kf0scale/ka,2.);
                        }
                        if (usetermerr)
                        {
                            if (addlinear)
                                kf0err += fabs(agahe*kf0scale/agah);
                            else
                                kf0err += pow(agahe*kf0scale/agah,2.);
                        }

                        kf0errsum[centbin] += kf0err;
                        kf0errsum[ncent] += kf0err;

                        TPad* pad = (TPad*) canall->cd(1);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiall,f0,kf0);

                        //subtract the contribution
                        AddToHisto(h1dphiall,kf0,-1.0);
                        AddToHisto(h1dphiall_nok,f0,-1.0);

                        //---------- initial cosine contribution ----------

                        double c0scale = 2.*PI*bgsbhs;
                        c0scale *= nevt;
                        c0scale /= (sigma*sigma);
                        c0sum[centbin] += c0scale;
                        c0sum[ncent] += c0scale;
                        TF1* c0 = new TF1("c0","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        c0->SetParameter(0,c0scale/ntrg);

                        double kc0scale = 2.*PI*bgsbhs;
                        kc0scale *= nevt;
                        kc0scale *= kb;
                        kc0sum[centbin] += kc0scale;
                        kc0sum[ncent] += kc0scale;
                        TF1* kc0 = new TF1("kc0","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        kc0->SetParameter(0,kc0scale/ntrg);

                        double kc0err = 0.;
                        if (usekerr)
                        {
                            if (addlinear)
                                kc0err += fabs(kbe*kc0scale/kb);
                            else
                                kc0err += pow(kbe*kc0scale/kb,2.);
                        }
                        if (usetermerr)
                        {
                            if (addlinear)
                                kc0err += fabs(bgsbhse*kc0scale/bgsbhs);
                            else
                                kc0err += pow(bgsbhse*kc0scale/bgsbhs,2.);
                        }
                        
                        kc0errsum[centbin] += kc0err;
                        kc0errsum[ncent] += kc0err;

                        pad = (TPad*) canall->cd(2);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiall,c0,kc0);

                        AddToHisto(h1dphiall,kc0,-1.0);
                        AddToHisto(h1dphiall_nok,c0,-1.0);

                        //write result to file
                        output->WriteTObject(h1dphiall);

                        //---------- draw the result ----------

                        TLine* line1 = new TLine(0.,0.,1.,0.);
                        line1->SetLineStyle(kDashed);

                        pad = (TPad*) canall->cd(3);
                        pad->Clear();
                        h1dphiall->Draw("E1");
                        line1->Draw("same");

                        pad = (TPad*) canall->cd(4);
                        pad->Clear();
                        h1dphiall_nok->Draw("E1");
                        {
                            double min = h1dphiall_nok->GetMinimum();
                            if (min > 0.)
                            {
                                double max = h1dphiall_nok->GetMaximum();
                                h1dphiall_nok->GetYaxis()->SetRangeUser(0.,max);
                                pad->Update();
                            }
                        }
                        line1->Draw("same");

                        canall->Update();
                        canall->SaveAs(Form("subtractbg/samepad/cent/%s.pdf",oss.str().c_str()));

                        //-------------------- in-plane case --------------------
                        h1trgcent = nullptr;
                        ingt->GetObject(Form("cent_adj%i_%s_ptt%i_dpsi%i",iadj,GAT::spec[itype],iptt,IN),h1trgcent);
                        ntrg = GetNTrg(h1trgcent,icent);
                        ntrgin[centbin] += ntrg;
                        ntrgin[ncent] += ntrg;

                        //get signal
                        oss.str("");
                        oss<<"dphi";
                        if (iadj)
                            oss<<"_adj";
                        oss<<"_"<<GAT::spec[itype]<<"_mix0_ptt"<<iptt<<"_cent"<<icent<<"_dpsi"<<IN<<"_xi"<<ixi;

                        TH1F* h1dphiin = nullptr;
                        ingt->GetObject(oss.str().c_str(),h1dphiin);
                        h1dphiin->SetTitle(";#Delta#phi/#pi;1/N_{t}  dN_{h}/d#Delta#phi");
                        h1dphiin->Scale(1./(h1dphiin->GetBinWidth(1)*PI));

                        TH1F* h1dphiin_nok = (TH1F*) h1dphiin->Clone(Form("%s_nok",oss.str().c_str()));

                        if (centbin != curcentbin)
                        {
                            ostringstream osssum;
                            osssum<<"dphi";
                            if (iadj)
                                osssum<<"_adj";
                            osssum<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_sumcent"<<centbin<<"_dpsi"<<IN<<"_xi"<<ixi;

                            cout<<"Cloning "<<h1dphiin->GetName()<<" to "<<osssum.str()<<endl;

                            h1sumin[centbin] = (TH1F*) h1dphiin->Clone(osssum.str().c_str());

                            if (centbin == 0)
                            {
                                ostringstream osssum;
                                osssum<<"dphi";
                                if (iadj)
                                    osssum<<"_adj";
                                osssum<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_sumcent"<<ncent<<"_dpsi"<<IN<<"_xi"<<ixi;

                                cout<<"Cloning "<<h1dphiin->GetName()<<" to "<<osssum.str()<<endl;

                                h1sumin[ncent] = (TH1F*) h1dphiin->Clone(osssum.str().c_str());
                            }
                            else
                            {
                                h1sumin[ncent]->Add(h1dphiin);
                            }
                        }
                        else
                        {
                            h1sumin[centbin]->Add(h1dphiin);
                            h1sumin[ncent]->Add(h1dphiin);
                        }

                        //normalize
                        h1dphiin->Scale(1./ntrg);
                        h1dphiin_nok->Scale(1./ntrg);

                        //---------- initial flat contribution ----------

                        f0scale /= 2.;
                        f0->SetParameter(0,f0scale/ntrg);

                        kf0scale /= 2.;
                        kf0->SetParameter(0,kf0scale/ntrg);

                        pad = (TPad*) canio->cd(1);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiin,f0,kf0);

                        //subtract the contribution
                        AddToHisto(h1dphiin,kf0,-1.0);
                        AddToHisto(h1dphiin_nok,f0,-1.0);

                        //---------- additional flat contribution ----------

                        double f1scale = 4.*ahbgs;
                        f1scale *= nevt;
                        f1sum[centbin] += f1scale;
                        f1sum[ncent] += f1scale;
                        TF1* f1 = new TF1("f1","[0]",0.,1.);
                        f1->SetParameter(0,f1scale/ntrg);

                        double kf1scale = khg*f1scale;
                        kf1sum[centbin] += kf1scale;
                        kf1sum[ncent] += kf1scale;
                        TF1* kf1 = new TF1("kf1","[0]",0.,1.);
                        kf1->SetParameter(0,kf1scale/ntrg);

                        double kf1err = 0.;
                        if (usekerr)
                        {
                            if (addlinear)
                                kf1err += fabs(khge*kf1scale/khg);
                            else
                                kf1err += pow(khge*kf1scale/khg,2.);
                        }
                        if (usetermerr)
                        {
                            if (addlinear)
                                kf1err += fabs(ahbgse*kf1scale/ahbgs);
                            else
                                kf1err += pow(ahbgse*kf1scale/ahbgs,2.);
                        }
                        
                        kf1errsum[centbin] += kf1err;
                        kf1errsum[ncent] += kf1err;
                        
                        pad = (TPad*) canio->cd(2);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiin,f1,kf1);

                        //subtract the contribution
                        AddToHisto(h1dphiin,kf1,-1.0);
                        AddToHisto(h1dphiin_nok,f1,-1.0);

                        //---------- additional cosine contribution ----------

                        double c1scale = 4.*agbhs;
                        c1scale *= nevt;
                        c1sum[centbin] += c1scale;
                        c1sum[ncent] += c1scale;
                        TF1* c1 = new TF1("c1","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        c1->SetParameter(0,c1scale/ntrg);

                        double kc1scale = kgh*c1scale;
                        kc1sum[centbin] += kc1scale;
                        kc1sum[ncent] += kc1scale;
                        TF1* kc1 = new TF1("kc1","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        kc1->SetParameter(0,kc1scale/ntrg);

                        double kc1err = 0.;
                        if (usekerr)
                        {
                            if (addlinear)
                                kc1err += fabs(kghe*kc1scale/kgh);
                            else
                                kc1err += pow(kghe*kc1scale/kgh,2.);
                        }
                        if (usetermerr)
                        {
                            if (addlinear)
                                kc1err += fabs(agbhse*kc1scale/agbhs);
                            else
                                kc1err += pow(agbhse*kc1scale/agbhs,2.);
                        }
                        
                        kc1errsum[centbin] += kc1err;
                        kc1errsum[ncent] += kc1err;

                        pad = (TPad*) canio->cd(3);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiin,c1,kc1);

                        AddToHisto(h1dphiin,kc1,-1.0);
                        AddToHisto(h1dphiin_nok,c1,-1.0);

                        //---------- initial cosine contribution ----------

                        c0scale /= 2.0;
                        c0->SetParameter(0,c0scale/ntrg);

                        kc0scale /= 2.0;
                        kc0->SetParameter(0,kc0scale/ntrg);

                        pad = (TPad*) canio->cd(4);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiin,c0,kc0);

                        AddToHisto(h1dphiin,kc0,-1.0);
                        AddToHisto(h1dphiin_nok,c0,-1.0);

                        //write result to file
                        output->WriteTObject(h1dphiin);

                        //---------- draw the result ----------

                        pad = (TPad*) canio->cd(5);
                        pad->Clear();
                        h1dphiin->Draw("E1");
                        line1->Draw("same");

                        pad = (TPad*) canio->cd(6);
                        pad->Clear();
                        h1dphiin_nok->Draw("E1");
                        {
                            double min = h1dphiin_nok->GetMinimum();
                            if (min > 0.)
                            {
                                double max = h1dphiin_nok->GetMaximum();
                                h1dphiin_nok->GetYaxis()->SetRangeUser(0.,max);
                                pad->Update();
                            }
                        }
                        line1->Draw("same");

                        canio->Update();
                        canio->SaveAs(Form("subtractbg/samepad/cent/%s.pdf",oss.str().c_str()));

                        //-------------------- out-of-plane case --------------------
                        h1trgcent = nullptr;
                        ingt->GetObject(Form("cent_adj%i_%s_ptt%i_dpsi%i",iadj,GAT::spec[itype],iptt,OUT),h1trgcent);
                        ntrg = GetNTrg(h1trgcent,icent);
                        ntrgout[centbin] += ntrg;
                        ntrgout[ncent] += ntrg;

                        //get signal
                        oss.str("");
                        oss<<"dphi";
                        if (iadj)
                            oss<<"_adj";
                        oss<<"_"<<GAT::spec[itype]<<"_mix0_ptt"<<iptt<<"_cent"<<icent<<"_dpsi"<<OUT<<"_xi"<<ixi;

                        TH1F* h1dphiout = nullptr;
                        ingt->GetObject(oss.str().c_str(),h1dphiout);
                        h1dphiout->SetTitle(";#Delta#phi/#pi;1/N_{t}  dN_{h}/d#Delta#phi");
                        h1dphiout->Scale(1./(h1dphiout->GetBinWidth(1)*PI));

                        TH1F* h1dphiout_nok = (TH1F*) h1dphiout->Clone(Form("%s_nok",oss.str().c_str()));

                        if (centbin != curcentbin)
                        {
                            ostringstream osssum;
                            osssum<<"dphi";
                            if (iadj)
                                osssum<<"_adj";
                            osssum<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_sumcent"<<centbin<<"_dpsi"<<OUT<<"_xi"<<ixi;

                            cout<<"Cloning "<<h1dphiout->GetName()<<" to "<<osssum.str()<<endl;
                            h1sumout[centbin] = (TH1F*) h1dphiout->Clone(osssum.str().c_str());

                            if (centbin == 0)
                            {
                                ostringstream osssum;
                                osssum<<"dphi";
                                if (iadj)
                                    osssum<<"_adj";
                                osssum<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_sumcent"<<ncent<<"_dpsi"<<OUT<<"_xi"<<ixi;

                                cout<<"Cloning "<<h1dphiout->GetName()<<" to "<<osssum.str()<<endl;
                                h1sumout[ncent] = (TH1F*) h1dphiout->Clone(osssum.str().c_str());
                            }
                            else
                            {
                                h1sumout[ncent]->Add(h1dphiout);
                            }
                        }
                        else
                        {
                            h1sumout[centbin]->Add(h1dphiout);
                            h1sumout[ncent]->Add(h1dphiout);
                        }

                        //normalize
                        h1dphiout->Scale(1./ntrg);
                        h1dphiout_nok->Scale(1./ntrg);

                        //---------- initial flat contribution ----------

                        f0->SetParameter(0,f0scale/ntrg);
                        kf0->SetParameter(0,kf0scale/ntrg);

                        pad = (TPad*) canio->cd(1);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiout,f0,kf0);

                        //subtract the contribution
                        AddToHisto(h1dphiout,kf0,-1.0);
                        AddToHisto(h1dphiout_nok,f0,-1.0);

                        //---------- additional flat contribution ----------

                        f1->SetParameter(0,-f1scale/ntrg);
                        kf1->SetParameter(0,-kf1scale/ntrg);

                        pad = (TPad*) canio->cd(2);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiout,f1,kf1);

                        //subtract the contribution
                        AddToHisto(h1dphiout,kf1,-1.0);
                        AddToHisto(h1dphiout_nok,f1,-1.0);

                        //---------- additional cosine contribution ----------

                        c1->SetParameter(0,-c1scale/ntrg);
                        kc1->SetParameter(0,-kc1scale/ntrg);

                        pad = (TPad*) canio->cd(3);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiout,c1,kc1);

                        AddToHisto(h1dphiout,kc1,-1.0);
                        AddToHisto(h1dphiout_nok,c1,-1.0);

                        //---------- initial cosine contribution ----------

                        c0->SetParameter(0,c0scale/ntrg);
                        kc0->SetParameter(0,kc0scale/ntrg);

                        pad = (TPad*) canio->cd(4);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiout,c0,kc0);

                        AddToHisto(h1dphiout,kc0,-1.0);
                        AddToHisto(h1dphiout_nok,c0,-1.0);

                        //write result to file
                        output->WriteTObject(h1dphiout);

                        //---------- draw the result ----------

                        pad = (TPad*) canio->cd(5);
                        pad->Clear();
                        h1dphiout->Draw("E1");
                        line1->Draw("same");

                        pad = (TPad*) canio->cd(6);
                        pad->Clear();
                        h1dphiout_nok->Draw("E1");
                        {
                            double min = h1dphiout_nok->GetMinimum();
                            if (min > 0.)
                            {
                                double max = h1dphiout_nok->GetMaximum();
                                h1dphiout_nok->GetYaxis()->SetRangeUser(0.,max);
                                pad->Update();
                            }
                        }
                        line1->Draw("same");

                        canio->Update();
                        canio->SaveAs(Form("subtractbg/samepad/cent/%s.pdf",oss.str().c_str()));

                        curcentbin = centbin;
                    }
                    
                    //---------- Draw the centrality merged results ----------
                    for (int icent = 0; icent <= ncent; icent++)
                    {
                        //-------------------- all-plane case --------------------

                        double ntrg = ntrgall[icent];
                        ofsntrgall<<sufg.str()<<" sumcent"<<icent<<" "<<ntrg<<endl;


                        //could move all these to appropriate loops so set only once
                        info.iadj = iadj;
                        info.itype = itype;
                        info.iptt = iptt;
                        info.ixi = ixi;
                        info.icent = icent;
                        
                        info.idpsi = ALL;
                        info.ntrg = ntrg;

                        h1sumall[icent]->Scale(1./ntrg);

                        TH1F* h1sumall_nok = (TH1F*) h1sumall[icent]->Clone(Form("%s_nok",h1sumall[icent]->GetName()));

                        //---------- initial flat contribution ----------

                        TF1* f0 = new TF1("f0","[0]",0.,1.);
                        f0->SetParameter(0,f0sum[icent]/ntrg);

                        TF1* kf0 = new TF1("kf0","[0]",0.,1.);
                        kf0->SetParameter(0,kf0sum[icent]/ntrg);

                        TPad* pad = (TPad*) canall->cd(1);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumall[icent],f0,kf0);
                        

                        {
                            double kf0err = kf0errsum[icent];
                            if (!addlinear)
                                kf0err = sqrt(kf0err);
                            kf0err/=ntrg;

                            info.flate[0] = kf0err;
                            info.flate[1] = 0.;
                        
                            SmallDrawer(h1sumall[icent],kf0,"(1/n_{t}) 4#pi#kappa_{a}#LTa_{t}#GT#LTa_{h}#GT","(a)",Form("subtractbg/%s_a.pdf",h1sumall[icent]->GetName()),kf0err);
                        }

                        //subtract the contribution
                        AddToHisto(h1sumall[icent],kf0,-1.0);
                        AddToHisto(h1sumall_nok,f0,-1.0);

                        //---------- initial cosine contribution ----------

                        TF1* c0 = new TF1("c0","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        c0->SetParameter(0,c0sum[icent]/ntrg);

                        TF1* kc0 = new TF1("kc0","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        kc0->SetParameter(0,kc0sum[icent]/ntrg);

                        pad = (TPad*) canall->cd(2);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumall[icent],c0,kc0);


                        {
                            double kc0err = kc0errsum[icent];
                            if (!addlinear)
                                kc0err = sqrt(kc0err);
                            kc0err/=ntrg;

                            info.cose[0] = kc0err;
                            info.cose[1] = 0.;
                            tree.Fill();
                        
                            SmallDrawer(h1sumall[icent],kc0,"(1/n_{t}) 2#pi#kappa_{b}#LTb_{t}#sigma#GT#LTb_{h}#sigma#GT cos2#Delta#phi","(b)",Form("subtractbg/%s_b.pdf",h1sumall[icent]->GetName()),kc0err);
                        }
                        
                        AddToHisto(h1sumall[icent],kc0,-1.0);
                        AddToHisto(h1sumall_nok,c0,-1.0);

                        //write result to file
                        output->WriteTObject(h1sumall[icent]);

                        //---------- draw the result ----------

                        TLine* line1 = new TLine(0.,0.,1.,0.);
                        line1->SetLineStyle(kDashed);

                        pad = (TPad*) canall->cd(3);
                        pad->Clear();
                        h1sumall[icent]->Draw("E1");
                        line1->Draw("same");

                        pad = (TPad*) canall->cd(4);
                        pad->Clear();
                        h1sumall_nok->Draw("E1");
                        {
                            double min = h1sumall_nok->GetMinimum();
                            if (min > 0.)
                            {
                                double max = h1sumall_nok->GetMaximum();
                                h1sumall_nok->GetYaxis()->SetRangeUser(0.,max);
                                pad->Update();
                            }
                        }
                        line1->Draw("same");

                        canall->Update();
                        canall->SaveAs(Form("subtractbg/samepad/%s.pdf",h1sumall[icent]->GetName()));


                        double fe_all = kf0errsum[icent];
                        double ce_all = kc0errsum[icent];

                        if (!addlinear)
                        {
                            fe_all = sqrt(fe_all);
                            ce_all = sqrt(ce_all);
                        }

                        fe_all/=ntrg;
                        ce_all/=ntrg;

                        BottomDrawer(h1sumall[icent],"(c)",Form("subtractbg/%s_c.pdf",h1sumall[icent]->GetName()),fe_all,ce_all);

                        TF1* kdf0 = new TF1("kdf0","[0]",0.,1.);
                        kdf0->SetParameter(0,(kf0sum[icent]-f0sum[icent])/ntrg);

                        TF1* kdc0 = new TF1("kdc0","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        kdc0->SetParameter(0,(kc0sum[icent]-c0sum[icent])/ntrg);


                        list<pair<TF1*,string>> kdmapall {
                            {kdf0,"(#kappa_{a}-1) 4#pi#LTa_{t}#GT#LTa_{h}#GT/n_{t}"},
                            {kdc0,"(#kappa_{b}-1) 2#pi#LTb_{t}#sigma#GT#LTb_{h}#sigma#GT cos(2#Delta#phi)/n_{t}"}
                        };
                                

                        BigDrawer(h1sumall[icent],kdmapall,Form("subtractbg/%s.pdf",h1sumall[icent]->GetName()));

                        //-------------------- in-plane case --------------------

                        ntrg = ntrgin[icent];
                        ofsntrgin<<sufg.str()<<" sumcent"<<icent<<" "<<ntrg<<endl;

                        info.idpsi = IN;
                        info.ntrg = ntrg;

                        h1sumin[icent]->Scale(1./ntrg);

                        TH1F* h1sumin_nok = (TH1F*) h1sumin[icent]->Clone(Form("%s_nok",h1sumin[icent]->GetName()));

                        //---------- initial flat contribution ----------

                        f0sum[icent] /= 2.;
                        f0->SetParameter(0,f0sum[icent]/ntrg);

                        kf0sum[icent] /= 2.;
                        kf0->SetParameter(0,kf0sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(1);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumin[icent],f0,kf0);
                        
                        {
                            double kf0err = kf0errsum[icent];
                            if (!addlinear)
                                kf0err = sqrt(kf0err);
                            kf0err/=(2.*ntrg);

                            info.flate[0] = kf0err;

                            SmallDrawer(h1sumin[icent],kf0,"(1/n_{t}) 2#pi#kappa_{a}#LTa_{t}#GT#LTa_{h}#GT","(a)",Form("subtractbg/%s_a.pdf",h1sumin[icent]->GetName()),kf0err);
                        }

                        //subtract the contribution
                        AddToHisto(h1sumin[icent],kf0,-1.0);
                        AddToHisto(h1sumin_nok,f0,-1.0);

                        //---------- additional flat contribution ----------

                        TF1* f1 = new TF1("f1","[0]",0.,1.);
                        f1->SetParameter(0,f1sum[icent]/ntrg);

                        TF1* kf1 = new TF1("kf1","[0]",0.,1.);
                        kf1->SetParameter(0,kf1sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(2);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumin[icent],f1,kf1);


                        {
                            double kf1err = kf1errsum[icent];
                            if (!addlinear)
                                kf1err = sqrt(kf1err);
                            kf1err/=ntrg;

                            info.flate[1] = kf1err;
                        
                            SmallDrawer(h1sumin[icent],kf1,"(1/n_{t}) 4#kappa_{ht}#LTa_{h}#GT#LTb_{t}#sigma#GT","(b)",Form("subtractbg/%s_b.pdf",h1sumin[icent]->GetName()),kf1err);
                        }

                        //subtract the contribution
                        AddToHisto(h1sumin[icent],kf1,-1.0);
                        AddToHisto(h1sumin_nok,f1,-1.0);

                        //---------- additional cosine contribution ----------

                        TF1* c1 = new TF1("c1","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        c1->SetParameter(0,c1sum[icent]/ntrg);

                        TF1* kc1 = new TF1("kc1","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        kc1->SetParameter(0,kc1sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(3);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumin[icent],c1,kc1);


                        {
                            double kc1err = kc1errsum[icent];
                            if (!addlinear)
                                kc1err = sqrt(kc1err);
                            kc1err/=ntrg;

                            info.cose[1] = kc1err;
                        
                            SmallDrawer(h1sumin[icent],kc1,"(1/n_{t}) 4#kappa_{th}#LTa_{t}#GT#LTb_{h}#sigma#GT cos2#Delta#phi","(c)",Form("subtractbg/%s_c.pdf",h1sumin[icent]->GetName()),kc1err);
                        }

                        AddToHisto(h1sumin[icent],kc1,-1.0);
                        AddToHisto(h1sumin_nok,c1,-1.0);

                        //---------- initial cosine contribution ----------

                        c0sum[icent] /= 2.0;
                        c0->SetParameter(0,c0sum[icent]/ntrg);

                        kc0sum[icent] /= 2.0;
                        kc0->SetParameter(0,kc0sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(4);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumin[icent],c0,kc0);


                        {
                            double kc0err = kc0errsum[icent];
                            if (!addlinear)
                                kc0err = sqrt(kc0err);
                            kc0err/=(2.*ntrg);

                            info.cose[0] = kc0err;
                            tree.Fill();

                            SmallDrawer(h1sumin[icent],kc0,"(1/n_{t}) #pi#kappa_{b}#LTb_{t}#sigma#GT#LTb_{h}#sigma#GT cos2#Delta#phi","(d)",Form("subtractbg/%s_d.pdf",h1sumin[icent]->GetName()),kc0err);
                        }                       

                        AddToHisto(h1sumin[icent],kc0,-1.0);
                        AddToHisto(h1sumin_nok,c0,-1.0);

                        //write result to file
                        output->WriteTObject(h1sumin[icent]);

                        //---------- draw the result ----------

                        pad = (TPad*) canio->cd(5);
                        pad->Clear();
                        h1sumin[icent]->Draw("E1");
                        line1->Draw("same");

                        pad = (TPad*) canio->cd(6);
                        pad->Clear();
                        h1sumin_nok->Draw("E1");
                        {
                            double min = h1sumin_nok->GetMinimum();
                            if (min > 0.)
                            {
                                double max = h1sumin_nok->GetMaximum();
                                h1sumin_nok->GetYaxis()->SetRangeUser(0.,max);
                                pad->Update();
                            }
                        }
                        line1->Draw("same");

                        canio->Update();
                        canio->SaveAs(Form("subtractbg/samepad/%s.pdf",h1sumin[icent]->GetName()));

                        //4 is from the 2 (now squared since taking sqrt) that divides the initial contributions which are originally scaled for all-plane
                        double fe_in = (kf0errsum[icent]/2 + kf1errsum[icent])/ntrg;
                        if (!addlinear)
                            fe_in = sqrt(kf0errsum[icent]/4 + kf1errsum[icent])/ntrg;
                        double ce_in = (kc0errsum[icent]/2 + kc1errsum[icent])/ntrg;
                        if (!addlinear)
                            ce_in = sqrt(kc0errsum[icent]/4 + kc1errsum[icent])/ntrg;
                        
                        BottomDrawer(h1sumin[icent],"(e)",Form("subtractbg/%s_e.pdf",h1sumin[icent]->GetName()),fe_in,ce_in);

                        kdf0->SetParameter(0,(kf0sum[icent]-f0sum[icent])/ntrg);
                        kdc0->SetParameter(0,(kc0sum[icent]-c0sum[icent])/ntrg);

                        TF1* kdf1 = new TF1("kdf1","[0]",0.,1.);
                        kdf1->SetParameter(0,(kf1sum[icent]-f1sum[icent])/ntrg);

                        TF1* kdc1 = new TF1("kdc1","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        kdc1->SetParameter(0,(kc1sum[icent]-c1sum[icent])/ntrg);

                        list<pair<TF1*,string>> kdmapin {
                            {kdf0,"(#kappa_{a}-1) 2#pi#LTa_{t}#GT#LTa_{h}#GT/n_{t}"},
                            {kdc0,"(#kappa_{b}-1) #pi#LTb_{t}#sigma#GT#LTb_{h}#sigma#GT cos(2#Delta#phi)/n_{t}"},
                            {kdf1,"(#kappa_{ht}-1) 4#LTa_{h}#GT#LTb_{t}#sigma#GT/n_{t}"},
                            {kdc1,"(#kappa_{th}-1) 4#LTa_{t}#GT#LTb_{h}#sigma#GT cos(2#Delta#phi)/n_{t}"}
                        };
                                

                        BigDrawer(h1sumin[icent],kdmapin,Form("subtractbg/%s.pdf",h1sumin[icent]->GetName()));

                        //-------------------- out-of-plane case --------------------

                        ntrg = ntrgout[icent];
                        ofsntrgout<<sufg.str()<<" sumcent"<<icent<<" "<<ntrg<<endl;

                        info.idpsi = OUT;
                        info.ntrg = ntrg;

                        h1sumout[icent]->Scale(1./ntrg);

                        TH1F* h1sumout_nok = (TH1F*) h1sumout[icent]->Clone(Form("%s_nok",h1sumout[icent]->GetName()));

                        //---------- initial flat contribution ----------

                        f0->SetParameter(0,f0sum[icent]/ntrg);
                        kf0->SetParameter(0,kf0sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(1);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumout[icent],f0,kf0);


                        {
                            double kf0err = kf0errsum[icent];
                            if (!addlinear)
                                kf0err = sqrt(kf0err);
                            kf0err/=(2.*ntrg);

                            info.flate[0] = kf0err;

                            SmallDrawer(h1sumout[icent],kf0,"(1/n_{t}) 2#pi#kappa_{a}#LTa_{t}#GT#LTa_{h}#GT","(a)",Form("subtractbg/%s_a.pdf",h1sumout[icent]->GetName()),kf0err);
                        }

                        //subtract the contribution
                        AddToHisto(h1sumout[icent],kf0,-1.0);
                        AddToHisto(h1sumout_nok,f0,-1.0);

                        //---------- additional flat contribution ----------

                        f1->SetParameter(0,-f1sum[icent]/ntrg);
                        kf1->SetParameter(0,-kf1sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(2);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumout[icent],f1,kf1);


                        {
                            double kf1err = kf1errsum[icent];
                            if (!addlinear)
                                kf1err = sqrt(kf1err);
                            kf1err/=ntrg;

                            info.flate[1] = kf1err;
                        
                            SmallDrawer(h1sumout[icent],kf1,"#minus (1/n_{t}) 4#kappa_{ht}#LTa_{h}#GT#LTb_{t}#sigma#GT","(b)",Form("subtractbg/%s_b.pdf",h1sumout[icent]->GetName()),kf1err);
                        }

                        //subtract the contribution
                        AddToHisto(h1sumout[icent],kf1,-1.0);
                        AddToHisto(h1sumout_nok,f1,-1.0);

                        //---------- additional cosine contribution ----------

                        c1->SetParameter(0,-c1sum[icent]/ntrg);
                        kc1->SetParameter(0,-kc1sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(3);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumout[icent],c1,kc1);


                        {
                            double kc1err = kc1errsum[icent];
                            if (!addlinear)
                                kc1err = sqrt(kc1err);
                            kc1err/=ntrg;

                            info.cose[1] = kc1err;
                        
                            SmallDrawer(h1sumout[icent],kc1,"#minus (1/n_{t}) 4#kappa_{th}#LTa_{t}#GT#LTb_{h}#sigma#GT cos2#Delta#phi","(c)",Form("subtractbg/%s_c.pdf",h1sumout[icent]->GetName()),kc1err);
                        }

                        AddToHisto(h1sumout[icent],kc1,-1.0);
                        AddToHisto(h1sumout_nok,c1,-1.0);

                        //---------- initial cosine contribution ----------

                        c0->SetParameter(0,c0sum[icent]/ntrg);
                        kc0->SetParameter(0,kc0sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(4);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumout[icent],c0,kc0);

                        {
                            double kc0err = kc0errsum[icent];
                            if (!addlinear)
                                kc0err = sqrt(kc0err);
                            kc0err/=(2.*ntrg);

                            info.cose[0] = kc0err;
                            tree.Fill();                            

                            SmallDrawer(h1sumout[icent],kc0,"(1/n_{t}) #pi#kappa_{b}#LTb_{t}#sigma#GT#LTb_{h}#sigma#GT cos2#Delta#phi","(d)",Form("subtractbg/%s_d.pdf",h1sumout[icent]->GetName()),kc0err);
                        }

                        AddToHisto(h1sumout[icent],kc0,-1.0);
                        AddToHisto(h1sumout_nok,c0,-1.0);

                        //write result to file
                        output->WriteTObject(h1sumout[icent]);

                        //---------- draw the result ----------

                        pad = (TPad*) canio->cd(5);
                        pad->Clear();
                        h1sumout[icent]->Draw("E1");
                        line1->Draw("same");

                        pad = (TPad*) canio->cd(6);
                        pad->Clear();
                        h1sumout_nok->Draw("E1");
                        {
                            double min = h1sumout_nok->GetMinimum();
                            min = min - fabs(min)*.1;
                            if (min > 0.)
                            {
                                double max = h1sumout_nok->GetMaximum();
                                max = max + fabs(max)*.1;
                                h1sumout_nok->GetYaxis()->SetRangeUser(0.,max);
                                pad->Update();
                            }
                        }
                        line1->Draw("same");

                        canio->Update();
                        canio->SaveAs(Form("subtractbg/samepad/%s.pdf",h1sumout[icent]->GetName()));

                        double fe_out = (kf0errsum[icent]/2 + kf1errsum[icent])/ntrg;
                        if (!addlinear)
                            fe_out = sqrt(kf0errsum[icent]/4 + kf1errsum[icent])/ntrg;
                        double ce_out = (kc0errsum[icent]/2 + kc1errsum[icent])/ntrg;
                        if (!addlinear)
                            ce_out = sqrt(kc0errsum[icent]/4 + kc1errsum[icent])/ntrg;
                        
                        BottomDrawer(h1sumout[icent],"(e)",Form("subtractbg/%s_e.pdf",h1sumout[icent]->GetName()),fe_out,ce_out);

                        kdf1->SetParameter(0,-(kf1sum[icent]-f1sum[icent])/ntrg);
                        kdc1->SetParameter(0,-(kc1sum[icent]-c1sum[icent])/ntrg);

                        list<pair<TF1*,string>> kdmapout {
                            {kdf0,"(#kappa_{a}-1) 2#pi#LTa_{t}#GT#LTa_{h}#GT/n_{t}"},
                            {kdc0,"(#kappa_{b}-1) #pi#LTb_{t}#sigma#GT#LTb_{h}#sigma#GT cos(2#Delta#phi)/n_{t}"},
                            {kdf1,"-(#kappa_{ht}-1) 4#LTa_{h}#GT#LTb_{t}#sigma#GT/n_{t}"},
                            {kdc1,"-(#kappa_{th}-1) 4#LTa_{t}#GT#LTb_{h}#sigma#GT cos(2#Delta#phi)/n_{t}"}
                        };
                                

                        BigDrawer(h1sumout[icent],kdmapout,Form("subtractbg/%s.pdf",h1sumout[icent]->GetName()));
                    }
                }
            }
        }
    }

    cout<<"Writing TTree "<<tree.GetName()<<" to "<<output->GetName()<<endl;
    output->WriteTObject(&tree);

    cout<<"Output written to "<<output->GetName()<<endl;
    output->Close();
    cout<<"Done!"<<endl;
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

double GetNTrg(TH1D* h1cent, int icent)
{
    double ntrg = 0.;
    for (int jcent = GAT::limcent[icent+1]; jcent > GAT::limcent[icent]; jcent--)
    {
        int tmpbin = h1cent->FindBin(jcent);
        ntrg += h1cent->GetBinContent(tmpbin);
    }
    return ntrg;
}

vector<double> GetNevt(TFile* file)
{

    TH1D* h1cent = nullptr;
    file->GetObject("bz0cent_py",h1cent);

    vector<double> vnevt;
    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        int b0 = h1cent->FindBin(GAT::limcent[icent]+1);
        int b1 = h1cent->FindBin(GAT::limcent[icent+1]);
        vnevt.push_back(h1cent->Integral(b0,b1));
    }
    return vnevt;
}

void OpenTheDrawer(TPad* pad, TH1* h1, TF1* f1, TF1* f2 = nullptr)
{
    pad->cd();
    h1->DrawCopy("E1");

    double y0 = h1->GetMinimum();
    double y1 = h1->GetMaximum();

    double f1min = f1->GetMinimum(0.,1.);
    y0 = min(y0,f1min);

    double f1max = f1->GetMaximum(0.,1.);
    y1 = max(y1,f1max);

    if (f2)
    {
        double f2min = f2->GetMinimum(0.,1.);
        y0 = min(y0,f2min);

        double f2max = f2->GetMaximum(0.,1.);
        y1 = max(y1,f2max);
    }

    y0 = y0 - fabs(y0)*.1;
    y1 = y1 + fabs(y1)*.1;
    h1->GetYaxis()->SetRangeUser(y0,y1);

    h1->DrawCopy("E1");
    
    f1->SetLineWidth(1);
    f1->SetLineColor(kRed);
    f1->Draw("same");
    
    if (f2)
    {
        f2->SetLineWidth(1);
        f2->SetLineColor(kBlue);
        f2->Draw("same");
    }

    return;
}

void SmallDrawer(TH1F* h1, TF1* f1, const char* fname, const char* label, const char* filename, double err)
{

    TCanvas can ("can","",2000,2000);
    can.SetLeftMargin(.175);

    h1->Draw("E1");

    h1->GetYaxis()->SetTitleOffset(1.7);
    h1->GetXaxis()->SetLabelOffset(.015);

    h1->GetXaxis()->SetNdivisions(505);


    double hmin = h1->GetBinContent(h1->GetMinimumBin()) - h1->GetBinError(h1->GetMinimumBin());
    double hmax = h1->GetBinContent(h1->GetMaximumBin()) + h1->GetBinError(h1->GetMaximumBin());
    
    double y0 = min(hmin,f1->GetMinimum(0.,1.));
    double y1 = max(hmax,f1->GetMaximum(0.,1.));


    y0 = y0 - (y1-y0)*.15;
    y1 = y1 + (y1-y0)*.25;
    
    h1->GetYaxis()->SetRangeUser(y0,y1);
    can.Update();

    if (err != -999.)
    {
        double par = f1->GetParameter(0);

        TF1* f1min = (TF1*) f1->Clone("f1min");
        f1min->SetParameter(0,par-err);

        TF1* f1max = (TF1*) f1->Clone("f1max");
        f1max->SetParameter(0,par+err);
        
        const Int_t n = 24;
        Double_t x[n], ymin[n], ymax[n];
        for (int i=0;i<n;i++) {
            x[i] = i/(n-1.);

            ymax[i] = f1max->Eval(x[i]);
            ymin[i] = f1min->Eval(x[i]);

        }
        TGraph *grmin = new TGraph(n,x,ymin);
        TGraph *grmax = new TGraph(n,x,ymax);
        TGraph *grshade = new TGraph(2*n);
        for (int i=0;i<n;i++) {
            grshade->SetPoint(i,x[i],ymax[i]);
            grshade->SetPoint(n+i,x[n-i-1],ymin[n-i-1]);
        }
        grshade->SetFillStyle(3013);
        grshade->SetFillColor(17);

        
        grshade->Draw("f");

        grmin->SetLineColor(17);
        grmin->SetLineWidth(1);
        grmin->Draw("l");

        grmax->SetLineColor(17);
        grmax->SetLineWidth(1);
        grmax->Draw("l");

    }

    TLine line(0.,0.,1.,0.);
    line.SetLineColor(1);
    line.SetLineWidth(1);
    line.SetLineStyle(9);
    if (y0 < 0 && y1 > 0)
        line.Draw();

    f1->SetLineColor(TColor::GetColorPalette(254*5/10));
    f1->SetLineWidth(14);
    f1->Draw("same");


    TLegend* leg;
    double legx0 = 0.425, legx1 = 0.89;
    if (strstr(fname,"cos") != NULL)
    {
        legx0 = .25;
        

    }


    leg = new TLegend(legx0,0.825,legx1,0.89);
    double legw = legx1 - legx0;
    leg->SetMargin(.08/legw);
    leg->SetBorderSize(0);
    leg->AddEntry(f1,fname,"l");
    leg->SetTextSize(0.05);
    leg->SetTextAlign(22);
    leg->Draw();


    TText txt (0.85,.775,label);
    txt.SetTextAlign(22);
    txt.SetNDC();
    txt.SetTextSize(0.045);
    txt.Draw();

    h1->Draw("E1 same");
    can.Update();
    can.SaveAs(filename);

}

void BottomDrawer(TH1F* h1, const char* label, const char* filename, double fe, double ce)
{

    TCanvas can ("can","",2000,2000);
    can.SetLeftMargin(.175);

    h1->GetYaxis()->SetTitleOffset(1.7);
    h1->GetXaxis()->SetLabelOffset(.015);

    h1->GetXaxis()->SetNdivisions(505);
    

    double hmin = h1->GetBinContent(h1->GetMinimumBin()) - h1->GetBinError(h1->GetMinimumBin());
    double hmax = h1->GetBinContent(h1->GetMaximumBin()) + h1->GetBinError(h1->GetMaximumBin());
    
    double y0 = min(hmin,0.);
    double y1 = max(hmax,0.);


    y0 = y0 - (y1-y0)*.15;
    y1 = y1 + (y1-y0)*.25;
    
    h1->GetYaxis()->SetRangeUser(y0,y1);

    h1->Draw("E1");

    TLine line(0.,0.,1.,0.);
    line.SetLineColor(1);
    line.SetLineWidth(1);
    line.SetLineStyle(9);
    line.Draw();

    if (fe != -999. && ce != -999.)
    {


        TF1* fp = new TF1(Form("%s_syserr1",h1->GetName()),"sqrt(pow([0],2.)+pow([1]*cos(x*2.*TMath::Pi()),2.))",0.,1.);       
        fp->SetParameter(0,fe);
        fp->SetParameter(1,ce);

        TF1* fm = new TF1(Form("%s_syserr0",h1->GetName()),"-sqrt(pow([0],2.)+pow([1]*cos(x*2.*TMath::Pi()),2.))",0.,1.);
        fm->SetParameter(0,fe);
        fm->SetParameter(1,ce);

        fp->SetLineColor(17);
        fp->Draw("same");
        fm->SetLineColor(17);
        fm->Draw("same");

        fp->Write();
        fm->Write();

        const Int_t n = 24;
        Double_t x[n], ym[n], yp[n];
        for (int i=0;i<n;i++) {
            x[i] = i/(n-1.);

            yp[i] = fp->Eval(x[i]);
            ym[i] = fm->Eval(x[i]);

        }

        TGraph *grshade = new TGraph(2*n);
        for (int i=0;i<n;i++) {
            grshade->SetPoint(i,x[i],yp[i]);
            grshade->SetPoint(n+i,x[n-i-1],ym[n-i-1]);
        }
        grshade->SetFillStyle(3013);
        grshade->SetFillColor(17);
        
        grshade->Draw("f");

    }

    
    TText txt (0.85,.85,label);
    
    txt.SetTextAlign(22);
    txt.SetNDC();
    txt.SetTextSize(0.04);
    txt.Draw();

    h1->Draw("E1 same");
    can.Update();
    can.SaveAs(filename);
}

void BigDrawer(TH1F* h1, const list<pair<TF1*,string>>& kdmap, const char* filename)
{

    TCanvas can ("can","",2000,2750);
    can.SetLeftMargin(.15);


    h1->Draw("E1");
    h1->SetMarkerSize(5);
    h1->GetXaxis()->SetTitleSize(.045);
    h1->GetXaxis()->SetLabelSize(.04);
    h1->GetYaxis()->SetTitleSize(.045);
    h1->GetYaxis()->SetLabelSize(.04);

    h1->GetXaxis()->SetLabelOffset(.005);
    h1->GetXaxis()->SetTitleOffset(.8);

    TLine line(0.,0.,1.,0.);
    line.SetLineWidth(4);
    line.Draw();

    TLegend leg;
    leg.SetX2(0.9);
    leg.SetY2(0.9);
    if (kdmap.size() == 2)
    {
        leg.SetX1(.21);
        leg.SetY1(.77);
    }
    else
    {
        leg.SetX1(.21);
        leg.SetY1(.64);

    }
    leg.SetTextSize(0.0425);

    double y0 = 0.;
    double y1 = 0.;
    int i = -1;
    for (const auto& kd : kdmap)
    {
        kd.first->SetLineColor(colors[++i]);
        kd.first->SetLineStyle(styles[i]);
        kd.first->SetLineWidth(10);
        kd.first->Draw("same");
        leg.AddEntry(kd.first,kd.second.c_str(),"l");

        y0 = min(y0,kd.first->GetMinimum());
        y1 = max(y1,kd.first->GetMaximum());
    }


    {
        double hmin = h1->GetBinContent(h1->GetMinimumBin()) - h1->GetBinError(h1->GetMinimumBin());
        double hmax = h1->GetBinContent(h1->GetMaximumBin()) + h1->GetBinError(h1->GetMaximumBin());

        y0 = min(hmin,min(y0,0.));
        y1 = max(hmax,max(y1,0.));
    }
    

    y0 = y0 - (y1-y0)*.05;
    if (kdmap.size() == 2)
        y1 = y1 + (y1-y0)*.29;
    else
        y1 = y1 + (y1-y0)*.58;

    h1->GetYaxis()->SetRangeUser(y0,y1);
    h1->GetYaxis()->SetMaxDigits(2);
    
    h1->Draw("E1 same");
    leg.Draw();

    can.Update();


    can.SaveAs(filename);
}

int GetCentBin(int ic)
{
    if (ic < 0 || ic >= GAT::ncent)
    {
        cerr<<"invalid value of "<<ic<<" passed to GetCentBin()"<<endl;
        return -1;
    }

    for (int icent = 0; icent < ncent; icent++)
    {

        if (GAT::limcent[ic] >= limcent[icent] && GAT::limcent[ic+1] <= limcent[icent+1])
            return icent;
    }

    cerr<<"WARNING: Could not find centbin for value "<<ic<<".  Returning -2"<<endl;
    return -2;
}

void AddToHisto(TH1F* h1, TF1* f1, double wt = 1.)
{
    vector<int> emptybins;
    for (int ix = 1; ix < h1->GetNbinsX(); ix++)
    {
        if (h1->GetBinContent(ix) == 0.)
        {
            emptybins.push_back(ix);
        }
    }

    h1->Add(f1,wt);

    for (const int& ix : emptybins)
    {
        h1->SetBinContent(ix,0);
    }
}
        
