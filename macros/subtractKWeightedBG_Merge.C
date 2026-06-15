
//main drawing functions are SmallDrawer and BottomDrawer
//leaving old drawing -- can cause confusion when reading code though

const int ncent = 2;
const int limcent[ncent+1] = {0,20,60};


enum RPO {IN,OUT,ALL};

const int nadj = 2;

const float PI = (float) TMath::Pi();

const EColor colors[4] = {kRed,kBlue,kGreen,kViolet};
const short styles[4] = {1,9,7,3};

namespace GAT {
    #include "gamAnaTree.h"
}


int ReadFile(const char* inname, TFile*& file);
double GetNTrg(TH1D* h1cent, int icent);
vector<double> GetNevt(TFile* file);
void OpenTheDrawer(TPad* pad, TH1* h1, TF1* f1, TF1* f2 = nullptr);
void SmallDrawer(TH1F* h1, TF1* f1, const char* fname, const char* label, const char* filename);
void BottomDrawer(TH1F* h1, const char* label, const char* filename);
void BigDrawer(TH1F* h1, const list<pair<TF1*,string>>& kdmap, const char* filename);

TFile* output;

void subtractKWeightedBG_Merge(
    const char* namegt = "gamtree_13466_acceptcorrect_merge.root",
    const char* namek = "kComponents.root",
    const char* namekmerge = "kComponents_Merge.root",
    const char* nameout = "weighted_bg_subtracted_merge.root"
    )
{
    cout<<"WARNING: please make sure input files are properly linked to the latest and greatest."<<endl;
    
    ofstream ofsntrgall("ntrgall_merge.txt");
    ofstream ofsntrgin("ntrgin_merge.txt");
    ofstream ofsntrgout("ntrgout_merge.txt");
    
    //style optimized for smalldrawer canvas
    //comment out most if interested in other plots
    gStyle->SetLineScalePS(1);
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0.);
    gStyle->SetEndErrorSize(10);
    gStyle->SetMarkerStyle(kFullCircle);
    gStyle->SetMarkerSize(1);
    gStyle->SetHistLineWidth(1);
    gStyle->SetMarkerColor(kBlack);
    gStyle->SetHistLineColor(kBlack);
    gStyle->SetTitleOffset(1.4,"Y");
    gROOT->ForceStyle();

    TFile* ingt = nullptr;
    if (ReadFile(namegt,ingt)<0) return;

    TFile* ink = nullptr;
    if (ReadFile(namek,ink)<0) return;

    TFile* inkmerge = nullptr;
    if (ReadFile(namekmerge,inkmerge)<0) return;

    TGraphErrors* grsigma = nullptr;
    inkmerge->GetObject("npartsigma",grsigma);

    output = new TFile(nameout,"RECREATE");
    output->WriteTObject(grsigma);

    vector<double> vnevt = GetNevt(ink);

    TCanvas* canio = new TCanvas("canio","",1200,1350);
    canio->Divide(2,3);

    TCanvas* canall = new TCanvas("canall","",600,1350);
    canall->Divide(1,4);

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;

            //saving too many histograms. skip these

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                //suffix of trigger related objects (histos, tgraphs, etc.)
                ostringstream sufg;
                sufg<<"adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt;

                //--------------- get ag component ----------
                TGraphErrors* grag = nullptr;
                inkmerge->GetObject(Form("ag_Np_%s",sufg.str().c_str()),grag);

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

                    double c0sum[ncent+1] = {};
                    double kc0sum[ncent+1] = {};

                    double f1sum[ncent+1] = {};
                    double kf1sum[ncent+1] = {};

                    double c1sum[ncent+1] = {};
                    double kc1sum[ncent+1] = {};


                    double ntrgall[ncent+1] = {};
                    double ntrgin[ncent+1] = {};
                    double ntrgout[ncent+1] = {};


                    //----- get remaining contributions and k corrections -----

                    //--- contributions ---
                    TGraphErrors* grbgs = nullptr;
                    inkmerge->GetObject(Form("bgs_Np_%s",sufh.str().c_str()),grbgs);

                    TGraphErrors* grah = nullptr;
                    inkmerge->GetObject(Form("ah_Np_%s",sufh.str().c_str()),grah);

                    TGraphErrors* grbhs = nullptr;
                    inkmerge->GetObject(Form("bhs_Np_%s",sufh.str().c_str()),grbhs);

                    TGraphErrors* grbgsbhs = nullptr;
                    inkmerge->GetObject(Form("bgsbhs_Np_%s",sufh.str().c_str()),grbgsbhs);

                    //--- k corrections ---

                    TGraphErrors* grka = nullptr;
                    ink->GetObject(Form("kamerge_Np_%s",sufh.str().c_str()),grka);

                    TGraphErrors* grkb = nullptr;
                    ink->GetObject(Form("kbmerge_Np_%s",sufh.str().c_str()),grkb);

                    TGraphErrors* grkgh = nullptr;
                    ink->GetObject(Form("kghmerge_Np_%s",sufh.str().c_str()),grkgh);

                    TGraphErrors* grkhg = nullptr;
                    ink->GetObject(Form("khgmerge_Np_%s",sufh.str().c_str()),grkhg);

                    //------------------------------------------------------------

                    for (int icent = 0; icent < ncent; icent++)
                    {
                        //this section does not pertain to _Merge.C analysis


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

                        //---------- get corrections for this centrality ----------
                        double ka = grka->GetY()[icent];

                        double kb = grkb->GetY()[icent];

                        double kgh = grkgh->GetY()[icent];

                        double khg = grkhg->GetY()[icent];


                        //-------------------- all-plane case --------------------
                        TH1D* h1trgcent = nullptr;
                        ingt->GetObject(Form("cent_adj%i_%s_ptt%i_dpsi%i",iadj,GAT::spec[itype],iptt,ALL),h1trgcent);
                        double ntrg = GetNTrg(h1trgcent,icent);
                        ntrgall[icent] = ntrg;
                        ntrgall[ncent] += ntrg;

                        //get signal
                        ostringstream oss;
                        oss<<"dphi";
                        if (iadj)
                            oss<<"_adj";
                        oss<<"_"<<GAT::spec[itype]<<"_mix0_ptt"<<iptt<<"_cent"<<icent<<"_dpsi"<<ALL<<"_xi"<<ixi;

                        TH1F* h1dphiall = nullptr;
                        ingt->GetObject(oss.str().c_str(),h1dphiall);
                        h1dphiall->SetTitle(";#Delta#phi/#pi;1/N_{#gamma}  dN_{h}/d#Delta#phi");
                        h1dphiall->Scale(1./(h1dphiall->GetBinWidth(1)*PI));

                        TH1F* h1dphiall_nok = (TH1F*) h1dphiall->Clone(Form("%s_nok",oss.str().c_str()));


                        {
                            
                            ostringstream osssum;
                            osssum<<"dphi";
                            if (iadj)
                                osssum<<"_adj";
                            osssum<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_sumcent"<<icent<<"_dpsi"<<ALL<<"_xi"<<ixi;

                            cout<<"Cloning "<<h1dphiall->GetName()<<" to "<<osssum.str()<<endl;
                            h1sumall[icent] = (TH1F*) h1dphiall->Clone(osssum.str().c_str());

                            if (icent == 0)
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
                        

                        //normalize
                        h1dphiall->Scale(1./ntrg);
                        h1dphiall_nok->Scale(1./ntrg);

                        //---------- initial flat contribution ----------


                        double f0scale = 4.*PI*ag*ah;
                        f0scale*=nevt;
                        f0sum[icent] += f0scale;
                        f0sum[ncent] += f0scale;
                        TF1* f0 = new TF1("f0","[0]",0.,1.);
                        f0->SetParameter(0,f0scale/ntrg);

                        double kf0scale = ka*f0scale;
                        kf0sum[icent] += kf0scale;
                        kf0sum[ncent] += kf0scale;
                        TF1* kf0 = new TF1("kf0","[0]",0.,1.);
                        kf0->SetParameter(0,kf0scale/ntrg);

                        TPad* pad = (TPad*) canall->cd(1);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiall,f0,kf0);

                        //subtract the contribution
                        h1dphiall->Add(kf0,-1.0);
                        h1dphiall_nok->Add(f0,-1.0);

                        //---------- initial cosine contribution ----------

                        double c0scale = 2.*PI*bgsbhs;
                        c0scale *= nevt;
                        c0scale /= (sigma*sigma);
                        c0sum[icent] += c0scale;
                        c0sum[ncent] += c0scale;
                        TF1* c0 = new TF1("c0","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        c0->SetParameter(0,c0scale/ntrg);

                        double kc0scale = 2.*PI*bgsbhs;
                        kc0scale *= nevt;
                        kc0scale *= kb;
                        kc0sum[icent] += kc0scale;
                        kc0sum[ncent] += kc0scale;
                        TF1* kc0 = new TF1("kc0","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        kc0->SetParameter(0,kc0scale/ntrg);

                        pad = (TPad*) canall->cd(2);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiall,c0,kc0);

                        h1dphiall->Add(kc0,-1.0);
                        h1dphiall_nok->Add(c0,-1.0);

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
                        canall->SaveAs(Form("subtractbgmerge/samepad/cent/%s.eps",oss.str().c_str()));

                        //-------------------- in-plane case --------------------
                        h1trgcent = nullptr;
                        ingt->GetObject(Form("cent_adj%i_%s_ptt%i_dpsi%i",iadj,GAT::spec[itype],iptt,IN),h1trgcent);
                        ntrg = GetNTrg(h1trgcent,icent);
                        ntrgin[icent] = ntrg;
                        ntrgin[ncent] += ntrg;

                        //get signal
                        oss.str("");
                        oss<<"dphi";
                        if (iadj)
                            oss<<"_adj";
                        oss<<"_"<<GAT::spec[itype]<<"_mix0_ptt"<<iptt<<"_cent"<<icent<<"_dpsi"<<IN<<"_xi"<<ixi;

                        TH1F* h1dphiin = nullptr;
                        ingt->GetObject(oss.str().c_str(),h1dphiin);
                        h1dphiin->SetTitle(";#Delta#phi/#pi;1/N_{#gamma}  dN_{h}/d#Delta#phi");
                        h1dphiin->Scale(1./(h1dphiin->GetBinWidth(1)*PI));

                        TH1F* h1dphiin_nok = (TH1F*) h1dphiin->Clone(Form("%s_nok",oss.str().c_str()));


                        {
                            ostringstream osssum;
                            osssum<<"dphi";
                            if (iadj)
                                osssum<<"_adj";
                            osssum<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_sumcent"<<icent<<"_dpsi"<<IN<<"_xi"<<ixi;

                            cout<<"Cloning "<<h1dphiin->GetName()<<" to "<<osssum.str()<<endl;

                            h1sumin[icent] = (TH1F*) h1dphiin->Clone(osssum.str().c_str());

                            if (icent == 0)
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
                        h1dphiin->Add(kf0,-1.0);
                        h1dphiin_nok->Add(f0,-1.0);

                        //---------- additional flat contribution ----------

                        double f1scale = 4.*ah*bgs;
                        f1scale *= nevt;
                        f1sum[icent] += f1scale;
                        f1sum[ncent] += f1scale;
                        TF1* f1 = new TF1("f1","[0]",0.,1.);
                        f1->SetParameter(0,f1scale/ntrg);

                        double kf1scale = khg*f1scale;
                        kf1sum[icent] += kf1scale;
                        kf1sum[ncent] += kf1scale;
                        TF1* kf1 = new TF1("kf1","[0]",0.,1.);
                        kf1->SetParameter(0,kf1scale/ntrg);

                        pad = (TPad*) canio->cd(2);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiin,f1,kf1);

                        //subtract the contribution
                        h1dphiin->Add(kf1,-1.0);
                        h1dphiin_nok->Add(f1,-1.0);

                        //---------- additional cosine contribution ----------

                        double c1scale = 4.*ag*bhs;
                        c1scale *= nevt;
                        c1sum[icent] += c1scale;
                        c1sum[ncent] += c1scale;
                        TF1* c1 = new TF1("c1","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        c1->SetParameter(0,c1scale/ntrg);

                        double kc1scale = kgh*c1scale;
                        kc1sum[icent] += kc1scale;
                        kc1sum[ncent] += kc1scale;
                        TF1* kc1 = new TF1("kc1","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        kc1->SetParameter(0,kc1scale/ntrg);

                        pad = (TPad*) canio->cd(3);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiin,c1,kc1);

                        h1dphiin->Add(kc1,-1.0);
                        h1dphiin_nok->Add(c1,-1.0);

                        //---------- initial cosine contribution ----------

                        c0scale /= 2.0;
                        c0->SetParameter(0,c0scale/ntrg);

                        kc0scale /= 2.0;
                        kc0->SetParameter(0,kc0scale/ntrg);

                        pad = (TPad*) canio->cd(4);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiin,c0,kc0);

                        h1dphiin->Add(kc0,-1.0);
                        h1dphiin_nok->Add(c0,-1.0);

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
                        canio->SaveAs(Form("subtractbgmerge/samepad/cent/%s.eps",oss.str().c_str()));

                        //-------------------- out-of-plane case --------------------
                        h1trgcent = nullptr;
                        ingt->GetObject(Form("cent_adj%i_%s_ptt%i_dpsi%i",iadj,GAT::spec[itype],iptt,OUT),h1trgcent);
                        ntrg = GetNTrg(h1trgcent,icent);
                        ntrgout[icent] = ntrg;
                        ntrgout[ncent] += ntrg;

                        //get signal
                        oss.str("");
                        oss<<"dphi";
                        if (iadj)
                            oss<<"_adj";
                        oss<<"_"<<GAT::spec[itype]<<"_mix0_ptt"<<iptt<<"_cent"<<icent<<"_dpsi"<<OUT<<"_xi"<<ixi;

                        TH1F* h1dphiout = nullptr;
                        ingt->GetObject(oss.str().c_str(),h1dphiout);
                        h1dphiout->SetTitle(";#Delta#phi/#pi;1/N_{#gamma}  dN_{h}/d#Delta#phi");
                        h1dphiout->Scale(1./(h1dphiout->GetBinWidth(1)*PI));

                        TH1F* h1dphiout_nok = (TH1F*) h1dphiout->Clone(Form("%s_nok",oss.str().c_str()));

                        
                        {
                            ostringstream osssum;
                            osssum<<"dphi";
                            if (iadj)
                                osssum<<"_adj";
                            osssum<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_sumcent"<<icent<<"_dpsi"<<OUT<<"_xi"<<ixi;

                            cout<<"Cloning "<<h1dphiout->GetName()<<" to "<<osssum.str()<<endl;
                            h1sumout[icent] = (TH1F*) h1dphiout->Clone(osssum.str().c_str());

                            if (icent == 0)
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
                        h1dphiout->Add(kf0,-1.0);
                        h1dphiout_nok->Add(f0,-1.0);

                        //---------- additional flat contribution ----------

                        f1->SetParameter(0,-f1scale/ntrg);
                        kf1->SetParameter(0,-kf1scale/ntrg);

                        pad = (TPad*) canio->cd(2);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiout,f1,kf1);

                        //subtract the contribution
                        h1dphiout->Add(kf1,-1.0);
                        h1dphiout_nok->Add(f1,-1.0);

                        //---------- additional cosine contribution ----------

                        c1->SetParameter(0,-c1scale/ntrg);
                        kc1->SetParameter(0,-kc1scale/ntrg);

                        pad = (TPad*) canio->cd(3);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiout,c1,kc1);

                        h1dphiout->Add(kc1,-1.0);
                        h1dphiout_nok->Add(c1,-1.0);

                        //---------- initial cosine contribution ----------

                        c0->SetParameter(0,c0scale/ntrg);
                        kc0->SetParameter(0,kc0scale/ntrg);

                        pad = (TPad*) canio->cd(4);
                        pad->Clear();
                        OpenTheDrawer(pad,h1dphiout,c0,kc0);

                        h1dphiout->Add(kc0,-1.0);
                        h1dphiout_nok->Add(c0,-1.0);

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
                        canio->SaveAs(Form("subtractbgmerge/samepad/cent/%s.eps",oss.str().c_str()));

                    }

                    //---------- Draw the centrality merged results ----------
                    for (int icent = 0; icent <= ncent; icent++)
                    {
                        //-------------------- all-plane case --------------------

                        double ntrg = ntrgall[icent];
                        ofsntrgall<<sufg.str()<<" sumcent"<<icent<<" "<<ntrg<<endl;

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

                        SmallDrawer(h1sumall[icent],kf0,"4#pi#kappa_{a}<a_{#gamma}><a_{h}>/N_{#gamma}","(a)",Form("subtractbgmerge/%s_a.eps",h1sumall[icent]->GetName()));

                        //subtract the contribution
                        h1sumall[icent]->Add(kf0,-1.0);
                        h1sumall_nok->Add(f0,-1.0);

                        //---------- initial cosine contribution ----------

                        TF1* c0 = new TF1("c0","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        c0->SetParameter(0,c0sum[icent]/ntrg);

                        TF1* kc0 = new TF1("kc0","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        kc0->SetParameter(0,kc0sum[icent]/ntrg);

                        pad = (TPad*) canall->cd(2);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumall[icent],c0,kc0);

                        SmallDrawer(h1sumall[icent],kc0,"2#pi#kappa_{b}<b_{#gamma}#sigma><b_{h}#sigma>/N_{#gamma}","(b)",Form("subtractbgmerge/%s_b.eps",h1sumall[icent]->GetName()));

                        h1sumall[icent]->Add(kc0,-1.0);
                        h1sumall_nok->Add(c0,-1.0);

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
                        canall->SaveAs(Form("subtractbgmerge/samepad/%s.eps",h1sumall[icent]->GetName()));

                        BottomDrawer(h1sumall[icent],"(c)",Form("subtractbgmerge/%s_c.eps",h1sumall[icent]->GetName()));

                        TF1* kdf0 = new TF1("kdf0","[0]",0.,1.);
                        kdf0->SetParameter(0,(kf0sum[icent]-f0sum[icent])/ntrg);

                        TF1* kdc0 = new TF1("kdc0","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        kdc0->SetParameter(0,(kc0sum[icent]-c0sum[icent])/ntrg);


                        list<pair<TF1*,string>> kdmapall {
                            {kdf0,"4#pi(#kappa_{a}-1)<a_{#gamma}><a_{h}>/N_{#gamma}"},
                            {kdc0,"2#pi(#kappa_{b}-1)<b_{#gamma}#sigma><b_{h}#sigma>/N_{#gamma}"}
                        };
                                

                        BigDrawer(h1sumall[icent],kdmapall,Form("subtractbgmerge/%s.eps",h1sumall[icent]->GetName()));

                        //-------------------- in-plane case --------------------

                        ntrg = ntrgin[icent];
                        ofsntrgin<<sufg.str()<<" sumcent"<<icent<<" "<<ntrg<<endl;

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

                        SmallDrawer(h1sumin[icent],kf0,"2#pi#kappa_{a}<a_{#gamma}><a_{h}>/N_{#gamma}","(a)",Form("subtractbgmerge/%s_a.eps",h1sumin[icent]->GetName()));

                        //subtract the contribution
                        h1sumin[icent]->Add(kf0,-1.0);
                        h1sumin_nok->Add(f0,-1.0);

                        //---------- additional flat contribution ----------

                        TF1* f1 = new TF1("f1","[0]",0.,1.);
                        f1->SetParameter(0,f1sum[icent]/ntrg);

                        TF1* kf1 = new TF1("kf1","[0]",0.,1.);
                        kf1->SetParameter(0,kf1sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(2);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumin[icent],f1,kf1);

                        SmallDrawer(h1sumin[icent],kf1,"4#kappa_{h#gamma}<a_{h}><b_{#gamma}#sigma>/N_{#gamma}","(b)",Form("subtractbgmerge/%s_b.eps",h1sumin[icent]->GetName()));

                        //subtract the contribution
                        h1sumin[icent]->Add(kf1,-1.0);
                        h1sumin_nok->Add(f1,-1.0);

                        //---------- additional cosine contribution ----------

                        TF1* c1 = new TF1("c1","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        c1->SetParameter(0,c1sum[icent]/ntrg);

                        TF1* kc1 = new TF1("kc1","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        kc1->SetParameter(0,kc1sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(3);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumin[icent],c1,kc1);

                        SmallDrawer(h1sumin[icent],kc1,"4#kappa_{#gammah}<a_{#gamma}><b_{h}#sigma>/N_{#gamma}","(c)",Form("subtractbgmerge/%s_c.eps",h1sumin[icent]->GetName()));

                        h1sumin[icent]->Add(kc1,-1.0);
                        h1sumin_nok->Add(c1,-1.0);

                        //---------- initial cosine contribution ----------

                        c0sum[icent] /= 2.0;
                        c0->SetParameter(0,c0sum[icent]/ntrg);

                        kc0sum[icent] /= 2.0;
                        kc0->SetParameter(0,kc0sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(4);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumin[icent],c0,kc0);

                        SmallDrawer(h1sumin[icent],kc0,"#pi#kappa_{b}<b_{#gamma}#sigma><b_{h}#sigma>/N_{#gamma}","(d)",Form("subtractbgmerge/%s_d.eps",h1sumin[icent]->GetName()));

                        h1sumin[icent]->Add(kc0,-1.0);
                        h1sumin_nok->Add(c0,-1.0);

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
                        canio->SaveAs(Form("subtractbgmerge/samepad/%s.eps",h1sumin[icent]->GetName()));

                        BottomDrawer(h1sumin[icent],"(e)",Form("subtractbgmerge/%s_e.eps",h1sumin[icent]->GetName()));

                        kdf0->SetParameter(0,(kf0sum[icent]-f0sum[icent])/ntrg);
                        kdc0->SetParameter(0,(kc0sum[icent]-c0sum[icent])/ntrg);

                        TF1* kdf1 = new TF1("kdf1","[0]",0.,1.);
                        kdf1->SetParameter(0,(kf1sum[icent]-f1sum[icent])/ntrg);

                        TF1* kdc1 = new TF1("kdc1","[0]*cos(x*2.*TMath::Pi())",0.,1.);
                        kdc1->SetParameter(0,(kc1sum[icent]-c1sum[icent])/ntrg);

                        list<pair<TF1*,string>> kdmapin {
                            {kdf0,"(#kappa_{a}-1)2#pi<a_{#gamma}><a_{h}>/N_{#gamma}"},
                            {kdc0,"(#kappa_{b}-1)#pi<b_{#gamma}#sigma><b_{h}#sigma>/N_{#gamma}"},
                            {kdf1,"(#kappa_{h#gamma}-1)4<a_{h}><b_{#gamma}#sigma>/N_{#gamma}"},
                            {kdc1,"(#kappa_{#gammah}-1)4<a_{#gamma}><b_{h}#sigma>/N_{#gamma}"}
                        };
                                

                        BigDrawer(h1sumin[icent],kdmapin,Form("subtractbgmerge/%s.eps",h1sumin[icent]->GetName()));

                        //-------------------- out-of-plane case --------------------

                        ntrg = ntrgout[icent];
                        ofsntrgout<<sufg.str()<<" sumcent"<<icent<<" "<<ntrg<<endl;

                        h1sumout[icent]->Scale(1./ntrg);

                        TH1F* h1sumout_nok = (TH1F*) h1sumout[icent]->Clone(Form("%s_nok",h1sumout[icent]->GetName()));

                        //---------- initial flat contribution ----------

                        f0->SetParameter(0,f0sum[icent]/ntrg);
                        kf0->SetParameter(0,kf0sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(1);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumout[icent],f0,kf0);

                        SmallDrawer(h1sumout[icent],kf0,"2#pi#kappa_{a}<a_{#gamma}><a_{h}>/N_{#gamma}","(a)",Form("subtractbgmerge/%s_a.eps",h1sumout[icent]->GetName()));

                        //subtract the contribution
                        h1sumout[icent]->Add(kf0,-1.0);
                        h1sumout_nok->Add(f0,-1.0);

                        //---------- additional flat contribution ----------

                        f1->SetParameter(0,-f1sum[icent]/ntrg);
                        kf1->SetParameter(0,-kf1sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(2);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumout[icent],f1,kf1);

                        SmallDrawer(h1sumout[icent],kf1,"4#kappa_{h#gamma}<a_{h}><b_{#gamma}#sigma>/N_{#gamma}","(b)",Form("subtractbgmerge/%s_b.eps",h1sumout[icent]->GetName()));

                        //subtract the contribution
                        h1sumout[icent]->Add(kf1,-1.0);
                        h1sumout_nok->Add(f1,-1.0);

                        //---------- additional cosine contribution ----------

                        c1->SetParameter(0,-c1sum[icent]/ntrg);
                        kc1->SetParameter(0,-kc1sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(3);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumout[icent],c1,kc1);

                        SmallDrawer(h1sumout[icent],kc1,"4#kappa_{#gammah}<a_{#gamma}><b_{h}#sigma>/N_{#gamma}","(c)",Form("subtractbgmerge/%s_c.eps",h1sumout[icent]->GetName()));

                        h1sumout[icent]->Add(kc1,-1.0);
                        h1sumout_nok->Add(c1,-1.0);

                        //---------- initial cosine contribution ----------

                        c0->SetParameter(0,c0sum[icent]/ntrg);
                        kc0->SetParameter(0,kc0sum[icent]/ntrg);

                        pad = (TPad*) canio->cd(4);
                        pad->Clear();
                        OpenTheDrawer(pad,h1sumout[icent],c0,kc0);

                        SmallDrawer(h1sumout[icent],kc0,"#pi#kappa_{b}<b_{#gamma}#sigma><b_{h}#sigma>/N_{#gamma}","(d)",Form("subtractbgmerge/%s_d.eps",h1sumout[icent]->GetName()));

                        h1sumout[icent]->Add(kc0,-1.0);
                        h1sumout_nok->Add(c0,-1.0);

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
                        canio->SaveAs(Form("subtractbgmerge/samepad/%s.eps",h1sumout[icent]->GetName()));

                        BottomDrawer(h1sumout[icent],"(e)",Form("subtractbgmerge/%s_e.eps",h1sumout[icent]->GetName()));

                        kdf1->SetParameter(0,-(kf1sum[icent]-f1sum[icent])/ntrg);
                        kdc1->SetParameter(0,-(kc1sum[icent]-c1sum[icent])/ntrg);

                        list<pair<TF1*,string>> kdmapout {
                            {kdf0,"(#kappa_{a}-1)2#pi<a_{#gamma}><a_{h}>/N_{#gamma}"},
                            {kdc0,"(#kappa_{b}-1)#pi<b_{#gamma}#sigma><b_{h}#sigma>/N_{#gamma}"},
                            {kdf1,"-(#kappa_{h#gamma}-1)4<a_{h}><b_{#gamma}#sigma>/N_{#gamma}"},
                            {kdc1,"-(#kappa_{#gammah}-1)4<a_{#gamma}><b_{h}#sigma>/N_{#gamma}"}
                        };
                                

                        BigDrawer(h1sumout[icent],kdmapout,Form("subtractbgmerge/%s.eps",h1sumout[icent]->GetName()));
                    }
                }
            }
        }
    }

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
    for (int jcent = limcent[icent+1]; jcent > limcent[icent]; jcent--)
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
    for (int icent = 0; icent < ncent; icent++)
    {
        int b0 = h1cent->FindBin(limcent[icent]+1);
        int b1 = h1cent->FindBin(limcent[icent+1]);
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
    f2->SetLineColor(kRed);
    f1->Draw("same");

    if (f2)
    {
        f2->SetLineWidth(1);
        f2->SetLineColor(kBlue);
        f2->Draw("same");
    }

    return;
}

void SmallDrawer(TH1F* h1, TF1* f1, const char* fname, const char* label, const char* filename)
{
    TCanvas can ("can","",2500,2000);//2x3

    double y0 = min(h1->GetMinimum(),f1->GetMinimum(0.,1.));
    double y1 = max(h1->GetMaximum(),f1->GetMaximum(0.,1.));

    y0 = y0 - (y1-y0)*.25;
    y1 = y1 + (y1-y0)*.25;
    h1->GetYaxis()->SetRangeUser(y0,y1);

    h1->Draw("E1");

    f1->SetLineColor(kBlue);
    f1->SetLineWidth(10);
    f1->Draw("same");

    TLegend leg (0.55,0.825,0.9,0.9);
    leg.AddEntry(f1,fname,"l");
    leg.SetTextSize(0.038);
    leg.SetTextAlign(22);
    leg.Draw();

    TText txt (0.5,.8625,label);
    txt.SetTextAlign(22);
    txt.SetNDC();
    txt.SetTextSize(0.04);
    txt.Draw();

    h1->Draw("E1 same");
    can.Update();
    can.SaveAs(filename);
}

void BottomDrawer(TH1F* h1, const char* label, const char* filename)
{
    TCanvas can ("can","",2500,2000);//2x3

    double y0 = min(h1->GetMinimum(),0.);
    double y1 = max(h1->GetMaximum(),0.);

    y0 = y0 - (y1-y0)*.25;
    y1 = y1 + (y1-y0)*.25;
    h1->GetYaxis()->SetRangeUser(y0,y1);

    h1->Draw("E1");

    TLine line(0.,0.,1.,0.);
    line.SetLineWidth(4);
    line.SetLineStyle(2);
    line.Draw();

    TText txt (0.5,.8625,label);
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
    TCanvas can ("can","",2500,2000);//2x3

    double y0 = min(h1->GetMinimum(),0.);
    double y1 = max(h1->GetMaximum(),0.);

    h1->Draw("E1");

    TLine line(0.,0.,1.,0.);
    line.SetLineWidth(4);
    line.Draw();

    TLegend leg;
    leg.SetX2(0.9);
    leg.SetY2(0.9);
    if (kdmap.size() == 2)
    {
        leg.SetX1(.525);
        leg.SetY1(.75);
    }
    else
    {
        leg.SetX1(.15);
        leg.SetY1(.75);
        leg.SetNColumns(2);
    }
    leg.SetTextSize(0.038);

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

    y0 = y0 - (y1-y0)*.25;
    y1 = y1 + (y1-y0)*.4;
    h1->GetYaxis()->SetRangeUser(y0,y1);
    
    leg.Draw(); 

    h1->Draw("E1 same");
    can.Update();
    can.SaveAs(filename);
}


