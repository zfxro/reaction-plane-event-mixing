//NOTES:
//ka and kamerge errors seem too large. others k's do as well. most pronounced in highest ptt, iptt=3. please check k error calculation

//reads acceptance corrected gamtree results
//for each trigger type, ptt, cent, xi
//     - calculate (or rather isolate) each k component
//     - correct each k component
//     - write resulting TGraphs to file

//use limit for convergence instead of limit of 10. (USE RELATIVE LIMIT)

//switch between tspline3 and tspline
//from tspline to tspline3 search and replace
//Eval(npart) with Eval(npart,0,"S")
//Eval(npartlo) with Eval(npartlo,0,"S")
//Eval(nparthi) with Eval(nparthi,0,"S")
//from tspline3 to tspline search and replace
//,0,"S") with )

const int nadj = 2;

const int npartmin = 2;
const int npartmax = 390;

//actual npart limits for our distributions
const double npartlo = 11.;
const double nparthi = 389.;

const float PI = (float) TMath::Pi();

const int ncentmerge = 2;
const int limcentmerge[ncentmerge+1] = {0,20,60};

// dropping this and moving to pt based lims
// keys must match dir names for new ks
map <string, pair <double,double>> mapklim = {
    {"kanew", make_pair(.996,1.094)},
    {"kbnew", make_pair(1.35,4.6)},
    {"kghnew", make_pair(0.94,1.16)},
    {"khgnew", make_pair(0.94,1.16)}
};

map <int, pair <double,double>> kalim = {
    {0, make_pair(.998,1.078)},
    {1, make_pair(.998,1.076)},
    {2, make_pair(.996,1.082)},
    {3, make_pair(.996,1.096)}
};

map <int, pair <double,double>> kblim = {
    {0, make_pair(1.3,3.9)},
    {1, make_pair(1.3,3.9)},
    {2, make_pair(1.3,4.4)},
    {3, make_pair(.6,6.2)}
};

map <int, pair <double,double>> kghlim = {
    {0, make_pair(.975,1.13)},
    {1, make_pair(.975,1.13)},
    {2, make_pair(.975,1.17)},
    {3, make_pair(.94,1.24)}
};

map <int, pair <double,double>> khglim = {
    {0, make_pair(.975,1.13)},
    {1, make_pair(.955,1.095)},
    {2, make_pair(.965,1.19)},
    {3, make_pair(.85,2.15)}
};

TRandom* rndm;

TGraphErrors* tgrsig;

namespace GAT {
    #include "gamAnaTree.h"
}

int GetNpartDists(TH1F* hnpart[GAT::ncent], vector<double>& vnpart, vector<double>& vnparte);
int GetNpartDists_Merge(TH1F* hnpart[ncentmerge], vector<double>& vnpart, vector<double>& vnparte);
int GetRPReso(vector<double>& vsig);
int Getng(TH1F* h1trgcent, vector<double>& vng, vector<double>& vnge, int icent, double nevt);

int MakeTGraphErrors(TGraphErrors*& tgr, string dir, string name, string xname, string yname, const int size, const double* x, const double* y, const double* ex = 0, const double* ey = 0, TCanvas* can = NULL);

double GetError(TH1F* h1, TF1* f1);
void DrawError(TF1* f1, double err, int ipar);
void DrawErrorExtra(TF1* f1, double err, int ipar);

void CalcK(TH1F* hnpart[GAT::ncent], TGraphErrors* tgra, TGraphErrors* tgrb, TGraphErrors* tgrae, TGraphErrors* tgrbe, vector<double>& vk, vector<double>& vke, vector<double>& num, vector<double>& nume, vector<double>& den, vector<double>& dene, TGraphErrors* tgrad = nullptr, TGraphErrors* tgrbd = nullptr, TGraphErrors* tgrade = nullptr, TGraphErrors* tgrbde = nullptr, bool usespline = true, bool flip = false);

void CalcK_Merge(TH1F* hnpart[ncentmerge], TGraphErrors* tgra, TGraphErrors* tgrb, TGraphErrors* tgrae, TGraphErrors* tgrbe, vector<double>& vk, vector<double>& vke, vector<double>& vnum, vector<double>& vnume, vector<double>& vden, vector<double>& vdene, TGraphErrors* tgrad = nullptr, TGraphErrors* tgrbd = nullptr, TGraphErrors* tgrade = nullptr, TGraphErrors* tgrbde = nullptr, bool usespline = true, bool flip = false);

double SetNewTGraphPoint(TGraphErrors* tgra, TGraphErrors* tgranew, int ipt);
double SetFirstTGraphPoint(TGraphErrors* tgr, TGraphErrors* tgrnew);
double SetLastTGraphPoint(TGraphErrors* tgr, TGraphErrors* tgrnew);
void SetPointWithRatio(TGraphErrors* tgr, TGraphErrors* tgrnew, int ipt, double rat);

vector<double> CalcKMC(TGraphErrors* tgra, TGraphErrors* tgrb, const vector<double>& vkold, const vector<double>& vkeold, const vector<double>& vnpart, TH1F* hnpart[GAT::ncent], TMultiGraph* mgrk, TMultiGraph* mgra, TMultiGraph* mgrb, TGraphErrors* tgrad = nullptr, TGraphErrors* tgrbd = nullptr, TMultiGraph* mgrad = nullptr, TMultiGraph* mgrbd = nullptr);

void calcKComponents(
    bool quicksavetrigs = false,
    const char* outname = "kComponents.root",
    const char* inname = "gamtree_acceptcorrect.root",
    const char* inname_gamhist = "gamAna_histos.root",
    const char* inname_trigs = "gamtree_triggers.root"
    )
{
    rndm = new TRandom3();//just use default. seed will always be set same

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(15);
    gStyle->SetPalette(kBird);

    gStyle->SetLineScalePS(1);
    //if change these settings please change in subtractKWeightedBG.C as well
    gStyle->SetTitleOffset(1.4,"Y");
    gStyle->SetLabelSize(.045,"Y");
    gStyle->SetTitleSize(.05,"Y");
    gStyle->SetLabelSize(.045,"X");
    gStyle->SetTitleSize(.05,"X");
    gROOT->ForceStyle();

    //Get Npart distributions
    TH1F* hnpart[GAT::ncent];
    vector<double> vnpart;
    vector<double> vnparte;
    if (GetNpartDists(hnpart,vnpart,vnparte)<0) return;

    //Get Npart merge distributions
    TH1F* hnpartmerge[ncentmerge];
    vector<double> vnpartmerge;
    vector<double> vnpartmergee;
    if (GetNpartDists_Merge(hnpartmerge,vnpartmerge,vnpartmergee)<0) return;

    // Set up centrality vector for TGraphs
    vector<double> vcent;
    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        double tmpcent = GAT::limcent[icent+1] + GAT::limcent[icent]+1;
        tmpcent/=2.;
        vcent.push_back(tmpcent);
        cout<<"Centrality midpoint of bin "<<icent<<" is "<<tmpcent<<endl;
    }

    //Fill vector of all npart values. Used for TGraphing products as function of Npart
    vector<double> vnpartall;
    for (int npart = npartmin; npart <= npartmax; npart++)
    {
        vnpartall.push_back(npart);
    }

    TFile* output = new TFile(outname,"RECREATE");
    if (!output || output->IsZombie())
    {
        cerr<<"Could not open "<<outname<<endl;
        return;
    }

    TCanvas* cmgr = new TCanvas("cmgr","",2200,2000);
    TCanvas* can = new TCanvas("c1","",2200,1440);

    cmgr->SetLeftMargin(.15);
    cmgr->SetRightMargin(.15);
    cmgr->SetTopMargin(.05);
    cmgr->SetBottomMargin(.125);

    //Get RP Reso
    vector<double> vsig;
    if (GetRPReso(vsig)<0) return;
    output->cd();//closing rpreso input does not reopen output? must cd..
    tgrsig = new TGraphErrors(vsig.size(),&vnpart[0],&vsig[0],&vnparte[0],0);
    tgrsig->SetTitle("npartsigma;N_{part};#sigma");
    tgrsig->SetName("npartsigma");
    cout<<"npartreso has been copied. Printing contents..."<<endl;
    tgrsig->Print();
    tgrsig->Write();
    tgrsig->Draw("AC");
    can->SaveAs("npartsigma.pdf");

    //Get centrality distribution from gamAna histo file
    cout<<"Reading "<<inname_gamhist<<endl;
    TFile* input_gamhist = new TFile(inname_gamhist,"READ");
    if (!input_gamhist || input_gamhist->IsZombie())
    {
        cerr<<"Could not open "<<inname_gamhist<<endl;
        return;
    }
    //this 2d histo is used for cent dist
    TH2F* h2bz0cent = (TH2F*) gROOT->FindObject("bz0cent");
    assert(h2bz0cent);
    TH1F* h1cent = (TH1F*) h2bz0cent->ProjectionY();
    h1cent->SetDirectory(output);
    input_gamhist->Close();
    vector<double> vnevt;
    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        int b0 = h1cent->FindBin(GAT::limcent[icent]+1);
        int b1 = h1cent->FindBin(GAT::limcent[icent+1]);
        vnevt.push_back(h1cent->Integral(b0,b1));
    }

    TFile* input;
    if (quicksavetrigs)
    {
        cout<<"quicksavetrigs set to true. Will not be reading file "<<inname<<endl;
    }
    else
    {
        cout<<"Reading "<<inname<<endl;
        input = new TFile(inname,"READ");
        if (!input || input->IsZombie())
        {
            cerr<<"Could not open "<<inname<<endl;
            return;
        }
    }

    cout<<"Reading "<<inname_trigs<<endl;
    TFile* input_trigs = new TFile(inname_trigs,"READ");
    if (!input_trigs || input_trigs->IsZombie())
    {
        cerr<<"Could not open "<<inname_trigs<<endl;
        return;
    }

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                TLatex* txtmom = new TLatex (0.5,.85,Form("%.1f #leq p_{T} < %.1f",GAT::limptt[iptt],GAT::limptt[iptt+1]));
                txtmom->SetTextAlign(22);
                txtmom->SetNDC();
                txtmom->SetTextSize(0.05);                  

                ostringstream sufg;
                sufg<<"adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt;

                printf("\nTrigger adj%i %s ptt%i dpsi%i\n",iadj,GAT::spec[itype],iptt,GAT::ndpsi);

                input_trigs->cd();

                TH1F* h1trgcent = (TH1F*) gROOT->FindObject(Form("cent_adj%i_%s_ptt%i_dpsi%i",iadj,GAT::spec[itype],iptt,GAT::ndpsi));
                if (!h1trgcent)
                {
                    cerr<<"Could not find ";
                    fprintf(stderr,"cent_adj%i_%s_ptt%i_dpsi%i\n",iadj,GAT::spec[itype],iptt,GAT::ndpsi);
                    return;
                }

                vector<double> vng;
                vector<double> vnge;

                for (int icent = 0; icent < GAT::ncent; icent++)
                {
                    if (Getng(h1trgcent,vng,vnge,icent,vnevt[icent]) < 0)
                        return;
                }

                output->cd();

                TGraph* tgrng = new TGraph(vng.size(),&vcent[0],&vng[0]);
                tgrng->SetTitle(Form("ng_adj%i_%s_ptt%i;centrality;#LTn_{#gamma}#GT",iadj,GAT::spec[itype],iptt));
                tgrng->Write(Form("ng_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt));

                //plot as function of Npart
                TGraphErrors* tgrng_Np = new TGraphErrors(vng.size(),&vnpart[0],&vng[0],&vnparte[0],&vnge[0]);
                tgrng_Np->SetTitle(Form("ng_Np_adj%i_%s_ptt%i;N_{part};#LTn_{#gamma}#GT",iadj,GAT::spec[itype],iptt));
                tgrng_Np->SetName(Form("ng_Np_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt));             
                tgrng_Np->Write();
                tgrng_Np->Draw("AC");
                can->SaveAs(Form("ng/ng_Np_adj%i_%s_ptt%i.pdf",iadj,GAT::spec[itype],iptt));

                //plot a_g too. why not?
                vector<double> vag (vng);
                transform(vag.begin(),vag.end(),vag.begin(), [](double x) { return x / (2.*PI); } );
                vector<double> vage (vnge);
                transform(vage.begin(),vage.end(),vage.begin(), [](double x) { return x / (2.*PI); } );
                // tgrag_Np->SetTitle(Form("ag_Np_adj%i_%s_ptt%i;N_{part};<a_{#gamma}>",iadj,GAT::spec[itype],iptt));
                TGraphErrors* tgrag_Np = nullptr;
                MakeTGraphErrors(tgrag_Np,"ag","ag_Np_"+sufg.str(),"N_{part}","#LTa_{#gamma}#GT",vag.size(),&vnpart[0],&vag[0],&vnparte[0],&vage[0],can);

                TGraphErrors* tgrage = nullptr;
                MakeTGraphErrors(tgrage,"age","age_Np_"+sufg.str(),"N_{part}","#LTa_{#gamma}#GT error",vag.size(),&vnpart[0],&vage[0],&vnparte[0],0,can);

                if (quicksavetrigs)
                {
                    cout<<"WARNING. Saving trigger counts only."<<endl;
                    continue;
                }
                else
                {
                    cout<<"quicksavetrigs bool set to false"<<endl;
                }

                input->cd();

                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    TLatex* txtxi = new TLatex (0.5,.775,Form("%.1f #leq #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]));
                    txtxi->SetTextAlign(22);
                    txtxi->SetNDC();
                    txtxi->SetTextSize(0.05);               

                    ostringstream suf;
                    suf<<"adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_xi"<<ixi;

                    vector<double> vah;//<a_h>
                    vector<double> vahe;//error on <a_h>

                    vector<double> vagah;//<a_g><a_h>
                    vector<double> vagahe;//error on <a_g><a_h>

                    vector<double> vbgs;//<b_g #sigma>
                    vector<double> vbgse;//error on <b_g #sigma>

                    vector<double> vbg;//<b_g>
                    vector<double> vbge;//error on <b_g>

                    vector<double> vvg;//<v_g>
                    vector<double> vvge;//error on <v_g>

                    vector<double> vbhs;//<b_h #sigma>
                    vector<double> vbhse;//error on <b_h #sigma>

                    vector<double> vbh;//<b_h>
                    vector<double> vbhe;//error on <b_h>

                    vector<double> vvh;//<b_h>
                    vector<double> vvhe;//error on <b_h>

                    vector<double> vbgsbhs;//<b_g #sigma><b_h #sigma>
                    vector<double> vbgsbhse;//error on <b_g #sigma><b_h #sigma>

                    vector<double> vahbgs;//<a_h><b_g #sigma>
                    vector<double> vahbgse;//error on <a_h><b_g #sigma>

                    vector<double> vagbhs;//<a_g><b_h #sigma>
                    vector<double> vagbhse;//error on <a_g><b_h #sigma>

                    for (int icent = 0; icent < GAT::ncent; icent++)
                    {
                        // =============================================
                        // ========== Get E_all E_in and E_out =========
                        // =============================================

                        input->cd();

                        // oss_Eall<<"_mix"<<1;//mixed

                        // ========== E_in ==========

                        ostringstream oss_Ein;
                        oss_Ein<<"dphi";
                        if (iadj)
                            oss_Ein<<"_adj";
                        oss_Ein<<"_"<<GAT::spec[itype];
                        oss_Ein<<"_mix"<<1;//mixed
                        oss_Ein<<"_ptt"<<iptt;
                        oss_Ein<<"_cent"<<icent;
                        oss_Ein<<"_dpsi"<<0;//idpsi=0 is in-plane
                        oss_Ein<<"_xi"<<ixi;

                        TH1F* h1Ein = (TH1F*) gROOT->FindObject(oss_Ein.str().c_str());
                        assert(h1Ein);

                        //scale by nevt in cent range!
                        h1Ein->Scale(1./vnevt[icent]);

                        //scale by binwidth (times pi)
                        h1Ein->Scale(1./(h1Ein->GetBinWidth(1)*PI));

                        // ========== E_out ==========

                        ostringstream oss_Eout;
                        oss_Eout<<"dphi";
                        if (iadj)
                            oss_Eout<<"_adj";
                        oss_Eout<<"_"<<GAT::spec[itype];
                        oss_Eout<<"_mix"<<1;//mixed
                        oss_Eout<<"_ptt"<<iptt;
                        oss_Eout<<"_cent"<<icent;
                        oss_Eout<<"_dpsi"<<1;//idpsi=0 is out-plane
                        oss_Eout<<"_xi"<<ixi;

                        TH1F* h1Eout = (TH1F*) gROOT->FindObject(oss_Eout.str().c_str());
                        assert(h1Eout);

                        //scale by nevt in cent range!
                        h1Eout->Scale(1./vnevt[icent]);

                        //scale by binwidth (times pi)
                        h1Eout->Scale(1./(h1Eout->GetBinWidth(1)*PI));

                        // ========== E_all ==========

                        TH1F* h1Eall = (TH1F*) h1Ein->Clone(Form("Eall_adj%i_%s_ptt%i_cent%i_xi%i",iadj,GAT::spec[itype],iptt,icent,ixi));
                        h1Eall->Add(h1Eout,1.0);
                        h1Eall->SetDirectory(output);

                        //==================== <a_h> ====================

                        double error_Eall = -9999.;
                        double integral_Eall = h1Eall->IntegralAndError(1,h1Eall->GetNbinsX(),error_Eall,"width")*PI;//width option multiplies by bin width.  the pi has been dropped and must be included manually
                        error_Eall*=PI;

                        //note: chatter scaled by PI to match function parameters
                        cout<<h1Eall->GetName()<<" integral = "<<integral_Eall/PI<<" error = "<<error_Eall/PI<<endl;

                        // note: N_g = n_g * N_evt
                        double ag = vag[icent];
                        double age = vage[icent];

                        //integral of Eall = pi * (4 pi) <a_g><a_h>
                        double ah = integral_Eall/(PI *(4.*PI)*ag);

                        double ahe = -9999.;
                        {
                            double denom = PI*(4.*PI)*ag;
                            double e0 = error_Eall/denom;
                            double e1 = -integral_Eall*PI*(4.*PI)*age;
                            e1 /= denom*denom;
                            ahe = sqrt(e0*e0 + e1*e1);
                        }                                         

                        vah.push_back(ah);
                        vahe.push_back(ahe);

                        //save full term...
                        double agah = integral_Eall/(PI*(4.*PI));
                        double agahe = error_Eall/(PI*(4.*PI));

                        vagah.push_back(agah);
                        vagahe.push_back(agahe);

                        //==================== <b_g s> <b_h s> ====================

                        cout<<"Fitting "<<h1Eall->GetName()<<endl;

                        h1Eall->SetTitle(";#Delta#phi/#pi;n_{t}Y^{#psi}_{mix}");

                        h1Eall->SetMarkerStyle(kFullCircle);
                        h1Eall->SetMarkerColor(kBlue);
                        h1Eall->SetMarkerSize(6);

                        TF1* f1all = new TF1(Form("fit_%s",h1Eall->GetName()),"[0]+[1]*cos(x*2.*TMath::Pi())",0.,1.);
                        f1all->FixParameter(0,integral_Eall/PI);
                        f1all->SetParName(0,"4#pi#LTa_{#gamma}>#LTa_{h}#GT ");
                        f1all->SetParName(1,"2#pi#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT ");
                        f1all->SetLineWidth(4);

                        h1Eall->Draw("E1");
                        can->SetLeftMargin(.15);
                        can->SetTopMargin(.15);
                        can->SetRightMargin(.1);
                        can->SetBottomMargin(.125);
                        h1Eall->GetYaxis()->SetTitleOffset(1.3);                        
                        h1Eall->GetYaxis()->SetMaxDigits(2);

                        h1Eall->GetYaxis()->SetTitleSize(.06);
                        h1Eall->GetYaxis()->SetLabelSize(.05);
                        h1Eall->GetXaxis()->SetTitleSize(.06);
                        h1Eall->GetXaxis()->SetLabelSize(.05);
                        h1Eall->GetXaxis()->SetTitleOffset(1.0);

                        h1Eall->Fit(f1all);

                        f1all->SetParError(0,error_Eall/PI);

                        //introducing a new error calculation
                        double coserr_all = GetError(h1Eall,f1all);             
                        DrawError(f1all,coserr_all,1);

                        f1all->SetParError(1,coserr_all);

                        TLine* tl = new TLine(0.,f1all->GetParameter(0),1.,f1all->GetParameter(0));
                        tl->SetLineStyle(9);
                        tl->SetLineWidth(3);
                        tl->Draw("same");
                        can->Modified();
                        can->Update();

                        TPaveText* pt = new TPaveText(.305,.85,.90,.99,"NDC");
                        pt->SetTextSize(.055);
                        pt->SetBorderSize(1);
                        pt->SetTextAlign(32);
                        pt->AddText(Form("4#pi#LTa_{#gamma}#GT#LTa_{h}#GT = %.3e #pm %.1e",f1all->GetParameter(0),f1all->GetParError(0)));
                        pt->AddText(Form("2#pi#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT = %.3e #pm %.2e",f1all->GetParameter(1),f1all->GetParError(1)));
                        pt->SetFillColorAlpha(0,0);
                        pt->Draw();

                        TLatex txtcent (0.15,.9375,Form("(%i-%i%%)",GAT::limcent[icent],GAT::limcent[icent+1]));
                        txtcent.SetTextAlign(12);
                        txtcent.SetNDC();
                        txtcent.SetTextSize(0.055);
                        txtcent.Draw();

                        //bring histo to front
                        h1Eall->Draw("E1 same");

                        can->SaveAs(Form("bgsbhs/bgsbhs_adj%i_%s_ptt%i_cent%i_xi%i.pdf",iadj,GAT::spec[itype],iptt,icent,ixi));

                        double bgsbhs = -9999., bgsbhse = -9999.;
                        bgsbhs = f1all->GetParameter(1)/(2.*PI);
                        bgsbhse = f1all->GetParError(1)/(2.*PI);

                        vbgsbhs.push_back(bgsbhs);
                        vbgsbhse.push_back(bgsbhse);

                        //==================== <b_g s> ====================

                        TH1F* h1Ediff = (TH1F*) h1Ein->Clone(Form("Ediff_adj%i_%s_ptt%i_cent%i_xi%i",iadj,GAT::spec[itype],iptt,icent,ixi));
                        h1Ediff->Add(h1Eout,-1.0);
                        h1Ediff->SetDirectory(output);

                        double error_Ediff = -9999.;
                        double integral_Ediff = h1Ediff->IntegralAndError(1,h1Ediff->GetNbinsX(),error_Ediff,"width")*PI;
                        error_Ediff*=PI;

                        cout<<h1Ediff->GetName()<<" integral = "<<integral_Ediff/PI<<" error = "<<error_Ediff/PI<<endl;

                        //WARING: this pi should not be here no? -> it's from integral
                        double ahbgs = integral_Ediff/(PI*8.);
                        double ahbgse = error_Ediff/(PI*8.);

                        vahbgs.push_back(ahbgs);
                        vahbgse.push_back(ahbgse);

                        double bgs = integral_Ediff/(PI*8.*ah);

                        double bgse = -9999.;
                        {
                            double denom = 8.*PI*ah;
                            double e0 = error_Ediff/denom;
                            double e1 = -integral_Ediff*8.*PI*ahe;
                            e1 /= denom*denom;
                            bgse = sqrt(e0*e0 + e1*e1);
                        }

                        vbgs.push_back(bgs);
                        vbgse.push_back(bgse);

                        double bg = bgs/vsig[icent];
                        double bge = bgse/vsig[icent];

                        vbg.push_back(bg);
                        vbge.push_back(bge);

                        double vg = bg/(2.*ag);
                        double vge = -9999.;
                        {
                            double e0 = bge/(2.*ag);
                            double e1 = -(bg*age)/(2.*ag*ag);
                            vge = sqrt(e0*e0 + e1*e1);
                        }

                        vvg.push_back(vg);
                        vvge.push_back(vge);

                        //==================== <b_h s> ====================

                        cout<<"Fitting "<<h1Ediff->GetName()<<endl;

                        h1Ediff->SetTitle(";#Delta#phi/#pi;n^{+}_{t}Y^{#psi+}_{mix} - n^{-}_{t}Y^{#psi-}_{mix}");

                        h1Ediff->SetMarkerStyle(kFullCircle);
                        h1Ediff->SetMarkerColor(kBlue);
                        h1Ediff->SetMarkerSize(6);

                        TF1* f1diff = new TF1(Form("bhscos_adj%i_%s_ptt%i_cent%i_xi%i",iadj,GAT::spec[itype],iptt,icent,ixi),"[0]+[1]*cos(x*2.*TMath::Pi())",0.,1.);
                        f1diff->FixParameter(0,8.*ahbgs);
                        f1diff->SetParName(0,"8#LTa_{h}#GT#LTb_{#gamma}#sigma#GT ");
                        f1diff->SetParName(1,"8#LTa_{#gamma}#GT#LTb_{h}#sigma#GT ");
                        f1diff->SetLineWidth(4);

                        h1Ediff->Draw("E1");
                        can->SetLeftMargin(.15);
                        can->SetTopMargin(.15);
                        can->SetRightMargin(.1);
                        can->SetBottomMargin(.125);
                        h1Ediff->GetYaxis()->SetTitleOffset(1.3);
                        h1Ediff->GetYaxis()->SetMaxDigits(2);

                        h1Ediff->GetYaxis()->SetTitleSize(.06);
                        h1Ediff->GetYaxis()->SetLabelSize(.05);
                        h1Ediff->GetXaxis()->SetTitleSize(.06);
                        h1Ediff->GetXaxis()->SetLabelSize(.05);
                        h1Ediff->GetXaxis()->SetTitleOffset(1.0);

                        h1Ediff->Fit(f1diff);

                        f1diff->SetParError(0,error_Ediff/PI);

                        //introducing a new error calculation
                        double coserr_diff = GetError(h1Ediff,f1diff);
                        DrawError(f1diff,coserr_diff,1);

                        f1diff->SetParError(1,coserr_diff);

                        tl = new TLine(0.,8.*ahbgs,1.,8.*ahbgs);
                        tl->SetLineStyle(9);
                        tl->SetLineWidth(3);
                        tl->Draw("same");
                        can->Modified();
                        can->Update();

                        pt->Clear();
                        pt->AddText(Form("8#LTa_{h}#GT#LTb_{#gamma}#sigma#GT = %.3e #pm %.1e",f1diff->GetParameter(0),f1diff->GetParError(0)));
                        pt->AddText(Form("8#LTa_{#gamma}#GT#LTb_{h}#sigma#GT = %.3e #pm %.2e",f1diff->GetParameter(1),f1diff->GetParError(1)));
                        pt->Draw();

                        txtcent.Draw();

                        //bring histo to front                  
                        h1Ediff->Draw("E1 same");

                        can->SaveAs(Form("bhscos/bhscos_adj%i_%s_ptt%i_cent%i_xi%i.pdf",iadj,GAT::spec[itype],iptt,icent,ixi));

                        double bhs = -9999., bhse = -9999.;
                        double agbhs = -9999., agbhse = -9999.;
                        {
                            double denom = 8.*ag;

                            double par = f1diff->GetParameter(1);
                            double parerr = f1diff->GetParError(1);

                            agbhs = par/8.;
                            agbhse = parerr/8.;

                            bhs = par/denom;

                            double e0 = parerr/denom;
                            double e1 = -(par*8.*age)/(denom*denom);

                            bhse = sqrt(e0*e0 + e1*e1);
                        }

                        vagbhs.push_back(agbhs);
                        vagbhse.push_back(agbhse);

                        vbhs.push_back(bhs);
                        vbhse.push_back(bhse);

                        double bh = bhs/vsig[icent];
                        double bhe = bhse/vsig[icent];

                        vbh.push_back(bh);
                        vbhe.push_back(bhe);

                        double vh = bh/(2.*ah);
                        double vhe = -9999.;
                        {
                            double e0 = bhe/(2.*ah);
                            double e1 = -(bh*ahe)/(2.*ah*ah);
                            vhe = sqrt(e0*e0 + e1*e1);
                        }

                        vvh.push_back(vh);
                        vvhe.push_back(vhe);
                    }

                    output->cd();

                    //==================== <a_h> ====================

                    TGraph* tgrah = new TGraph(vah.size(),&vcent[0],&vah[0]);
                    tgrah->SetTitle(Form("ah_adj%i_%s_ptt%i_xi%i;centrality;<a_h}#GT",iadj,GAT::spec[itype],iptt,ixi));
                    tgrah->Write(Form("ah_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));

                    TGraphErrors* tgrah_Np = nullptr;
                    MakeTGraphErrors(tgrah_Np,"ah","ah_Np_"+suf.str(),"N_{part}","#LTa_{h}#GT",vah.size(),&vnpart[0],&vah[0],&vnparte[0],&vahe[0],can);

                    TGraphErrors* tgrahe = nullptr;
                    MakeTGraphErrors(tgrahe,"ahe","ahe_Np_"+suf.str(),"N_{part}","#LTa_{h}#GT error",vah.size(),&vnpart[0],&vahe[0],&vnparte[0],0,can);

                    //==================== <a_g><a_h> ====================

                    TGraphErrors* tgragah = nullptr;
                    MakeTGraphErrors(tgragah,"agah","agah_Np_"+suf.str(),"N_{part}","#LTa_{#gamma}#GT#LTa_{h}#GT",vagah.size(),&vnpart[0],&vagah[0],&vnparte[0],&vagahe[0],can);

                    TGraphErrors* tgragahe = nullptr;
                    MakeTGraphErrors(tgragahe,"agahe","agahe_Np_"+suf.str(),"N_{part}","#LTa_{h}#GT error",vagah.size(),&vnpart[0],&vagahe[0],&vnparte[0],0,can);                   

                    //==================== <b_g sigma> ====================

                    TGraph* tgrbgs = new TGraph(vbgs.size(),&vcent[0],&vbgs[0]);
                    tgrbgs->SetTitle(Form("bgs_adj%i_%s_ptt%i_xi%i;centrality;#LTb_{#gamma} #sigma#GT",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbgs->SetName(Form("bgs_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbgs->Write();

                    TGraphErrors* tgrbgs_Np = nullptr;
                    MakeTGraphErrors(tgrbgs_Np,"bgs","bgs_Np_"+suf.str(),"N_{part}","#LTb_{#gamma} #sigma#GT",vbgs.size(),&vnpart[0],&vbgs[0],&vnparte[0],&vbgse[0],can);

                    TGraphErrors* tgrbgse = nullptr;
                    MakeTGraphErrors(tgrbgse,"bgse","bgse_Np_"+suf.str(),"N_{part}","#LTb_{#gamma} #sigma#GT error",vbgse.size(),&vnpart[0],&vbgse[0],&vnparte[0],0,can);

                    //==================== <b_g> ====================

                    TGraph* tgrbg = new TGraph(vbg.size(),&vcent[0],&vbg[0]);
                    tgrbg->SetTitle(Form("bg_adj%i_%s_ptt%i_xi%i;centrality;#LTb_{#gamma} #sigma#GT/#LT#sigma#GT",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbg->SetName(Form("bg_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbg->Write();

                    TGraphErrors* tgrbg_Np = nullptr;
                    MakeTGraphErrors(tgrbg_Np,"bg","bg_Np_"+suf.str(),"N_{part}","#LTb_{#gamma}#GT",vbg.size(),&vnpart[0],&vbg[0],&vnparte[0],&vbge[0],can);

                    TGraphErrors* tgrbge = nullptr;
                    MakeTGraphErrors(tgrbge,"bge","bge_Np_"+suf.str(),"N_{part}","#LTb_{#gamma}#GT error",vbge.size(),&vnpart[0],&vbge[0],&vnparte[0],0,can);

                    //==================== <v_g> ====================

                    TGraph* tgrvg = new TGraph(vvg.size(),&vcent[0],&vvg[0]);
                    tgrvg->SetTitle(Form("vg_adj%i_%s_ptt%i_xi%i;centrality;#frac{#LTb_{#gamma}#sigma#GT}{#LT#sigma#GT#LTa_{#gamma}#GT} #approx v_{#gamma}",iadj,GAT::spec[itype],iptt,ixi));
                    tgrvg->SetName(Form("vg_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));
                    tgrvg->Write();

                    TGraphErrors* tgrvg_Np = nullptr;
                    MakeTGraphErrors(tgrvg_Np,"vg","vg_Np_"+suf.str(),"N_{part}","#LTv_{#gamma}#GT",vvg.size(),&vnpart[0],&vvg[0],&vnparte[0],&vvge[0],can);

                    //==================== <b_h sigma> ====================

                    TGraph* tgrbhs = new TGraph(vbhs.size(),&vcent[0],&vbhs[0]);
                    tgrbhs->SetTitle(Form("bhs_adj%i_%s_ptt%i_xi%i;centrality;#LTb_{h} #sigma#GT",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbhs->Write(Form("bhs_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));

                    TGraphErrors* tgrbhs_Np = nullptr;
                    MakeTGraphErrors(tgrbhs_Np,"bhs","bhs_Np_"+suf.str(),"N_{part}","#LTb_{h} #sigma#GT",vbhs.size(),&vnpart[0],&vbhs[0],&vnparte[0],&vbhse[0],can);

                    TGraphErrors* tgrbhse = nullptr;
                    MakeTGraphErrors(tgrbhse,"bhse","bhse_Np_"+suf.str(),"N_{part}","#LTb_{h} #sigma#GT error",vbhse.size(),&vnpart[0],&vbhse[0],&vnparte[0],0,can);

                    //==================== <b_h> ====================

                    TGraph* tgrbh = new TGraph(vbh.size(),&vcent[0],&vbh[0]);
                    tgrbh->SetTitle(Form("bh_adj%i_%s_ptt%i_xi%i;centrality;#LTb_{h} #sigma#GT/#LT#sigma#GT",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbh->Write(Form("bh_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));

                    TGraphErrors* tgrbh_Np = nullptr;
                    MakeTGraphErrors(tgrbh_Np,"bh","bh_Np_"+suf.str(),"N_{part}","#LTb_{h}#GT",vbh.size(),&vnpart[0],&vbh[0],&vnparte[0],&vbhe[0],can);

                    TGraphErrors* tgrbhe = nullptr;
                    MakeTGraphErrors(tgrbhe,"bhe","bhe_Np_"+suf.str(),"N_{part}","#LTb_{h}#GT error",vbhe.size(),&vnpart[0],&vbhe[0],&vnparte[0],0,can);

                    //==================== <v_h> ====================

                    TGraph* tgrvh = new TGraph(vvh.size(),&vcent[0],&vvh[0]);
                    tgrvh->SetTitle(Form("vh_adj%i_%s_ptt%i_xi%i;centrality;#frac{#LTb_{h}#sigma#GT}{#LT#sigma#GT#LTa_{h}#GT} #approx v_{h}",iadj,GAT::spec[itype],iptt,ixi));
                    tgrvh->Write(Form("vh_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));

                    TGraphErrors* tgrvh_Np = nullptr;
                    MakeTGraphErrors(tgrvh_Np,"vh","vh_Np_"+suf.str(),"N_{part}","#LTv_{h}#GT",vvh.size(),&vnpart[0],&vvh[0],&vnparte[0],&vvhe[0],can);

                    //==================== <b_g sigma> <b_h sigma> ====================

                    TGraph* tgrbgsbhs = new TGraph(vbgsbhs.size(),&vcent[0],&vbgsbhs[0]);
                    tgrbgsbhs->SetTitle(Form("bgsbhs_adj%i_%s_ptt%i_xi%i;centrality;#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbgsbhs->Write(Form("bgsbhs_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));

                    TGraphErrors* tgrbgsbhs_Np = nullptr;
                    MakeTGraphErrors(tgrbgsbhs_Np,"bgsbhs","bgsbhs_Np_"+suf.str(),"N_{part}","#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT",vbgsbhs.size(),&vnpart[0],&vbgsbhs[0],&vnparte[0],&vbgsbhse[0],can);

                    TGraphErrors* tgrbgsbhse = nullptr;
                    MakeTGraphErrors(tgrbgsbhse,"bgsbhse","bgsbhse_Np_"+suf.str(),"N_{part}","#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT error",vbgsbhse.size(),&vnpart[0],&vbgsbhse[0],&vnparte[0],0,can);

                    //==================== <a_h> <b_g sigma> ====================

                    TGraphErrors* tgrahbgs = nullptr;
                    MakeTGraphErrors(tgrahbgs,"ahbgs","ahbgs_Np_"+suf.str(),"N_{part}","#LTa_{h}#GT#LTb_{#gamma}#sigma#GT",vahbgs.size(),&vnpart[0],&vahbgs[0],&vnparte[0],&vahbgse[0],can);

                    TGraphErrors* tgrahbgse = nullptr;
                    MakeTGraphErrors(tgrahbgse,"ahbgse","ahbgse_Np_"+suf.str(),"N_{part}","#LTa_{h}#GT#LTb_{#gamma}#sigma#GT error",vahbgse.size(),&vnpart[0],&vahbgse[0],&vnparte[0],0,can);

                    //==================== <a_g> <b_h sigma> ====================

                    TGraphErrors* tgragbhs = nullptr;
                    MakeTGraphErrors(tgragbhs,"agbhs","agbhs_Np_"+suf.str(),"N_{part}","#LTa_{#gamma}#GT#LTb_{h}#sigma#GT",vagbhs.size(),&vnpart[0],&vagbhs[0],&vnparte[0],&vagbhse[0],can);

                    TGraphErrors* tgragbhse = nullptr;
                    MakeTGraphErrors(tgragbhse,"agbhse","agbhse_Np_"+suf.str(),"N_{part}","#LTa_{#gamma}#GT#LTb_{h}#sigma#GT error",vagbhse.size(),&vnpart[0],&vagbhse[0],&vnparte[0],0,can);

                    //==================================================
                    // Calculate K Corrections with Obtained Components
                    //==================================================

                    //-------------------- k_a --------------------

                    cout<<"Calculating "<<"ka_Np_"<<suf.str()<<endl;

                    vector<double> vka, vkae, vkanum, vkanume, vkaden, vkadene;
                    CalcK(hnpart,tgrag_Np,tgrah_Np,tgrage,tgrahe,vka,vkae,vkanum,vkanume,vkaden,vkadene);

                    TGraphErrors* tgrka = nullptr;
                    MakeTGraphErrors(tgrka,"ka","ka_Np_"+suf.str(),"N_{part}","#frac{#LTa_{#gamma} a_{h}#GT}{#LTa_{#gamma}#GT#LTa_{h}#GT}",vka.size(),&vnpart[0],&vka[0],&vnparte[0],&vkae[0],can);

                    //----- make ka with monte carlo errors -----

                    TMultiGraph* mgrka = new TMultiGraph();
                    TMultiGraph* mgrag = new TMultiGraph();
                    TMultiGraph* mgrah = new TMultiGraph();
                    vector<double> vkaenew = CalcKMC(tgrag_Np,tgrah_Np,vka,vkae,vnpart,hnpart,mgrka,mgrag,mgrah);

                    TGraphErrors* tgrkanew = nullptr;
                    MakeTGraphErrors(tgrkanew,"kanew","kanew_Np_"+suf.str(),"N_{part}","#frac{#LTa_{#gamma} a_{h}#GT}{#LTa_{#gamma}#GT#LTa_{h}#GT}",vka.size(),&vnpart[0],&vka[0],&vnparte[0],&vkaenew[0],can);

                    cmgr->cd();

                    mgrka->Add(tgrkanew,"p");
                    mgrka->SetTitle(";N_{part};#kappa_{a}");
                    mgrka->Draw("a");
                    if (itype == GAT::PI0)
                        mgrka->GetYaxis()->SetRangeUser(kalim.at(iptt).first,kalim.at(iptt).second);
                    mgrka->GetYaxis()->SetLabelOffset(.01);
                    mgrka->GetXaxis()->SetLabelOffset(.01);
                    mgrka->GetYaxis()->SetTitleOffset(1.5);
                    mgrka->GetXaxis()->SetTitleOffset(1.1);
                    txtmom->Draw();
                    txtxi->Draw();
                    cmgr->SaveAs(Form("kanew/multigraph_ka_%s.pdf",suf.str().c_str()));

                    mgrag->SetTitle(";N_{part};a_{#gamma}");
                    mgrag->Draw("a");
                    cmgr->SaveAs(Form("kanew/multigraph_ag_%s.pdf",suf.str().c_str()));

                    mgrah->SetTitle(";N_{part};a_{h}");
                    mgrah->Draw("a");
                    cmgr->SaveAs(Form("kanew/multigraph_ah_%s.pdf",suf.str().c_str()));
                    //---------- end ka monte carlo ----------

                    TGraphErrors* tgrkanum = nullptr;
                    MakeTGraphErrors(tgrkanum,"kanum","kanum_Np_"+suf.str(),"N_{part}","#LTa_{#gamma} a_{h}#GT",vkanum.size(),&vnpart[0],&vkanum[0],&vnparte[0],&vkanume[0],can);

                    TGraphErrors* tgrkaden = nullptr;
                    MakeTGraphErrors(tgrkaden,"kaden","kaden_Np_"+suf.str(),"N_{part}","#LTa_{#gamma}#GT#LTa_{h}#GT",vkaden.size(),&vnpart[0],&vkaden[0],&vnparte[0],&vkadene[0],can);

                    //-------------------- k_a merge --------------------

                    cout<<"Calculating "<<"kamerge_Np_"<<suf.str()<<endl;

                    vector<double> vkamerge, vkamergee, vkamergenum, vkamergenume, vkamergeden, vkamergedene;
                    CalcK_Merge(hnpartmerge,tgrag_Np,tgrah_Np,tgrage,tgrahe,vkamerge,vkamergee,vkamergenum,vkamergenume,vkamergeden,vkamergedene);

                    TGraphErrors* tgrkamerge = nullptr;
                    MakeTGraphErrors(tgrkamerge,"kamerge","kamerge_Np_"+suf.str(),"N_{part}","#frac{#LTa_{#gamma} a_{h}#GT}{#LTa_{#gamma}#GT#LTa_{h}#GT}",vkamerge.size(),&vnpartmerge[0],&vkamerge[0],&vnpartmergee[0],&vkamergee[0],can);

                    TGraphErrors* tgrkamergenum = nullptr;
                    MakeTGraphErrors(tgrkamergenum,"kamergenum","kamergenum_Np_"+suf.str(),"N_{part}","#LTa_{#gamma} a_{h}#GT",vkamergenum.size(),&vnpartmerge[0],&vkamergenum[0],&vnpartmergee[0],&vkamergenume[0],can);

                    TGraphErrors* tgrkamergeden = nullptr;
                    MakeTGraphErrors(tgrkamergeden,"kamergeden","kamergeden_Np_"+suf.str(),"N_{part}","#LTa_{#gamma}#GT#LTa_{h}#GT",vkamergeden.size(),&vnpartmerge[0],&vkamergeden[0],&vnpartmergee[0],&vkamergedene[0],can);              

                    //-------------------- k_b --------------------

                    cout<<"Calculating "<<"kb_Np_"<<suf.str()<<endl;

                    //before calculating k, create new functions in denominator
                    TGraphErrors* tgrbgs_new = new TGraphErrors(GAT::ncent+2);
                    tgrbgs_new->SetPoint(0,nparthi,tgrbg_Np->Eval(nparthi)*tgrsig->Eval(nparthi));
                    tgrbgs_new->SetPointError(0,tgrbg_Np->GetErrorX(0),tgrbg_Np->GetErrorY(0)*tgrsig->Eval(nparthi));
                    for (int icent = 0; icent < GAT::ncent; icent++)
                    {
                        tgrbgs_new->SetPoint(icent+1,vnpart[icent],vbgs[icent]);
                        tgrbgs_new->SetPointError(icent+1,vnparte[icent],vbgse[icent]);
                    }
                    tgrbgs_new->SetPoint(GAT::ncent+1,npartlo,tgrbg_Np->Eval(npartlo)*tgrsig->Eval(npartlo));
                    tgrbgs_new->SetPointError(GAT::ncent+1,tgrbg_Np->GetErrorX(GAT::ncent+1),tgrbg_Np->GetErrorY(GAT::ncent+1)*tgrsig->Eval(npartlo));

                    string name_bgs = "bgs_new_"+suf.str();
                    tgrbgs_new->SetTitle(Form("%s;N_{part};#LTb_{#gamma} #sigma#GT",name_bgs.c_str()));
                    tgrbgs_new->SetName(name_bgs.c_str());
                    cout<<"Made TgraphErrors "<<tgrbgs_new->GetName()<<endl;
                    tgrbgs_new->Write();
                    tgrbgs_new->Draw("AL");
                    can->SaveAs(Form("bgsnew/%s.pdf",name_bgs.c_str()));

                    TGraphErrors* tgrbhs_new = new TGraphErrors(GAT::ncent+2);
                    tgrbhs_new->SetPoint(0,nparthi,tgrbh_Np->Eval(nparthi)*tgrsig->Eval(nparthi));
                    tgrbhs_new->SetPointError(0,tgrbh_Np->GetErrorX(0),tgrbh_Np->GetErrorY(0)*tgrsig->Eval(nparthi));
                    for (int icent = 0; icent < GAT::ncent; icent++)
                    {
                        tgrbhs_new->SetPoint(icent+1,vnpart[icent],vbhs[icent]);
                        tgrbhs_new->SetPointError(icent+1,vnparte[icent],vbhse[icent]);
                    }
                    tgrbhs_new->SetPoint(GAT::ncent+1,npartlo,tgrbh_Np->Eval(npartlo)*tgrsig->Eval(npartlo));
                    tgrbhs_new->SetPointError(GAT::ncent+1,tgrbh_Np->GetErrorX(GAT::ncent),tgrbh_Np->GetErrorY(GAT::ncent)*tgrsig->Eval(npartlo));

                    string name_bhs = "bhs_new_"+suf.str();
                    tgrbhs_new->SetTitle(Form("%s;N_{part};#LTb_{#gamma} #sigma#GT",name_bhs.c_str()));
                    tgrbhs_new->SetName(name_bhs.c_str());
                    cout<<"Made TgraphErrors "<<tgrbhs_new->GetName()<<endl;
                    tgrbhs_new->Write();
                    tgrbhs_new->Draw("AL");
                    can->SaveAs(Form("bhsnew/%s.pdf",name_bhs.c_str()));

                    vector<double> vkb, vkbe, vkbnum, vkbnume, vkbden, vkbdene;
                    CalcK(hnpart,tgrbg_Np,tgrbh_Np,tgrbge,tgrbhe,vkb,vkbe,vkbnum,vkbnume,vkbden,vkbdene,tgrbgs_new,tgrbhs_new,tgrbgse,tgrbhse);

                    TGraphErrors* tgrkb = nullptr;
                    MakeTGraphErrors(tgrkb,"kb","kb_Np_"+suf.str(),"N_{part}","#frac{#LTb_{#gamma} b_{h}#GT}{#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT}",vkb.size(),&vnpart[0],&vkb[0],&vnparte[0],&vkbe[0],can);

                    //----- make kb with monte carlo errors -----
                    TMultiGraph* mgrkb = new TMultiGraph();
                    TMultiGraph* mgrbg = new TMultiGraph();
                    TMultiGraph* mgrbh = new TMultiGraph();
                    TMultiGraph* mgrbgs = new TMultiGraph();
                    TMultiGraph* mgrbhs = new TMultiGraph();
                    vector<double> vkbenew = CalcKMC(tgrbg_Np,tgrbh_Np,vkb,vkbe,vnpart,hnpart,mgrkb,mgrbg,mgrbh,tgrbgs_Np,tgrbhs_Np,mgrbgs,mgrbhs);

                    TGraphErrors* tgrkbnew = nullptr;
                    MakeTGraphErrors(tgrkbnew,"kbnew","kbnew_Np_"+suf.str(),"N_{part}","#frac{#LTb_{#gamma} b_{h}#GT}{#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT}",vkb.size(),&vnpart[0],&vkb[0],&vnparte[0],&vkbenew[0],can);

                    cmgr->cd();

                    mgrkb->Add(tgrkbnew,"p");
                    mgrkb->SetTitle(";N_{part};#kappa_{b}/#LT#sigma#GT^{2}");
                    mgrkb->Draw("a");
                    if (itype == GAT::PI0)
                        mgrkb->GetYaxis()->SetRangeUser(kblim.at(iptt).first,kblim.at(iptt).second);
                    mgrkb->GetYaxis()->SetLabelOffset(.01);
                    mgrkb->GetXaxis()->SetLabelOffset(.01);
                    mgrkb->GetYaxis()->SetTitleOffset(1.5);
                    mgrkb->GetXaxis()->SetTitleOffset(1.1);

                    txtmom->Draw();
                    txtxi->Draw();
                    cmgr->SaveAs(Form("kbnew/multigraph_kb_%s.pdf",suf.str().c_str()));

                    mgrbg->SetTitle(";N_{part};b_{#gamma}");
                    mgrbg->Draw("a");
                    cmgr->SaveAs(Form("kbnew/multigraph_bg_%s.pdf",suf.str().c_str()));

                    mgrbh->SetTitle(";N_{part};b_{h}");
                    mgrbh->Draw("a");
                    cmgr->SaveAs(Form("kbnew/multigraph_bh_%s.pdf",suf.str().c_str()));

                    mgrbgs->SetTitle(";N_{part};b_{#gamma} #sigma");
                    mgrbgs->Draw("a");
                    cmgr->SaveAs(Form("kbnew/multigraph_bgs_%s.pdf",suf.str().c_str()));

                    mgrbhs->SetTitle(";N_{part};b_{h} #sigma");
                    mgrbhs->Draw("a");
                    cmgr->SaveAs(Form("kbnew/multigraph_bhs_%s.pdf",suf.str().c_str()));
                    //---------- end kb monte carlo ----------

                    TGraphErrors* tgrkbnum = nullptr;
                    MakeTGraphErrors(tgrkbnum,"kbnum","kbnum_Np_"+suf.str(),"N_{part}","#LTb_{#gamma} b_{h}#GT",vkbnum.size(),&vnpart[0],&vkbnum[0],&vnparte[0],&vkbnume[0],can);

                    TGraphErrors* tgrkbden = nullptr;
                    MakeTGraphErrors(tgrkbden,"kbden","kbden_Np_"+suf.str(),"N_{part}","#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT",vkbden.size(),&vnpart[0],&vkbden[0],&vnparte[0],&vkbdene[0],can);

                    //-------------------- k_b merge --------------------

                    cout<<"Calculating "<<"kbmerge_Np_"<<suf.str()<<endl;

                    vector<double> vkbmerge, vkbmergee, vkbmergenum, vkbmergenume, vkbmergeden, vkbmergedene;
                    CalcK_Merge(hnpartmerge,tgrbg_Np,tgrbh_Np,tgrbge,tgrbhe,vkbmerge,vkbmergee,vkbmergenum,vkbmergenume,vkbmergeden,vkbmergedene,tgrbgs_Np,tgrbhs_Np,tgrbgse,tgrbhse);

                    TGraphErrors* tgrkbmerge = nullptr;
                    MakeTGraphErrors(tgrkbmerge,"kbmerge","kbmerge_Np_"+suf.str(),"N_{part}","#frac{#LTb_{#gamma} b_{h}#GT}{#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT}",vkbmerge.size(),&vnpartmerge[0],&vkbmerge[0],&vnpartmergee[0],&vkbmergee[0],can);

                    TGraphErrors* tgrkbmergenum = nullptr;
                    MakeTGraphErrors(tgrkbmergenum,"kbmergenum","kbmergenum_Np_"+suf.str(),"N_{part}","#LTb_{#gamma} b_{h}#GT",vkbmergenum.size(),&vnpartmerge[0],&vkbmergenum[0],&vnpartmergee[0],&vkbmergenume[0],can);

                    TGraphErrors* tgrkbmergeden = nullptr;
                    MakeTGraphErrors(tgrkbmergeden,"kbmergeden","kbmergeden_Np_"+suf.str(),"N_{part}","#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT",vkbmergeden.size(),&vnpartmerge[0],&vkbmergeden[0],&vnpartmergee[0],&vkbmergedene[0],can);

                    //-------------------- k_b2 --------------------

                    //this kb includes sigma in num and denom

                    cout<<"Calculating "<<"kb2_Np_"<<suf.str()<<endl;

                    vector<double> vkb2, vkb2e, vkb2num, vkb2nume, vkb2den, vkb2dene;

                    for (int i = 0; i < vkb.size(); i++)
                    {
                        double kb2 = vkb[i];
                        kb2 *= vsig[i]*vsig[i];
                        vkb2.push_back(kb2);

                        double kb2e = vkbe[i];
                        kb2e *= vsig[i]*vsig[i];
                        vkb2e.push_back(kb2e);

                        double kb2num = vkbnum[i];
                        kb2num *= vsig[i]*vsig[i];
                        vkb2num.push_back(kb2num);

                        double kb2nume = vkbnume[i];
                        kb2nume *= vsig[i]*vsig[i];
                        vkb2nume.push_back(kb2nume);

                        double kb2den = vkbden[i];
                        vkb2den.push_back(kb2den);

                        double kb2dene = vkbdene[i];
                        vkb2dene.push_back(kb2dene);
                    }

                    TGraphErrors* tgrkb2 = nullptr;
                    MakeTGraphErrors(tgrkb2,"kb2","kb2_Np_"+suf.str(),"N_{part}","#frac{#LTb_{#gamma} b_{h}#GT#LT#sigma#GT^{2}}{#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT}",vkb2.size(),&vnpart[0],&vkb2[0],&vnparte[0],&vkb2e[0],can);

                    TGraphErrors* tgrkb2num = nullptr;
                    MakeTGraphErrors(tgrkb2num,"kb2num","kb2num_Np_"+suf.str(),"N_{part}","#LTb_{#gamma} b_{h}#GT#LT#sigma#GT^{2}",vkb2num.size(),&vnpart[0],&vkb2num[0],&vnparte[0],&vkb2nume[0],can);

                    TGraphErrors* tgrkb2den = nullptr;
                    MakeTGraphErrors(tgrkb2den,"kb2den","kb2den_Np_"+suf.str(),"N_{part}","#LTb_{#gamma}#sigma#GT#LTb_{h}#sigma#GT",vkb2den.size(),&vnpart[0],&vkb2den[0],&vnparte[0],&vkb2dene[0],can);

                    //-------------------- k_bg --------------------

                    cout<<"Calculating "<<"kbgs_Np_"<<suf.str()<<endl;
                    vector<double> vkbgs, vkbgse, vkbgsnum, vkbgsnume, vkbgsden, vkbgsdene;
                    CalcK(hnpart,tgrbg_Np,tgrsig,tgrbge,nullptr,vkbgs,vkbgse,vkbgsnum,vkbgsnume,vkbgsden,vkbgsdene,nullptr,nullptr,nullptr,nullptr,true,true);

                    TGraphErrors* tgrkbgs = nullptr;
                    MakeTGraphErrors(tgrkbgs,"kbgs","kbgs_Np_"+suf.str(),"N_{part}","#frac{#LTb_{#gamma}#GT#LT#sigma#GT}{#LTb_{#gamma} #sigma#GT}",vkbgs.size(),&vnpart[0],&vkbgs[0],&vnparte[0],&vkbgse[0],can);

                    TGraphErrors* tgrkbgsnum = nullptr;
                    MakeTGraphErrors(tgrkbgsnum,"kbgsnum","kbgsnum_Np_"+suf.str(),"N_{part}","#LTb_{#gamma}#GT#LT#sigma#GT",vkbgsnum.size(),&vnpart[0],&vkbgsnum[0],&vnparte[0],&vkbgsnume[0],can);

                    TGraphErrors* tgrkbgsden = nullptr;
                    MakeTGraphErrors(tgrkbgsden,"kbgsden","kbgsden_Np_"+suf.str(),"N_{part}","#LTb_{#gamma} #sigma#GT",vkbgsden.size(),&vnpart[0],&vkbgsden[0],&vnparte[0],&vkbgsdene[0],can);

                    //-------------------- k_bh --------------------

                    cout<<"Calculating "<<"kbhs_Np_"<<suf.str()<<endl;
                    vector<double> vkbhs, vkbhse, vkbhsnum, vkbhsnume, vkbhsden, vkbhsdene;
                    CalcK(hnpart,tgrbh_Np,tgrsig,tgrbhe,nullptr,vkbhs,vkbhse,vkbhsnum,vkbhsnume,vkbhsden,vkbhsdene,nullptr,nullptr,nullptr,nullptr,true,true);

                    TGraphErrors* tgrkbhs = nullptr;
                    MakeTGraphErrors(tgrkbhs,"kbhs","kbhs_Np_"+suf.str(),"N_{part}","#frac{#LTb_{h}#GT#LT#sigma#GT}{#LTb_{h} #sigma#GT}",vkbhs.size(),&vnpart[0],&vkbhs[0],&vnparte[0],&vkbhse[0],can);

                    TGraphErrors* tgrkbhsnum = nullptr;
                    MakeTGraphErrors(tgrkbhsnum,"kbhsnum","kbhsnum_Np_"+suf.str(),"N_{part}","#LTb_{h}#GT#LT#sigma#GT",vkbhsnum.size(),&vnpart[0],&vkbhsnum[0],&vnparte[0],&vkbhsnume[0],can);

                    TGraphErrors* tgrkbhsden = nullptr;
                    MakeTGraphErrors(tgrkbhsden,"kbhsden","kbhsden_Np_"+suf.str(),"N_{part}","#LTb_{h} #sigma#GT",vkbhsden.size(),&vnpart[0],&vkbhsden[0],&vnparte[0],&vkbhsdene[0],can);

                    //-------------------- k_gh --------------------

                    cout<<"Calculating "<<"kgh_Np_"<<suf.str()<<endl;
                    vector<double> vkgh, vkghe, vkghnum, vkghnume, vkghden, vkghdene;
                    CalcK(hnpart,tgrag_Np,tgrbhs_Np,tgrage,tgrbhse,vkgh,vkghe,vkghnum,vkghnume,vkghden,vkghdene);

                    TGraphErrors* tgrkgh = nullptr;
                    MakeTGraphErrors(tgrkgh,"kgh","kgh_Np_"+suf.str(),"N_{part}","#frac{#LTa_{#gamma} b_{h} #sigma#GT}{#LTa_{#gamma}#GT#LTb_{h}#sigma#GT}",vkgh.size(),&vnpart[0],&vkgh[0],&vnparte[0],&vkghe[0],can);

                    //----- make kgh with monte carlo errors -----

                    TMultiGraph* mgrkgh = new TMultiGraph();
                    TMultiGraph* mgrkgh_ag = new TMultiGraph();
                    TMultiGraph* mgrkgh_bhs = new TMultiGraph();
                    vector<double> vkghenew = CalcKMC(tgrag_Np,tgrbhs_Np,vkgh,vkghe,vnpart,hnpart,mgrkgh,mgrkgh_ag,mgrkgh_bhs);

                    TGraphErrors* tgrkghnew = nullptr;
                    MakeTGraphErrors(tgrkghnew,"kghnew","kghnew_Np_"+suf.str(),"N_{part}","#frac{#LTa_{#gamma} b_{h} #sigma#GT}{#LTa_{#gamma}#GT#LTb_{h}#sigma#GT}",vkgh.size(),&vnpart[0],&vkgh[0],&vnparte[0],&vkghenew[0],can);

                    cmgr->cd();

                    mgrkgh->Add(tgrkghnew,"p");
                    mgrkgh->SetTitle(";N_{part};#kappa_{#gamma h}");
                    mgrkgh->Draw("a");
                    if (itype == GAT::PI0)
                        mgrkgh->GetYaxis()->SetRangeUser(kghlim.at(iptt).first,kghlim.at(iptt).second);
                    mgrkgh->GetYaxis()->SetLabelOffset(.01);
                    mgrkgh->GetXaxis()->SetLabelOffset(.01);
                    mgrkgh->GetYaxis()->SetTitleOffset(1.5);
                    mgrkgh->GetXaxis()->SetTitleOffset(1.1);

                    txtmom->Draw();
                    txtxi->Draw();
                    cmgr->SaveAs(Form("kghnew/multigraph_kgh_%s.pdf",suf.str().c_str()));

                    mgrkgh_ag->SetTitle(";N_{part};a_{#gamma}");
                    mgrkgh_ag->Draw("a");
                    cmgr->SaveAs(Form("kghnew/multigraph_ag_%s.pdf",suf.str().c_str()));

                    mgrkgh_bhs->SetTitle(";N_{part};b_{h} #sigma");
                    mgrkgh_bhs->Draw("a");
                    cmgr->SaveAs(Form("kghnew/multigraph_bhs_%s.pdf",suf.str().c_str()));
                    //---------- end kgh monte carlo ----------

                    TGraphErrors* tgrkghnum = nullptr;
                    MakeTGraphErrors(tgrkghnum,"kghnum","kghnum_Np_"+suf.str(),"N_{part}","#LTa_{#gamma} b_{h} #sigma#GT",vkghnum.size(),&vnpart[0],&vkghnum[0],&vnparte[0],&vkghnume[0],can);

                    TGraphErrors* tgrkghden = nullptr;
                    MakeTGraphErrors(tgrkghden,"kghden","kghden_Np_"+suf.str(),"N_{part}","#LTa_{#gamma}#GT#LTb_{h}#sigma#GT",vkghden.size(),&vnpart[0],&vkghden[0],&vnparte[0],&vkghdene[0],can);

                    //-------------------- k_gh merge --------------------

                    cout<<"Calculating "<<"kghmerge_Np_"<<suf.str()<<endl;
                    vector<double> vkghmerge, vkghmergee, vkghmergenum, vkghmergenume, vkghmergeden, vkghmergedene;
                    CalcK_Merge(hnpartmerge,tgrag_Np,tgrbhs_Np,tgrage,tgrbhse,vkghmerge,vkghmergee,vkghmergenum,vkghmergenume,vkghmergeden,vkghmergedene);

                    TGraphErrors* tgrkghmerge = nullptr;
                    MakeTGraphErrors(tgrkghmerge,"kghmerge","kghmerge_Np_"+suf.str(),"N_{part}","#frac{#LTa_{#gamma} b_{h} #sigma#GT}{#LTa_{#gamma}#GT#LTb_{h}#sigma#GT}",vkghmerge.size(),&vnpartmerge[0],&vkghmerge[0],&vnpartmergee[0],&vkghmergee[0],can);

                    TGraphErrors* tgrkghmergenum = nullptr;
                    MakeTGraphErrors(tgrkghmergenum,"kghmergenum","kghmergenum_Np_"+suf.str(),"N_{part}","#LTa_{#gamma} b_{h} #sigma#GT",vkghmergenum.size(),&vnpartmerge[0],&vkghmergenum[0],&vnpartmergee[0],&vkghmergenume[0],can);

                    TGraphErrors* tgrkghmergeden = nullptr;
                    MakeTGraphErrors(tgrkghmergeden,"kghmergeden","kghmergeden_Np_"+suf.str(),"N_{part}","#LTa_{#gamma}#GT#LTb_{h}#sigma#GT",vkghmergeden.size(),&vnpartmerge[0],&vkghmergeden[0],&vnpartmergee[0],&vkghmergedene[0],can);

                    //-------------------- k_hg --------------------

                    cout<<"Calculating "<<"khg_Np_"<<suf.str()<<endl;
                    vector<double> vkhg, vkhge, vkhgnum, vkhgnume, vkhgden, vkhgdene;
                    CalcK(hnpart,tgrah_Np,tgrbgs_Np,tgrahe,tgrbgse,vkhg,vkhge,vkhgnum,vkhgnume,vkhgden,vkhgdene);

                    TGraphErrors* tgrkhg = nullptr;
                    MakeTGraphErrors(tgrkhg,"khg","khg_Np_"+suf.str(),"N_{part}","#frac{#LTa_{h} b_{#gamma} #sigma#GT}{#LTa_{h}#GT#LTb_{#gamma}#sigma#GT}",vkhg.size(),&vnpart[0],&vkhg[0],&vnparte[0],&vkhge[0],can);

                    //----- make khg with monte carlo errors -----
                    TMultiGraph* mgrkhg = new TMultiGraph();
                    TMultiGraph* mgrkhg_ah = new TMultiGraph();
                    TMultiGraph* mgrkhg_bgs = new TMultiGraph();
                    vector<double> vkhgenew = CalcKMC(tgrah_Np,tgrbgs_Np,vkhg,vkhge,vnpart,hnpart,mgrkhg,mgrkhg_ah,mgrkhg_bgs);

                    TGraphErrors* tgrkhgnew = nullptr;
                    MakeTGraphErrors(tgrkhgnew,"khgnew","khgnew_Np_"+suf.str(),"N_{part}","#frac{#LTa_{h} b_{#gamma} #sigma#GT}{#LTa_{h}#GT#LTb_{#gamma}#sigma#GT}",vkhg.size(),&vnpart[0],&vkhg[0],&vnparte[0],&vkhgenew[0],can);

                    cmgr->cd();

                    mgrkhg->Add(tgrkhgnew,"p");
                    mgrkhg->SetTitle(";N_{part};#kappa_{h #gamma}");
                    mgrkhg->Draw("a");
                    if (itype == GAT::PI0)
                        mgrkhg->GetYaxis()->SetRangeUser(khglim.at(iptt).first,khglim.at(iptt).second);
                    mgrkhg->GetYaxis()->SetLabelOffset(.01);
                    mgrkhg->GetXaxis()->SetLabelOffset(.01);
                    mgrkhg->GetYaxis()->SetTitleOffset(1.5);
                    mgrkhg->GetXaxis()->SetTitleOffset(1.1);

                    txtmom->Draw();
                    txtxi->Draw();
                    cmgr->SaveAs(Form("khgnew/multigraph_khg_%s.pdf",suf.str().c_str()));

                    mgrkhg_ah->SetTitle(";N_{part};a_{h}");
                    mgrkhg_ah->Draw("a");
                    cmgr->SaveAs(Form("khgnew/multigraph_ah_%s.pdf",suf.str().c_str()));

                    mgrkhg_bgs->SetTitle(";N_{part};b_{#gamma} #sigma");
                    mgrkhg_bgs->Draw("a");
                    cmgr->SaveAs(Form("khgnew/multigraph_bgs_%s.pdf",suf.str().c_str()));
                    //---------- end khg monte carlo ----------

                    TGraphErrors* tgrkhgnum = nullptr;
                    MakeTGraphErrors(tgrkhgnum,"khgnum","khgnum_Np_"+suf.str(),"N_{part}","#LTa_{h} b_{#gamma} #sigma#GT",vkhgnum.size(),&vnpart[0],&vkhgnum[0],&vnparte[0],&vkhgnume[0],can);

                    TGraphErrors* tgrkhgden = nullptr;
                    MakeTGraphErrors(tgrkhgden,"khgden","khgden_Np_"+suf.str(),"N_{part}","#LTa_{h}#GT#LTb_{#gamma}#sigma#GT",vkhgden.size(),&vnpart[0],&vkhgden[0],&vnparte[0],&vkhgdene[0],can);

                    //-------------------- k_hg merge --------------------

                    cout<<"Calculating "<<"khgmerge_Np_"<<suf.str()<<endl;
                    vector<double> vkhgmerge, vkhgmergee, vkhgmergenum, vkhgmergenume, vkhgmergeden, vkhgmergedene;
                    CalcK_Merge(hnpartmerge,tgrah_Np,tgrbgs_Np,tgrahe,tgrbgse,vkhgmerge,vkhgmergee,vkhgmergenum,vkhgmergenume,vkhgmergeden,vkhgmergedene);

                    TGraphErrors* tgrkhgmerge = nullptr;
                    MakeTGraphErrors(tgrkhgmerge,"khgmerge","khgmerge_Np_"+suf.str(),"N_{part}","#frac{#LTa_{h} b_{#gamma} #sigma#GT}{#LTa_{h}#GT#LTb_{#gamma}#sigma#GT}",vkhgmerge.size(),&vnpartmerge[0],&vkhgmerge[0],&vnpartmergee[0],&vkhgmergee[0],can);

                    TGraphErrors* tgrkhgmergenum = nullptr;
                    MakeTGraphErrors(tgrkhgmergenum,"khgmergenum","khgmergenum_Np_"+suf.str(),"N_{part}","#LTa_{h} b_{#gamma} #sigma#GT",vkhgmergenum.size(),&vnpartmerge[0],&vkhgmergenum[0],&vnpartmergee[0],&vkhgmergenume[0],can);

                    TGraphErrors* tgrkhgmergeden = nullptr;
                    MakeTGraphErrors(tgrkhgmergeden,"khgmergeden","khgmergeden_Np_"+suf.str(),"N_{part}","#LTa_{h}#GT#LTb_{#gamma}#sigma#GT",vkhgmergeden.size(),&vnpartmerge[0],&vkhgmergeden[0],&vnpartmergee[0],&vkhgmergedene[0],can);

                }
            }
        }
    }
    cout<<"Writing output to "<<output->GetName()<<endl;
    output->Write();
    cout<<"Finished"<<endl;
}

int Getng(TH1F* h1trgcent, vector<double>& vng, vector<double>& vnge, int icent, double nevt)
{
    double ng = 0.;
    double err = 0.;
    for (int jcent = GAT::limcent[icent+1]; jcent > GAT::limcent[icent]; jcent--)
    {
        int bin = h1trgcent->FindBin(jcent);
        ng += h1trgcent->GetBinContent(bin);
        err += pow(h1trgcent->GetBinError(bin),2.);
    }
    ng/=nevt;    
    err = sqrt(err);
    err/=nevt;

    vng.push_back(ng);
    vnge.push_back(err);

    //check
    if (ng != vng[icent])
    {
        cerr<<"ng != vng[icent]"<<endl;
        cerr<<ng<<" != "<<vng[icent]<<endl;
        return -1;
    }

    return 0;
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

int GetNpartDists_Merge(TH1F* hnpart[ncentmerge], vector<double>& vnpart, vector<double>& vnparte)
{
    const char* centdistsname = "CentDists_2bins.root";
    TFile* fcentdists = new TFile(centdistsname,"READ");
    if (!fcentdists || fcentdists->IsZombie())
    {
        cerr<<"Could not open "<<centdistsname<<endl;
        return -1;
    }

    for (int icent = 0; icent < ncentmerge; icent++)
    {
        hnpart[icent] = (TH1F*) gROOT->FindObject(Form("AuAu_Npart_%i_%i",limcentmerge[icent],limcentmerge[icent+1]));
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

int GetRPReso(vector<double>& vsig)
{
    const char* inname = "rpreso.root";
    TFile* input = new TFile(inname,"READ");
    if (!input || input->IsZombie())
    {
        cerr<<"Could not open "<<inname<<endl;
        return -1;
    }

    TGraphErrors* tgr = (TGraphErrors*) gROOT->FindObject("npartreso");
    assert(tgr);
    cout<<"Found reaction plane resolution tgraph "<<tgr->GetName()<<". Printing contents..."<<endl;
    tgr->Print();

    double* yary = tgr->GetY();

    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        vsig.push_back(yary[icent]);
    }

    input->Close();
    delete input;

    return 0;
}

// vector<double> CorrectVector(vector<double>& vog, TGraph* tgr, TH1F* hnpart[GAT::ncent], vector<double>& vnpart, int cnt)

vector<double> CalcKMC(TGraphErrors* tgra, TGraphErrors* tgrb, const vector<double>& vkold, const vector<double>& vkeold, const vector<double>& vnpart, TH1F* hnpart[GAT::ncent], TMultiGraph* mgrk, TMultiGraph* mgra, TMultiGraph* mgrb, TGraphErrors* tgrad = nullptr, TGraphErrors* tgrbd = nullptr, TMultiGraph* mgrad = nullptr, TMultiGraph* mgrbd = nullptr)
{
    //storing new error on k in this vector to be returned
    vector<double> vkerr(GAT::ncent,0.);

    const int NMONTE = 100;
    int ntrial = 0;
    for (int imonte = 0; imonte < NMONTE; imonte++)
    {
        //let's create new functions
        TGraphErrors* tgranew = new TGraphErrors(GAT::ncent+2);
        TGraphErrors* tgrbnew = new TGraphErrors(GAT::ncent+2);

        //new functions for denominator if given
        TGraphErrors* tgradnew = nullptr;
        if (tgrad != nullptr) tgradnew = new TGraphErrors(GAT::ncent+2);
        TGraphErrors* tgrbdnew = nullptr;
        if (tgrbd != nullptr) tgrbdnew = new TGraphErrors(GAT::ncent+2);

        //set first point to evaluated value at npartmin
        {
            double ynewa = SetFirstTGraphPoint(tgra,tgranew);
            if (tgrad != nullptr)
                tgradnew->SetPoint(0,nparthi,ynewa*tgrsig->Eval(nparthi));

            double ynewb = SetFirstTGraphPoint(tgrb,tgrbnew);
            if (tgrbd != nullptr)
                tgrbdnew->SetPoint(0,nparthi,ynewb*tgrsig->Eval(nparthi));
        }
        {
            double ynewa = SetLastTGraphPoint(tgra,tgranew);
            if (tgrad != nullptr)
                tgradnew->SetPoint(GAT::ncent+1,npartlo,ynewa*tgrsig->Eval(npartlo));

            double ynewb = SetLastTGraphPoint(tgrb,tgrbnew);
            if (tgrbd != nullptr)
                tgrbdnew->SetPoint(GAT::ncent+1,npartlo,ynewb*tgrsig->Eval(npartlo));

        }

        for (int icent = 1; icent <= GAT::ncent; icent++)
        {
            double yrata = SetNewTGraphPoint(tgra,tgranew,icent);
            double yratb = SetNewTGraphPoint(tgrb,tgrbnew,icent);
            //not using these ratios now but will be useful when
            //denominator differs from numerator

            if (tgradnew != nullptr)
            {
                SetPointWithRatio(tgrad, tgradnew, icent, yrata);

            }

            if (tgrbdnew != nullptr)
            {
                SetPointWithRatio(tgrbd, tgrbdnew, icent, yratb);

            }
        }

        //new tgraphs have been made. find new k

        //function requires these vars but most wont be used
        vector<double> vk, vke, vknum, vknume, vkden, vkdene;

        //passing nullptr for tgrae and tgrbe since
        //we care not about the CalcK method of estimated the error on k
        CalcK(hnpart,tgranew,tgrbnew,nullptr,nullptr,vk,vke,vknum,vknume,vkden,vkdene,tgradnew,tgrbdnew);

        bool skip = false;
        for (int icent = 0; icent < GAT::ncent; icent++)
        {
            if (fabs(vk[icent]-vkold[icent]) > 3.0*2.0*vkeold[icent])
                skip = true;
        }
        if (skip) continue;

        bool shalldraw = true;
        for (int icent = 0; icent < GAT::ncent; icent++)
        {                   
            //use initial calc of k as mean in variance calc
            vkerr[icent] += pow(vk[icent]-vkold[icent],2.);

        }

        //attemtping to isolate extreme cases
        if (shalldraw)
        {
            TGraph* grk = new TGraph((int)vk.size(),vnpart.data(),vk.data());
            grk->SetLineColor(TColor::GetColorPalette(254*imonte/(NMONTE-1)));
            mgrk->Add(grk,"lx");

            tgranew->SetLineColor(TColor::GetColorPalette(254*imonte/(NMONTE-1)));
            mgra->Add(tgranew,"lx");

            tgrbnew->SetLineColor(TColor::GetColorPalette(254*imonte/(NMONTE-1)));
            mgrb->Add(tgrbnew,"lx");

            if (mgrad != nullptr)
            {
                tgradnew->SetLineColor(TColor::GetColorPalette(254*imonte/(NMONTE-1)));
                mgrad->Add(tgradnew,"lx");
            }

            if (mgrbd != nullptr)
            {
                tgrbdnew->SetLineColor(TColor::GetColorPalette(254*imonte/(NMONTE-1)));
                mgrbd->Add(tgrbdnew,"lx");
            }
        }
    }
    //monte carlo concluded
    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        //average to get variance. sqrt to get sigma
        vkerr[icent] = sqrt(vkerr[icent]/NMONTE);
    }

    //add original graphs last so theyre plotted on top
    mgra->Add(tgra,"p");
    mgrb->Add(tgrb,"p");
    if (mgrad != nullptr) mgrad->Add(tgrad,"p");
    if (mgrbd != nullptr) mgrbd->Add(tgrbd,"p");

    return vkerr;
}

//this method was manually copied to acceptCorrect.C

void CalcK(TH1F* hnpart[GAT::ncent], TGraphErrors* tgra, TGraphErrors* tgrb, TGraphErrors* tgrae, TGraphErrors* tgrbe, vector<double>& vk, vector<double>& vke, vector<double>& vnum, vector<double>& vnume, vector<double>& vden, vector<double>& vdene, TGraphErrors* tgrad = nullptr, TGraphErrors* tgrbd = nullptr, TGraphErrors* tgrade = nullptr, TGraphErrors* tgrbde = nullptr, bool usespline = true, bool flip = false)
{
    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        //have to comment out. too much chatter
        // cout<<"---------- icent = "<<icent<<" ----------"<<endl;

        int nbinx = hnpart[icent]->GetNbinsX();

        double num = 0.;
        double nume = 0.;
        double dene = 0.;

        //must recalculate exp values
        //differ to measured values
        //causing less than 1 k values
        double expa = 0.;
        double expb = 0.;

        double expae = 0.;
        double expbe = 0.;

        for (int ibin = nbinx; ibin >= 1; ibin--)
        {
            double w = hnpart[icent]->GetBinContent(ibin);
            if (w <= 0.0) continue;

            double we = hnpart[icent]->GetBinError(ibin);

            double npart = hnpart[icent]->GetBinCenter(ibin);

            double a{0.}, ae{0.}, b{0.}, be{0.};

            if (usespline)//first order spline
            {
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

                if (tgrae != nullptr)
                {
                    ae = tgrae->Eval(npart);
                    if (ae < 0.) ae = 0.;
                }

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

                if (tgrbe != nullptr)
                {
                    be = tgrbe->Eval(npart);
                    if (be < 0.) be = 0.;
                }
            }
            else
            {
                cerr<<"THIS IS NO ELSE!"<<endl;
                return;
            }

            num += w*a*b;

            if (tgrad && tgrbd)
            {
                a = tgrad->Eval(npart);
                if (a < 0.)
                {
                    ostringstream oss;
                    oss<<tgrad->GetName();
                    string str = oss.str();
                    if (str[0] == 'a')
                    {
                        cout<<tgrad->GetName()<<" evaluated to "<<a<<". Setting to 0."<<endl;
                        a = 0;
                    }
                }
                if (tgrade != nullptr)
                    ae = tgrade->Eval(npart);
                else
                    ae = -9999.;//if tgrade not specified error is not calculated. set to absurd val

                b = tgrbd->Eval(npart);
                if (b < 0.)
                {
                    ostringstream oss;
                    oss<<tgrbd->GetName();
                    string str = oss.str();
                    if (str[0] == 'a')
                    {
                        cout<<tgrbd->GetName()<<" evaluated to "<<b<<". Setting to 0."<<endl;
                        b = 0;
                    }
                }
                if (tgrbde != nullptr)
                    be = tgrbde->Eval(npart);
                else
                    be = -9999.;//if tgrbde not specified error is not calculated. set to absurd val
            }

            expa += w*a;

            expb += w*b;

        }

        nume = pow(tgra->GetEY()[icent]*tgrb->GetY()[icent],2.);
        nume += pow(tgra->GetY()[icent]*tgrb->GetEY()[icent],2.);

        nume = sqrt(nume);

        double den = expa*expb;
        if (tgrad && tgrbd)
        {
            dene = pow(tgrad->GetEY()[icent]*tgrbd->GetY()[icent],2.);
            dene += pow(tgrad->GetY()[icent]*tgrbd->GetEY()[icent],2.);
            dene = sqrt(dene);
        }
        else dene = nume;

        double k;
        if (flip)
        {
            vnum.push_back(den);
            vnume.push_back(dene);
            vden.push_back(num);
            vdene.push_back(nume);

            k = den/num;
            cout<<"denom and num have been swapped. following output does not apply"<<endl;
        }
        else
        {
            vnum.push_back(num);
            vnume.push_back(nume);
            vden.push_back(den);
            vdene.push_back(dene);

            k = num/den;
        }

        vk.push_back(k);        

        double ke = k*sqrt(pow(nume/num,2)+pow(dene/den,2))/2.;//average between perfectly correlated and uncorrelated
        vke.push_back(ke);

        //what about the error in the denominator
        //but what about covariance? same term (tho weighted) is used in above loop
        //all weighted errors are actually based off this error and neighboring ones
        //forget it. it is just an approximation

        //have to comment out. too much chatter
    }

}

//only changes are ncent and hnpart
void CalcK_Merge(TH1F* hnpart[ncentmerge], TGraphErrors* tgra, TGraphErrors* tgrb, TGraphErrors* tgrae, TGraphErrors* tgrbe, vector<double>& vk, vector<double>& vke, vector<double>& vnum, vector<double>& vnume, vector<double>& vden, vector<double>& vdene, TGraphErrors* tgrad = nullptr, TGraphErrors* tgrbd = nullptr, TGraphErrors* tgrade = nullptr, TGraphErrors* tgrbde = nullptr, bool usespline = true, bool flip = false)
{
    cout<<"CalcK_Merge():"<<endl;

    for (int icent = 0; icent < ncentmerge; icent++)
    {
        //have to comment out. too much chatter
        // cout<<"---------- icent = "<<icent<<" ----------"<<endl;

        int nbinx = hnpart[icent]->GetNbinsX();

        double num = 0.;
        double nume = 0.;
        double dene = 0.;

        //must recalculate exp values
        //differ to measured values
        //causing less than 1 k values
        double expa = 0.;
        double expb = 0.;

        double expae = 0.;
        double expbe = 0.;

        for (int ibin = 1; ibin <= nbinx; ibin++)
        {
            double w = hnpart[icent]->GetBinContent(ibin);
            if (w <= 0.0) continue;

            double we = hnpart[icent]->GetBinError(ibin);

            double npart = hnpart[icent]->GetBinCenter(ibin);

            double a, ae, b, be;

            if (usespline)//first order spline
            {
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

                ae = tgrae->Eval(npart);
                if (ae < 0.) ae = 0.;

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

                if (tgrbe != nullptr)
                {
                    be = tgrbe->Eval(npart);
                    if (be < 0.) be = 0.;
                }
                else be = 0.;
            }
            else
            {
                cerr<<"THIS IS NO ELSE!"<<endl;
                return;
            }

            num += w*a*b;

            if (tgrad && tgrbd)
            {
                a = tgrad->Eval(npart);
                if (a < 0.)
                {
                    ostringstream oss;
                    oss<<tgrad->GetName();
                    string str = oss.str();
                    if (str[0] == 'a')
                    {
                        cout<<tgrad->GetName()<<" evaluated to "<<a<<". Setting to 0."<<endl;
                        a = 0;
                    }
                }
                ae = tgrade->Eval(npart);

                b = tgrbd->Eval(npart);
                if (b < 0.)
                {
                    ostringstream oss;
                    oss<<tgrbd->GetName();
                    string str = oss.str();
                    if (str[0] == 'a')
                    {
                        cout<<tgrbd->GetName()<<" evaluated to "<<b<<". Setting to 0."<<endl;
                        b = 0;
                    }
                }
                be = tgrbde->Eval(npart);
            }

            expa += w*a;

            expb += w*b;

        }

        nume = pow(tgra->GetEY()[icent]*tgrb->GetY()[icent],2.);
        nume += pow(tgra->GetY()[icent]*tgrb->GetEY()[icent],2.);

        nume = sqrt(nume);

        double den = expa*expb;
        if (tgrad && tgrbd)
        {
            dene = pow(tgrad->GetEY()[icent]*tgrbd->GetY()[icent],2.);
            dene += pow(tgrad->GetY()[icent]*tgrbd->GetEY()[icent],2.);
            dene = sqrt(dene);
        }
        else dene = nume;

        double k;
        if (flip)
        {
            vnum.push_back(den);
            vnume.push_back(dene);
            vden.push_back(num);
            vdene.push_back(nume);

            k = den/num;
            cout<<"denom and num have been swapped. following output does not apply"<<endl;
        }
        else
        {
            vnum.push_back(num);
            vnume.push_back(nume);
            vden.push_back(den);
            vdene.push_back(dene);

            k = num/den;
        }

        vk.push_back(k);        

        double ke = k*sqrt(pow(nume/num,2)+pow(dene/den,2))/2.;//average between perfectly correlated and uncorrelated
        vke.push_back(ke);

        //what about the error in the denominator
        //but what about covariance? same term (tho weighted) is used in above loop
        //all weighted errors are actually based off this error and neighboring ones
        //forget it. it is just an approximation

        //have to comment out. too much chatter
    }

}

//      else

//      cout<<"---------- icent = "<<icent<<" ----------"<<endl;

//          else

// vector<double> CalcK(TH1F* hnpart[GAT::ncent], TGraphErrors* tgra, TGraphErrors* tgrb, const vector<double>& va, const vector<double>& vb, bool usespline = true)

int MakeTGraphErrors(TGraphErrors*& tgr, string dir, string name, string xname, string yname, const int size, const double* x, const double* y, const double* ex = 0, const double* ey = 0, TCanvas* can = NULL)
{
    tgr = new TGraphErrors(size,x,y,ex,ey);

    ostringstream oss;
    oss<<name<<";"<<xname<<";"<<yname;
    tgr->SetTitle(oss.str().c_str());

    tgr->SetName(name.c_str());
    cout<<"Made TGraphErrors "<<tgr->GetName()<<endl;
    tgr->Write();

    if (can)
    {
        can->cd();

        oss.clear();
        oss.str("");
        oss<<dir<<"/"<<name<<".pdf";

        tgr->Draw("AL");

        // I meant to set limits for the multigraphs... not these...
        // knew limits only set for k*new for pi0s

        can->SaveAs(oss.str().c_str());
    }
    return 0;
}

double GetError(TH1F* h1, TF1* f1)
{
    double num = 0., den = 0.;

    int nx = h1->GetNbinsX();

    for (int ix = 1; ix <= nx; ix++)
    {
        double hval = h1->GetBinContent(ix);
        double x = h1->GetBinCenter(ix);
        double fval = f1->Eval(x);
        double err = h1->GetBinError(ix);

        num += pow((hval-fval)/err,2.);
        den += pow(1./err,2.);
    }

    return sqrt(num/den);
}

void DrawError(TF1* f1, double err, int ipar)
{
    double par = f1->GetParameter(ipar);

    TF1* f1min = (TF1*) f1->Clone("f1min");
    f1min->SetParameter(ipar,par-err);

    TF1* f1max = (TF1*) f1->Clone("f1max");
    f1max->SetParameter(ipar,par+err);

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
    grshade->SetFillColor(16);

    grshade->Draw("f");

    grmin->SetLineColor(16);
    grmin->SetLineWidth(1);
    grmin->Draw("l");

    grmax->SetLineColor(16);
    grmax->SetLineWidth(1);
    grmax->Draw("l");

    f1->Draw("same");
}

void DrawErrorExtra(TF1* f1, double err, int ipar)
{
    double par = f1->GetParameter(ipar);

    TF1* f1min = (TF1*) f1->Clone("f1min");
    f1min->SetParameter(ipar,par-err);

    TF1* f1max = (TF1*) f1->Clone("f1max");
    f1max->SetParameter(ipar,par+err);

    const Int_t n = 24;
    Double_t x[n], ymin[n], ymax[n];
    for (int i=0;i<n;i++) {
        x[i] = i/(n-1.);

        ymax[i] = f1max->Eval(x[i]);
        ymin[i] = f1min->Eval(x[i]);

    }
    TGraph *grmin = new TGraph(n,x,ymin);
    TGraph *grmax = new TGraph(n,x,ymax);

    grmin->SetLineColor(14);
    grmin->SetLineWidth(1);
    grmin->Draw("l");

    grmax->SetLineColor(14);
    grmax->SetLineWidth(1);
    grmax->Draw("l");

    f1->Draw("same");
}

double SetNewTGraphPoint(TGraphErrors* tgr, TGraphErrors* tgrnew, int ipt)
{
    double x, y;
    tgr->GetPoint(ipt-1,x,y);
    double ex = tgr->GetErrorX(ipt-1);
    double ey = tgr->GetErrorY(ipt-1);

    //get new y value using error as sigma
    double ynew = rndm->Gaus(y,ey);
    tgrnew->SetPoint(ipt,x,ynew);
    tgrnew->SetPointError(ipt,ex,ey);

    return ynew/y;

}

double SetLastTGraphPoint(TGraphErrors* tgr, TGraphErrors* tgrnew)
{
    double x = npartlo;
    double y = tgr->Eval(x);
    double ex = tgr->GetErrorX(GAT::ncent-1);
    double ey = tgr->GetErrorY(GAT::ncent-1);

    double ynew = rndm->Gaus(y,ey);
    tgrnew->SetPoint(GAT::ncent+1,x,ynew);
    tgrnew->SetPointError(GAT::ncent+1,ex,ey);

    return ynew;
}

double SetFirstTGraphPoint(TGraphErrors* tgr, TGraphErrors* tgrnew)
{
    double x = nparthi;
    double y = tgr->Eval(x);
    double ex = tgr->GetErrorX(0);
    double ey = tgr->GetErrorY(0);

    double ynew = rndm->Gaus(y,ey);
    tgrnew->SetPoint(0,x,ynew);
    tgrnew->SetPointError(0,ex,ey);

    return ynew;
}

void SetPointWithRatio(TGraphErrors* tgr, TGraphErrors* tgrnew, int ipt, double rat)
{
    double x, y;
    tgr->GetPoint(ipt-1,x,y);
    double ex = tgr->GetErrorX(ipt-1);
    double ey = tgr->GetErrorY(ipt-1);

    tgrnew->SetPoint(ipt,x,y*rat);
    tgrnew->SetPointError(ipt,ex,ey);
}
