//altered from xi to pta

//reads acceptance corrected gamtree results
//for each trigger type, ptt, cent, xi
//     - calculate (or rather isolate) each k component
//     - correct each k component
//     - write resulting TGraphs to file

//use limit for convergence instead of limit of 10. (USE RELATIVE LIMIT)

const int nadj = 2;

const int npartmin = 2;
const int npartmax = 390;

const float PI = (float) TMath::Pi();

namespace GAT {
    #include "gamAnaTree.h"
}

int GetNpartDists(TH1F* hnpart[GAT::ncent], vector<double>& vnpart, vector<double>& vnparte);
int GetRPReso(vector<double>& vsig);
vector<double> CorrectVector(vector<double>& vog, TGraph* tgr, TH1F* hnpart[GAT::ncent], vector<double>& vnpart, int cnt=0);
int Getng(TH1F* h1trgcent, vector<double>& vng, vector<double>& vnge, int icent, double nevt);

int MakeTGraphErrors(TGraphErrors*& tgr, string dir, string name, string xname, string yname, const int size, const double* x, const double* y, const double* ex = 0, const double* ey = 0, TCanvas* can = NULL);

void CalcK(TH1F* hnpart[GAT::ncent], TGraphErrors* tgra, TGraphErrors* tgrb, TGraphErrors* tgrae, TGraphErrors* tgrbe, vector<double>& vk, vector<double>& vke, vector<double>& num, vector<double>& nume, vector<double>& den, vector<double>& dene, TGraphErrors* tgrad = nullptr, TGraphErrors* tgrbd = nullptr, TGraphErrors* tgrade = nullptr, TGraphErrors* tgrbde = nullptr, bool usespline = true);

void calcKComponents_pta(
    const char* inname = "gamtree_acceptcorrect.root",
    const char* inname_gamhist = "gamAna_histos.root",
    const char* outname = "kComponents_pta.root"
    )
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(12);
    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(15);

    //Get Npart distributions
    TH1F* hnpart[GAT::ncent];
    vector<double> vnpart;
    vector<double> vnparte;
    if (GetNpartDists(hnpart,vnpart,vnparte)<0) return;

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

    TCanvas* can = new TCanvas("c1","",2000,2000);

    //Get RP Reso
    vector<double> vsig;
    if (GetRPReso(vsig)<0) return;
    output->cd();//closing rpreso input does not reopen output? must cd..
    TGraphErrors* tgrsig = new TGraphErrors(vsig.size(),&vnpart[0],&vsig[0],&vnparte[0],0);
    tgrsig->SetTitle("npartsigma;N_{part};#sigma");
    tgrsig->SetName("npartsigma");
    tgrsig->Write();
    tgrsig->Draw("AC");
    can->SaveAs("npartsigma.png");

    //Get centrality distribution from gamAna histo file
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

    TFile* input = new TFile(inname,"READ");
    if (!input || input->IsZombie())
    {
        cerr<<"Could not open "<<inname<<endl;
        return;
    }    

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            cout<<"WARNING! Only looking at Pi0!"<<endl;
            if (itype != GAT::PI0) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                ostringstream sufg;
                sufg<<"adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt;

                printf("\nTrigger adj%i %s ptt%i dpsi%i\n",iadj,GAT::spec[itype],iptt,GAT::ndpsi);

                input->cd();

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
                tgrng->SetTitle(Form("ng_adj%i_%s_ptt%i;centrality;<n_{#gamma}>",iadj,GAT::spec[itype],iptt));
                tgrng->Write(Form("ng_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt));

                //plot as function of Npart
                TGraphErrors* tgrng_Np = new TGraphErrors(vng.size(),&vnpart[0],&vng[0],&vnparte[0],&vnge[0]);
                tgrng_Np->SetTitle(Form("ng_Np_adj%i_%s_ptt%i;N_{part};<n_{#gamma}>",iadj,GAT::spec[itype],iptt));
                tgrng_Np->SetName(Form("ng_Np_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt));             
                tgrng_Np->Write();
                tgrng_Np->Draw("AC");
                can->SaveAs(Form("ng/ng_Np_adj%i_%s_ptt%i.png",iadj,GAT::spec[itype],iptt));

                //plot a_g too. why not?
                vector<double> vag (vng);
                transform(vag.begin(),vag.end(),vag.begin(), [](double x) { return x / (2.*PI); } );
                vector<double> vage (vnge);
                transform(vage.begin(),vage.end(),vage.begin(), [](double x) { return x / (2.*PI); } );
                // tgrag_Np->SetTitle(Form("ag_Np_adj%i_%s_ptt%i;N_{part};<a_{#gamma}>",iadj,GAT::spec[itype],iptt));
                TGraphErrors* tgrag_Np = nullptr;
                MakeTGraphErrors(tgrag_Np,"ag","ag_Np_"+sufg.str(),"N_{part}","<a_{#gamma}>",vag.size(),&vnpart[0],&vag[0],&vnparte[0],&vage[0],can);

                TGraphErrors* tgrage = nullptr;
                MakeTGraphErrors(tgrage,"age","age_Np_"+sufg.str(),"N_{part}","<a_{#gamma}> error",vag.size(),&vnpart[0],&vage[0],&vnparte[0],0,can);

                // tgrng_Np_C->SetTitle(Form("ng_Np_C_adj%i_%s_ptt%i;N_{part};<n_{#gamma}>",iadj,GAT::spec[itype],iptt));

                for (int ipta = 0; ipta < GAT::npta; ipta++)
                {

                    ostringstream suf;
                    suf<<"adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_pta"<<ipta;

                    vector<double> vah;//<a_h>
                    vector<double> vahe;//error on <a_h>

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
                        oss_Ein<<"_pta"<<ipta;

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
                        oss_Eout<<"_pta"<<ipta;

                        TH1F* h1Eout = (TH1F*) gROOT->FindObject(oss_Eout.str().c_str());
                        assert(h1Eout);

                        //scale by nevt in cent range!
                        h1Eout->Scale(1./vnevt[icent]);

                        //scale by binwidth (times pi)
                        h1Eout->Scale(1./(h1Eout->GetBinWidth(1)*PI));

                        // ========== E_all ==========

                        TH1F* h1Eall = (TH1F*) h1Ein->Clone(Form("Eall_adj%i_%s_ptt%i_cent%i_pta%i",iadj,GAT::spec[itype],iptt,icent,ipta));
                        h1Eall->Add(h1Eout,1.0);
                        h1Eall->SetDirectory(output);

                        //==================== <a_h> ====================

                        double error_Eall = -9999.;
                        double integral_Eall = h1Eall->IntegralAndError(1,h1Eall->GetNbinsX(),error_Eall,"width")*PI;//width option multiplies by bin width.  the pi has been dropped and must be included manually
                        error_Eall*=PI;

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

                        //==================== <b_g s> <b_h s> ====================

                        cout<<"Fitting "<<h1Eall->GetName()<<endl;

                        h1Eall->SetTitle(";#Delta#phi/#pi;n_{#gamma}Y^{#psi}_{mix}");

                        h1Eall->SetMarkerStyle(kFullCircle);
                        h1Eall->SetMarkerColor(kBlue);
                        h1Eall->SetMarkerSize(6);

                        TF1* f1all = new TF1(Form("fit_%s",h1Eall->GetName()),"[0]+[1]*cos(x*2.*TMath::Pi())",0.,1.);
                        f1all->FixParameter(0,integral_Eall/PI);
                        f1all->SetParName(0,"4#pi<a_{#gamma}><a_{h}> ");
                        f1all->SetParName(1,"2#pi<b_{#gamma}#sigma><b_{h}#sigma> ");
                        f1all->SetLineWidth(4);

                        h1Eall->Draw("E1");
                        can->SetLeftMargin(.15);
                        h1Eall->GetYaxis()->SetTitleOffset(2.0);
                        h1Eall->GetYaxis()->SetMaxDigits(2);

                        h1Eall->Fit(f1all,"");

                        TLine* tl = new TLine(0.,f1all->GetParameter(0),1.,f1all->GetParameter(0));
                        tl->SetLineStyle(9);
                        tl->SetLineWidth(3);
                        tl->Draw("same");
                        can->Modified();
                        can->Update();

                        TPaveStats *ps = (TPaveStats*)can->GetPrimitive("stats");
                        ps->Draw();
                        ps->SetX1NDC(.25);
                        ps->SetX2NDC(.80);
                        ps->SetY1NDC(.80);
                        ps->SetY2NDC(.925);

                        can->SaveAs(Form("bgsbhs/bgsbhs_adj%i_%s_ptt%i_cent%i_pta%i.png",iadj,GAT::spec[itype],iptt,icent,ipta));

                        double bgsbhs = -9999., bgsbhse = -9999.;
                        bgsbhs = f1all->GetParameter(1)/(2.*PI);
                        bgsbhse = f1all->GetParError(1)/(2.*PI);

                        vbgsbhs.push_back(bgsbhs);
                        vbgsbhse.push_back(bgsbhse);

                        //==================== <b_g s> ====================

                        TH1F* h1Ediff = (TH1F*) h1Ein->Clone(Form("Ediff_adj%i_%s_ptt%i_cent%i_pta%i",iadj,GAT::spec[itype],iptt,icent,ipta));
                        h1Ediff->Add(h1Eout,-1.0);
                        h1Ediff->SetDirectory(output);

                        double error_Ediff {-9999.};
                        double integral_Ediff = h1Ediff->IntegralAndError(1,h1Ediff->GetNbinsX(),error_Ediff,"width")*PI;
                        error_Ediff*=PI;

                        //WARING: this pi should not be here no? -> it's from integral
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

                        h1Ediff->SetTitle(";#Delta#phi/#pi;n^{+}_{#gamma}Y^{#psi+}_{mix} - n^{-}_{#gamma}Y^{#psi-}_{mix}");

                        h1Ediff->SetMarkerStyle(kFullCircle);
                        h1Ediff->SetMarkerColor(kBlue);
                        h1Ediff->SetMarkerSize(6);

                        TF1* f1 = new TF1(Form("bhscos_adj%i_%s_ptt%i_cent%i_pta%i",iadj,GAT::spec[itype],iptt,icent,ipta),"[0]+[1]*cos(x*2.*TMath::Pi())",0.,1.);
                        f1->FixParameter(0,8*ah*bgs);
                        f1->SetParName(0,"8<a_{h}><b_{#gamma}#sigma> ");
                        f1->SetParName(1,"8<a_{#gamma}><b_{h}#sigma> ");
                        f1->SetLineWidth(4);

                        h1Ediff->Draw("E1");
                        can->SetLeftMargin(.15);
                        h1Ediff->GetYaxis()->SetTitleOffset(2.0);
                        h1Ediff->GetYaxis()->SetMaxDigits(2);

                        h1Ediff->Fit(f1,"E");

                        tl = new TLine(0.,8*ah*bgs,1.,8*ah*bgs);
                        tl->SetLineStyle(9);
                        tl->SetLineWidth(3);
                        tl->Draw("same");
                        can->Modified();
                        can->Update();

                        ps = (TPaveStats*)can->GetPrimitive("stats");
                        ps->Draw();
                        ps->SetX1NDC(.25);
                        ps->SetX2NDC(.80);
                        ps->SetY1NDC(.80);
                        ps->SetY2NDC(.925);

                        can->SaveAs(Form("bhscos/bhscos_adj%i_%s_ptt%i_cent%i_pta%i.png",iadj,GAT::spec[itype],iptt,icent,ipta));

                        double bhs = -9999., bhse = -9999.;
                        {
                            double denom = 8.*ag;

                            double par = f1->GetParameter(1);
                            double parerr = f1->GetParError(1);

                            bhs = par/denom;

                            double e0 = parerr/denom;
                            double e1 = -(par*8.*age)/(denom*denom);

                            bhse = sqrt(e0*e0 + e1*e1);
                        }

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
                    tgrah->SetTitle(Form("ah_adj%i_%s_ptt%i_pta%i;centrality;<a_h}>",iadj,GAT::spec[itype],iptt,ipta));
                    tgrah->Write(Form("ah_adj%i_%s_ptt%i_pta%i",iadj,GAT::spec[itype],iptt,ipta));

                    TGraphErrors* tgrah_Np = nullptr;
                    MakeTGraphErrors(tgrah_Np,"ah","ah_Np_"+suf.str(),"N_{part}","<a_{h}>",vah.size(),&vnpart[0],&vah[0],&vnparte[0],&vahe[0],can);

                    TGraphErrors* tgrahe = nullptr;
                    MakeTGraphErrors(tgrahe,"ahe","ahe_Np_"+suf.str(),"N_{part}","<a_{h}> error",vah.size(),&vnpart[0],&vahe[0],&vnparte[0],0,can);

                    //==================================================

                    //==================== <b_g sigma> ====================

                    TGraph* tgrbgs = new TGraph(vbgs.size(),&vcent[0],&vbgs[0]);
                    tgrbgs->SetTitle(Form("bgs_adj%i_%s_ptt%i_pta%i;centrality;<b_{#gamma} #sigma>",iadj,GAT::spec[itype],iptt,ipta));
                    tgrbgs->SetName(Form("bgs_adj%i_%s_ptt%i_pta%i",iadj,GAT::spec[itype],iptt,ipta));
                    tgrbgs->Write();

                    TGraphErrors* tgrbgs_Np = nullptr;
                    MakeTGraphErrors(tgrbgs_Np,"bgs","bgs_Np_"+suf.str(),"N_{part}","<b_{#gamma} #sigma>",vbgs.size(),&vnpart[0],&vbgs[0],&vnparte[0],&vbgse[0],can);

                    TGraphErrors* tgrbgse = nullptr;
                    MakeTGraphErrors(tgrbgse,"bgse","bgse_Np_"+suf.str(),"N_{part}","<b_{#gamma} #sigma> error",vbgse.size(),&vnpart[0],&vbgse[0],&vnparte[0],0,can);

                    //==================================================

                    //==================== <b_g> ====================

                    TGraph* tgrbg = new TGraph(vbg.size(),&vcent[0],&vbg[0]);
                    tgrbg->SetTitle(Form("bg_adj%i_%s_ptt%i_pta%i;centrality;<b_{#gamma} #sigma>/<#sigma>",iadj,GAT::spec[itype],iptt,ipta));
                    tgrbg->SetName(Form("bg_adj%i_%s_ptt%i_pta%i",iadj,GAT::spec[itype],iptt,ipta));
                    tgrbg->Write();

                    TGraphErrors* tgrbg_Np = nullptr;
                    MakeTGraphErrors(tgrbg_Np,"bg","bg_Np_"+suf.str(),"N_{part}","<b_{#gamma}>",vbg.size(),&vnpart[0],&vbg[0],&vnparte[0],&vbge[0],can);

                    TGraphErrors* tgrbge = nullptr;
                    MakeTGraphErrors(tgrbge,"bge","bge_Np_"+suf.str(),"N_{part}","<b_{#gamma}> error",vbge.size(),&vnpart[0],&vbge[0],&vnparte[0],0,can);

                    //==================================================

                    //==================== <v_g> ====================

                    TGraph* tgrvg = new TGraph(vvg.size(),&vcent[0],&vvg[0]);
                    tgrvg->SetTitle(Form("vg_adj%i_%s_ptt%i_pta%i;centrality;#frac{<b_{#gamma}#sigma>}{<#sigma><a_{#gamma}>} #approx v_{#gamma}",iadj,GAT::spec[itype],iptt,ipta));
                    tgrvg->SetName(Form("vg_adj%i_%s_ptt%i_pta%i",iadj,GAT::spec[itype],iptt,ipta));
                    tgrvg->Write();

                    TGraphErrors* tgrvg_Np = nullptr;
                    MakeTGraphErrors(tgrvg_Np,"vg","vg_Np_"+suf.str(),"N_{part}","<v_{#gamma}>",vvg.size(),&vnpart[0],&vvg[0],&vnparte[0],&vvge[0],can);

                    //==================================================

                    //==================== <b_h sigma> ====================

                    TGraph* tgrbhs = new TGraph(vbhs.size(),&vcent[0],&vbhs[0]);
                    tgrbhs->SetTitle(Form("bhs_adj%i_%s_ptt%i_pta%i;centrality;<b_{h} #sigma>",iadj,GAT::spec[itype],iptt,ipta));
                    tgrbhs->Write(Form("bhs_adj%i_%s_ptt%i_pta%i",iadj,GAT::spec[itype],iptt,ipta));

                    TGraphErrors* tgrbhs_Np = nullptr;
                    MakeTGraphErrors(tgrbhs_Np,"bhs","bhs_Np_"+suf.str(),"N_{part}","<b_{h} #sigma>",vbhs.size(),&vnpart[0],&vbhs[0],&vnparte[0],&vbhse[0],can);

                    TGraphErrors* tgrbhse = nullptr;
                    MakeTGraphErrors(tgrbhse,"bhse","bhse_Np_"+suf.str(),"N_{part}","<b_{h} #sigma> error",vbhse.size(),&vnpart[0],&vbhse[0],&vnparte[0],0,can);

                    //==================================================

                    //==================== <b_h> ====================

                    TGraph* tgrbh = new TGraph(vbh.size(),&vcent[0],&vbh[0]);
                    tgrbh->SetTitle(Form("bh_adj%i_%s_ptt%i_pta%i;centrality;<b_{h} #sigma>/<#sigma>",iadj,GAT::spec[itype],iptt,ipta));
                    tgrbh->Write(Form("bh_adj%i_%s_ptt%i_pta%i",iadj,GAT::spec[itype],iptt,ipta));

                    TGraphErrors* tgrbh_Np = nullptr;
                    MakeTGraphErrors(tgrbh_Np,"bh","bh_Np_"+suf.str(),"N_{part}","<b_{h}>",vbh.size(),&vnpart[0],&vbh[0],&vnparte[0],&vbhe[0],can);

                    TGraphErrors* tgrbhe = nullptr;
                    MakeTGraphErrors(tgrbhe,"bhe","bhe_Np_"+suf.str(),"N_{part}","<b_{h}> error",vbhe.size(),&vnpart[0],&vbhe[0],&vnparte[0],0,can);

                    //==================================================

                    //==================== <v_h> ====================

                    TGraph* tgrvh = new TGraph(vvh.size(),&vcent[0],&vvh[0]);
                    tgrvh->SetTitle(Form("vh_adj%i_%s_ptt%i_pta%i;centrality;#frac{<b_{h}#sigma>}{<#sigma><a_{h}>} #approx v_{h}",iadj,GAT::spec[itype],iptt,ipta));
                    tgrvh->Write(Form("vh_adj%i_%s_ptt%i_pta%i",iadj,GAT::spec[itype],iptt,ipta));

                    TGraphErrors* tgrvh_Np = nullptr;
                    MakeTGraphErrors(tgrvh_Np,"vh","vh_Np_"+suf.str(),"N_{part}","<v_{h}>",vvh.size(),&vnpart[0],&vvh[0],&vnparte[0],&vvhe[0],can);

                    //==================================================

                    //==================== <b_g sigma> <b_h sigma> ====================

                    TGraph* tgrbgsbhs = new TGraph(vbgsbhs.size(),&vcent[0],&vbgsbhs[0]);
                    tgrbgsbhs->SetTitle(Form("bgsbhs_adj%i_%s_ptt%i_pta%i;centrality;<b_{#gamma}#sigma><b_{h}#sigma>",iadj,GAT::spec[itype],iptt,ipta));
                    tgrbgsbhs->Write(Form("bgsbhs_adj%i_%s_ptt%i_pta%i",iadj,GAT::spec[itype],iptt,ipta));

                    TGraphErrors* tgrbgsbhs_Np = nullptr;
                    MakeTGraphErrors(tgrbgsbhs_Np,"bgsbhs","bgsbhs_Np_"+suf.str(),"N_{part}","<b_{#gamma}#sigma><b_{h}#sigma>",vbgsbhs.size(),&vnpart[0],&vbgsbhs[0],&vnparte[0],&vbgsbhse[0],can);

                    TGraphErrors* tgrbgsbhse = nullptr;
                    MakeTGraphErrors(tgrbgsbhse,"bgsbhse","bgsbhse_Np_"+suf.str(),"N_{part}","<b_{#gamma}#sigma><b_{h}#sigma> error",vbgsbhse.size(),&vnpart[0],&vbgsbhse[0],&vnparte[0],0,can);

                    //==================================================

                    //==================================================
                    // Calculate K Corrections with Obtained Components
                    //==================================================

                    //-------------------- k_a --------------------

                    cout<<"Calculating "<<"ka_Np_"<<suf.str()<<endl;

                    vector<double> vka, vkae, vkanum, vkanume, vkaden, vkadene;
                    CalcK(hnpart,tgrag_Np,tgrah_Np,tgrage,tgrahe,vka,vkae,vkanum,vkanume,vkaden,vkadene);

                    TGraphErrors* tgrka = nullptr;
                    MakeTGraphErrors(tgrka,"ka","ka_Np_"+suf.str(),"N_{part}","#frac{<a_{#gamma} a_{h}>}{<a_{#gamma}><a_{h}>}",vka.size(),&vnpart[0],&vka[0],&vnparte[0],&vkae[0],can);

                    TGraphErrors* tgrkanum = nullptr;
                    MakeTGraphErrors(tgrkanum,"kanum","kanum_Np_"+suf.str(),"N_{part}","<a_{#gamma} a_{h}>",vkanum.size(),&vnpart[0],&vkanum[0],&vnparte[0],&vkanume[0],can);

                    TGraphErrors* tgrkaden = nullptr;
                    MakeTGraphErrors(tgrkaden,"kaden","kaden_Np_"+suf.str(),"N_{part}","<a_{#gamma}><a_{h}>",vkaden.size(),&vnpart[0],&vkaden[0],&vnparte[0],&vkadene[0],can);              

                    //-------------------- k_b --------------------

                    cout<<"Calculating "<<"kb_Np_"<<suf.str()<<endl;

                    vector<double> vkb, vkbe, vkbnum, vkbnume, vkbden, vkbdene;
                    CalcK(hnpart,tgrbg_Np,tgrbh_Np,tgrbge,tgrbhe,vkb,vkbe,vkbnum,vkbnume,vkbden,vkbdene,tgrbgs_Np,tgrbhs_Np,tgrbgse,tgrbhse);

                    TGraphErrors* tgrkb = nullptr;
                    MakeTGraphErrors(tgrkb,"kb","kb_Np_"+suf.str(),"N_{part}","#frac{<b_{#gamma} b_{h}>}{<b_{#gamma}#sigma><b_{h}#sigma>}",vkb.size(),&vnpart[0],&vkb[0],&vnparte[0],&vkbe[0],can);

                    TGraphErrors* tgrkbnum = nullptr;
                    MakeTGraphErrors(tgrkbnum,"kbnum","kbnum_Np_"+suf.str(),"N_{part}","<b_{#gamma} b_{h}>",vkbnum.size(),&vnpart[0],&vkbnum[0],&vnparte[0],&vkbnume[0],can);

                    TGraphErrors* tgrkbden = nullptr;
                    MakeTGraphErrors(tgrkbden,"kbden","kbden_Np_"+suf.str(),"N_{part}","<b_{#gamma}#sigma><b_{h}#sigma>",vkbden.size(),&vnpart[0],&vkbden[0],&vnparte[0],&vkbdene[0],can);

                    //-------------------- k_gh --------------------

                    cout<<"Calculating "<<"kgh_Np_"<<suf.str()<<endl;
                    vector<double> vkgh, vkghe, vkghnum, vkghnume, vkghden, vkghdene;
                    CalcK(hnpart,tgrag_Np,tgrbhs_Np,tgrage,tgrbhse,vkgh,vkghe,vkghnum,vkghnume,vkghden,vkghdene);

                    TGraphErrors* tgrkgh = nullptr;
                    MakeTGraphErrors(tgrkgh,"kgh","kgh_Np_"+suf.str(),"N_{part}","#frac{<a_{#gamma} b_{h} #sigma>}{<a_{#gamma}><b_{h}#sigma>}",vkgh.size(),&vnpart[0],&vkgh[0],&vnparte[0],&vkghe[0],can);

                    TGraphErrors* tgrkghnum = nullptr;
                    MakeTGraphErrors(tgrkghnum,"kghnum","kghnum_Np_"+suf.str(),"N_{part}","<a_{#gamma} b_{h} #sigma>",vkghnum.size(),&vnpart[0],&vkghnum[0],&vnparte[0],&vkghnume[0],can);

                    TGraphErrors* tgrkghden = nullptr;
                    MakeTGraphErrors(tgrkghden,"kghden","kghden_Np_"+suf.str(),"N_{part}","<a_{#gamma}><b_{h}#sigma>",vkghden.size(),&vnpart[0],&vkghden[0],&vnparte[0],&vkghdene[0],can);

                    //-------------------- k_hg --------------------

                    cout<<"Calculating "<<"khg_Np_"<<suf.str()<<endl;
                    vector<double> vkhg, vkhge, vkhgnum, vkhgnume, vkhgden, vkhgdene;
                    CalcK(hnpart,tgrag_Np,tgrbhs_Np,tgrage,tgrbhse,vkhg,vkhge,vkhgnum,vkhgnume,vkhgden,vkhgdene);

                    TGraphErrors* tgrkhg = nullptr;
                    MakeTGraphErrors(tgrkhg,"khg","khg_Np_"+suf.str(),"N_{part}","#frac{<a_{h} b_{#gamma} #sigma>}{<a_{h}><b_{#gamma}#sigma>}",vkhg.size(),&vnpart[0],&vkhg[0],&vnparte[0],&vkhge[0],can);

                    TGraphErrors* tgrkhgnum = nullptr;
                    MakeTGraphErrors(tgrkhgnum,"khgnum","khgnum_Np_"+suf.str(),"N_{part}","<a_{h} b_{#gamma} #sigma>",vkhgnum.size(),&vnpart[0],&vkhgnum[0],&vnparte[0],&vkhgnume[0],can);

                    TGraphErrors* tgrkhgden = nullptr;
                    MakeTGraphErrors(tgrkhgden,"khgden","khgden_Np_"+suf.str(),"N_{part}","<a_{h}><b_{#gamma}#sigma>",vkhgden.size(),&vnpart[0],&vkhgden[0],&vnparte[0],&vkhgdene[0],can);

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
    assert(("Could not find npartreso",tgr));
    cout<<"Found reaction plane resolution tgraph "<<tgr->GetName()<<endl;
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

vector<double> CorrectVector(vector<double>& vog, TGraph* tgr, TH1F* hnpart[GAT::ncent], vector<double>& vnpart, int cnt)
{
    cout<<"Iteration "<<cnt<<endl;
    if (cnt == 100)
    {
        cerr<<"no convergence after 100 attempts!"<<endl;
        return vector<double> {};
    }

    vector<double> vnew;
    double reldiffmax = 0.0;
    for (int icent = 0; icent < GAT::ncent; icent++)
    {

        double expval = 0.0;
        int nbinsx = hnpart[icent]->GetNbinsX();
        for (int ibin = 1; ibin <= nbinsx; ibin++)
        {
            double w = hnpart[icent]->GetBinContent(ibin);
            if (w == 0.0) continue;

            double npart = hnpart[icent]->GetBinCenter(ibin);

            double eval = tgr->Eval(npart);//TSpline
            expval += eval*w;
        }

        double reldiff = fabs(vog[icent]-expval)/expval;
        reldiffmax = max(reldiffmax,reldiff);

        double newval = (tgr->GetY())[icent]*(vog[icent]/expval);

        vnew.push_back(newval);
    }

    cout<<"reldiffmax = "<<reldiffmax<<endl;
    if (reldiffmax < 0.001)
    {
        cout<<"Finished correcting vector."<<endl;
        return vnew;
    }
    cout<<"replacing tgraph with new vector..."<<endl;

    delete tgr;
    tgr = new TGraph(vnew.size(),&vnpart[0],&vnew[0]);

    return CorrectVector(vog,tgr,hnpart,vnpart,++cnt);
}

//this method was manually copied to acceptCorrect.C

void CalcK(TH1F* hnpart[GAT::ncent], TGraphErrors* tgra, TGraphErrors* tgrb, TGraphErrors* tgrae, TGraphErrors* tgrbe, vector<double>& vk, vector<double>& vke, vector<double>& vnum, vector<double>& vnume, vector<double>& vden, vector<double>& vdene, TGraphErrors* tgrad = nullptr, TGraphErrors* tgrbd = nullptr, TGraphErrors* tgrade = nullptr, TGraphErrors* tgrbde = nullptr, bool usespline = true)
{
    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        cout<<"---------- icent = "<<icent<<" ----------"<<endl;
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

                be = tgrbe->Eval(npart);
                if (be < 0.) be = 0.;
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
            cout<<expa<<endl;

            expb += w*b;
            cout<<expb<<endl;

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

        vnum.push_back(num);
        vnume.push_back(nume);
        vden.push_back(den);
        vdene.push_back(dene);

        double k = num/den;
        vk.push_back(k);        

        double ke = k*sqrt(pow(nume/num,2)+pow(dene/den,2))/2.;//average between perfectly correlated and uncorrelated
        vke.push_back(ke);

        //what about the error in the denominator
        //but what about covariance? same term (tho weighted) is used in above loop
        //all weighted errors are actually based off this error and neighboring ones
        //forget it. it is just an approximation

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
    tgr->Print();
    tgr->Write();

    if (can)
    {
        oss.clear();
        oss.str("");
        oss<<dir<<"/"<<name<<".png";

        tgr->Draw("AL");
        can->SaveAs(oss.str().c_str());
    }
    return 0;
}
