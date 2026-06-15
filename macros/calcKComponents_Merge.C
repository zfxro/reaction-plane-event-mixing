//I am working to adapt this macro now that gamanatree properly does the merging at beginning.
//OK. this macro reads dphi_*_centmerge_* histos from def acceptCorrect output. just read dphi_*_cent_* histos from acceptCorrect_merge output
//just changed input to *_merge.root and changed dphi name from centmerge to cent. done

//copied from calcKComponents.C
//this version is more true to its name. calculates just components not k values themselves
//this is for 0-20 and 20-60 centralities

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

const int ncent = 2;
const int limcent[ncent+1] = {0,20,60};

namespace GAT {
    #include "gamAnaTree.h"
}

int GetNpartDists(TH1F* hnpart[ncent], vector<double>& vnpart, vector<double>& vnparte);

int GetRPReso(vector<double>& vsig);
int Getng(TH1F* h1trgcent, vector<double>& vng, vector<double>& vnge, int icent, double nevt);

int MakeTGraphErrors(TGraphErrors*& tgr, string dir, string name, string xname, string yname, const int size, const double* x, const double* y, const double* ex = 0, const double* ey = 0, TCanvas* can = NULL);

void calcKComponents_Merge(
    bool quicksavetrigs = false,
    const char* outname = "kComponents_Merge.root",
    const char* inname = "gamtree_acceptcorrect_merge.root",
    const char* inname_gamhist = "gamAna_histos.root",
    const char* inname_trigs = "gamtree_triggers.root"
    )
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(12);
    gStyle->SetErrorX(0);
    gStyle->SetEndErrorSize(15);   

    //Get Npart merge distributions
    TH1F* hnpart[ncent];
    vector<double> vnpart;
    vector<double> vnparte;
    if (GetNpartDists(hnpart,vnpart,vnparte)<0) return;

    // Set up centrality vector for TGraphs
    vector<double> vcent;
    for (int icent = 0; icent < ncent; icent++)
    {
        double tmpcent = limcent[icent+1] + limcent[icent]+1;
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
    can->SaveAs("npartsigmamerge.eps");

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
    for (int icent = 0; icent < ncent; icent++)
    {
        int b0 = h1cent->FindBin(limcent[icent]+1);
        int b1 = h1cent->FindBin(limcent[icent+1]);
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
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
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

                for (int icent = 0; icent < ncent; icent++)
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
                can->SaveAs(Form("ngmerge/ng_Np_adj%i_%s_ptt%i.eps",iadj,GAT::spec[itype],iptt));

                //plot a_g too. why not?
                vector<double> vag (vng);
                transform(vag.begin(),vag.end(),vag.begin(), [](double x) { return x / (2.*PI); } );
                vector<double> vage (vnge);
                transform(vage.begin(),vage.end(),vage.begin(), [](double x) { return x / (2.*PI); } );
                // tgrag_Np->SetTitle(Form("ag_Np_adj%i_%s_ptt%i;N_{part};<a_{#gamma}>",iadj,GAT::spec[itype],iptt));
                TGraphErrors* tgrag_Np = nullptr;
                MakeTGraphErrors(tgrag_Np,"agmerge","ag_Np_"+sufg.str(),"N_{part}","<a_{#gamma}>",vag.size(),&vnpart[0],&vag[0],&vnparte[0],&vage[0],can);

                TGraphErrors* tgrage = nullptr;
                MakeTGraphErrors(tgrage,"agemerge","age_Np_"+sufg.str(),"N_{part}","<a_{#gamma}> error",vag.size(),&vnpart[0],&vage[0],&vnparte[0],0,can);

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

                // tgrng_Np_C->SetTitle(Form("ng_Np_C_adj%i_%s_ptt%i;N_{part};<n_{#gamma}>",iadj,GAT::spec[itype],iptt));

                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {

                    ostringstream suf;
                    suf<<"adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt<<"_xi"<<ixi;

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

                    for (int icent = 0; icent < ncent; icent++)
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

                        can->SaveAs(Form("bgsbhsmerge/bgsbhs_adj%i_%s_ptt%i_cent%i_xi%i.eps",iadj,GAT::spec[itype],iptt,icent,ixi));

                        double bgsbhs = -9999., bgsbhse = -9999.;
                        bgsbhs = f1all->GetParameter(1)/(2.*PI);
                        bgsbhse = f1all->GetParError(1)/(2.*PI);

                        vbgsbhs.push_back(bgsbhs);
                        vbgsbhse.push_back(bgsbhse);

                        //==================== <b_g s> ====================

                        TH1F* h1Ediff = (TH1F*) h1Ein->Clone(Form("Ediff_adj%i_%s_ptt%i_cent%i_xi%i",iadj,GAT::spec[itype],iptt,icent,ixi));
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

                        TF1* f1 = new TF1(Form("bhscos_adj%i_%s_ptt%i_cent%i_xi%i",iadj,GAT::spec[itype],iptt,icent,ixi),"[0]+[1]*cos(x*2.*TMath::Pi())",0.,1.);
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

                        can->SaveAs(Form("bhscosmerge/bhscos_adj%i_%s_ptt%i_cent%i_xi%i.eps",iadj,GAT::spec[itype],iptt,icent,ixi));

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
                    tgrah->SetTitle(Form("ah_adj%i_%s_ptt%i_xi%i;centrality;<a_h}>",iadj,GAT::spec[itype],iptt,ixi));
                    tgrah->Write(Form("ah_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));

                    TGraphErrors* tgrah_Np = nullptr;
                    MakeTGraphErrors(tgrah_Np,"ahmerge","ah_Np_"+suf.str(),"N_{part}","<a_{h}>",vah.size(),&vnpart[0],&vah[0],&vnparte[0],&vahe[0],can);

                    TGraphErrors* tgrahe = nullptr;
                    MakeTGraphErrors(tgrahe,"ahemerge","ahe_Np_"+suf.str(),"N_{part}","<a_{h}> error",vah.size(),&vnpart[0],&vahe[0],&vnparte[0],0,can);

                    //==================================================

                    //==================== <b_g sigma> ====================

                    TGraph* tgrbgs = new TGraph(vbgs.size(),&vcent[0],&vbgs[0]);
                    tgrbgs->SetTitle(Form("bgs_adj%i_%s_ptt%i_xi%i;centrality;<b_{#gamma} #sigma>",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbgs->SetName(Form("bgs_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbgs->Write();

                    TGraphErrors* tgrbgs_Np = nullptr;
                    MakeTGraphErrors(tgrbgs_Np,"bgsmerge","bgs_Np_"+suf.str(),"N_{part}","<b_{#gamma} #sigma>",vbgs.size(),&vnpart[0],&vbgs[0],&vnparte[0],&vbgse[0],can);

                    TGraphErrors* tgrbgse = nullptr;
                    MakeTGraphErrors(tgrbgse,"bgsemerge","bgse_Np_"+suf.str(),"N_{part}","<b_{#gamma} #sigma> error",vbgse.size(),&vnpart[0],&vbgse[0],&vnparte[0],0,can);

                    //==================================================

                    //==================== <b_g> ====================

                    TGraph* tgrbg = new TGraph(vbg.size(),&vcent[0],&vbg[0]);
                    tgrbg->SetTitle(Form("bg_adj%i_%s_ptt%i_xi%i;centrality;<b_{#gamma} #sigma>/<#sigma>",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbg->SetName(Form("bg_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbg->Write();

                    TGraphErrors* tgrbg_Np = nullptr;
                    MakeTGraphErrors(tgrbg_Np,"bgmerge","bg_Np_"+suf.str(),"N_{part}","<b_{#gamma}>",vbg.size(),&vnpart[0],&vbg[0],&vnparte[0],&vbge[0],can);

                    TGraphErrors* tgrbge = nullptr;
                    MakeTGraphErrors(tgrbge,"bgemerge","bge_Np_"+suf.str(),"N_{part}","<b_{#gamma}> error",vbge.size(),&vnpart[0],&vbge[0],&vnparte[0],0,can);

                    //==================================================

                    //==================== <v_g> ====================

                    TGraph* tgrvg = new TGraph(vvg.size(),&vcent[0],&vvg[0]);
                    tgrvg->SetTitle(Form("vg_adj%i_%s_ptt%i_xi%i;centrality;#frac{<b_{#gamma}#sigma>}{<#sigma><a_{#gamma}>} #approx v_{#gamma}",iadj,GAT::spec[itype],iptt,ixi));
                    tgrvg->SetName(Form("vg_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));
                    tgrvg->Write();

                    TGraphErrors* tgrvg_Np = nullptr;
                    MakeTGraphErrors(tgrvg_Np,"vgmerge","vg_Np_"+suf.str(),"N_{part}","<v_{#gamma}>",vvg.size(),&vnpart[0],&vvg[0],&vnparte[0],&vvge[0],can);

                    //==================================================

                    //==================== <b_h sigma> ====================

                    TGraph* tgrbhs = new TGraph(vbhs.size(),&vcent[0],&vbhs[0]);
                    tgrbhs->SetTitle(Form("bhs_adj%i_%s_ptt%i_xi%i;centrality;<b_{h} #sigma>",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbhs->Write(Form("bhs_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));

                    TGraphErrors* tgrbhs_Np = nullptr;
                    MakeTGraphErrors(tgrbhs_Np,"bhsmerge","bhs_Np_"+suf.str(),"N_{part}","<b_{h} #sigma>",vbhs.size(),&vnpart[0],&vbhs[0],&vnparte[0],&vbhse[0],can);

                    TGraphErrors* tgrbhse = nullptr;
                    MakeTGraphErrors(tgrbhse,"bhsemerge","bhse_Np_"+suf.str(),"N_{part}","<b_{h} #sigma> error",vbhse.size(),&vnpart[0],&vbhse[0],&vnparte[0],0,can);

                    //==================================================

                    //==================== <b_h> ====================

                    TGraph* tgrbh = new TGraph(vbh.size(),&vcent[0],&vbh[0]);
                    tgrbh->SetTitle(Form("bh_adj%i_%s_ptt%i_xi%i;centrality;<b_{h} #sigma>/<#sigma>",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbh->Write(Form("bh_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));

                    TGraphErrors* tgrbh_Np = nullptr;
                    MakeTGraphErrors(tgrbh_Np,"bhmerge","bh_Np_"+suf.str(),"N_{part}","<b_{h}>",vbh.size(),&vnpart[0],&vbh[0],&vnparte[0],&vbhe[0],can);

                    TGraphErrors* tgrbhe = nullptr;
                    MakeTGraphErrors(tgrbhe,"bhemerge","bhe_Np_"+suf.str(),"N_{part}","<b_{h}> error",vbhe.size(),&vnpart[0],&vbhe[0],&vnparte[0],0,can);

                    //==================================================

                    //==================== <v_h> ====================

                    TGraph* tgrvh = new TGraph(vvh.size(),&vcent[0],&vvh[0]);
                    tgrvh->SetTitle(Form("vh_adj%i_%s_ptt%i_xi%i;centrality;#frac{<b_{h}#sigma>}{<#sigma><a_{h}>} #approx v_{h}",iadj,GAT::spec[itype],iptt,ixi));
                    tgrvh->Write(Form("vh_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));

                    TGraphErrors* tgrvh_Np = nullptr;
                    MakeTGraphErrors(tgrvh_Np,"vhmerge","vh_Np_"+suf.str(),"N_{part}","<v_{h}>",vvh.size(),&vnpart[0],&vvh[0],&vnparte[0],&vvhe[0],can);

                    //==================================================

                    //==================== <b_g sigma> <b_h sigma> ====================

                    TGraph* tgrbgsbhs = new TGraph(vbgsbhs.size(),&vcent[0],&vbgsbhs[0]);
                    tgrbgsbhs->SetTitle(Form("bgsbhs_adj%i_%s_ptt%i_xi%i;centrality;<b_{#gamma}#sigma><b_{h}#sigma>",iadj,GAT::spec[itype],iptt,ixi));
                    tgrbgsbhs->Write(Form("bgsbhs_adj%i_%s_ptt%i_xi%i",iadj,GAT::spec[itype],iptt,ixi));

                    TGraphErrors* tgrbgsbhs_Np = nullptr;
                    MakeTGraphErrors(tgrbgsbhs_Np,"bgsbhsmerge","bgsbhs_Np_"+suf.str(),"N_{part}","<b_{#gamma}#sigma><b_{h}#sigma>",vbgsbhs.size(),&vnpart[0],&vbgsbhs[0],&vnparte[0],&vbgsbhse[0],can);

                    TGraphErrors* tgrbgsbhse = nullptr;
                    MakeTGraphErrors(tgrbgsbhse,"bgsbhsemerge","bgsbhse_Np_"+suf.str(),"N_{part}","<b_{#gamma}#sigma><b_{h}#sigma> error",vbgsbhse.size(),&vnpart[0],&vbgsbhse[0],&vnparte[0],0,can);

                    //==================================================

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
    for (int jcent = limcent[icent+1]; jcent > limcent[icent]; jcent--)
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

int GetNpartDists(TH1F* hnpart[ncent], vector<double>& vnpart, vector<double>& vnparte)
{
    const char* centdistsname = "CentDists_2bins.root";
    TFile* fcentdists = new TFile(centdistsname,"READ");
    if (!fcentdists || fcentdists->IsZombie())
    {
        cerr<<"Could not open "<<centdistsname<<endl;
        return -1;
    }

    for (int icent = 0; icent < ncent; icent++)
    {
        hnpart[icent] = (TH1F*) gROOT->FindObject(Form("AuAu_Npart_%i_%i",limcent[icent],limcent[icent+1]));
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

    TGraphErrors* tgr = (TGraphErrors*) gROOT->FindObject("centresomerge");
    assert(("Could not find centresomerge",tgr));
    cout<<"Found reaction plane resolution tgraph "<<tgr->GetName()<<endl;
    tgr->Print();

    double* yary = tgr->GetY();

    for (int icent = 0; icent < ncent; icent++)
    {
        vsig.push_back(yary[icent]);
    }

    input->Close();
    delete input;

    return 0;
}

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
        oss<<dir<<"/"<<name<<".eps";

        tgr->Draw("AL");
        can->SaveAs(oss.str().c_str());
    }
    return 0;
}
