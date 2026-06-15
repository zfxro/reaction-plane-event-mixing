//reads acceptance corrected gamtree results
//for each trigger type, ptt, cent, xi
//     - calculate (or rather isolate) each k component
//     - correct each k component
//     - write resulting TGraphs to file

//use limit for convergence instead of limit of 10. (USE RELATIVE LIMIT)

const int nadj = 2;

const int npartmin = 2;
const int npartmax = 390;

const int ncent = 2;
const int limcent[ncent+1] = {0,20,60};

const int markstyle[] = {20,21,33,47,34};

const float PI = (float) TMath::Pi();

namespace GAT {
    #include "gamAnaTree.h"
}

int GetNpartDists(TH1F* hnpart[GAT::ncent]);
int MergeNpartDists(TH1F* hnpartold[GAT::ncent], TH1F* hnpart[ncent+1], vector<double>& vnpart, vector<double>& vnparte);
int GetCentBin(int ic);

int MakeTGraphErrors(TGraphErrors*& tgr, string dir, string name, string xname, string yname, const int size, const double* x, const double* y, const double* ex = 0, const double* ey = 0, TCanvas* can = NULL);

void CalcK(TH1F* hnpart[ncent+1], TGraphErrors* tgra, TGraphErrors* tgrb, vector<double>& vk, vector<double>& vke, TGraphErrors* tgrae = nullptr, TGraphErrors* tgrbe = nullptr, TGraphErrors* tgrad = nullptr, TGraphErrors* tgrbd = nullptr, bool usespline = true);

void calcSigmaK(
    const char* inname = "kComponents.root",
    const char* outname = "sigmak_justdrawing.root"
    )
{
    gStyle->SetOptStat(0);
    gStyle->SetMarkerSize(2);
    gStyle->SetEndErrorSize(15);
    gROOT->ForceStyle();

    //Get Npart distributions
    TH1F* hnpartold[GAT::ncent];
    if (GetNpartDists(hnpartold)<0) return;

    //Get Npart distributions
    TH1F* hnpart[ncent] = {};
    vector<double> vnpart;
    vector<double> vnparte;
    if (MergeNpartDists(hnpartold,hnpart,vnpart,vnparte)<0) return;

    TFile* input = new TFile(inname,"READ");
    if (!input || input->IsZombie())
    {
        cerr<<"Could not open "<<outname<<endl;
        return;
    }

    TFile* output = new TFile(outname,"RECREATE");
    if (!output || output->IsZombie())
    {
        cerr<<"Could not open "<<outname<<endl;
        return;
    }

    TGraphErrors* tgrsig = nullptr;
    input->GetObject("npartsigma",tgrsig);
    tgrsig->Write();

    //calc average sigma in new cent bins
    //result will be written but not needed for rest of macro
    TGraphErrors* tgrexpsig = new TGraphErrors();
    for (int icent = 0; icent < ncent+1; icent++)
    {
        cout<<"Calculating expectation value of sigma in cent bin "<<icent<<endl;

        int nbinx = hnpart[icent]->GetNbinsX();

        double expsig = 0.;
        double wtot = 0.;

        for (int ibin = 1; ibin <= nbinx; ibin++)
        {
            double w = hnpart[icent]->GetBinContent(ibin);
            if (w <= 0.0) continue;

            wtot+=w;

            double npart = hnpart[icent]->GetBinCenter(ibin);

            double sig = tgrsig->Eval(npart);

            expsig += sig*w;
        }

        if (wtot!=1.0)
        {
            cerr<<"Ooops.  Weight total = "<<setprecision(10)<<wtot<<endl;
        }

        tgrexpsig->SetPoint(icent,icent,expsig);
    }

    tgrexpsig->SetTitle("expsig;icent;<#sigma>");
    tgrexpsig->SetName("expsig");
    output->WriteTObject(tgrexpsig);

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < GAT::ntype+1; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                ostringstream sufg;
                sufg<<"adj"<<iadj<<"_"<<GAT::spec[itype]<<"_ptt"<<iptt;

                TGraphErrors* tgrng = nullptr;
                input->GetObject(Form("ng_Np_adj%i_%s_ptt%i",iadj,GAT::spec[itype],iptt),tgrng);
                tgrng->Write();

                //error on ng not saved to input. must make here
                //it seems im passing a directory name (2nd arg) for no reason (canvas is null so not saving pic)
                TGraphErrors* tgrnge = nullptr;
                MakeTGraphErrors(tgrnge,"nge","nge_Np_"+sufg.str(),"<N_{part}>","<n_{t}> error",tgrng->GetN(),tgrng->GetX(),tgrng->GetEY(),tgrng->GetEX());

                vector<double> vk, vke;         
                CalcK(hnpart,tgrng,tgrsig,vk,vke,tgrnge);

                TGraphErrors* tgrk = nullptr;
                MakeTGraphErrors(tgrk,"ksigma","ksigma_Np_"+sufg.str(),"<N_{part}>","#frac{<n_{t} #sigma>}{<n_{t}><#sigma>}",vk.size(),&vnpart[0],&vk[0],&vnparte[0],&vke[0]);

                output->WriteTObject(tgrk);

                //draw tgrng and tgrsig same
                //see root cern transpad.C
                TCanvas* can = new TCanvas("can","",1200,1200);
                TPad* pad1 = new TPad("pad1","",0,0,1,1);
                TPad* pad2 = new TPad("pad2","",0,0,1,1);
                pad2->SetFillStyle(4000);
                pad1->Draw();
                pad1->cd();

                tgrng->SetTitle(";<N_{part}>;n_{t}");
                tgrng->SetMarkerStyle(markstyle[0]);
                tgrng->SetMarkerColor(TColor::GetColorPalette(440/3*0));
                tgrng->SetLineColor(TColor::GetColorPalette(440/3*0));
                tgrng->SetLineWidth(4);
                tgrng->GetYaxis()->SetLabelColor(TColor::GetColorPalette(440/3*0));
                tgrng->GetYaxis()->SetTitleColor(TColor::GetColorPalette(440/3*0));
                tgrng->GetYaxis()->SetAxisColor(TColor::GetColorPalette(440/3*0));
                tgrng->GetYaxis()->SetTitleOffset(1.2);
                double labsize = tgrng->GetYaxis()->GetLabelSize();
                short labfont = tgrng->GetYaxis()->GetLabelFont();
                double titsize = tgrng->GetYaxis()->GetTitleSize();
                short titfont = tgrng->GetYaxis()->GetTitleFont();
                tgrng->Draw("ap");

                double xmin = 0.;//or 0?
                double xmax = npartmax;
                tgrng->GetXaxis()->SetLimits(xmin,xmax);

                can->cd();
                double ymin = .425;
                double ymax = .9;
                double dy = (ymax-ymin)/.8;
                double dx = (xmax-xmin)/.8;
                pad2->Range(xmin-.1*dx,ymin-.1*dy,xmax+.1*dx,ymax+.1*dy);
                pad2->Draw();
                pad2->cd();

                tgrsig->SetTitle(";<N_{part}>;#sigma");
                tgrsig->SetMarkerStyle(markstyle[1]);
                tgrsig->SetMarkerColor(TColor::GetColorPalette(440/3*1));
                tgrsig->SetLineColor(TColor::GetColorPalette(440/3*1));
                tgrsig->GetXaxis()->SetLimits(xmin,xmax);
                tgrsig->SetLineWidth(4);
                tgrsig->Draw("p same");

                TGaxis* axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,508,"+L");
                axis->SetTitle("#sigma");
                axis->SetLabelSize(labsize);
                axis->SetLabelColor(TColor::GetColorPalette(440/3*1));
                axis->SetLabelFont(labfont);
                axis->SetTitleSize(titsize);
                axis->SetTitleOffset(1.2);
                axis->SetTitleColor(TColor::GetColorPalette(440/3*1));
                axis->SetTitleFont(titfont);
                axis->SetLineColor(TColor::GetColorPalette(440/3*1));           
                axis->Draw();

                TLegend leg (.33,.175,.73,.275);
                leg.SetTextSize(.04);
                TLegendEntry* legentry = leg.AddEntry(tgrng,Form("n_{t}  p_{T}#in [%1.1f,%1.1f)",GAT::limptt[iptt],GAT::limptt[iptt+1]),"lp");
                legentry->SetTextColor(tgrng->GetYaxis()->GetTitleColor());
                legentry = leg.AddEntry(tgrsig,"#sigma","lp");
                legentry->SetTextColor(axis->GetLabelColor());
                leg.Draw();

                can->SaveAs(Form("sigmak/ng_sigma_adj%i_%s_ptt%i.pdf",iadj,GAT::spec[itype],iptt));

                delete can;
                can = new TCanvas("can","",1200,1200);
                tgrk->SetMarkerStyle(markstyle[0]);
                // tgrk->SetTitle(";<N_{part}>;#frac{<n_{t} #sigma>}{<n_{t}><#sigma>}");
                tgrk->SetTitle(";<N_{part}>;<n_{t} #sigma>/<n_{t}><#sigma>");
                tgrk->Draw("AP");
                tgrk->GetHistogram()->GetXaxis()->SetLimits(xmin,xmax);
                tgrk->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
                can->Update();
                TLine line (xmin,1.0,xmax,1.0);
                line.SetLineStyle(7);
                line.Draw();
                TLegend leg2 (.45,.525,.85,.6);
                leg2.SetTextAlign(22);
                leg2.AddEntry(tgrk,Form("#kappa_{#sigma}  p_{T}#in [%1.1f,%1.1f)",GAT::limptt[iptt],GAT::limptt[iptt+1]),"lp");
                leg2.SetTextSize(.04);
                leg2.Draw();
                can->SaveAs(Form("sigmak/ksigma_adj%i_%s_ptt%i.pdf",iadj,GAT::spec[itype],iptt));
                delete can;

            }
        }
    }

    cout<<"Writing output to "<<output->GetName()<<endl;
    output->Close();
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

int GetNpartDists(TH1F* hnpart[GAT::ncent])
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
    }
    fcentdists->Close();
    return 0;
}

//rewriting this method
int MergeNpartDists(TH1F* hnpartold[GAT::ncent], TH1F* hnpart[ncent+1], vector<double>& vnpart, vector<double>& vnparte)
{
    //0-60%
    hnpart[ncent] = (TH1F*) hnpartold[0]->Clone(Form("AuAu_Npart_%i_%i",limcent[0],limcent[ncent]));    

    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        int centbin = GetCentBin(icent);
        double wt = GAT::limcent[icent+1] - GAT::limcent[icent];

        if (!hnpart[centbin])
        {
            hnpart[centbin] = (TH1F*) hnpartold[icent]->Clone(Form("AuAu_Npart_%i_%i",limcent[centbin],limcent[centbin+1]));
            hnpart[centbin]->Scale(wt);
        }
        else
            hnpart[centbin]->Add(hnpartold[icent],wt);

        if (icent == 0)
        {
            hnpart[ncent] = (TH1F*) hnpartold[0]->Clone(Form("AuAu_Npart_%i_%i",limcent[0],limcent[ncent]));
            hnpart[ncent]->Scale(wt);
        }
        else
            hnpart[ncent]->Add(hnpartold[icent],wt);

    }

    for (int icent = 0; icent < ncent+1; icent++)
    {
        double scl {-9999.};
        if (icent == ncent)
            scl = limcent[ncent] - limcent[0];
        else
            scl = limcent[icent+1] - limcent[icent];

        hnpart[icent]->Scale(1./scl);

        vnpart.push_back(hnpart[icent]->GetMean());     
        vnparte.push_back(hnpart[icent]->GetRMS());
    }

    return 0;
}

//      else

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

//this method was manually copied to acceptCorrect.C

void CalcK(TH1F* hnpart[ncent+1], TGraphErrors* tgra, TGraphErrors* tgrb, vector<double>& vk, vector<double>& vke, TGraphErrors* tgrae = nullptr, TGraphErrors* tgrbe = nullptr, TGraphErrors* tgrad = nullptr, TGraphErrors* tgrbd = nullptr, bool usespline = true)
{
    for (int icent = 0; icent < ncent+1; icent++)
    {
        cout<<"---------- icent = "<<icent<<" ----------"<<endl;
        int nbinx = hnpart[icent]->GetNbinsX();

        double k = 0.;
        double ke = 0.;

        //must recalculate exp values
        //differ to measured values
        //causing less than 1 k values
        double expa = 0.;
        double expb = 0.;

        for (int ibin = 1; ibin <= nbinx; ibin++)
        {
            double w = hnpart[icent]->GetBinContent(ibin);
            if (w <= 0.0) continue;

            double we = hnpart[icent]->GetBinError(ibin);

            double npart = hnpart[icent]->GetBinCenter(ibin);

            double a = -999.;
            double ae = -999.;
            double b = -999.;
            double be = -999.;

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

                if (tgrae)
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

                if (tgrbe)
                    be = tgrbe->Eval(npart);
                if (be < 0.) be = 0.;
            }
            else
            {
                cerr<<"THIS IS NO ELSE!"<<endl;
                return;
            }

            k += w*a*b;

            ke += pow(w*ae*b,2.) + pow(w*a*be,2.);

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
            }

            expa += w*a;

            expb += w*b;

        }

        cout<<"k = "<<k<<"/("<<expa<<"*"<<expb<<") = ";
        k /= (expa*expb);
        cout<<k<<endl;
        vk.push_back(k);

        ke = sqrt(ke);

        ke /= (expa*expb);
        vke.push_back(ke);

        //what about the error in the denominator
        //but what about covariance? same term (tho weighted) is used in above loop
        //all weighted errors are actually based off this error and neighboring ones
        //forget it. it is just an approximation

        cout<<icent<<" "<<k<<" +- "<<ke<<endl;
    }

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
        oss<<dir<<"/"<<name<<".pdf";

        tgr->Draw("AL");
        can->SaveAs(oss.str().c_str());
    }
    return 0;
}
