#include "gamAnaTree.h"

ofstream ofs;
ofstream ofs_other;
bool watcher= false;
int ipt;
int icent;
const float nsigma = 3.0;
bool drawmany = false;
TH2F* h2ptcent;

Double_t background(Double_t *x, Double_t *par) {
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t gaussian(Double_t *x,Double_t *par)
{
    Double_t arg = 0;
    if (par[2] != 0) arg = (x[0] - par[1])/par[2];
    Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}

Double_t fitFunction(Double_t *x, Double_t *par) {
    return background(x,par) + gaussian(x,&par[3]);
}

void Fitting(TH1F *histo, TCanvas* can, int loc = 1)
{
    TVirtualPad* tvp = can->cd(loc);

    double amp = histo->GetMaximum();
    double mu = histo->GetBinCenter(histo->GetMaximumBin());
    double sig = .01;

    cout<<"amp = "<<amp<<"; mu = "<<mu<<"; sig = "<<sig<<endl;

    int ibin = histo->FindBin(mu-3.*sig);
    int maxbin = histo->FindBin(mu+3.*sig);

    cout<<"ibin = "<<ibin<<"; maxbin = "<<maxbin<<endl;

    TH1F* hclone = (TH1F*)histo->Clone("clone");

    if (watcher)
    {
        hclone->Draw();
        cout<<"Histogram clone has been drawn.  Press return to continue."<<endl;
        tvp->Update();
        gSystem->ProcessEvents();
        cin.ignore();
    }

    for (;ibin<=maxbin;ibin++)
    {
        hclone->SetBinContent(ibin,0);
        hclone->SetBinError(ibin,0);
    }

    if (watcher)
    {
        hclone->Draw();
        cout<<"pi0 peak bins set to 0.  Press return to continue."<<endl;
        tvp->Update();
        gSystem->ProcessEvents();
        cin.ignore();
    }

    double polyXi = .05;
    double polyXf = .3;
    TF1* iniPoly = new TF1("iniPoly",background,polyXi,polyXf,3);
    iniPoly->SetParameters(amp/100.,10000.,-10000.);
    hclone->Fit("iniPoly","R");

    if (watcher)
    {
        cout<<"Background has been fit with 2nd order polynomial.  Press return to continue."<<endl;
        tvp->Update();
        gSystem->ProcessEvents();
        cin.ignore();
    }

    double polyPars[3];
    iniPoly->GetParameters(polyPars);
    amp -= iniPoly->Eval(mu);
    TF1 *fitFcn = new TF1("fitFcn",fitFunction,.05,.3,6);
    fitFcn->SetParNames("p0","p1","p2","Constant","Mean","Sigma");
    fitFcn->SetLineColor(kMagenta);
    fitFcn->SetParameters(polyPars[0],polyPars[1],polyPars[2],amp,mu,sig);

    histo->Fit("fitFcn","R","ep");

    if (watcher)
    {
        cout<<"Background and peak have been fit.  Press return to continue."<<endl;
        tvp->Update();
        gSystem->ProcessEvents();
        cin.ignore();
    }

    TF1 *backFcn = new TF1("backFcn",background,0,3,3);
    backFcn->SetLineColor(kRed);

    TF1 *signalFcn = new TF1("signalFcn",gaussian,0,3,3);
    signalFcn->SetLineColor(kBlue);
    Double_t par[6];

    fitFcn->GetParameters(par);

    backFcn->SetParameters(par);
    backFcn->Draw("same");

    signalFcn->SetParameters(&par[3]);
    signalFcn->Draw("same");

    TLegend *legend=new TLegend(0.70,0.75,.99,.95);
    legend->SetTextFont(72);
    legend->SetTextSize(0.04);
    legend->SetFillColor(0);
    legend->AddEntry(histo,"Data","lpe");
    legend->AddEntry(backFcn,"Background fit","l");
    legend->AddEntry(signalFcn,"Signal fit","l");
    legend->AddEntry(fitFcn,"Global Fit","l");
    legend->Draw();

    double integral = par[3]*sqrt(2*TMath::Pi())*par[5];

    TPaveStats *ps = new TPaveStats();

    tvp->Update();
    ps = (TPaveStats*)tvp->GetPrimitive("stats");
    ps->SetName("mystats");

    ps->SetX1NDC(.70);
    ps->SetX2NDC(.99);
    ps->SetY1NDC(.35);
    ps->SetY2NDC(.75);

    TList *list = ps->GetListOfLines();

    if (watcher)
    {
        cout<<"Fitting method has ended.  Press return to continue."<<endl;
        tvp->Update();
        gSystem->ProcessEvents();
        cin.ignore();
    }

    float mean = fitFcn->GetParameter("Mean");
    float sigma = fitFcn->GetParameter("Sigma");
    mean = .138;
    sigma = .01;
    float xlo = mean - nsigma*sigma;
    float xhi = mean + nsigma*sigma;
    float bgintgrl = par[0]*(xhi-xlo) + (1./2.)*par[1]*(xhi*xhi - xlo*xlo) + (1./3.)*par[2]*(xhi*xhi*xhi - xlo*xlo*xlo);
    float ngam = h2ptcent->GetBinContent(ipt+1,icent+1);
    ofs<<ipt<<" "<<icent<<" "<<bgintgrl<<" "<<ngam<<" "<<bgintgrl/ngam<<endl;
    float sgnlintgrl = signalFcn->Integral(xlo,xhi);
    ofs_other<<ipt<<" "<<icent<<" "<<bgintgrl<<" "<<ngam<<" "<<bgintgrl/ngam<<" "<<integral<<" "<<sgnlintgrl<<" "<<sgnlintgrl/integral<<endl;

    TLatex *tltot = new TLatex(0,0,Form("#gamma total = %.0f",ngam));
    tltot->SetTextFont(72);
    tltot->SetTextSize(0.04);
    list->Add(tltot);

    TLatex *tlmis = new TLatex(0,0,Form("#gamma mistagged = %.0f",bgintgrl));
    tlmis->SetTextFont(72);
    tlmis->SetTextSize(0.04);
    list->Add(tlmis);

    TLatex *tlm = new TLatex(0,0,Form("m = %.5f",bgintgrl/ngam));
    tlm->SetTextFont(72);
    tlm->SetTextSize(0.04);
    list->Add(tlm);

    ps->Draw();
    histo->SetStats(0);
    tvp->Modified();

    return;
}


void drawPi0Plots_ptcent(
    const char* input
    )
{
    cout<<"WARNING!!! nsigma = "<<nsigma<<"  Make sure this agrees with gamAna.cc"<<endl;

    cout<<"Using ncent from gamAnaTree. do not redefine. ncent = "<<ncent<<endl;

    gStyle->SetOptStat(0);
    gStyle->SetStatFont(72);
    gStyle->SetStatFontSize(0.02);
    gStyle->SetTitleFontSize(0.035);
    gStyle->SetOptFit(0011);
    gStyle->SetStatW(.25);
    gStyle->SetStatX(.99);
    gStyle->SetStatY(.75);
    gStyle->SetPadRightMargin(0.01);
    gStyle->SetPadBottomMargin(0.075);
    gStyle->SetPadTopMargin(0.05);

    ofs.open("pi0bg_ptcent.txt");
    ofs_other.open("pi0bg_ptcent_other.txt");

    TFile in(input);

    h2ptcent = (TH2F*) gROOT->FindObject("ptcent_rebin");

    for (ipt = 0; ipt < nptt; ipt++)
    {
        TCanvas *can = new TCanvas(Form("can%i",ipt),"",1000,600);
        if (drawmany)
            can->Divide(2,2,0.001,0.001);
        if (ncent > 4)
        {
            cout<<__LINE__<<" please fix canvas division. ncent too large."<<endl;
            return;
        }

        for (icent = 0; icent < ncent+1; icent++)
        {
            TH1F* h1invm = (TH1F*) gROOT->FindObject(Form("invm_rebin_cent%i_pt%i",icent,ipt));

            int locent;
            int hicent;
            if (icent == ncent)
            {
                locent = limcent[0];
                hicent = limcent[ncent];
            }
            else
            {
                locent = limcent[icent];
                hicent = limcent[icent+1];
            }


            h1invm->SetTitle(Form("%.1f < p_{T} < %.1f   %d-%d%%;m_{pair};N_{pair}",limptt[ipt],limptt[ipt+1],locent,hicent));

            h1invm->Rebin(4);
            float binw = h1invm->GetBinWidth(1);
            h1invm->Scale(1./binw);

            if (drawmany)
                Fitting(h1invm,can,icent+1);
            else
            {
                Fitting(h1invm,can);
                can->SaveAs(Form("pi0fit_pt%i_cent%i.pdf",ipt,icent));
            }

        }

        if (drawmany)
            can->SaveAs(Form("pi0fit_pt%i.pdf",ipt));

    }
    cout<<"\nAll Done!\n";
}
