//Construct decay photon integrated away side yield (IASY) using:
//
//     - parent pt of measured decay photons (from simulation, inclusive)
//     - iasy of pi0 as function of trigger pt
//
//The pi0 iasy vs ptt graphs are fit with pol1 and pol2.
//These functions, along with just using bincontents (and extending first and last bin),
//     are used to develop decay iasy.

const int ncent = 2;
const int limcent[ncent+1][2] = {
    {0,20},
    {20,60},
    {0,60}
};

namespace GAT {
#include "gamAnaTree.h"
}

TFile* OpenFile(const char* name);

int GetPtBin(double pt)
{
    for (int ipt = 0; ipt < GAT::nptt; ipt++)
    {
        if (pt <= GAT::limptt[ipt+1])
            return ipt;
    }
    return GAT::nptt - 1;
}

void makeDecayIASY(
    const char* inname_iasy = "integawaysideyield.root",
    const char* inname_parentpt = "parentptana.root",
    const char* outname = "iasy_dec.root"
    )
{
    cout<<"Making decay photon integrated away side yield..."<<endl;

    TFile* input_iasy = OpenFile(inname_iasy);
    TFile* input_parentpt = OpenFile(inname_parentpt);

    TFile* output = new TFile(outname,"RECREATE");

    TCanvas can ("can","",1500,800);
    gStyle->SetOptFit(0);
    gStyle->SetPalette(kBird);
    gStyle->SetEndErrorSize(7);

    TH1D* h1parentpt[ncent+1][GAT::nptt] = {};

    for (int icent = 0; icent <= ncent; icent++)
    {
        for (int iptt = 0; iptt < GAT::nptt; iptt++)
        {
            TH1D* h1tagged;
            input_parentpt->GetObject(Form("tagged_cent%i_pt%i",icent,iptt),h1tagged);

            TH1D* h1missed;
            input_parentpt->GetObject(Form("missed_cent%i_pt%i",icent,iptt),h1missed);

            h1parentpt[icent][iptt] = (TH1D*) h1tagged->Clone(Form("parentpt_cent%i_pt%i",icent,iptt));
            h1parentpt[icent][iptt]->Add(h1missed);
        }
    }
    for (int idpsi = 0; idpsi <= 2; idpsi++)
    {
        for (int ixi = 0; ixi < GAT::nxi; ixi++)
        {
            for (int icent = 0; icent <= ncent; icent++)
            {
                TGraphErrors* iasy_stat;
                input_iasy->GetObject(Form("ipttyield_adj0_pi0_cent%i_dpsi%i_xi%i",icent,idpsi,ixi),iasy_stat);

                TGraphErrors* iasy_sys;
                input_iasy->GetObject(Form("ipttyield_syserr_adj0_pi0_cent%i_dpsi%i_xi%i",icent,idpsi,ixi),iasy_sys);

                TGraphErrors* iasy = iasy_stat;

                iasy->SetMarkerStyle(kFullCircle);
                iasy->SetMarkerSize(2.0);
                iasy->SetLineWidth(1);
                iasy->SetMarkerColor(TColor::GetColorPalette(0));
                iasy->SetLineColor(TColor::GetColorPalette(0));
                iasy->SetTitle(";#pi^{0} p_{T};1/N_{t} dN/d#xi");
                iasy->Draw("AP");

                iasy->GetYaxis()->SetTitleOffset(1.4);
                iasy->GetXaxis()->SetTitleOffset(1.3);

                //----- set ranges -----
                double xmin = 3.;
                double xmax = 30.;
                iasy->GetXaxis()->SetLimits(xmin,xmax);

                double ymin = iasy->GetYaxis()->GetXmin();
                double ymax = iasy->GetYaxis()->GetXmax();

                //----- fit with straight line -----
                iasy->Fit("pol1","0");
                TF1* tmp = iasy->GetFunction("pol1");
                TF1* pol1 = (TF1*) tmp->Clone();
                pol1->SetLineColor(TColor::GetColorPalette(80));

                //draw fit function
                pol1->SetLineWidth(1);
                pol1->Draw("same");

                //----- set y range to include fit -----
                ymin = min(ymin,pol1->GetMinimum(xmin,xmax));
                ymax = max(ymax,pol1->GetMaximum(xmin,xmax));

                iasy->GetYaxis()->SetRangeUser(ymin,ymax);

                //---------- Confidence Interval (1 sigma) ----------
                TGraphErrors* grint = new TGraphErrors(GAT::nptt);
                for (int i=0; i<GAT::nptt; i++)
                    grint->SetPoint(i, iasy->GetX()[i], 0);
                (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,.68);

                //----- set minx max graphs -----
                TGraph* grmin = new TGraphErrors(GAT::nptt);
                TGraph* grmax = new TGraphErrors(GAT::nptt);
                for (int i=0; i<GAT::nptt; i++) {
                    grmin->SetPoint(i,grint->GetX()[i],grint->GetY()[i]-grint->GetEY()[i]);
                    grmax->SetPoint(i,grint->GetX()[i],grint->GetY()[i]+grint->GetEY()[i]);
                }

                //----- extend min max graphs over entire pt range by using eval then draw -----
                TF1* mineval = new TF1("mineval",[&](double* x,double* p) {return grmin->Eval(x[0]);},xmin,xmax,0);
                mineval->SetLineWidth(1);
                mineval->Draw("same");
                mineval->SetLineColor(TColor::GetColorPalette(120));

                TF1* maxeval = new TF1("maxeval",[&](double* x,double* p) {return grmax->Eval(x[0]);},xmin,xmax,0);
                maxeval->SetLineWidth(1);
                maxeval->Draw("same");
                maxeval->SetLineColor(TColor::GetColorPalette(120));

                //----- confidence interval band -----
                int nx = 100;
                TGraph* grshade = new TGraph(2*nx);
                for (int ix=0; ix<nx; ix++)
                {
                    double x = xmin + ix*(xmax-xmin)/(nx-1.);
                    grshade->SetPoint(ix,x,grmax->Eval(x));
                    x = xmax - ix*(xmax-xmin)/(nx-1.);
                    grshade->SetPoint(nx+ix,x,grmin->Eval(x));
                }

                grshade->SetFillStyle(1001);
                grshade->SetFillColorAlpha(TColor::GetColorPalette(100),.1);
                grshade->Draw("fsame");

                //----- step function -----
                TF1* stepfun = new TF1("stepfun",[&](double* x,double* p) {return iasy->GetY()[GetPtBin(x[0])];},xmin,xmax,0);
                stepfun->SetLineWidth(1);
                stepfun->Draw("same");
                stepfun->SetLineColor(TColor::GetColorPalette(200));

                //----- draw legend -----
                TLegend* leg = new TLegend(0.7,0.8375,0.9,0.9);
                leg->SetTextSize(.035);
                leg->SetTextAlign(22);
                leg->AddEntry(iasy,Form("%.1f #leq  #xi < %.1f",GAT::limxi[ixi],GAT::limxi[ixi+1]),"P");
                leg->Draw();

                can.SaveAs(Form("pics/iasy_ptt_pi0_cent%i_dpsi%i_xi%i.pdf",icent,idpsi,ixi));

                //------------------------------------------------------------
                // decay yield using pol1 fit
                TGraphErrors* decY_pol1 = new TGraphErrors(GAT::nptt);

                decY_pol1->SetTitle(iasy->GetTitle());

                // dec yield using step function
                TGraphErrors* decY_step = new TGraphErrors(GAT::nptt);

                // combine step and pol1
                TGraphErrors* decY_combo = new TGraphErrors(GAT::nptt);

                TGraphErrors* decY_sys = new TGraphErrors(GAT::nptt);
                decY_sys->SetTitle(iasy->GetTitle());

                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    double totparent = 0.;
                    double sumyield = 0.;
                    double sumyieldmin = 0.;
                    double sumyieldmax = 0.;

                    double sumyieldstep = 0.;

                    double sumstaterr = 0.;
                    double sumsyserr = 0.;

                    for (int iptp = 1; iptp <= h1parentpt[icent][iptt]->GetNbinsX(); iptp++)
                    {
                        double ptp = h1parentpt[icent][iptt]->GetBinCenter(iptp);
                        if (ptp < GAT::limptt[iptt]-3. || ptp > 30.) continue;

                        double nparent = h1parentpt[icent][iptt]->GetBinContent(iptp);
                        totparent += nparent;

                        double yield = pol1->Eval(ptp);
                        sumyield += yield*nparent;

                        double yieldmin = mineval->Eval(ptp);
                        sumyieldmin += yieldmin*nparent;

                        double yieldmax = maxeval->Eval(ptp);
                        sumyieldmax += yieldmax*nparent;

                        double yieldstep = stepfun->Eval(ptp);
                        sumyieldstep += yieldstep*nparent;

                        double staterr = iasy_stat->GetEY()[GetPtBin(ptp)];
                        sumstaterr += staterr*nparent;

                        double syserr = iasy_sys->GetEY()[GetPtBin(ptp)];
                        sumsyserr += syserr*nparent;
                    }

                    double avgyield = sumyield/totparent;
                    double avgyieldmin = sumyieldmin/totparent;
                    double avgyieldmax = sumyieldmax/totparent;

                    decY_pol1->SetPoint(iptt,(GAT::limptt[iptt]+GAT::limptt[iptt+1])/2.-.2,avgyield);
                    double avgerr = ( abs(avgyieldmax-avgyield) + abs(avgyieldmin-avgyield) )/2.;
                    decY_pol1->SetPointError(iptt,0,avgerr);

                    // step function
                    double avgyieldstep = sumyieldstep/totparent;
                    decY_step->SetPoint(iptt,(GAT::limptt[iptt]+GAT::limptt[iptt+1])/2.+.2,avgyieldstep);

                    double avgstaterr = sumstaterr/totparent;
                    decY_step->SetPointError(iptt,0,avgstaterr);

                    double comboyield = (avgyield + avgyieldstep)/2.;
                    decY_combo->SetPoint(iptt,(GAT::limptt[iptt]+GAT::limptt[iptt+1])/2.,comboyield);
                    double comboerr = 0.;
                    comboerr += fabs(comboyield - avgyield - avgerr);
                    comboerr += fabs(comboyield - avgyield + avgerr);
                    comboerr += fabs(comboyield - avgyieldstep - avgstaterr);
                    comboerr += fabs(comboyield - avgyieldstep + avgstaterr);
                    comboerr /= 4.;
                    decY_combo->SetPointError(iptt,0,comboerr);

                    decY_sys->SetPoint(iptt,(GAT::limptt[iptt]+GAT::limptt[iptt+1])/2.,comboyield);

                    double avgsyserr = sumsyserr/totparent;

                    double syserr = sqrt( avgsyserr*avgsyserr + pow((avgyield-avgyieldstep)/2.,2.) );
                    decY_sys->SetPointError(iptt,.06,syserr);
                }

                TMultiGraph* mgr = new TMultiGraph();

                decY_pol1->SetMarkerStyle(kFullCircle);
                decY_pol1->SetMarkerSize(1.5);
                decY_pol1->SetLineWidth(2);
                decY_pol1->SetMarkerColor(TColor::GetColorPalette(80));
                decY_pol1->SetLineColor(TColor::GetColorPalette(80));

                decY_step->SetMarkerStyle(kFullCircle);
                decY_step->SetMarkerSize(1.5);
                decY_step->SetLineWidth(2);
                decY_step->SetMarkerColor(TColor::GetColorPalette(200));
                decY_step->SetLineColor(TColor::GetColorPalette(200));

                decY_combo->SetMarkerStyle(kFullCircle);
                decY_combo->SetMarkerSize(1.5);
                decY_combo->SetLineWidth(2);
                decY_combo->SetMarkerColor(TColor::GetColorPalette(140));
                decY_combo->SetLineColor(TColor::GetColorPalette(140));

                decY_sys->SetMarkerSize(0);
                decY_sys->SetLineWidth(0);
                decY_sys->SetFillColorAlpha(TColor::GetColorPalette(200),.50);

                mgr->Add(decY_sys,"2");
                decY_pol1->SetLineWidth(1);
                mgr->Add(decY_pol1,"P");
                decY_step->SetLineWidth(1);
                mgr->Add(decY_step,"P");
                decY_combo->SetLineWidth(1);
                mgr->Add(decY_combo,"P");

                decY_combo->Write(Form("ipttyield_adj2_dec_cent%i_dpsi%i_xi%i",icent,idpsi,ixi));
                decY_sys->Write(Form("ipttyield_syserr_adj2_dec_cent%i_dpsi%i_xi%i",icent,idpsi,ixi));

                mgr->SetTitle(";#gamma_{dec} p_{T};1/N_{t} dN/d#xi");
                mgr->Draw("A");

                mgr->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
                mgr->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);

                TLatex* txtall = new TLatex(0.5,.85,Form("Au+Au %i-%i%% 200GeV       %.1f #leq #xi < %.1f",limcent[icent][0],limcent[icent][1],GAT::limxi[ixi],GAT::limxi[ixi+1]));
                txtall->SetTextSize(.035);
                txtall->SetTextAlign(22);
                txtall->SetNDC();
                txtall->Draw();

                can.Update();
                can.SaveAs(Form("pics/iasy_ptt_dec_cent%i_dpsi%i_xi%i.pdf",icent,idpsi,ixi));

            }
        }
    }

    cout<<"Writing to "<<output->GetName()<<endl;
    output->Write();
    cout<<"Closing "<<output->GetName()<<endl;
    output->Close();

}

TFile* OpenFile(const char* name)
{
    cout<<"Reading file "<<name<<endl;

    TFile* file = new TFile(name,"READ");

    if (!file)
    {
        cerr<<"Could not open "<<name<<endl;
    }
    if (file->IsZombie())
    {
        cerr<<"File "<<name<<" is zombie!"<<endl;
    }

    return file;
}
