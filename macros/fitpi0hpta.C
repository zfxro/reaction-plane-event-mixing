//normalize pi0-h away side yield vs ptt curves
//fit the pta combined result

enum RPO {IN,OUT,ALL};

const float PI = TMath::Pi();
const int nadj = 2;
const int ntype = 3;

const int ncent = 8;
const int nsumcent = 2;

const int markstyle[] = {20,21,33,47,34};
const int markcol[] = {kBlue,kRed,kGreen,kMagenta,kCyan,kYellow};

namespace GAT {
#include "gamAnaTree.h"
}

TGraphErrors* tgripttyield[nadj][ntype][GAT::npta][nsumcent][GAT::ndpsi+1] = {};

void fitpi0hpta(
    const char* inname = "integawaysideyield_pta.root",
    const char* outname = "fitpi0hpta.root"
    )
{

    TFile* input = new TFile(inname,"READ");
    if (!input || input->IsZombie())
    {
        cerr<<"Error opening "<<inname<<endl;
        return;
    }

    //get objects
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            //doesnt exist
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype != GAT::PI0) continue;

            for (int ipta = 0; ipta < GAT::npta; ipta++)
            {
                for (int icent = 0; icent < nsumcent; icent++)
                {

                    for (int idpsi = 0; idpsi < GAT::ndpsi+1; idpsi++)
                    {

                        input->GetObject(Form("ipttyield_adj%i_%s_pta%i_cent%i_dpsi%i",iadj,GAT::spec[itype],ipta,icent,idpsi),tgripttyield[iadj][itype][ipta][icent][idpsi]);
                    }
                }
            }
        }
    }

    //normalize curves
    for (int icent = 0; icent < nsumcent; icent++)
    {
        for (int ipta = 0; ipta < GAT::npta; ipta++)
        {
            TGraphErrors* tgr = tgripttyield[0][GAT::PI0][ipta][icent][ALL];

            //get error weighted average
            double sum = 0.;
            double wt = 0.;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                double x = tgr->GetX()[iptt];
                double y = tgr->GetY()[iptt];
                double ey = tgr->GetEY()[iptt];

                sum += y/(ey*ey);
                wt += 1/(ey*ey);
            }

            double mean = sum/wt;
            double err = sqrt(1./wt);

            //normalize by weighted mean
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                tgr->GetY()[iptt] /= mean;
                tgr->GetEY()[iptt] /= mean;
            }
        }
    }

    for (int icent = 0; icent < nsumcent; icent++)
    {
        TCanvas can ("can","",1000,1000);
        can.SetLogy();

        TMultiGraph* mgr = new TMultiGraph();
        TLegend leg (0.3,0.775,0.9,0.9);
        leg.SetTextSize(.035);
        leg.SetTextAlign(22);
        leg.SetNColumns(2);
        for (int ipta = GAT::npta-1; ipta >= 0; ipta--)
        {
            TGraphErrors* tgr = tgripttyield[0][GAT::PI0][ipta][icent][ALL];

            tgr->SetMarkerColor(TColor::GetColorPalette(220*ipta/(GAT::npta-1)));
            tgr->SetLineColor(TColor::GetColorPalette(220*ipta/(GAT::npta-1)));

            tgr->SetMarkerSize(2.);
            tgr->SetLineWidth(3);
            tgr->SetMarkerStyle(markstyle[ipta]);

            mgr->Add(tgr,"LP");

            mgr->SetTitle(tgr->GetTitle());
        }

        for (int ipta = 0; ipta < GAT::npta; ipta++)
        {
            TGraphErrors* tgr = tgripttyield[0][GAT::PI0][ipta][icent][ALL];
            ostringstream oss;
            leg.AddEntry(tgr,Form("%.1f #leq p_{Ta} < %.1f",GAT::limpta[ipta],GAT::limpta[ipta+1]),"lp");
        }

        mgr->GetHistogram()->GetYaxis()->SetRangeUser(1.e-1,1.e2);
        mgr->Draw("A");

        mgr->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
        mgr->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
        mgr->GetHistogram()->GetXaxis()->SetLimits(5.,15.);
        mgr->GetHistogram()->GetXaxis()->SetNdivisions(-5);

        leg.Draw();

        double xmin = mgr->GetHistogram()->GetXaxis()->GetXmin();
        double xmax = mgr->GetHistogram()->GetXaxis()->GetXmax();
        TLine* line = new TLine(xmin,0.0,xmax,0.0);
        line->SetLineStyle(7);
        line->Draw();

        mgr->Fit("expo");

        can.Modified();
        can.SaveAs(Form("fitpi0hpta_cent%i.png",icent));
    }
}
