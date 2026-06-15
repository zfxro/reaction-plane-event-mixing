//draw reaction plane resolution plots

const float PI = TMath::Pi();
const int nadj = 2;
const int nmix = 2;
const int ntype = 3;

const int markstyle[] = {20,24,21,25,22,26,23,32,34,47};

namespace GAT {
#include "gamAnaTree.h"
}

TGraph* tgrsig;
TH1D* hpsidiff[GAT::ncent];
TGraph* tgrpsidiff[GAT::ncent];

int GetHistos(const char* name);

void drawRPReso(const char* inname = "rpreso.root")
{
    GetHistos(inname);

    TCanvas can ("can","",1000,1000);
    gStyle->SetPalette(kBird);

    tgrsig->SetMarkerSize(2);
    tgrsig->SetMarkerStyle(kFullCircle);
    tgrsig->SetMarkerColor(kBlue);
    tgrsig->SetLineColor(kBlue);

    tgrsig->Draw("ACP");
    tgrsig->SetTitle("");

    TAxis* xaxis = tgrsig->GetHistogram()->GetXaxis();
    xaxis->SetTitle("Centrality (%)");
    xaxis->SetTitleOffset(1.1);
    xaxis->SetLimits(GAT::limcent[0],GAT::limcent[GAT::ncent]);
    xaxis->SetNdivisions(200+GAT::ncent);

    TAxis* yaxis = tgrsig->GetHistogram()->GetYaxis();
    yaxis->SetTitle("#sigma");
    yaxis->SetTitleOffset(1.5);

    can.Update();
    can.SaveAs("resopics/sigma.pdf");

    TMultiGraph* mgr_[2];
    TLegend* leg[2];
    for (int i = 0; i < 2; i++)
    {
        mgr_[i] = new TMultiGraph();
        // leg[i]->SetMargin(.3);

        leg[i] = new TLegend(0.7,0.65,.89,.89);
        leg[i]->SetBorderSize(0);
    }
    
    TMultiGraph* mgr = new TMultiGraph();
    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        TGraph* tgr = tgrpsidiff[icent];
        
        tgr->SetMarkerStyle(markstyle[icent]);

        tgr->SetMarkerColor(TColor::GetColorPalette(240*icent/(GAT::ncent-1)));
        tgr->SetLineColor(TColor::GetColorPalette(240*icent/(GAT::ncent-1)));
        
        tgr->SetMarkerSize(0);
        tgr->SetLineWidth(5);

        mgr->Add(tgr,"CP");

        ostringstream oss;
        oss<<GAT::limcent[icent]<<"-"<<GAT::limcent[icent+1]<<"%"<<endl;
        if (icent < 4)
        {
            mgr_[0]->Add(tgr,"CP");
            leg[0]->AddEntry(tgr,oss.str().c_str(),"lp");
        }
        else
        {
            mgr_[1]->Add(tgr,"CP");
            leg[1]->AddEntry(tgr,oss.str().c_str(),"lp");
        }
    }


    mgr->SetTitle(";#psi_{2,FVTXN} #minus #psi_{2,FVTXS} (1/#pi);n_{evt}");
    mgr->Draw("A");

    mgr->GetXaxis()->SetLimits(-1.,1.);
    mgr->GetYaxis()->SetMaxDigits(2);
    mgr->GetXaxis()->SetTitleOffset(1.1);
    can.SaveAs("resopics/psidiff.pdf");

    for (int i = 0; i < 2; i++)
    {
        mgr_[i]->SetTitle(";#psi_{2,FVTXN} #minus #psi_{2,FVTXS} (1/#pi);n_{evt}");
        mgr_[i]->Draw("A");

        leg[i]->Draw();

        mgr_[i]->GetXaxis()->SetLimits(-1.,1.);
        mgr_[i]->GetYaxis()->SetMaxDigits(2);
        mgr_[i]->GetXaxis()->SetTitleOffset(1.1);
        can.SaveAs(Form("resopics/psidiff_%i.pdf",i));
    }

    return;
}


int GetHistos(const char* name)
{
    TFile* input = new TFile(name,"READ");

    tgrsig = (TGraph*) gROOT->FindObject("centreso");
    cout<<"\nPrinting centreso"<<endl;
    tgrsig->Print();

    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        hpsidiff[icent] = (TH1D*) gROOT->FindObject(Form("psidiff_%i",icent));
        hpsidiff[icent]->Sumw2();
        hpsidiff[icent]->Scale(1./hpsidiff[icent]->GetEntries());
        tgrpsidiff[icent] = new TGraph(hpsidiff[icent]);
        // tgrpsidiff[icent]->Print();
    }


    return 0;
}

