//draw bgsbhs the 2 ways in which it is found

enum RPO {IN,OUT,ALL};

const int nadj = 2;

const float PI = (float) TMath::Pi();

namespace GAT {
#include "../src/gamAnaTree/gamAnaTree.h"
}

void drawKComponents_bgsbhs(
    const string inname = "kComponents.root"
    )
{
    TFile* input = new TFile(inname.c_str(),"READ");
    if (!input || input->IsZombie())
    {
        cerr<<"Could not open "<<inname<<endl;
        return;
    }
    else
    {
        cout<<"Reading "<<inname<<endl;
    }

    TCanvas* can = new TCanvas("can","",2000,2000);
    TPad* ptop = new TPad("ptop","",0,.3,1,1);
    ptop->Draw();
    TPad* pbot = new TPad("pbot","",0,0,1,.3);
    pbot->Draw();

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

                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    ostringstream sufh;
                    sufh<<sufg.str()<<"_xi"<<ixi;

                    TGraphErrors* grbgsbhs = nullptr;
                    input->GetObject(Form("bgsbhs_Np_%s",sufh.str().c_str()),grbgsbhs);
                    if (grbgsbhs == nullptr)
                    {
                        printf("Could not find bgsbhs_Np_%s. Exiting.\n",sufh.str().c_str());
                        return;
                    }
                    else
                    {
                        printf("Found bgsbhs_Np_%s.\n",sufh.str().c_str());
                    }
                                    

                    TGraphErrors* grbgs = nullptr;
                    input->GetObject(Form("bgs_Np_%s",sufh.str().c_str()),grbgs);
                    if (grbgs == nullptr)
                    {
                        printf("Could not find bgs_Np_%s. Exiting.\n",sufh.str().c_str());
                        return;
                    }
                    else
                    {
                        printf("Found bgs_Np_%s.\n",sufh.str().c_str());
                    }
                    
                    TGraphErrors* grbhs = nullptr;
                    input->GetObject(Form("bhs_Np_%s",sufh.str().c_str()),grbhs);
                    if (grbhs == nullptr)
                    {
                        printf("Could not find bhs_Np_%s. Exiting.\n",sufh.str().c_str());
                        return;
                    }
                    else
                    {
                        printf("Found bhs_Np_%s.\n",sufh.str().c_str());
                    }

                    TGraphErrors* grprod = new TGraphErrors();
                    TGraphErrors* grrat = new TGraphErrors();
                    for (int i = 0; i < grbgs->GetN(); i++)
                    {
                        double np = grbgs->GetX()[i];

                        double prod = grbgs->GetY()[i] * grbhs->GetY()[i];

                        grprod->SetPoint(i,np,prod);

                        double e0 = grbgs->GetY()[i] * grbhs->GetEY()[i];
                        double e1 = grbgs->GetEY()[i] * grbhs->GetY()[i];
                        double err = sqrt(e0*e0 + e1*e1);

                        grprod->SetPointError(i,grbgs->GetEX()[i],err);

                        grrat->SetPoint(i,np,prod/grbgsbhs->GetY()[i]);

                        e0 = grprod->GetEY()[i]/grbgsbhs->GetY()[i];
                        e1 = grbgsbhs->GetEY()[i]*prod/grbgsbhs->GetY()[i]/grbgsbhs->GetY()[i];
                        err = sqrt(e0*e0 + e1*e1);

                        grrat->SetPointError(i,grbgs->GetEX()[i],err);
                    }

                    TMultiGraph* mgr = new TMultiGraph();

                    grprod->SetLineColor(kRed);
                    grprod->SetLineWidth(4);                
                    mgr->Add(grprod,"l");
                    
                    grbgsbhs->SetLineColor(kBlue);
                    grbgsbhs->SetLineWidth(4);
                    mgr->Add(grbgsbhs,"l");

                    ptop->cd();
                    mgr->Draw("a");

                    can->Update();

                    pbot->cd();

                    grrat->SetLineColor(kMagenta);
                    grrat->SetLineWidth(4);;
                    grrat->Draw("ALP");

                    TLine* tl = new TLine(0.,1.,400.,1.);
                    tl->SetLineStyle(7);
                    tl->Draw("same");

                    can->Update();
                    
                    can->SaveAs(Form("bgsbhs_study/%s.png",grbgsbhs->GetName()));

                }
            }
        }
    }
}
