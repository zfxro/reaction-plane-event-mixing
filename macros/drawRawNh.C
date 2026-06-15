//draw acceptance correction (result of acceptCorrect.C)

const float PI = TMath::Pi();
const int nadj = 2;
const int nmix = 2;
const int ntype = 3;

const int markstyle[] = {20,24,21,25,22,26,23,32,34,47};

int GetRawNhHistos(const char* name);

namespace GAT {
    #include "gamAnaTree.h"
}

TH1F* h1had[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi+1][GAT::nxi] = {};
TGraphErrors* tgrhad[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi+1] = {};

void drawRawNh(
    const char* inname = "gamtree_13466_rawnh.root",
    const char* outname = "tgraphs_rawnh.root"
    )
{

    GetRawNhHistos(inname);

    TFile* output = new TFile(outname,"RECREATE");

    TCanvas* can = new TCanvas("can","",1000,1500);
    const int npad = 4;
    TPad* pads[npad] = {};

    pads[0] = new TPad("pad0","",0.,.5,.5,1.);
    pads[0]->SetTopMargin(.12);
    pads[0]->SetBottomMargin(0.);
    pads[0]->SetLeftMargin(.18);
    pads[0]->SetRightMargin(0.);
    pads[0]->SetTicks(0,1);
    pads[0]->Draw();

    pads[1] = new TPad("pad1","",.5,.5,1.,1.);
    pads[1]->SetTopMargin(.12);
    pads[1]->SetBottomMargin(0.);
    pads[1]->SetLeftMargin(0.);
    pads[1]->SetRightMargin(.18);
    pads[1]->Draw();

    pads[2] = new TPad("pad2","",0.,0.,.5,.5);
    pads[2]->SetTopMargin(0.);
    pads[2]->SetBottomMargin(.12);
    pads[2]->SetLeftMargin(.18);
    pads[2]->SetRightMargin(0.);
    pads[2]->SetTicks(1,1);
    pads[2]->Draw();

    pads[3] = new TPad("pad3","",.5,0.,1.,.5);
    pads[3]->SetTopMargin(0.);
    pads[3]->SetBottomMargin(.12);
    pads[3]->SetLeftMargin(0.);
    pads[3]->SetRightMargin(.18);
    pads[3]->SetTicks(1,0);
    pads[3]->Draw();

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;

            for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
            {
                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    pads[iptt]->cd();

                    TLegend* leg = nullptr;
                    if (iptt == 0)
                    {
                        leg = new TLegend(0.69,0.03,.96,.36);
                    }
                    if (iptt == 1)
                    {
                        leg = new TLegend(0.53,0.03,.80,.25);
                    }

                    TMultiGraph* multi = new TMultiGraph();

                    for (int icent = 0; icent < GAT::ncent; icent++)
                    {
                        TGraphErrors* tgr = tgrhad[iadj][itype][iptt][icent][idpsi];

                        gStyle->SetPalette(kBird);
                        tgr->SetMarkerColor(TColor::GetColorPalette(254*icent/(GAT::ncent-1)));
                        tgr->SetLineColor(TColor::GetColorPalette(254*icent/(GAT::ncent-1)));
                        tgr->SetMarkerSize(3.0);
                        tgr->SetLineWidth(3);
                        tgr->SetMarkerStyle(markstyle[icent]);

                        tgr->Write();

                        multi->Add(tgr,"LP");

                        if ((iptt == 0 && icent < 6) || (iptt == 1 && icent >= 6))
                        {
                            ostringstream oss;
                            oss<<GAT::limcent[icent]<<"-"<<GAT::limcent[icent+1]<<"%"<<endl;
                            leg->AddEntry(tgr,oss.str().c_str(),"lp");
                        }
                    }

                    multi->SetTitle(";#xi;dN_{h}/dN_{evt}");
                    multi->Draw("A");
                    multi->GetHistogram()->GetYaxis()->SetRangeUser(0.003,15);
                    multi->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
                    multi->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
                    multi->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
                    multi->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
                    if (leg)
                    {
                        leg->SetMargin(.4);
                        leg->SetTextSize(.05);
                        leg->Draw();
                    }

                    pads[iptt]->SetLogy();
                    pads[iptt]->Modified();
                    pads[iptt]->Update();

                    TLatex latex;
                    latex.SetTextSize(0.055);
                    latex.DrawLatex(.13,6,Form("#pi^{0} p_{T} #in [%.1f,%.1f)",GAT::limptt[iptt],GAT::limptt[iptt+1]));


                }

                can->SaveAs(Form("rawnhpics/rawnh_adj%i_%s_dpsi%i.png",iadj,GAT::spec[itype],idpsi));
            }
        }
    }

    cout<<"TGraphs have been written to "<<outname<<endl;
    output->Close();

    cout<<"Finished drawing acceptance*efficiency correction."<<endl;
    return;
}

int GetRawNhHistos(const char* name)
{
    TFile* input = new TFile(name,"READ");

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;

            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                for (int icent = 0; icent < GAT::ncent; icent++)
                {
                    for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                    {
                        vector<double> vxi;
                        vector<double> vhad;
                        vector<double> verr;

                        for (int ixi = 0; ixi < GAT::nxi; ixi++)
                        {
                            h1had[iadj][itype][iptt][icent][idpsi][ixi] = (TH1F*) gROOT->FindObject(Form("acc_adj%i_%s_ptt%i_cent%i_dpsi%i_xi%i",iadj,GAT::spec[itype],iptt,icent,idpsi,ixi));

                            TH1F* h1ptr = h1had[iadj][itype][iptt][icent][idpsi][ixi];

                            double nbin = h1ptr->GetNbinsX();

                            double error = -9999.;
                            double integral = h1ptr->IntegralAndError(1,nbin,error);

                            double avghad = integral;
                            error = error;

                            vhad.push_back(avghad);
                            verr.push_back(error);

                            vxi.push_back(ixi);
                        }

                        tgrhad[iadj][itype][iptt][icent][idpsi] = new TGraphErrors(GAT::nxi,&vxi[0],&vhad[0],0,&verr[0]);
                        tgrhad[iadj][itype][iptt][icent][idpsi]->SetName(Form("tgrhad_adj%i_%s_ptt%i_cent%i_dpsi%i",iadj,GAT::spec[itype],iptt,icent,idpsi));
                    }
                }
            }
        }
    }

    input->Close();
    delete input;

    return 0;
}
