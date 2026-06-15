//draw the act of acceptance
//one column three plots
//     raw result (no acceptance)
//     the acceptance
//     the ratio (the corrected result)

namespace GAT {
#include "../src/gamAnaTree/gamAnaTree.h"
}

const double PI = TMath::Pi();

const int nadj = 2;
const int ntype = 3;
const int nmix = 2;

TH1F* h1acc[nadj][ntype][GAT::nptt][GAT::ncent][GAT::ndpsi+1][GAT::nxi];
TH1F* h1dphi_acc[nadj][ntype][nmix][GAT::nptt][GAT::ncent][GAT::ndpsi+1][GAT::nxi];
TH1F* h1dphi_raw[nadj][ntype][nmix][GAT::nptt][GAT::ncent][GAT::ndpsi+1][GAT::nxi];

TFile* file_acc;
TFile* file_raw;

void GetHistosAcc(const char* inname);
void GetHistosRaw(const char* inname);
double GetNTrg(TH1D* h1cent, int icent);
void ScaleDphiHistos();
void DrawHistos(TCanvas* can, TH1F* dphiraw, TH1F* acc, TH1F* dphiacc);

void drawActOfAcceptance(
    const char* inname_acc = "gamtree_13466_acceptcorrect.root",
    const char* inname_raw = "gamtree_13466_psimix.root"
    )
{
    GetHistosAcc(inname_acc);   
    GetHistosRaw(inname_raw);

    ScaleDphiHistos();

    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);

    const int iadj = 0;
    const int itype = GAT::PI0;
    const int imix = 0;

    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
        {
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {       
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    TCanvas* can = new TCanvas("can","",2500,4500);
                    DrawHistos(can,h1dphi_raw[iadj][itype][imix][iptt][icent][idpsi][ixi],h1acc[iadj][itype][iptt][icent][idpsi][ixi],h1dphi_acc[iadj][itype][imix][iptt][icent][idpsi][ixi]);
                    can->SaveAs(Form("accact/accact_adj%i_%s_ptt%i_cent%i_dpsi%i_xi%i.pdf",iadj,GAT::spec[itype],iptt,icent,idpsi,ixi));
                    delete can;
                }
            }
        }
    }

    cout<<"drew the act of acceptance"<<endl;
}

void DrawHistos(TCanvas* can, TH1F* dphiraw, TH1F* acc, TH1F* dphiacc)
{
    
    TPad* top = new TPad("top","",0.,.65,1.,1.);
    top->SetTopMargin(1./7.);
    top->SetBottomMargin(0.);
    top->SetLeftMargin(.15);
    top->Draw();
    
    TPad* mid = new TPad("mid","",0.,.35,1.,.65);
    mid->SetTopMargin(0.);
    mid->SetBottomMargin(0.);
    mid->SetLeftMargin(.15);
    mid->Draw();
    
    TPad* bot = new TPad("mid","",0.,0.,1.,.35);
    bot->SetTopMargin(0.);
    bot->SetBottomMargin(1./7.);
    bot->SetLeftMargin(.15);
    bot->Draw();

    top->cd();
    dphiraw->SetTitle(";#Delta#phi (1/#pi);1/N_{#gamma}  dN_{h}/d#Delta#phi");
    dphiraw->SetMarkerColor(TColor::GetColorPalette(254*0/3));
    dphiraw->SetMarkerColor(TColor::GetColorPalette(254*0/3));
    dphiraw->SetMarkerStyle(20);
    dphiraw->SetMarkerSize(5);
    dphiraw->GetYaxis()->SetTitleSize(.045);
    dphiraw->GetYaxis()->SetLabelSize(.04);
    dphiraw->Draw();
    TLatex* tex = new TLatex();
    tex->SetTextAlign(22);
    tex->SetTextSize(.05);
    tex->DrawLatexNDC(.75,.75,"raw");

    mid->cd();
    acc->SetTitle(";#Delta#phi (1/#pi);#epsilon");
    acc->SetMarkerColor(TColor::GetColorPalette(254*1/3));
    acc->SetLineColor(TColor::GetColorPalette(254*1/3));
    acc->SetMarkerStyle(20);
    acc->SetMarkerSize(5);
    acc->GetYaxis()->SetTitleSize(.045*.35/.3);
    acc->GetYaxis()->SetLabelSize(.04*.35/.3);
    acc->GetYaxis()->SetRangeUser(.0325,.235);
    acc->Draw();
    tex->DrawLatexNDC(.75,.85,"acceptance");

    bot->cd();
    dphiacc->SetTitle(";#Delta#phi (1/#pi);1/N_{#gamma}  dN_{h}/d#Delta#phi");
    dphiacc->SetMarkerColor(TColor::GetColorPalette(254*2/3));
    dphiacc->SetLineColor(TColor::GetColorPalette(254*2/3));
    dphiacc->SetMarkerStyle(20);
    dphiacc->SetMarkerSize(5);
    dphiacc->GetYaxis()->SetTitleSize(.045);
    dphiacc->GetYaxis()->SetLabelSize(.04);
    dphiacc->GetYaxis()->SetRangeUser(2.28,2.82);
    dphiacc->GetXaxis()->SetTitleSize(.045);
    dphiacc->GetXaxis()->SetLabelSize(.04);
    dphiacc->Draw();
    tex->DrawLatexNDC(.75,.85,"corrected");
}

void GetHistosAcc(const char* inname)
{
    file_acc = new TFile(inname,"READ");

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            for (int icent = 0; icent < GAT::ncent; icent++)
            {
                for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                {
                    for (int iptt = 0; iptt < GAT::nptt; iptt++)
                    {       
                        for (int ixi = 0; ixi < GAT::nxi; ixi++)
                        {
                            ostringstream osssuf;
                            osssuf<<"_ptt"<<iptt;
                            osssuf<<"_cent"<<icent;
                            osssuf<<"_dpsi"<<idpsi;
                            osssuf<<"_xi"<<ixi;

                            ostringstream ossacc;
                            ossacc<<"acc";
                            ossacc<<"_adj"<<iadj;
                            ossacc<<"_"<<GAT::spec[itype];
                            ossacc<<osssuf.str();

                            h1acc[iadj][itype][iptt][icent][idpsi][ixi] = (TH1F*) gROOT->FindObject(ossacc.str().c_str());

                            for (int imix = 0; imix < nmix; imix++)
                            {
                                ostringstream ossdphi;
                                ossdphi<<"dphi";
                                if (iadj)
                                    ossdphi<<"_adj";
                                ossdphi<<"_"<<GAT::spec[itype];
                                ossdphi<<"_mix"<<imix;
                                ossdphi<<osssuf.str();
                                
                                h1dphi_acc[iadj][itype][imix][iptt][icent][idpsi][ixi] = (TH1F*) gROOT->FindObject(ossdphi.str().c_str());
                            }
                            
                        }
                    }
                }
            }
        }
    }                       
    return;
}

void GetHistosRaw(const char* inname)
{
    file_raw = new TFile(inname,"READ");

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            for (int imix = 0; imix < nmix; imix++)
            {
                for (int icent = 0; icent < GAT::ncent; icent++)
                {
                    for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                    {
                        for (int iptt = 0; iptt < GAT::nptt; iptt++)
                        {           
                            for (int ixi = 0; ixi < GAT::nxi; ixi++)
                            {
                                ostringstream oss;
                                oss<<"dphi";
                                if (iadj)
                                    oss<<"_adj";
                                oss<<"_"<<GAT::spec[itype];
                                oss<<"_mix"<<imix;
                                oss<<"_ptt"<<iptt;
                                oss<<"_cent"<<icent;
                                oss<<"_dpsi"<<idpsi;
                                oss<<"_xi"<<ixi;

                                if (idpsi == GAT::ndpsi)
                                {
                                    h1dphi_raw[iadj][itype][imix][iptt][icent][idpsi][ixi] = (TH1F*) h1dphi_raw[iadj][itype][imix][iptt][icent][0][ixi]->Clone(oss.str().c_str());

                                    h1dphi_raw[iadj][itype][imix][iptt][icent][idpsi][ixi]->Add(h1dphi_raw[iadj][itype][imix][iptt][icent][1][ixi]);
                                    continue;
                                }                           

                                h1dphi_raw[iadj][itype][imix][iptt][icent][idpsi][ixi] = (TH1F*) gROOT->FindObject(oss.str().c_str());
                                if (!h1dphi_raw[iadj][itype][imix][iptt][icent][idpsi][ixi])
                                {
                                    cerr<<"Could not find "<<oss.str()<<endl;
                                    return -1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

void ScaleDphiHistos()
{

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < GAT::ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            for (int icent = 0; icent < GAT::ncent; icent++)
            {
                for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
                {
                    for (int iptt = 0; iptt < GAT::nptt; iptt++)
                    {    

                        TH1D* h1trgcent = nullptr;
                        file_acc->GetObject(Form("cent_adj%i_%s_ptt%i_dpsi%i",iadj,GAT::spec[itype],iptt,idpsi),h1trgcent);
                        double ntrg = GetNTrg(h1trgcent,icent);

                        for (int ixi = 0; ixi < GAT::nxi; ixi++)
                        {
                            for (int imix = 0; imix < nmix; imix++)
                            {
                                double scale = ntrg*h1dphi_acc[iadj][itype][imix][iptt][icent][idpsi][ixi]->GetBinWidth(1)*PI*2.;
                                
                                h1dphi_acc[iadj][itype][imix][iptt][icent][idpsi][ixi]->Scale(1./scale);
                                h1dphi_raw[iadj][itype][imix][iptt][icent][idpsi][ixi]->Scale(1./scale);
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

double GetNTrg(TH1D* h1cent, int icent)
{
    double ntrg = 0.;
    for (int jcent = GAT::limcent[icent+1]; jcent > GAT::limcent[icent]; jcent--)
    {
        int tmpbin = h1cent->FindBin(jcent);
        ntrg += h1cent->GetBinContent(tmpbin);
    }
    return ntrg;
}
