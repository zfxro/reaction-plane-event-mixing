//produce inclusive photon histograms to be merged with main gamAnaTree output
//
//histos to save:
//cent_ptt_inc_dpsi?
//dphi_pi0_mix?_ptt?_cent?_dpsi?_xi?
//dphi_pi0_mix?_ptt?_cent?_dpsi?_pta?
//dphiacc_pi0_ptt?_cent?_dpsi?_xi?

const float PI = TMath::Pi();

namespace GAT {
     #include "gamAnaTree.h"
}

void MakeInc(const char* inname, const char* outname);
TFile* OpenFile(const char* name);

void makeIncGamTree(
    const char* inname_rp = "gamtree_13466_psimix.root",
    const char* outname_rp = "gamtree_13466_inc_psimix.root",
    const char* inname_norp = "gamtree_13466_nopsimix.root",
    const char* outname_norp = "gamtree_13466_inc_nopsimix.root"
    )
{
    MakeInc(inname_rp,outname_rp);
    MakeInc(inname_norp,outname_norp);
    cout<<"Done."<<endl;
}

void MakeInc(const char* inname, const char* outname)
{
    TFile* output = new TFile(outname,"RECREATE");

    TFile* input = OpenFile(inname);

    for (int idpsi = 0; idpsi < GAT::ndpsi; idpsi++)
    {
        TH1F* cent_pi0 = (TH1F*) gROOT->FindObject(Form("cent_ptt_pi0_dpsi%d",idpsi));
        cent_pi0->SetDirectory(output);

        TH1F* cent_adj_dir = (TH1F*) gROOT->FindObject(Form("cent_ptt_adj_dir_dpsi%d",idpsi));
        cent_adj_dir->SetDirectory(output);
        TH1F* cent_adj_dec = (TH1F*) gROOT->FindObject(Form("cent_ptt_adj_dec_dpsi%d",idpsi));
        cent_adj_dec->SetDirectory(output);

        TH1F* cent_dir = (TH1F*) gROOT->FindObject(Form("cent_ptt_dir_dpsi%d",idpsi));
        cent_dir->SetDirectory(output);
        TH1F* cent_dec = (TH1F*) gROOT->FindObject(Form("cent_ptt_dec_dpsi%d",idpsi));
        cent_dec->SetDirectory(output);

        TH1F* cent_inc = (TH1F*) cent_dir->Clone(Form("cent_ptt_inc_dpsi%d",idpsi));
        cent_inc->Add(cent_dec);
        cent_inc->SetDirectory(output);

        for (int icent = 0; icent < GAT::ncent; icent++)
        {
            for (int iptt = 0; iptt < GAT::nptt; iptt++)
            {
                for (int ixi = 0; ixi < GAT::nxi; ixi++)
                {
                    TH1F* dphiacc_pi0 = (TH1F*) gROOT->FindObject(Form("dphiacc_pi0_ptt%d_cent%d_dpsi%d_xi%d",iptt,icent,idpsi,ixi));
                    dphiacc_pi0->SetDirectory(output);

                    TH1F* dphiacc_adj_dir = (TH1F*) gROOT->FindObject(Form("dphiacc_adj_dir_ptt%d_cent%d_dpsi%d_xi%d",iptt,icent,idpsi,ixi));
                    dphiacc_adj_dir->SetDirectory(output);
                    TH1F* dphiacc_adj_dec = (TH1F*) gROOT->FindObject(Form("dphiacc_adj_dec_ptt%d_cent%d_dpsi%d_xi%d",iptt,icent,idpsi,ixi));
                    dphiacc_adj_dec->SetDirectory(output);

                    TH1F* dphiacc_dir = (TH1F*) gROOT->FindObject(Form("dphiacc_dir_ptt%d_cent%d_dpsi%d_xi%d",iptt,icent,idpsi,ixi));
                    dphiacc_dir->SetDirectory(output);
                    TH1F* dphiacc_dec = (TH1F*) gROOT->FindObject(Form("dphiacc_dec_ptt%d_cent%d_dpsi%d_xi%d",iptt,icent,idpsi,ixi));
                    dphiacc_dec->SetDirectory(output);

                    TH1F* dphiacc_inc = (TH1F*) dphiacc_dir->Clone(Form("dphiacc_inc_ptt%d_cent%d_dpsi%d_xi%d",iptt,icent,idpsi,ixi));
                    dphiacc_inc->Add(dphiacc_dec);
                    dphiacc_inc->SetDirectory(output);
                }
                for (int imix = 0; imix < 2; imix++)
                {
                    for (int ixi = 0; ixi < GAT::nxi; ixi++)
                    {
                        TH1F* dphixi_pi0 = (TH1F*) gROOT->FindObject(Form("dphi_pi0_mix%d_ptt%d_cent%d_dpsi%d_xi%d",imix,iptt,icent,idpsi,ixi));
                        dphixi_pi0->SetDirectory(output);

                        TH1F* dphixi_adj_dir = (TH1F*) gROOT->FindObject(Form("dphi_adj_dir_mix%d_ptt%d_cent%d_dpsi%d_xi%d",imix,iptt,icent,idpsi,ixi));
                        dphixi_adj_dir->SetDirectory(output);
                        TH1F* dphixi_adj_dec = (TH1F*) gROOT->FindObject(Form("dphi_adj_dec_mix%d_ptt%d_cent%d_dpsi%d_xi%d",imix,iptt,icent,idpsi,ixi));
                        dphixi_adj_dec->SetDirectory(output);

                        TH1F* dphixi_dir = (TH1F*) gROOT->FindObject(Form("dphi_dir_mix%d_ptt%d_cent%d_dpsi%d_xi%d",imix,iptt,icent,idpsi,ixi));
                        dphixi_dir->SetDirectory(output);
                        TH1F* dphixi_dec = (TH1F*) gROOT->FindObject(Form("dphi_dec_mix%d_ptt%d_cent%d_dpsi%d_xi%d",imix,iptt,icent,idpsi,ixi));
                        dphixi_dec->SetDirectory(output);

                        TH1F* dphixi_inc = (TH1F*) dphixi_dir->Clone(Form("dphi_inc_mix%d_ptt%d_cent%d_dpsi%d_xi%d",imix,iptt,icent,idpsi,ixi));
                        dphixi_inc->Add(dphixi_dec);
                        dphixi_inc->SetDirectory(output);
                    }
                    for (int ipta = 0; ipta < GAT::npta; ipta++)
                    {
                        TH1F* dphipta_pi0 = (TH1F*) gROOT->FindObject(Form("dphi_pi0_mix%d_ptt%d_cent%d_dpsi%d_pta%d",imix,iptt,icent,idpsi,ipta));
                        dphipta_pi0->SetDirectory(output);

                        TH1F* dphipta_adj_dir = (TH1F*) gROOT->FindObject(Form("dphi_adj_dir_mix%d_ptt%d_cent%d_dpsi%d_pta%d",imix,iptt,icent,idpsi,ipta));
                        dphipta_adj_dir->SetDirectory(output);
                        TH1F* dphipta_adj_dec = (TH1F*) gROOT->FindObject(Form("dphi_adj_dec_mix%d_ptt%d_cent%d_dpsi%d_pta%d",imix,iptt,icent,idpsi,ipta));
                        dphipta_adj_dec->SetDirectory(output);

                        TH1F* dphipta_dir = (TH1F*) gROOT->FindObject(Form("dphi_dir_mix%d_ptt%d_cent%d_dpsi%d_pta%d",imix,iptt,icent,idpsi,ipta));
                        dphipta_dir->SetDirectory(output);
                        TH1F* dphipta_dec = (TH1F*) gROOT->FindObject(Form("dphi_dec_mix%d_ptt%d_cent%d_dpsi%d_pta%d",imix,iptt,icent,idpsi,ipta));
                        dphipta_dec->SetDirectory(output);

                        TH1F* dphipta_inc = (TH1F*) dphipta_dir->Clone(Form("dphi_inc_mix%d_ptt%d_cent%d_dpsi%d_pta%d",imix,iptt,icent,idpsi,ipta));
                        dphipta_inc->Add(dphipta_dec);
                        dphipta_inc->SetDirectory(output);
                    }

                }
            }
        }
    }

    input->Close();

    cout<<"Writing to "<<output->GetName()<<endl;
    output->Write();
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
