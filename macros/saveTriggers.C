const float PI = TMath::Pi();

const int nadj = 2;
const int nmix = 2;
const int ntype = 4;

const double detah = 0.7;


namespace GAT {
     #include "gamAnaTree.h"
}

TH2F* h2centptt[nadj][ntype][GAT::ndpsi+1];

TGraphErrors* tgrng[nadj][ntype][GAT::nptt][GAT::ndpsi+1];
TGraphErrors* tgreff[nadj][ntype][GAT::nptt][GAT::ndpsi+1][GAT::nxi];    

void saveTriggers(
    const char* name_gamtree = "gamtree.root",
    const char* name_out = "gamtree_triggers.root"
    )
{

    cout<<"Counting Triggers"<<endl;

    TFile* output = new TFile(name_out,"RECREATE");
    cout<<"To be written to "<<name_out<<endl;

    //==================================================

    cout<<"Reading h2centptt's from "<<name_gamtree<<endl;
    TFile* gamtree = new TFile(name_gamtree,"READ");
    if (!gamtree || gamtree->IsZombie())
    {
        cerr<<"Could not open "<<name_gamtree<<endl;
        return;
    }

    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
            {
                ostringstream oss;
                oss<<"cent_ptt";
                if (iadj)
                    oss<<"_adj";
                oss<<"_"<<GAT::spec[itype];
                oss<<"_dpsi"<<idpsi;

                if (idpsi == GAT::ndpsi)
                {
                    h2centptt[iadj][itype][idpsi] = (TH2F*) h2centptt[iadj][itype][0]->Clone(oss.str().c_str());

                    h2centptt[iadj][itype][idpsi]->Add(h2centptt[iadj][itype][1]);
                    h2centptt[iadj][itype][idpsi]->SetDirectory(output);
                    continue;
                }       

                h2centptt[iadj][itype][idpsi] = (TH2F*) gROOT->FindObject(oss.str().c_str());
                if (!h2centptt[iadj][itype][idpsi])
                {
                    cerr<<"Could not find "<<oss.str()<<endl;
                    return -1;
                }
                h2centptt[iadj][itype][idpsi]->SetDirectory(output);
            }
        }
    }

    //==================================================

    cout<<"Taking centrality projections"<<endl;
    
    for (int iadj = 0; iadj < nadj; iadj++)
    {
        for (int itype = 0; itype < ntype; itype++)
        {
            if (itype == GAT::PI0 && iadj == 1) continue;
            if (itype == GAT::INC && iadj == 1) continue;

            for (int idpsi = 0; idpsi <= GAT::ndpsi; idpsi++)
            {
                TH1F* h1ptt_tmp = (TH1F*) h2centptt[iadj][itype][idpsi]->ProjectionY();

                for (int iptt = 0; iptt < GAT::nptt; iptt++)
                {
                    int pttbinlo = h1ptt_tmp->FindBin(GAT::limptt[iptt]);
                    int pttbinhi = h1ptt_tmp->FindBin(GAT::limptt[iptt+1])-1;
                    if (pttbinhi > h1ptt_tmp->GetNbinsX())
                    {
                        cerr<<"THIS SHOULD NEVER HAPPEN"<<endl;
                        pttbinhi = h1ptt_tmp->GetNbinsX();
                    }

                    TH1F* h1cent = (TH1F*) h2centptt[iadj][itype][idpsi]->ProjectionX(Form("cent_adj%i_%s_ptt%i_dpsi%i",iadj,GAT::spec[itype],iptt,idpsi),pttbinlo,pttbinhi);
                    h1cent->SetDirectory(output);
                }
            }
        }
    }

    cout<<"Writing to "<<name_out<<endl;
    output->Write();
    output->Close();
}
