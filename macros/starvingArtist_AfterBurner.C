namespace GAT {
#include "gamAnaTree.h"
}

const int ncent = GAT::ncent;
const int ntype = 2;
const int ndpsi = GAT::ndpsi;
const int nptt = GAT::nptt;
const int nxi = GAT::nxi;

TH1F* h1dphi[ntype][nptt][ncent+1][ndpsi+1][nxi];
TH1F* h1dphi_true[ntype][nptt][ncent+1][ndpsi+1][nxi];

double ntrg[ntype][nptt][ncent+1][ndpsi+1] = {};

const char* spec[ntype] = {"dir","dec"};

void starvingArtist_AfterBurner(
    const char* name_gamtree = "gamtree.root",
    const char* name_artist = "artist_out.root",
    const char* name_output = "artist_out_wavg.root"
    )
{
    cout<<"Running starvingArtist_AfterBurner.C"<<endl;

    TFile* gamtree = new TFile(name_gamtree,"READ");
    for (int itype = 0; itype < ntype; itype++)
    {
        for (int iptt = 0; iptt < nptt; iptt++)
        {
            for (int icent = 0; icent < ncent; icent++)
            {
                for (int idpsi = 0; idpsi < ndpsi; idpsi++)
                {
                    TH1F* h1tmp = (TH1F*) gROOT->FindObject(Form("ntrg_%s_pttb%i_centb%i_dpsib%i",spec[itype],iptt,icent,idpsi));

                    float tmpbc = h1tmp->GetBinContent(1);
                    ntrg[itype][iptt][icent][idpsi] = tmpbc;
                    ntrg[itype][iptt][ncent][idpsi] += tmpbc;
                    ntrg[itype][iptt][ncent][ndpsi] += tmpbc;
                    ntrg[itype][iptt][icent][ndpsi] += tmpbc;
                }
            }
        }
    }

    for (int itype = 0; itype < ntype; itype++)
    {
        for (int iptt = 0; iptt < nptt; iptt++)
        {
            for (int icent = 0; icent < ncent+1; icent++)
            {
                for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
                {
                    cout<<"ntrg["<<itype<<"]["<<iptt<<"]["<<icent<<"]["<<idpsi<<"] = "<<ntrg[itype][iptt][icent][idpsi]<<endl;
                }
            }
        }
    }

    gamtree->Close();


    TFile* output = new TFile(name_output,"RECREATE");

    TFile* artist = new TFile(name_artist,"READ");

    for (int iptt = 0; iptt < nptt; iptt++)
    {
        for (int icent = 0; icent < ncent+1; icent++)
        {
            for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
            {
                for (int ixi = 0; ixi < nxi; ixi++)
                {                   
                    h1dphi[GAT::DIR][iptt][icent][idpsi][ixi] = (TH1F*) gROOT->FindObject(Form("dphi_dir_mix2_ptt%i_cent%i_dpsi%i_xi%i",iptt,icent,idpsi,ixi));
                    h1dphi[GAT::DIR][iptt][icent][idpsi][ixi]->SetDirectory(0);

                    h1dphi[GAT::DEC][iptt][icent][idpsi][ixi] = (TH1F*) gROOT->FindObject(Form("dphi_dec_mix2_ptt%i_cent%i_dpsi%i_xi%i",iptt,icent,idpsi,ixi));
                    h1dphi[GAT::DEC][iptt][icent][idpsi][ixi]->SetDirectory(0);
                }
            }
        }
    }

    artist->Close();
    output->cd();

    const int nbin = h1dphi[0][0][0][0][0]->GetNbinsX();

    // Read mistag and tageff weighted averages (ksic applied)
    //========== MISTAG ==========
    double mistag[nptt][ncent] = {};
    const char* name_mistag = "mistag_wavg.txt";
    ifstream ifs_mistag(name_mistag,ifstream::in);
    if (!ifs_mistag)
    {
        cerr<<"Could not open "<<name_mistag<<endl;
        cout<<"Exiting..."<<endl;
        return;
    }
    string line;
    while (getline(ifs_mistag,line))
    {
        if (!isdigit(line.at(0))) continue;

        istringstream iss_line(line.data());
        int iptt, icent;
        double val;
        iss_line >> iptt >> icent >> val;

        mistag[iptt][icent] = val;

        cout<<iptt<<"  "<<icent<<"  "<<mistag[iptt][icent]<<endl;
    }

    //========== TAGEFF ==========
    double tageff[nptt][ncent] = {};
    const char* name_tageff = "tageff_wavg.txt";
    ifstream ifs_tageff(name_tageff,ifstream::in);
    if (!ifs_tageff)
    {
        cerr<<"Could not open "<<name_tageff<<endl;
        cout<<"Exiting..."<<endl;
        return;
    }
    line.clear();
    while (getline(ifs_tageff,line))
    {
        if (!isdigit(line.at(0))) continue;

        istringstream iss_line(line.data());
        int iptt, icent;
        double val;
        iss_line >> iptt >> icent >> val;

        tageff[iptt][icent] = val;

        cout<<iptt<<"  "<<icent<<"  "<<tageff[iptt][icent]<<endl;
    }

    //==================================================

    for (int iptt = 0; iptt < nptt; iptt++)
    {
        for (int icent = 0; icent < ncent+1; icent++)
        {
            cout<<"\n--------------------------------------------------"<<endl;                                                           
            cout<<"iptt = "<<iptt<<"  icent = "<<icent<<endl;
            cout<<"--------------------------------------------------"<<endl;
                                                                          
            double t = tageff[iptt][icent];
            double m = mistag[iptt][icent];

            cout<<"tageff = "<<t<<endl;
            cout<<"mistag = "<<m<<endl;
                                 
            double dirtdirm = t/(t-m);
            double dirtdecm = -(1.-t)/(t-m);

            cout<<"ndirt = "<<dirtdirm<<" ndirm  + "<<dirtdecm<<" ndecm"<<endl;

            double dectdecm = (1.-m)/(t-m);
            double dectdirm = -m/(t-m);

            cout<<"ndect = "<<dectdirm<<" ndirm  + "<<dectdecm<<" ndecm"<<endl;
                
            for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
            {
                double ndirm = ntrg[GAT::DIR][iptt][icent][idpsi];
                double ndecm = ntrg[GAT::DEC][iptt][icent][idpsi];

                double ndirt = dirtdirm*ndirm + dirtdecm*ndecm;
                double ndect = dectdirm*ndirm + dectdecm*ndecm;
                
                for (int ixi = 0; ixi < nxi; ixi++)
                {
                    TH1F* h1dirmeas = h1dphi[GAT::DIR][iptt][icent][idpsi][ixi];
                    TH1F* h1decmeas = h1dphi[GAT::DEC][iptt][icent][idpsi][ixi];
                    
                    h1dphi_true[GAT::DIR][iptt][icent][idpsi][ixi] = new TH1F(Form("dphi_dir_true_ptt%i_cent%i_dpsi%i_xi%i",iptt,icent,idpsi,ixi),Form("dphi_dir_true_ptt%i_cent%i_dpsi%i_xi%i",iptt,icent,idpsi,ixi),nbin,0.,1.);

                    TH1F* h1dirtrue = h1dphi_true[GAT::DIR][iptt][icent][idpsi][ixi];

                    h1dirtrue->Add(h1dirmeas,(dirtdirm*ndirm)/ndirt);
                    h1dirtrue->Add(h1decmeas,(dirtdecm*ndecm)/ndirt);

                    
                    h1dphi_true[GAT::DEC][iptt][icent][idpsi][ixi] = new TH1F(Form("dphi_dec_true_ptt%i_cent%i_dpsi%i_xi%i",iptt,icent,idpsi,ixi),Form("dphi_dec_true_ptt%i_cent%i_dpsi%i_xi%i",iptt,icent,idpsi,ixi),nbin,0.,1.);                  

                    TH1F* h1dectrue = h1dphi_true[GAT::DEC][iptt][icent][idpsi][ixi];

                    h1dectrue->Add(h1dirmeas,(dectdirm*ndirm)/ndect);
                    h1dectrue->Add(h1decmeas,(dectdecm*ndecm)/ndect);
                }
            }
        }
    }

    output->Write();
    output->Close();

    cout<<"Wrote to "<<output->GetName()<<endl;

    return;
                 
}
