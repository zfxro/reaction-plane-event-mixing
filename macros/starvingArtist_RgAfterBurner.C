namespace GAT {
#include "gamAnaTree.h"
}

const int ncent = GAT::ncent;
const int ntype = GAT::ntype;
const int ndpsi = GAT::ndpsi;
const int nptt = GAT::nptt;
const int nxi = GAT::nxi;

TFile* output;

double rgavg[ncent+1][ndpsi+1][nptt];

TH1F* h1dphi[ntype+1][nptt][ncent+1][ndpsi+1][nxi];

int ReadRgAvg(const char* loc);

void GetHistos(const char* loc, int itype, bool isadj);

void starvingArtist_RgAfterBurner(
    const char* name_artist_raw = "../raw/artist_out.root",
    const char* name_artist_adj = "../adj/artist_out.root",
    const char* name_rgavg = "RgAvg.txt",
    const char* name_output = "artist_out.root"
    )
{
    cout<<"Running starvingArtist_RgAfterBurner.C"<<endl;

    output = new TFile(name_output,"RECREATE");

    if (ReadRgAvg(name_rgavg) < 0) return;

    GetHistos(name_artist_raw,GAT::INC,0);
    GetHistos(name_artist_adj,GAT::DEC,1);

    for (int iptt = 0; iptt < nptt; iptt++)
    {
        for (int icent = 0; icent < ncent+1; icent++)
        {
            for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
            {
                for (int ixi = 0; ixi < nxi; ixi++)
                {
                    h1dphi[GAT::INC][iptt][icent][idpsi][ixi]->Print();
                    h1dphi[GAT::DEC][iptt][icent][idpsi][ixi]->Print();
                }
            }
        }
    }

    for (int iptt = 0; iptt < nptt; iptt++)
    {
        for (int icent = 0; icent < ncent+1; icent++)
        {
            for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
            {
                for (int ixi = 0; ixi < nxi; ixi++)
                {
                    TH1F* hinc = h1dphi[GAT::INC][iptt][icent][idpsi][ixi];
                    TH1F* hdec = h1dphi[GAT::DEC][iptt][icent][idpsi][ixi];
                    TH1F* hdir = h1dphi[GAT::DIR][iptt][icent][idpsi][ixi];

                    ostringstream ossh;
                    ossh<<"dphi_adj_dir_mix2_ptt"<<iptt<<"_cent"<<icent<<"_dpsi"<<idpsi<<"_xi"<<ixi;
                    hdir = (TH1F*) hinc->Clone(ossh.str().c_str());
                    hdir->SetDirectory(output);

                    hdir->Add(hdec,-1./rgavg[icent][idpsi][iptt]);
                }
            }
        }
    }

    cout<<"Writing to "<<output->GetName()<<endl;
    output->Write();
    output->Close();
    cout<<"Done!"<<endl;
}

void GetHistos(const char* loc, int itype, bool isadj)
{
    TFile* input = new TFile(loc,"READ");

    for (int iptt = 0; iptt < nptt; iptt++)
    {
        for (int icent = 0; icent < ncent+1; icent++)
        {
            for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
            {
                for (int ixi = 0; ixi < nxi; ixi++)
                {
                    ostringstream ossh;
                    ossh<<"dphi_";
                    if (isadj)
                        ossh<<"adj_";
                    ossh<<GAT::spec[itype]<<"_mix2_ptt"<<iptt<<"_cent"<<icent<<"_dpsi"<<idpsi<<"_xi"<<ixi;

                    h1dphi[itype][iptt][icent][idpsi][ixi] = NULL;
                    h1dphi[itype][iptt][icent][idpsi][ixi] = (TH1F*) gROOT->FindObject(ossh.str().c_str());
                    assert(h1dphi[itype][iptt][icent][idpsi][ixi]);

                    ostringstream osshnew;
                    osshnew<<"dphi_adj_"<<GAT::spec[itype]<<"_mix2_ptt"<<iptt<<"_cent"<<icent<<"_dpsi"<<idpsi<<"_xi"<<ixi;
                    h1dphi[itype][iptt][icent][idpsi][ixi]->SetName(osshnew.str().c_str());
                    
                    h1dphi[itype][iptt][icent][idpsi][ixi]->SetDirectory(output);
                }
            }
        }
    }

    input->Close();
    output->cd();
}

int ReadRgAvg(const char* loc)
{
    cout<<"Reading average Rg values."<<endl;
    
    ifstream ifs(loc, ifstream::in);
    if (!ifs.is_open())
    {
        cerr << "Error opening " << loc << endl;
        return -1;
    }

    string line;
    while (getline(ifs,line))
    {
        if (!isdigit(line.at(0))) continue;
        istringstream issline(line.data());

        int icent = -999;
        int idpsi = -999;
        int iptt = -999;
        double avg = 0.0;
        
        issline >> icent >> idpsi >> iptt >> avg;

        rgavg[icent][idpsi][iptt] = avg;
    }


    cout<<"Printing average Rg values. Should match "<<loc<<endl;
    for (int icent = 0; icent < ncent+1; icent++)
    {
        for (int idpsi = 0; idpsi < ndpsi+1; idpsi++)
        {
            for (int iptt = 0; iptt < nptt; iptt++)
            {
                cout<<icent<<" "<<idpsi<<" "<<iptt<<" "<<rgavg[icent][idpsi][iptt]<<endl;
            }
        }
    }

    return 0;
    
}
