//see backup/CalcRPReso.C for more complicated macro with inc photon cent weighting

#include "gamConst.h"

namespace GAT {
    #include "gamAnaTree.h"
}

int GetNpartDists(TH1F* hnpart[GAT::ncent], vector<double>& vnpart, vector<double>& vnparte);

const float PI = TMath::Pi();

const int ncentmerge = 2;
const int limcentmerge[ncentmerge+1] = {0,20,60};

void CalcRPReso(
    const char * inname="gamAna_histos.root",
    const char * outname="rpreso.root",
    bool isBBC = false
    )
{
    //Get Npart distributions. (really only need mean and error)
    TH1F* hnpart[GAT::ncent];
    vector<double> vnpart;//mean npart
    vector<double> vnparte;//sigma npart
    if (GetNpartDists(hnpart,vnpart,vnparte)<0) return;

    // Set up centrality vector
    vector<double> vcent;
    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        double tmpcent = GAT::limcent[icent+1] + GAT::limcent[icent]+1;
        tmpcent/=2.;
        vcent.push_back(tmpcent);
        cout<<"Centrality midpoint of bin "<<icent<<" is "<<tmpcent<<endl;
    }

    // Set up centrality vector for merged bins
    vector<double> vcentmerge;
    for (int icent = 0; icent < ncentmerge; icent++)
    {
        double tmpcent = limcentmerge[icent+1] + limcentmerge[icent]+1;
        tmpcent/=2.;
        vcentmerge.push_back(tmpcent);
        cout<<"Centrality midpoint of merged bin "<<icent<<" is "<<tmpcent<<endl;
    }

    TFile* outfile = new TFile(outname,"RECREATE");

    //open file
    TFile* infile = new TFile(inname);

    TH2F* h2centpsidiff = nullptr;
    if (isBBC)
        h2centpsidiff = (TH2F*) gROOT->FindObject("centpsidiff_bbcns");
    else
        h2centpsidiff = (TH2F*) gROOT->FindObject("centpsidiff_fvtxns");

    TH1D* h1cent = h2centpsidiff->ProjectionX("cent",1,h2centpsidiff->GetNbinsY());
    h1cent->SetDirectory(outfile);

    vector<double> vreso;

    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        int cbin0 = h1cent->FindBin(GAT::limcent[icent]+1);//do not include lower limit
        int cbin1 = h1cent->FindBin(GAT::limcent[icent+1]);//include upper limit

        TH1D* h1psidiff = h2centpsidiff->ProjectionY(Form("psidiff_%i",icent),cbin0,cbin1);
        h1psidiff->SetDirectory(outfile);

        int nx = h1psidiff->GetNbinsX();

        float cosine_sum = 0.0;

        for (int ix = 1; ix <= nx; ix++)
        {
            float bincenter = h1psidiff->GetBinCenter(ix);

            float bincontent = h1psidiff->GetBinContent(ix);

            float cosine_tmp = cos(2.*PI*bincenter);

            cosine_sum += cosine_tmp*bincontent;
        }

        //must get <Cos(2. * Psi_BBCN - Psi_BBCS)> so divide by number of entries
        //note: Pi is left out of above definition
        //Psi had been divided by Pi when filling histogram from infile

        float nentry = h1psidiff->GetEntries();
        //see https://arxiv.org/abs/nucl-ex/9805001
        //specifically Eqs. (14) & (15)
        float reso = sqrt(2.*cosine_sum/nentry);

        vreso.push_back(reso);

    }

    //==================================================
    // Repeat for merged centrality bins
    //==================================================
    vector<double> vresomerge;

    for (int icent = 0; icent < ncentmerge; icent++)
    {
        int cbin0 = h1cent->FindBin(limcentmerge[icent]+1);//do not include lower limit
        int cbin1 = h1cent->FindBin(limcentmerge[icent+1]);//include upper limit

        TH1D* h1psidiff = h2centpsidiff->ProjectionY(Form("psidiff_merge%i",icent),cbin0,cbin1);
        h1psidiff->SetDirectory(outfile);

        int nx = h1psidiff->GetNbinsX();

        float cosine_sum = 0.0;

        for (int ix = 1; ix <= nx; ix++)
        {
            float bincenter = h1psidiff->GetBinCenter(ix);

            float bincontent = h1psidiff->GetBinContent(ix);

            float cosine_tmp = cos(2.*PI*bincenter);

            cosine_sum += cosine_tmp*bincontent;
        }

        //must get <Cos(2. * Psi_BBCN - Psi_BBCS)> so divide by number of entries
        //note: Pi is left out of above definition
        //Psi had been divided by Pi when filling histogram from infile

        float nentry = h1psidiff->GetEntries();
        //see https://arxiv.org/abs/nucl-ex/9805001
        //specifically Eqs. (14) & (15)
        float reso = sqrt(2.*cosine_sum/nentry);

        vresomerge.push_back(reso);

    }

    infile->Close();
    outfile->cd();

    TGraph* tgrcentreso = new TGraph(vcent.size(),&vcent[0],&vreso[0]);
    tgrcentreso->SetTitle("centreso;centrality;#sigma");
    tgrcentreso->Write("centreso");

    TGraph* tgrcentresomerge = new TGraph(vcentmerge.size(),&vcentmerge[0],&vresomerge[0]);
    tgrcentresomerge->SetTitle("centresomerge;centrality;#sigma");
    tgrcentresomerge->Write("centresomerge");

    TGraphErrors* tgrnpartreso = new TGraphErrors(vnpart.size(),&vnpart[0],&vreso[0],&vnparte[0]);
    tgrnpartreso->SetTitle("npartreso;N_{part};#sigma");
    tgrnpartreso->Write("npartreso");

    outfile->Write();
    outfile->Close();
}


int GetNpartDists(TH1F* hnpart[GAT::ncent], vector<double>& vnpart, vector<double>& vnparte)
{
    const char* centdistsname = "CentDists_9bins.root";
    TFile* fcentdists = new TFile(centdistsname,"READ");
    if (!fcentdists || fcentdists->IsZombie())
    {
        cerr<<"Could not open "<<centdistsname<<endl;
        return -1;
    }

    for (int icent = 0; icent < GAT::ncent; icent++)
    {
        hnpart[icent] = (TH1F*) gROOT->FindObject(Form("AuAu_Npart_%i_%i",GAT::limcent[icent],GAT::limcent[icent+1]));
        if (!hnpart[icent])
        {
            return -1;
        }
        hnpart[icent]->SetDirectory(0);
        vnpart.push_back(hnpart[icent]->GetMean());
        vnparte.push_back(hnpart[icent]->GetRMS());
    }
    fcentdists->Close();
    return 0;
}
