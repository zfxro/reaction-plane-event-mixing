#include "gamConst.h"

const float PI = TMath::Pi();

void PlotRPResoCentrality(
    const char * inname="merged_histos.root"
    )
{
    TH1F* h1centreso = new TH1F("centreso","BBC #Psi_{2} Resolution;Centrality Bin;Resolution",GA::ncent_tree,0,GA::ncent_tree);

    //open file
    TFile* infile = new TFile(inname);

    for (int icent = 0; icent < GA::ncent_tree; icent++)
    {
        //obtain histo of interest
        TH1F* h1psidiff = (TH1F*) gROOT->FindObject(Form("psidiff_bbcns_cent%i",icent));

        //get nbins (needed for new histogram)
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
        float reso = sqrt(cosine_sum/nentry);

        h1centreso->SetBinContent(icent+1,reso);

    }
    infile->Close();
}
