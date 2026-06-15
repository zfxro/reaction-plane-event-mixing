//reads in root file with Psi_BBCN - Psi_BBCS
//plots a more intuitive representation
//calculates reaction plane resolution

const float PI = TMath::Pi();

void PlotRPReso(const char * inname =
        "gamana_150k_requireprecise.root"
    )
{
    //open file
    TFile* infile = new TFile(inname);

    //obtain histo of interest
    TH1F* h1psidiff = (TH1F*) gROOT->FindObject("psidiff_bbcns");

    //get nbins (needed for new histogram)
    int nx = h1psidiff->GetNbinsX();

    //book histograms that will be created
    TH1F* h1psidiff_adj = new TH1F("psidiff_adj","psidiff_adj;Psi_{BBCN} - Psi_{BBCS} (Adjusted);Events",nx/4,0.0,0.5);
    //hopefully nx is divisible by 4...

    float cosine_sum = 0.0;

    for (int ix = 1; ix <= nx; ix++)
    {
        float bincenter = h1psidiff->GetBinCenter(ix);
        float bincenter_new = fabs(bincenter);
        if (bincenter_new > .5) bincenter_new = 1. - bincenter_new;

        float bincontent = h1psidiff->GetBinContent(ix);

        h1psidiff_adj->Fill(bincenter_new, bincontent);

        float cosine_tmp = cos(2.*PI*bincenter);

        cosine_sum += cosine_tmp*bincontent;
    }

    //must get <Cos(2. * Psi_BBCN - Psi_BBCS)> so divide by number of entries
    //note: Pi is left out of above definition
    //Psi had been divided by Pi when filling histogram from infile

    float nentry = h1psidiff->GetEntries();
    float reso = sqrt(cosine_sum/nentry);

    cout<<reso<<endl;

    h1psidiff_adj->Draw();
}
