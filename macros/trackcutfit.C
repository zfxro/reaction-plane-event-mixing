//used to determine track cuts
//prints mean and sigma of cut variables to file for "good" tracks
//reads output from trackCutDetermination.cc

#include <TMath.h>

const int ncut = 7;
const int npt = 9;
const int ncent = 5;

const string cutstr[ncut] =
{
    "pc3dphi",
    "pc3dz",
    "emcdphi",
    "emcdz",
    "chi2ndf",
    "dcaz",
    "dca2d"
};

void trackcutfit()
{
    ostringstream oss;
    oss<<"mkdir";
    for (int igood = 0; igood < 2; igood++)
    {
        oss<<" good"<<igood;
        oss<<" good"<<igood<<"/cent";
        oss<<" good"<<igood<<"/pt";
        for (int icent = 0; icent < ncent; icent++)
        {
            oss<<" good"<<igood<<"/cent/cent"<<icent;
        }
        for (int ipt = 0; ipt < npt; ipt++)
        {
            oss<<" good"<<igood<<"/pt/pt"<<ipt;
        }
    }
    system(oss.str().c_str());

    gStyle->SetTitleFontSize(0.03);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(111);
    gStyle->SetStatH(.15);
    gStyle->SetStatW(.2);
    gStyle->SetStatX(.9);
    gStyle->SetStatY(.9);


    TF1* f1 = new TF1("f1","[Constant]*TMath::Landau(x,[MPV],[Sigma],false)",0.,4.);

    TFile* infile = new TFile("trackcuts.root","READ");

    TH1F* histo[ncut][ncent][npt][2];
    TF1* fit[ncut][ncent][npt][2];

    for (int icut = 0; icut < ncut; icut++)
    {
        for (int icent = 0; icent < ncent; icent++)
        {
            for (int ipt = 0; ipt < npt; ipt++)
            {
                for (int isgood = 0; isgood < 2; isgood++)
                {
                    histo[icut][icent][ipt][isgood] = (TH1F*) infile->Get(Form("%s_cent%i_pt%i_goodtrk%i",cutstr[icut].c_str(),icent,ipt,isgood));
                    printf("%s_cent%i_pt%i_goodtrk%i\n",cutstr[icut].c_str(),icent,ipt,isgood);
                    if (cutstr[icut] == "chi2ndf")
                    {
                        fit[icut][icent][ipt][isgood] = new TF1(Form("fit_%s_cent%i_pt%i_goodtrk%i",cutstr[icut].c_str(),icent,ipt,isgood),"[Constant]*TMath::Landau(x,[MPV],[Sigma],false)",0.,4.);
                        histo[icut][icent][ipt][isgood]->GetXaxis()->SetRangeUser(0.,4.);
                    }
                    else
                    {
                        int firstbin = histo[icut][icent][ipt][isgood]->GetXaxis()->GetFirst();
                        int lastbin = histo[icut][icent][ipt][isgood]->GetXaxis()->GetLast();
                        float xmin = histo[icut][icent][ipt][isgood]->GetBinCenter(firstbin);
                        float xmax = histo[icut][icent][ipt][isgood]->GetBinCenter(lastbin);
                        fit[icut][icent][ipt][isgood] = new TF1(Form("fit_%s_cent%i_pt%i_goodtrk%i",cutstr[icut].c_str(),icent,ipt,isgood),"[Constant]*exp(-0.5*((x-[Mean])/[Sigma])*((x-[Mean])/[Sigma]))",xmin,xmax);
                    }
                    float maxval = histo[icut][icent][ipt][isgood]->GetMaximum();
                    float maxbin = histo[icut][icent][ipt][isgood]->GetMaximumBin();
                    float maxbincenter = histo[icut][icent][ipt][isgood]->GetBinCenter(maxbin);
                    float rms = histo[icut][icent][ipt][isgood]->GetRMS();
                    fit[icut][icent][ipt][isgood]->SetParameters(maxval,maxbincenter,rms*.1);
                    fit[icut][icent][ipt][isgood]->SetParLimits(2,0.,rms*1.);


                    histo[icut][icent][ipt][isgood]->Fit(Form("fit_%s_cent%i_pt%i_goodtrk%i",cutstr[icut].c_str(),icent,ipt,isgood),"RB");
                }
            }
        }
    }

    //=================================================================
    // Save fit parameters to file.  Will read later to use for cuts
    //=================================================================
    ofstream* ofs = new ofstream("fitpar.txt");
    for (int icut = 0; icut < ncut; icut++)
    {
        if (cutstr[icut] == "chi2ndf") continue;
        for (int icent = 0; icent < ncent; icent++)
        {
            for (int ipt = 0; ipt < npt; ipt++)
            {
                float mean = fit[icut][icent][ipt][1]->GetParameter("Mean");
                float sigma = fit[icut][icent][ipt][1]->GetParameter("Sigma");
                *ofs<<cutstr[icut]<<" "<<icent<<" "<<ipt<<" "<<mean<<" "<<sigma<<endl;
            }
        }
    }

    //=================================================================
    // Draw some stuff
    //=================================================================

    TCanvas* canpt[ncut][npt][2];

    for (int icut = 0; icut < ncut; icut++)
    {
        for (int ipt = 0; ipt < npt; ipt++)
        {
            for (int isgood = 0; isgood < 2; isgood++)
            {
                canpt[icut][ipt][isgood] = new TCanvas(Form("can_%s_pt%i_goodtrk%i",cutstr[icut].c_str(),ipt,isgood),"",1600,900);
                canpt[icut][ipt][isgood]->Divide(3,2);
                for (int icent = 0; icent < ncent; icent++)
                {
                    TPad* tp = (TPad*)canpt[icut][ipt][isgood]->cd(icent+1);
                    if (cutstr[icut] != "chi2ndf")
                    {
                        tp->SetLogy();
                        float mean = fit[icut][icent][ipt][isgood]->GetParameter("Mean");
                        float sigma = fit[icut][icent][ipt][isgood]->GetParameter("Sigma");
                        float nsig = 10.;
                        float xmin = mean - nsig*sigma;
                        float xmax = mean + nsig*sigma;
                        histo[icut][icent][ipt][isgood]->GetXaxis()->SetRangeUser(xmin,xmax);
                    }

                    histo[icut][icent][ipt][isgood]->Draw();
                    fit[icut][icent][ipt][isgood]->Draw("same");
                }
                canpt[icut][ipt][isgood]->SaveAs(Form("good%i/pt/pt%i/%s_pt%i_goodtrk%i.png",isgood,ipt,cutstr[icut].c_str(),ipt,isgood));
            }
        }
    }


    TCanvas* cancent[ncut][ncent][2];

    for (int icut = 0; icut < ncut; icut++)
    {
        for (int icent = 0; icent < ncent; icent++)
        {
            for (int isgood = 0; isgood < 2; isgood++)
            {
                cancent[icut][icent][isgood] = new TCanvas(Form("can_%s_cent%i_goodtrk%i",cutstr[icut].c_str(),icent,isgood),"",1600,900);
                cancent[icut][icent][isgood]->Divide(3,3);

                for (int ipt = 0; ipt < npt; ipt++)
                {
                    TPad* tp = (TPad*)cancent[icut][icent][isgood]->cd(ipt+1);
                    if (cutstr[icut] != "chi2ndf")
                    {
                        tp->SetLogy();
                        float mean = fit[icut][icent][ipt][isgood]->GetParameter("Mean");
                        float sigma = fit[icut][icent][ipt][isgood]->GetParameter("Sigma");
                        float nsig = 15.;
                        float xmin = mean - nsig*sigma;
                        float xmax = mean + nsig*sigma;
                        histo[icut][icent][ipt][isgood]->GetXaxis()->SetRangeUser(xmin,xmax);
                    }

                    histo[icut][icent][ipt][isgood]->Draw();
                    fit[icut][icent][ipt][isgood]->Draw("same");
                }
                cancent[icut][icent][isgood]->SaveAs(Form("good%i/cent/cent%i/%s_cent%i_goodtrk%i.png",isgood,icent,cutstr[icut].c_str(),icent,isgood));
            }
        }
    }
}
