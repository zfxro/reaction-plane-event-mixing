void findlownmixrun()
{
    const int ntype = 3;
    const char* spec[ntype] = {"dir","dec","pi0"};

    const int ncent = 2;
    const int limcent[ncent+1] = {0,30,100};

    const int nptt = 3;
    const float limptt[nptt] = {4.,6.,10.};

    const char* filelistname = "filelist.txt";
    string filename;
    ifstream filelist(filelistname,ifstream::in);
    while (getline(filelist,filename))
    {
        TFile* file = TFile::Open(filename.c_str());

        float nmix_pttcent[ntype][nptt+1][ncent+1];
        memset(nmix_pttcent,0,sizeof(nmix_pttcent));

        for (int itype = 0; itype < ntype; itype++)
        {
            TH2F* h2tmp = (TH2F*) gROOT->FindObject(Form("nmixevt_centptt_%s",spec[itype]));
            if (h2tmp == NULL)
            {
                printf("nmixevt_centptt_%s not found!",spec[itype]);
                return;
            }
            for (int iptt = 0; iptt < nptt; iptt++)
            {
                for (int icent = 0; icent < ncent; icent++)
                {
                    float tmpbc = h2tmp->GetBinContent(icent+1,iptt+1);
                    nmix_pttcent[itype][iptt][icent] = tmpbc;
                    nmix_pttcent[itype][iptt][ncent] += tmpbc;
                    nmix_pttcent[itype][nptt][icent] += tmpbc;
                    nmix_pttcent[itype][nptt][ncent] += tmpbc;
                }
            }
        }

        float ntrg_pttcent[ntype][nptt+1][ncent+1];
        memset(ntrg_pttcent,0,sizeof(ntrg_pttcent));

        for (int itype = 0; itype < ntype; itype++)
        {
            TH2F* h2tmp = (TH2F*) gROOT->FindObject(Form("ctr_centptt_%s",spec[itype]));
            if (h2tmp == NULL)
            {
                printf("ntrgevt_centptt_%s not found!",spec[itype]);
                return;
            }
            for (int iptt = 0; iptt < nptt; iptt++)
            {
                for (int icent = 0; icent < ncent; icent++)
                {
                    float tmpbc = h2tmp->GetBinContent(icent+1,iptt+1);
                    ntrg_pttcent[itype][iptt][icent] = tmpbc;
                    ntrg_pttcent[itype][iptt][ncent] += tmpbc;
                    ntrg_pttcent[itype][nptt][icent] += tmpbc;
                    ntrg_pttcent[itype][nptt][ncent] += tmpbc;
                }
            }
        }

        //==================================================
        for (int itype = 0; itype < ntype; itype++)
        {
            for (int iptt = 0; iptt < nptt; iptt++)
            {
                for (int icent = 0; icent < ncent; icent++)
                {
                    if (ntrg_pttcent[itype][iptt][icent] != 0)
                    {
                        if (nmix_pttcent[itype][iptt][icent] < 5)
                        {
                            printf("File %s; Trigger: type %i, ptt %i, cent %i, ntrg %i; Nmix = %i;",
                                   file.c_str(),itype,iptt,icent,ntrg_pttcent[itype][iptt][icent],nmix_pttcent[itype][iptt][icent]);
                        }
                    }
                }
            }
        }
        file->Close();
        delete file;
    }
}
