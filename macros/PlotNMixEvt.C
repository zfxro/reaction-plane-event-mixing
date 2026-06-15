//reads root file.  plots number of mixed events per tree bin

#include "gamConst.h"

void PlotNMixEvt(int iprocess = -1)
{
    const int verbosity = 0;

    const char* inlistname = "inlist.txt";
    const string outdir = "./";

    string inname;
    if (iprocess >= 0)
    {
        int filectr = 0;
        ifstream inlist(inlistname,ifstream::in);
        if (inlist.is_open())
        {
            while (getline(inlist,inname))
            {
                if (filectr==iprocess) break;
                filectr++;
            }
            inlist.close();
        }
        else
        {
            cerr<<"Cannot open "<<inlistname<<endl;
            return;
        }
    }
    else
    {
        inname = "datafile.root";
    }

    //==================================================
    // Open input file
    //==================================================
    TFile* infile = TFile::Open(inname.c_str());
    if (!infile)
    {
        cout<<"Line "<<__LINE__<<": Unable to open "<<inname <<endl;
        cout<<"Exiting."<<endl;

        return;
    }
    if (infile->IsZombie())
    {
        cout<<"Line "<<__LINE__<<": "<<inname<<" is zombie!"<<endl;
        cout<<"Exiting."<<endl;

        return;
    }
    else if (verbosity > 0)
    {
        cout<<"Input opened successfully"<<endl;
    }
    //==================================================
    // Input file has been opened
    //==================================================

    ostringstream outoss;
    if (iprocess >= 0)
    {
        int dotpos = inname.find(".");
        string str1 = inname.substr(dotpos-6,6);
        outoss<<outdir<<str1<<"_analyzed.root";
    }
    else
    {
        outoss<<outdir<<"test_analyzed.root";
    }

    //==================================================
    // Create output file and initialize histograms
    //==================================================
    TFile* outfile = new TFile(outoss.str().c_str(),"RECREATE");
    if (!outfile)
    {
        cout<<"Line "<<__LINE__<<": Unable to open "<<outoss.str() <<endl;
        cout<<"Exiting."<<endl;

        return;
    }
    else if (verbosity > 0)
    {
        cout<<"Output opened successfully"<<endl;
    }

    //============= Output ready =======================

    TH1F* h1nmix = new TH1F("nmix","nmixevtperbin",1000,0.,1000.);

    int nmix_intcent[93][GA::npsi_tree][GA::nzv_tree];
    memset(nmix_intcent,0,sizeof(nmix_intcent));

    for (int icent_tree = 0; icent_tree < GA::ncent_tree; icent_tree++)
    {
        for (int ipsi_tree = 0; ipsi_tree < GA::npsi_tree; ipsi_tree++)
        {
            for (int izv_tree = 0; izv_tree < GA::nzv_tree; izv_tree++)
            {
                TTree* trktree = NULL;
                trktree = (TTree*) infile->Get(Form("trktree_C%i_P%i_Z%i",icent_tree,ipsi_tree,izv_tree));
                if (!trktree)
                {
                    cout<<"Line "<<__LINE__<<": Unable to find trktree in "<<inname <<endl;
                    cout<<"Exiting."<<endl;
                    return;
                }
                else if (verbosity > 1)
                    cout<<"Found trktree "<<trktree->GetName()<<endl;

                int nmix = trktree->GetEntries();

                if (verbosity > 1)
                    cout<<"N mix event = "<<nmix<<endl;

                h1nmix->Fill(nmix);

                float cent_trk;
                TBranch* centBranch_trk = trktree->GetBranch("Centrality");
                centBranch_trk->SetAddress(&cent_trk);

                for (unsigned int ievt = 0; ievt < nmix; ievt++)
                {
                    trktree->GetEntry(ievt);

                    nmix_intcent[(int)cent_trk-1][ipsi_tree][izv_tree]++;
                }
            }
        }
    }

    TH1F* h1nmix_intcent = new TH1F("nmix_intcent","nmixevtperbin integer centrality",1000,0.,1000.);

    for (int icent = 0; icent < 93; icent++)
    {
        for (int ipsi = 0; ipsi < GA::npsi_tree; ipsi++)
        {
            for (int izv = 0; izv < GA::nzv_tree; izv++)
            {
                h1nmix_intcent->Fill(nmix_intcent[icent][ipsi][izv]);
                if (nmix_intcent[icent][ipsi][izv] < 20)
                {
                    printf("cent = %i; ipsi = %i; izv = %i; nmix = %i\n",icent,ipsi,izv,nmix_intcent[icent][ipsi][izv]);
                }
            }
        }
    }


    if (verbosity > 0)
        cout<<"Closing input file..."<<endl;
    infile->Close();

    if (verbosity > 0)
        cout<<"Writing output file "<<outfile->GetName()<<endl;
    outfile->Write();

    if (verbosity > 0)
        cout<<"Closing output file..."<<endl;
    outfile->Close();

    if (verbosity > 0)
        cout<<"Successful completion!"<<endl;

    return;

}
