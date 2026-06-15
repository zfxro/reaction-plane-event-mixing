#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TTree.h>

using namespace std;

void saveobject(const char* objectname, const char* filelist = "list.txt", const char* outdir = "/phenix/scratch/zrowan/")
{
    ifstream ifs (filelist);
    string filename;
    while(getline(ifs,filename))
    {

        int dotpos = filename.find(".");
        string runno = filename.substr(dotpos-6,6);

        TFile* file = new TFile(filename.c_str(),"READ");
        TTree* oldtree = (TTree*)file->Get(objectname);
        if (oldtree == NULL) 
        {
            cout<<"Could not find object."<<endl;
            return;
        }

        ostringstream oss;
        oss<<outdir<<objectname<<"_"<<runno<<".root";
        cout<<oss.str()<<endl;

        TFile* outfile = new TFile(oss.str().c_str(),"RECREATE");
        TTree* newtree = oldtree->CloneTree();
        outfile->Write();
    }
}
