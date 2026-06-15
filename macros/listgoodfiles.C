#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TTree.h>

using namespace std;

void listgoodfiles(const char* inlistname = "list.txt", const char* outlistname = "goodlist.txt")
{
    cout<<"inlist = "<<inlistname<<"   outlist = "<<outlistname<<endl;
    ofstream ofs (outlistname,ofstream::out);
    ifstream ifs (inlistname);
    string filename;
    while(getline(ifs,filename))
    {
        cout<<"Opening "<<filename<<endl;
        TFile* file = new TFile(filename.c_str(),"READ");
        if (file->IsZombie()) 
        {
            cout<<"   Is zombie"<<endl;
            file->Close();
            delete file;
            continue;
        }
        else
        {
            cout<<"   Is good"<<endl;
            ofs<<filename<<endl;
            file->Close();
            delete file;
        }
    }
    ifs.close();
    ofs.close();
    cout<<"All done!"<<endl;
}
