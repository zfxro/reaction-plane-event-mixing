#ifndef __CENTRALTRACKDATA_H_
#define __CENTRALTRACKDATA_H_

#include <Rtypes.h>

class CentralTrackData
{
public:
    bool Pc3Ok();
    bool EmcOk();
    bool ZedOk(float);
    bool QualOk(int);

    unsigned char GetEmcSigChar() {return emcsig;};
    unsigned char GetPc3SigChar() {return pc3sig;};

    void SetEmcSig(float);
    void SetPc3Sig(float);

    bool EmcMatch(float);
    bool Pc3Match(float);

    bool CheckDEP(float);
    int SetDEP(float);
    char GetDEPChar() {return dep;};

    void Print();

    CentralTrackData();
    virtual ~CentralTrackData();

private:

    unsigned char pc3sig;
    unsigned char emcsig;
    char dep;

public:

    bool isn0neg;
    bool arm;

public:

    ClassDef(CentralTrackData,1)
};

#endif
