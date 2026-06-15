#ifndef __CENTRALARMTRACK_H_
#define __CENTRALARMTRACK_H_

#include "ChargedTrack.h"
#include "CentralTrackData.h"
#include "gamConst.h"
#include <Rtypes.h>

class CentralArmTrack: public ChargedTrack
{
public:

    GA::tracktype_t GetTrackType() {return GA::CNT;};

    unsigned char GetEmcSigChar() {return cndata.GetEmcSigChar();};
    unsigned char GetPc3SigChar() {return cndata.GetPc3SigChar();};

    unsigned char GetDEPChar() {return cndata.GetDEPChar();};

    float GetIsN0Neg() {return cndata.isn0neg;};

    void SetEmcSig(float emcsig) {cndata.SetEmcSig(emcsig);};
    void SetPc3Sig(float pc3sig) {cndata.SetPc3Sig(pc3sig);};
    int SetDEP(float dep) {return cndata.SetDEP(dep);};
    void SetIsN0Neg(bool isn0neg) {cndata.isn0neg = isn0neg;};

    bool Pc3Ok() {return cndata.Pc3Ok();};
    bool EmcOk() {return cndata.EmcOk();};
    bool ZedOk(float zed) {return cndata.ZedOk(zed);};
    bool QualOk(int quality) {return cndata.QualOk(quality);};

    bool EmcMatch(float nsig) {return cndata.EmcMatch(nsig);};
    bool Pc3Match(float nsig) {return cndata.Pc3Match(nsig);};

    CentralTrackData* GetCentralTrackData() {return &cndata;};

    void Print() {ChargedTrack::Print(); cndata.Print();};

    CentralArmTrack();
    virtual ~CentralArmTrack();

private:

    CentralTrackData cndata;

public:

    ClassDef(CentralArmTrack,1)
};

#endif
