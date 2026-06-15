#ifndef __TRACKCOLLECTIONBASE_H_
#define __TRACKCOLLECTIONBASE_H_

#include <Rtypes.h>
#include "gamConst.h"

class ChargedTrack;

class TrackCollectionBase
{
protected:

    bool mixable;

public:

    virtual GA::tracktype_t GetTrackType() = 0;

    TrackCollectionBase();
    TrackCollectionBase(bool);
    virtual ~TrackCollectionBase();

    void SetMixable(bool);
    bool IsMixable();

    virtual unsigned int GetEntries() = 0;
    virtual bool GetEntry(unsigned int, ChargedTrack**) = 0;
    virtual void AddEntry(ChargedTrack*) = 0;

    ClassDef(TrackCollectionBase,1)
};

#endif
