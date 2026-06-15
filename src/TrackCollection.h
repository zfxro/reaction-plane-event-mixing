#ifndef __TRACKCOLLECTION_H_
#define __TRACKCOLLECTION_H_

#include <vector>
#include <Rtypes.h>

#include "TrackCollectionBase.h"
#include "ChargedTrack.h"
#include "CentralArmTrack.h"
#include "gamConst.h"

template <typename TrackType>
class TrackCollection : public TrackCollectionBase
{
public:

    std::vector <TrackType> trk_v;

    TrackCollection();
    TrackCollection(bool);
    virtual ~TrackCollection();

    GA::tracktype_t GetTrackType();
    unsigned int GetEntries();
    bool GetEntry(unsigned int, ChargedTrack**);
    bool GetEntry(unsigned int, TrackType**);
    void AddEntry(ChargedTrack*);
    void AddEntry(TrackType*);

    ClassDef(TrackCollection,1)
};

#endif
