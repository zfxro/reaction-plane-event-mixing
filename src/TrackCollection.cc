#include "TrackCollection.h"

template <typename TrackType>
TrackCollection<TrackType>::TrackCollection()
    : TrackCollectionBase()
{
}

template <typename TrackType>
TrackCollection<TrackType>::TrackCollection(bool canmix)
    : TrackCollectionBase(canmix)
{
}

template <typename TrackType>
TrackCollection<TrackType>::~TrackCollection()
{
    trk_v.clear();
}

template <typename TrackType>
GA::tracktype_t TrackCollection<TrackType>::GetTrackType()
{
    TrackType tmptrk;
    return tmptrk.GetTrackType();
}

template <typename TrackType>
unsigned int TrackCollection<TrackType>::GetEntries()
{
    return trk_v.size();
}

template <typename TrackType>
bool TrackCollection<TrackType>::GetEntry(unsigned int itrk, ChargedTrack** trk)
{
    if (itrk >= trk_v.size())
    {
        return false;
    }

    *trk = &trk_v[itrk];
    return true;
}

template <typename TrackType>
bool TrackCollection<TrackType>::GetEntry(unsigned int itrk, TrackType** trk)
{
    if (itrk >= trk_v.size())
    {
        return false;
    }

    *trk = &trk_v[itrk];
    return true;
}

template <typename TrackType>
void TrackCollection<TrackType>::AddEntry(ChargedTrack* trk)
{
    TrackType* casttrk = dynamic_cast<TrackType*>(trk);
    trk_v.push_back(*casttrk);
}

template <typename TrackType>
void TrackCollection<TrackType>::AddEntry(TrackType* trk)
{
    trk_v.push_back(*trk);
}

template class TrackCollection<CentralArmTrack>;
