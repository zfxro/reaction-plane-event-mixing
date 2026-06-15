#include "TrackCollectionBase.h"

ClassImp(TrackCollectionBase);

TrackCollectionBase::TrackCollectionBase()
{
    mixable = false;
}

TrackCollectionBase::TrackCollectionBase(bool canmix)
{
    mixable = canmix;
}

TrackCollectionBase::~TrackCollectionBase()
{
}

void TrackCollectionBase::SetMixable(bool canmix)
{
    mixable = canmix;
}

bool TrackCollectionBase::IsMixable()
{
    return mixable;
}
