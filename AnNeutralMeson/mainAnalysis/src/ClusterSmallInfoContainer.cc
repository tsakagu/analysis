#include "ClusterSmallInfoContainer.h"

#include <iostream>
#include <stdexcept>
#include <bitset>

ClusterSmallInfoContainer::ClusterSmallInfoContainer()
{
  _clones = new TClonesArray("ClusterSmallInfo", 50); // 50 is the maximum per event (very rarely reached in pp collisions)
  _clones->SetOwner();
  _clones->SetName("ClusterSmallInfoContainer");

  for (int i = 0; i < 50; ++i)
  {
    // as tower numbers are fixed per event
    // construct towers once per run, and clear the towers for first use
    _clones->ConstructedAt(i, "C");
  }
}

ClusterSmallInfoContainer::~ClusterSmallInfoContainer()
{
  delete _clones;
}

void ClusterSmallInfoContainer::Reset()
{
  _size = 0;
  _live_trigger = 0;
  _scaled_trigger = 0;
  _bunchnumber = 0;
  //_total_E = 0;
  // _total_E_quiet = 0;
  // delete _clones;
  
  // clear content of ClusterSmallInfo in the container for the next event
  for (Int_t i = 0; i < _clones->GetEntriesFast(); ++i)
  {
    TObject* obj = _clones->UncheckedAt(i);

    if (obj == nullptr)
    {
      std::cout << __PRETTY_FUNCTION__ << " Fatal access error:"
                << " _clones->GetSize() = " << _clones->GetSize()
                << " _clones->GetEntriesFast() = " << _clones->GetEntriesFast()
                << " i = " << i << std::endl;
      _clones->Print();
    }

    assert(obj);
    // same as TClonesArray::Clear() but only clear but not to erase all towers
    obj->Clear();
    obj->ResetBit(kHasUUID);
    obj->ResetBit(kIsReferenced);
    obj->SetUniqueID(0);
  }
}

void ClusterSmallInfoContainer::identify(std::ostream &os) const
{
  os << "ClusterSmallInfoContainer: live trigger = " << std::bitset<64>(_live_trigger) << std::endl;
  os << "ClusterSmallInfoContainer: scaled trigger = " << std::bitset<64>(_scaled_trigger) << std::endl;
  os << "ClusterSmallInfoContainer: bunchnumber = " << _bunchnumber << std::endl;
  //os << "ClusterSmallInfoContainer: bad flag = " << _bad_flag << std::endl;
  os << "ClusterSmallInfoContainer: size = " << _size << std::endl;
  for (size_t i = 0; i < _size; i++) {
    const ClusterSmallInfo *csi = get_cluster_at(i);
    csi->identify(os);
  }
}

bool ClusterSmallInfoContainer::add_cluster(const float& eta,
                                            const float& phi,
                                            const float& ecore,
                                            const float& energy,
                                            const float& chi2)
{
  if (_size >= 50)
  {
    std::cout << "ClusterSmallInfoContainer::_size has already reached its maximum (50). Next cluster won't be filled." << std::endl;
    return false;
  }
  
  ClusterSmallInfo *csi = (ClusterSmallInfo*) _clones->ConstructedAt(_size);
  csi->set(eta, phi, ecore, energy, chi2);
  _size++;

  return true;
}

ClusterSmallInfo* ClusterSmallInfoContainer::get_cluster_at(int pos)
{
  return (ClusterSmallInfo*) _clones->At(pos);
}

const ClusterSmallInfo* ClusterSmallInfoContainer::get_cluster_at(int pos) const
{
  return (ClusterSmallInfo*) _clones->At(pos);
}

void ClusterSmallInfoContainer::compress()
{
  // To remove all the O(50) zero entries in the final TTree.
  _clones->ExpandCreateFast(_size);
}
