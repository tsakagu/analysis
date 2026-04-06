#ifndef __CLUSTER_SMALL_INFO_CONTAINER_H__
#define __CLUSTER_SMALL_INFO_CONTAINER_H__

#include "ClusterSmallInfo.h"
#include <phool/PHObject.h>
#include <TClonesArray.h>
#include <cstdint>
#include <iostream>

class ClusterSmallInfoContainer : public PHObject
{
 public:
  ClusterSmallInfoContainer();
  ~ClusterSmallInfoContainer() override;
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;

  void set_live_trigger(uint64_t live_trigger) { _live_trigger = live_trigger; }
  void set_scaled_trigger(uint64_t scaled_trigger) { _scaled_trigger = scaled_trigger; }
  void set_bunchnumber(float bunchnumber) { _bunchnumber = bunchnumber; }
  //void set_bad_flag(bool bad_flag) { _bad_flag = bad_flag; }
  //void set_total_E(float total_E) { _total_E = total_E; }
  //void set_total_E_quiet(float total_E_quiet) { _total_E_quiet = total_E_quiet; }
  bool add_cluster(const float&, const float&, const float&, const float&, const float&);
  
  ClusterSmallInfo* get_cluster_at(int);
  const ClusterSmallInfo* get_cluster_at(int) const;

  uint64_t get_live_trigger() const { return _live_trigger; }
  uint64_t get_scaled_trigger() const { return _scaled_trigger; }
  float get_bunch_number() const { return _bunchnumber; }
  //bool get_bad_flag() const { return _bad_flag; }
  //bool get_total_E() const { return _total_E; }
  //bool get_total_E_quiet() const { return _total_E_quiet; }
  size_t size() const { return _size; }
  void compress();
  
 protected:
  size_t _size;
  TClonesArray *_clones = nullptr;
  uint64_t _live_trigger;
  uint64_t _scaled_trigger;
  float _bunchnumber;
  //bool _bad_flag = false; // If total EMCal E < 0 (event selection)
  //float _total_E = 0;
  //float _total_E_quiet = 0;

 private:
  ClassDefOverride(ClusterSmallInfoContainer, 1);
};

#endif
