#ifndef __CLUSTER_SMALL_INFO_H__
#define __CLUSTER_SMALL_INFO_H__

#include <phool/PHObject.h>

class ClusterSmallInfo : public PHObject
{
 public:
  ClusterSmallInfo();
  ClusterSmallInfo(const float&, const float&, const float&, const float&, const float&);

  void set(const float&, const float&, const float&, const float&, const float&);

  void identify(std::ostream &os = std::cout) const override;

  void Reset() override;

  //! Clear is used by TClonesArray to reset the small cluster to initial state
  // without calling destructor constructor
  void Clear(Option_t* = "") override;

  float get_eta() const { return _eta; }
  float get_phi() const { return _phi; }
  float get_ecore() const { return _ecore; }
  float get_energy() const { return _energy; }
  float get_chi2() const { return _chi2; }
  
 private:
  float _eta;
  float _phi;
  float _ecore;
  float _energy;
  float _chi2;

  ClassDefOverride(ClusterSmallInfo, 1);
};

#endif
