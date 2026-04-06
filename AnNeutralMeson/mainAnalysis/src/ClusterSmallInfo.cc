#include "ClusterSmallInfo.h"

ClusterSmallInfo::ClusterSmallInfo() :
  _eta(0),
  _phi(0),
  _ecore(0),
  _energy(0),
  _chi2(0)
{
}

ClusterSmallInfo::ClusterSmallInfo(const float& eta,
                                   const float& phi,
                                   const float& ecore,
                                   const float& energy,
                                   const float& chi2) :
  _eta(eta),
  _phi(phi),
  _ecore(ecore),
  _energy(energy),
  _chi2(chi2)
{
}

void ClusterSmallInfo::Reset()
{
  _eta = 0;
  _phi = 0;
  _ecore = 0;
  _energy = 0;
  _chi2 = 0;
}

void ClusterSmallInfo::Clear(Option_t* /*unused*/)
{
  _eta = 0;
  _phi = 0;
  _ecore = 0;
  _energy = 0;
  _chi2 = 0;
}

void ClusterSmallInfo::identify(std::ostream &os) const
{
  os << "ClusterSmallInfo with properties: "
     << "eta = " << _eta << ", "
     << "phi = " << _phi << ", "
     << "ecore = " << _ecore << ", "
     << "energy = " << _energy << ", "
     << "chi2 = " << _chi2 << std::endl;
}

void ClusterSmallInfo::set(const float& eta,
                           const float& phi,
                           const float& ecore,
                           const float& energy,
                           const float& chi2)
{
  _eta = eta;
  _phi = phi;
  _ecore = ecore;
  _energy = energy;
  _chi2 = chi2;
}
