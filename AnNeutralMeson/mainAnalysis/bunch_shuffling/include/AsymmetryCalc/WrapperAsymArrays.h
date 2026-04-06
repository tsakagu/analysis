#ifndef __WRAPPER_ASYM_ARRAYS_H__
#define __WRAPPER_ASYM_ARRAYS_H__

#include "AsymmetryCalc/Constants.h"

#include <iostream>
#include <array>
#include <string>
#include <cstdlib>
#include <cstdint>

namespace AsymmetryCalc {
  class pt_asym_array {
  public:
    static constexpr int nDimensions = 6;
    
    static constexpr size_t dim_sizes[nDimensions] = {ASYM_CONSTANTS::nBeams, ASYM_CONSTANTS::nParticles, ASYM_CONSTANTS::nRegions, ASYM_CONSTANTS::nPtBins, ASYM_CONSTANTS::nDirections, ASYM_CONSTANTS::nPhiBins};
    
    static constexpr const char* dim_names[nDimensions] = {"beam", "particle", "region", "pt", "direction", "phi"};

    static constexpr size_t size_array = dim_sizes[0] * dim_sizes[1] * dim_sizes[2] * dim_sizes[3] * dim_sizes[4] * dim_sizes[5];
    
    using storage_type = std::array<double, size_array>;

    double& operator()(size_t iBeam,
                         size_t iParticle,
                         size_t iRegion,
                         size_t iPt,
                         size_t iDir,
                         size_t iPhi)
    {
      check_bounds(iBeam, iParticle, iRegion, iPt, iDir, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iPt, iDir, iPhi)];
    }

    const double& operator()(size_t iBeam,
                                 size_t iParticle,
                                 size_t iRegion,
                                 size_t iPt,
                                 size_t iDir,
                                 size_t iPhi) const
    {
      check_bounds(iBeam, iParticle, iRegion, iPt, iDir, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iPt, iDir, iPhi)];
    }

  private:

    storage_type data_{};

    void check_bounds(size_t iBeam,
                      size_t iParticle,
                      size_t iRegion,
                      size_t iPt,
                      size_t iDir,
                      size_t iPhi) const
    {
      size_t indices[nDimensions] = {iBeam, iParticle, iRegion, iPt, iDir, iPhi};
      for (int iDim = 0; iDim < nDimensions; iDim++)
      {
        if (indices[iDim] > dim_sizes[iDim])
        {
          std::cerr << "Error: " << dim_names[iDim] << " index out of bounds" << std::endl;
          std::exit(1);
        }
      }
    }
    
    static constexpr size_t index(size_t iBeam,
                                  size_t iParticle,
                                  size_t iRegion,
                                  size_t iPt,
                                  size_t iDir,
                                  size_t iPhi)
    {
      return (((((iBeam * ASYM_CONSTANTS::nParticles + iParticle)
                 * ASYM_CONSTANTS::nRegions + iRegion)
                * ASYM_CONSTANTS::nPtBins + iPt)
               * ASYM_CONSTANTS::nDirections + iDir)
              * ASYM_CONSTANTS::nPhiBins +iPhi);
    }
  };

  class pt_nodir_asym_array {
  public:
    static constexpr int nDimensions = 5;
    
    static constexpr size_t dim_sizes[nDimensions] = {ASYM_CONSTANTS::nBeams, ASYM_CONSTANTS::nParticles, ASYM_CONSTANTS::nRegions, ASYM_CONSTANTS::nPtBins, ASYM_CONSTANTS::nPhiBins};
    
    static constexpr const char* dim_names[nDimensions] = {"beam", "particle", "region", "pt", "phi"};

    static constexpr size_t size_array = dim_sizes[0] * dim_sizes[1] * dim_sizes[2] * dim_sizes[3] * dim_sizes[4];
    
    using storage_type = std::array<double, size_array>;

    double& operator()(size_t iBeam,
                         size_t iParticle,
                         size_t iRegion,
                         size_t iPt,
                         size_t iPhi)
    {
      check_bounds(iBeam, iParticle, iRegion, iPt, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iPt, iPhi)];
    }

    const double& operator()(size_t iBeam,
                             size_t iParticle,
                             size_t iRegion,
                             size_t iPt,
                             size_t iPhi) const
    {
      check_bounds(iBeam, iParticle, iRegion, iPt, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iPt, iPhi)];
    }

  private:

    storage_type data_{};

    void check_bounds(size_t iBeam,
                      size_t iParticle,
                      size_t iRegion,
                      size_t iPt,
                      size_t iPhi) const
    {
      size_t indices[nDimensions] = {iBeam, iParticle, iRegion, iPt, iPhi};
      for (int iDim = 0; iDim < nDimensions; iDim++)
      {
        if (indices[iDim] > dim_sizes[iDim])
        {
          std::cerr << "Error: " << dim_names[iDim] << " index out of bounds" << std::endl;
          std::exit(1);
        }
      }
    }
    
    static constexpr size_t index(size_t iBeam,
                                  size_t iParticle,
                                  size_t iRegion,
                                  size_t iPt,
                                  size_t iPhi)
    {
      return ((((iBeam * ASYM_CONSTANTS::nParticles + iParticle)
                * ASYM_CONSTANTS::nRegions + iRegion)
               * ASYM_CONSTANTS::nPtBins + iPt)
              * ASYM_CONSTANTS::nPhiBins +iPhi);
    }
  };

  class eta_asym_array {
  public:
    static constexpr int nDimensions = 5;
    
    static constexpr size_t dim_sizes[nDimensions] = {ASYM_CONSTANTS::nBeams, ASYM_CONSTANTS::nParticles, ASYM_CONSTANTS::nRegions, ASYM_CONSTANTS::nEtaBins, ASYM_CONSTANTS::nPhiBins};
    
    static constexpr const char* dim_names[nDimensions] = {"beam", "particle", "region", "eta", "phi"};

    static constexpr size_t size_array = dim_sizes[0] * dim_sizes[1] * dim_sizes[2] * dim_sizes[3] * dim_sizes[4];
    
    using storage_type = std::array<double, size_array>;

    double& operator()(size_t iBeam,
                           size_t iParticle,
                           size_t iRegion,
                           size_t iEta,
                           size_t iPhi)
    {
      check_bounds(iBeam, iParticle, iRegion, iEta, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iEta, iPhi)];
    }

    const double& operator()(size_t iBeam,
                                 size_t iParticle,
                                 size_t iRegion,
                                 size_t iEta,
                                 size_t iPhi) const
    {
      check_bounds(iBeam, iParticle, iRegion, iEta, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iEta, iPhi)];
    }

  private:

    storage_type data_{};

    void check_bounds(size_t iBeam,
                      size_t iParticle,
                      size_t iRegion,
                      size_t iEta,
                      size_t iPhi) const
    {
      size_t indices[nDimensions] = {iBeam, iParticle, iRegion, iEta, iPhi};
      for (int iDim = 0; iDim < nDimensions; iDim++)
      {
        if (indices[iDim] > dim_sizes[iDim])
        {
          std::cerr << "Error: " << dim_names[iDim] << " index out of bounds" << std::endl;
          std::exit(1);
        }
      }
    }
    
    static constexpr size_t index(size_t iBeam,
                                  size_t iParticle,
                                  size_t iRegion,
                                  size_t iEta,
                                  size_t iPhi)
    {
      return ((((iBeam * ASYM_CONSTANTS::nParticles + iParticle)
                * ASYM_CONSTANTS::nRegions + iRegion)
               * ASYM_CONSTANTS::nEtaBins + iEta)
              * ASYM_CONSTANTS::nPhiBins +iPhi);
    }
  };

  class xf_asym_array {
  public:
    static constexpr int nDimensions = 5;
    
    static constexpr size_t dim_sizes[nDimensions] = {ASYM_CONSTANTS::nBeams, ASYM_CONSTANTS::nParticles, ASYM_CONSTANTS::nRegions, ASYM_CONSTANTS::nXfBins, ASYM_CONSTANTS::nPhiBins};
    
    static constexpr const char* dim_names[nDimensions] = {"beam", "particle", "region", "xf", "phi"};

    static constexpr size_t size_array = dim_sizes[0] * dim_sizes[1] * dim_sizes[2] * dim_sizes[3] * dim_sizes[4];
    
    using storage_type = std::array<double, size_array>;

    double& operator()(size_t iBeam,
                           size_t iParticle,
                           size_t iRegion,
                           size_t iXf,
                           size_t iPhi)
    {
      check_bounds(iBeam, iParticle, iRegion, iXf, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iXf, iPhi)];
    }

    const double& operator()(size_t iBeam,
                                 size_t iParticle,
                                 size_t iRegion,
                                 size_t iXf,
                                 size_t iPhi) const
    {
      check_bounds(iBeam, iParticle, iRegion, iXf, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iXf, iPhi)];
    }

  private:

    storage_type data_{};

    void check_bounds(size_t iBeam,
                      size_t iParticle,
                      size_t iRegion,
                      size_t iXf,
                      size_t iPhi) const
    {
      size_t indices[nDimensions] = {iBeam, iParticle, iRegion, iXf, iPhi};
      for (int iDim = 0; iDim < nDimensions; iDim++)
      {
        if (indices[iDim] > dim_sizes[iDim])
        {
          std::cerr << "Error: " << dim_names[iDim] << " index out of bounds" << std::endl;
          std::exit(1);
        }
      }
    }
    
    static constexpr size_t index(size_t iBeam,
                                  size_t iParticle,
                                  size_t iRegion,
                                  size_t iXf,
                                  size_t iPhi)
    {
      return ((((iBeam * ASYM_CONSTANTS::nParticles + iParticle)
                * ASYM_CONSTANTS::nRegions + iRegion)
               * ASYM_CONSTANTS::nXfBins + iXf)
              * ASYM_CONSTANTS::nPhiBins +iPhi);
    }
  };
};

#endif
