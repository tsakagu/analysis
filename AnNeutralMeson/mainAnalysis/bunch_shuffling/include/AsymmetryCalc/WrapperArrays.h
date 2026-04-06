#ifndef __WRAPPER_ARRAYS_H__
#define __WRAPPER_ARRAYS_H__

#include "AsymmetryCalc/Constants.h"

#include <iostream>
#include <array>
#include <string>
#include <cstdlib>
#include <cstdint>

namespace AsymmetryCalc {
  class pt_yield_array {
  public:
    static constexpr int nDimensions = 7;
    
    static constexpr size_t dim_sizes[nDimensions] = {ASYM_CONSTANTS::nBeams, ASYM_CONSTANTS::nParticles, ASYM_CONSTANTS::nRegions, ASYM_CONSTANTS::nPtBins, ASYM_CONSTANTS::nDirections, ASYM_CONSTANTS::nSpins, ASYM_CONSTANTS::nPhiBins};
    
    static constexpr const char* dim_names[nDimensions] = {"beam", "particle", "region", "pt", "direction", "spin", "phi"};

    static constexpr size_t size_array = dim_sizes[0] * dim_sizes[1] * dim_sizes[2] * dim_sizes[3] * dim_sizes[4] * dim_sizes[5] * dim_sizes[6];
    
    using storage_type = std::array<uint32_t, size_array>;

    uint32_t& operator()(size_t iBeam,
                         size_t iParticle,
                         size_t iRegion,
                         size_t iPt,
                         size_t iDir,
                         size_t iSpin,
                         size_t iPhi)
    {
      check_bounds(iBeam, iParticle, iRegion, iPt, iDir, iSpin, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iPt, iDir, iSpin, iPhi)];
    }

    const uint32_t& operator()(size_t iBeam,
                                 size_t iParticle,
                                 size_t iRegion,
                                 size_t iPt,
                                 size_t iDir,
                                 size_t iSpin,
                                 size_t iPhi) const
    {
      check_bounds(iBeam, iParticle, iRegion, iPt, iDir, iSpin, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iPt, iDir, iSpin, iPhi)];
    }

    pt_yield_array& operator+=(const pt_yield_array& other)
    {
      for (size_t i = 0; i < size_array; ++i) {
        data_[i] += other.data_[i];
      }
      return *this;
    }

  private:

    storage_type data_{};

    void check_bounds(size_t iBeam,
                      size_t iParticle,
                      size_t iRegion,
                      size_t iPt,
                      size_t iDir,
                      size_t iSpin,
                      size_t iPhi) const
    {
      size_t indices[nDimensions] = {iBeam, iParticle, iRegion, iPt, iDir, iSpin, iPhi};
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
                                  size_t iSpin,
                                  size_t iPhi)
    {
      return ((((((iBeam * ASYM_CONSTANTS::nParticles + iParticle)
                  * ASYM_CONSTANTS::nRegions + iRegion)
                 * ASYM_CONSTANTS::nPtBins + iPt)
                * ASYM_CONSTANTS::nDirections + iDir)
               * ASYM_CONSTANTS::nSpins + iSpin)
              * ASYM_CONSTANTS::nPhiBins +iPhi);
    }
  };

  class eta_yield_array {
  public:
    static constexpr int nDimensions = 6;
    
    static constexpr size_t dim_sizes[nDimensions] = {ASYM_CONSTANTS::nBeams, ASYM_CONSTANTS::nParticles, ASYM_CONSTANTS::nRegions, ASYM_CONSTANTS::nEtaBins, ASYM_CONSTANTS::nSpins, ASYM_CONSTANTS::nPhiBins};
    
    static constexpr const char* dim_names[nDimensions] = {"beam", "particle", "region", "eta", "spin", "phi"};

    static constexpr size_t size_array = dim_sizes[0] * dim_sizes[1] * dim_sizes[2] * dim_sizes[3] * dim_sizes[4] * dim_sizes[5];
    
    using storage_type = std::array<uint32_t, size_array>;

    uint32_t& operator()(size_t iBeam,
                           size_t iParticle,
                           size_t iRegion,
                           size_t iEta,
                           size_t iSpin,
                           size_t iPhi)
    {
      check_bounds(iBeam, iParticle, iRegion, iEta, iSpin, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iEta, iSpin, iPhi)];
    }

    const uint32_t& operator()(size_t iBeam,
                                 size_t iParticle,
                                 size_t iRegion,
                                 size_t iEta,
                                 size_t iSpin,
                                 size_t iPhi) const
    {
      check_bounds(iBeam, iParticle, iRegion, iEta, iSpin, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iEta, iSpin, iPhi)];
    }

    eta_yield_array& operator+=(const eta_yield_array& other)
    {
      for (size_t i = 0; i < size_array; ++i) {
        data_[i] += other.data_[i];
      }
      return *this;
    }

  private:

    storage_type data_{};

    void check_bounds(size_t iBeam,
                      size_t iParticle,
                      size_t iRegion,
                      size_t iEta,
                      size_t iSpin,
                      size_t iPhi) const
    {
      size_t indices[nDimensions] = {iBeam, iParticle, iRegion, iEta, iSpin, iPhi};
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
                                  size_t iSpin,
                                  size_t iPhi)
    {
      return (((((iBeam * ASYM_CONSTANTS::nParticles + iParticle)
                  * ASYM_CONSTANTS::nRegions + iRegion)
                 * ASYM_CONSTANTS::nEtaBins + iEta)
                * ASYM_CONSTANTS::nSpins + iSpin)
               * ASYM_CONSTANTS::nPhiBins +iPhi);
    }
  };

  class xf_yield_array {
  public:
    static constexpr int nDimensions = 6;
    
    static constexpr size_t dim_sizes[nDimensions] = {ASYM_CONSTANTS::nBeams, ASYM_CONSTANTS::nParticles, ASYM_CONSTANTS::nRegions, ASYM_CONSTANTS::nXfBins, ASYM_CONSTANTS::nSpins, ASYM_CONSTANTS::nPhiBins};
    
    static constexpr const char* dim_names[nDimensions] = {"beam", "particle", "region", "xf", "spin", "phi"};

    static constexpr size_t size_array = dim_sizes[0] * dim_sizes[1] * dim_sizes[2] * dim_sizes[3] * dim_sizes[4] * dim_sizes[5];
    
    using storage_type = std::array<uint32_t, size_array>;

    uint32_t& operator()(size_t iBeam,
                           size_t iParticle,
                           size_t iRegion,
                           size_t iXf,
                           size_t iSpin,
                           size_t iPhi)
    {
      check_bounds(iBeam, iParticle, iRegion, iXf, iSpin, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iXf, iSpin, iPhi)];
    }

    const uint32_t& operator()(size_t iBeam,
                                 size_t iParticle,
                                 size_t iRegion,
                                 size_t iXf,
                                 size_t iSpin,
                                 size_t iPhi) const
    {
      check_bounds(iBeam, iParticle, iRegion, iXf, iSpin, iPhi);
      return data_[index(iBeam, iParticle, iRegion, iXf, iSpin, iPhi)];
    }

    xf_yield_array& operator+=(const xf_yield_array& other)
    {
      for (size_t i = 0; i < size_array; ++i) {
        data_[i] += other.data_[i];
      }
      return *this;
    }

  private:

    storage_type data_{};

    void check_bounds(size_t iBeam,
                      size_t iParticle,
                      size_t iRegion,
                      size_t iXf,
                      size_t iSpin,
                      size_t iPhi) const
    {
      size_t indices[nDimensions] = {iBeam, iParticle, iRegion, iXf, iSpin, iPhi};
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
                                  size_t iSpin,
                                  size_t iPhi)
    {
      return (((((iBeam * ASYM_CONSTANTS::nParticles + iParticle)
                  * ASYM_CONSTANTS::nRegions + iRegion)
                 * ASYM_CONSTANTS::nXfBins + iXf)
                * ASYM_CONSTANTS::nSpins + iSpin)
               * ASYM_CONSTANTS::nPhiBins +iPhi);
    }
  };
};

#endif
