#ifndef __WRAPPER_FIT_FINAL_ARRAYS_H__
#define __WRAPPER_FIT_FINAL_ARRAYS_H__

#include "AsymmetryCalc/Constants.h"

#include <iostream>
#include <array>
#include <string>
#include <cstdlib>
#include <cstdint>

namespace AsymmetryCalc {
  class pt_final_array {
  public:
    static constexpr int nDimensions = 3;
    
    static constexpr size_t dim_sizes[nDimensions] = {ASYM_CONSTANTS::nParticles, ASYM_CONSTANTS::nPtBins, ASYM_CONSTANTS::nDirections};
    
    static constexpr const char* dim_names[nDimensions] = {"particle", "pt", "direction"};

    static constexpr size_t size_array = dim_sizes[0] * dim_sizes[1] * dim_sizes[2];
    
    using storage_type = std::array<double, size_array>;

    double& operator()(size_t iParticle,
                       size_t iPt,
                       size_t iDir)
    {
      check_bounds(iParticle, iPt, iDir);
      return data_[index(iParticle, iPt, iDir)];
    }

    const double& operator()(size_t iParticle,
                             size_t iPt,
                             size_t iDir) const

    {
      check_bounds(iParticle, iPt, iDir);
      return data_[index(iParticle, iPt, iDir)];
    }

  private:

    storage_type data_{};

    void check_bounds(size_t iParticle,
                      size_t iPt,
                      size_t iDir) const
    {
      size_t indices[nDimensions] = {iParticle, iPt, iDir};
      for (int iDim = 0; iDim < nDimensions; iDim++)
      {
        if (indices[iDim] > dim_sizes[iDim])
        {
          std::cerr << "Error: " << dim_names[iDim] << " index out of bounds" << std::endl;
          std::exit(1);
        }
      }
    }
    
    static constexpr size_t index(size_t iParticle,
                                  size_t iPt,
                                  size_t iDir)
    {
      return (((iParticle)
               * ASYM_CONSTANTS::nPtBins + iPt)
              * ASYM_CONSTANTS::nDirections + iDir);
    }
  };

  class pt_nodir_final_array {
  public:
    static constexpr int nDimensions = 2;
    
    static constexpr size_t dim_sizes[nDimensions] = {ASYM_CONSTANTS::nParticles, ASYM_CONSTANTS::nPtBins};
    
    static constexpr const char* dim_names[nDimensions] = {"particle", "pt"};

    static constexpr size_t size_array = dim_sizes[0] * dim_sizes[1];
    
    using storage_type = std::array<double, size_array>;

    double& operator()(size_t iParticle,
                       size_t iPt)
    {
      check_bounds(iParticle, iPt);
      return data_[index(iParticle, iPt)];
    }

    const double& operator()(size_t iParticle,
                             size_t iPt) const

    {
      check_bounds(iParticle, iPt);
      return data_[index(iParticle, iPt)];
    }

  private:

    storage_type data_{};

    void check_bounds(size_t iParticle,
                      size_t iPt) const
    {
      size_t indices[nDimensions] = {iParticle, iPt};
      for (int iDim = 0; iDim < nDimensions; iDim++)
      {
        if (indices[iDim] > dim_sizes[iDim])
        {
          std::cerr << "Error: " << dim_names[iDim] << " index out of bounds" << std::endl;
          std::exit(1);
        }
      }
    }
    
    static constexpr size_t index(size_t iParticle,
                                  size_t iPt)
    {
      return ((iParticle)
               * ASYM_CONSTANTS::nPtBins + iPt);
    }
  };

  class eta_final_array {
  public:
    static constexpr int nDimensions = 2;
    
    static constexpr size_t dim_sizes[nDimensions] = {ASYM_CONSTANTS::nParticles, ASYM_CONSTANTS::nEtaBins};
    
    static constexpr const char* dim_names[nDimensions] = {"particle", "eta"};

    static constexpr size_t size_array = dim_sizes[0] * dim_sizes[1];
    
    using storage_type = std::array<double, size_array>;

    double& operator()(size_t iParticle,
                       size_t iEta)
    {
      check_bounds(iParticle, iEta);
      return data_[index(iParticle, iEta)];
    }

    const double& operator()(size_t iParticle,
                             size_t iEta) const
    {
      check_bounds(iParticle, iEta);
      return data_[index(iParticle, iEta)];
    }

  private:

    storage_type data_{};

    void check_bounds(size_t iParticle,
                      size_t iEta) const
    {
      size_t indices[nDimensions] = {iParticle, iEta};
      for (int iDim = 0; iDim < nDimensions; iDim++)
      {
        if (indices[iDim] > dim_sizes[iDim])
        {
          std::cerr << "Error: " << dim_names[iDim] << " index out of bounds" << std::endl;
          std::exit(1);
        }
      }
    }
    
    static constexpr size_t index(size_t iParticle,
                                  size_t iEta)
    {
      return ((iParticle)
              * ASYM_CONSTANTS::nEtaBins + iEta);
    }
  };

  class xf_final_array {
  public:
    static constexpr int nDimensions = 2;
    
    static constexpr size_t dim_sizes[nDimensions] = {ASYM_CONSTANTS::nParticles, ASYM_CONSTANTS::nXfBins};
    
    static constexpr const char* dim_names[nDimensions] = {"particle", "xf"};

    static constexpr size_t size_array = dim_sizes[0] * dim_sizes[1];
    
    using storage_type = std::array<double, size_array>;

    double& operator()(size_t iParticle,
                       size_t iXf)
    {
      check_bounds(iParticle, iXf);
      return data_[index(iParticle, iXf)];
    }

    const double& operator()(size_t iParticle,
                             size_t iXf) const
    {
      check_bounds(iParticle, iXf);
      return data_[index(iParticle, iXf)];
    }

  private:

    storage_type data_{};

    void check_bounds(size_t iParticle,
                      size_t iXf) const
    {
      size_t indices[nDimensions] = {iParticle, iXf};
      for (int iDim = 0; iDim < nDimensions; iDim++)
      {
        if (indices[iDim] > dim_sizes[iDim])
        {
          std::cerr << "Error: " << dim_names[iDim] << " index out of bounds" << std::endl;
          std::exit(1);
        }
      }
    }
    
    static constexpr size_t index(size_t iParticle,
                                  size_t iXf)
    {
      return ((iParticle)
              * ASYM_CONSTANTS::nXfBins + iXf);
    }
  };
};

#endif
