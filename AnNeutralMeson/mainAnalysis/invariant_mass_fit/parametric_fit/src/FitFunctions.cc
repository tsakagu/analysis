#include "CustomHistogramFit/FitFunctions.h"
#include <iostream>
#include <cmath>

namespace FitFunctions {
  
  // Implementation of the fitting functions

  double gauss(const double* x, const double* par) {
    double shape =
      par[0] * std::exp(-0.5 * std::pow((x[0] - par[1]) / par[2], 2));
    return shape;
  }
  
  double poly2(const double* x, const double* par) {
    // Second order polynomial parametrized in terms of Bernstein polynomials for the background
    double poly =
      par[0] * std::pow(1 - x[0], 2) +
      par[1] * 2 * x[0] * (1 - x[0]) +
      par[2] * std::pow(x[0], 2);
    return poly;
  }

  double poly3(const double* x, const double* par) {
    // Third order polynomial parametrized in terms of Bernstein polynomials for the background
    double poly =
      par[0] * std::pow(1 - x[0], 3) +
      par[1] * 3 * x[0] * std::pow(1 - x[0], 2) +
      par[2] * 3 * std::pow(x[0], 2) * (1 - x[0]) +
      par[3] * std::pow(x[0], 3);
    return poly;
  }
  
  double poly5(const double* x, const double* par) {
    // Fifth order polynomial parametrized in terms of Bernstein polynomials for the background
    double poly =
      par[0] * std::pow(1 - x[0], 5) +
      par[1] * 5 * x[0] * std::pow(1 - x[0], 4) +
      par[2] * 10 * std::pow(x[0], 2) * std::pow(1 - x[0], 3) +
      par[3] * 10 * std::pow(x[0], 3) * std::pow(1 - x[0], 2) +
      par[4] * 5 * std::pow(x[0], 4) * (1 - x[0]) +
      par[5] * std::pow(x[0], 5);
    return poly;
  }
  
  double argusRoofit(const double* x, const double* par)   {
    // Argus Background just as it is defined in Roofit
    double argusBG =
      par[0] * x[0] *
    std::pow(1 - std::pow(x[0] / par[1], 2), par[3]) *
      std::exp(par[2] * std::pow(1 - x[0] / par[1], 2));
    return argusBG;
  }

  double argusModified(const double* x, const double* par) {
    // Revisited formula for the Argus background
    // Observable
    double m = x[0]; // Input invariant mass

    // Parameters
    double M = par[0]; // Maximum of the background
    double m_offset = par[1]; // start of the shape
    double m_max = par[2]; // x-coord for which the maximum is reached
    double p = par[3]; // Defines how much the background varies as x changes. The smaller p the flatter.
    double q = par[4]; // Defines how the background increases at low x. The smaller q, the steeper increase

    // Shift invariant mass wrt the offset
    double m_shift = std::abs(m - m_offset);
    double m_max_shift = std::max(std::abs(m_max - m_offset), 1e-6);
    
    // Compute the background
    double argusBG =
      M * std::exp(q * std::log(m_shift / m_max_shift) +
                   q / p * (1 - std::pow(m_shift / m_max_shift, p))); 
    return argusBG;
  }

  // Modified exponential function
  double chat1(const double* x, const double* par) {
    // Observable
    double M = x[0];

    // Parameters
    double A = par[0]; // Normalization constant
    double Mth = par[1]; // Mass threshold
    double B = par[2]; // Control for the steepness of the initial increase
    double C = par[3]; // Control for the steepness of the exponential decrease

    double result =
      A * std::exp(B * std::log(M - Mth) - C * ( M - Mth));
    return result;
  }

  // Threshold function
  double chat2(const double* x, const double* par) {
    // Observable
    double M = x[0];

    // Parameters
    double A = par[0]; // Normalization constant
    double Mth = par[1]; // Mass threshold
    double B = par[2]; // Control for the steepness of the initial increase
    double C = par[3]; // Control for the steepness of the exponential decrease

    double result =
      A * (1 - std::exp(- (M - Mth) / B )) * std::exp(- C * (M - Mth));
    return result;
  }

  // Polynomial function with threshold
  double chat3(const double* x, const double* par) {
    // Observable
    double M = x[0];

    // Parameters
    double A = par[0]; // Normalization constant
    double Mth = par[1]; // Mass threshold
    double B = par[2]; // First shape parameter
    double C = par[3]; // Second shape parameter
    double D = par[4]; // Third shape parameter

    double result =
      A * std::pow(M - Mth, B) * (1 + C * (M - Mth) + D * std::pow(M - Mth, 2));
    return result;
  }
  
  // Table of models 
  const std::map<FitModel, modelStructure> models = {
    {FitModel::GAUSS, {gauss, "Gaussian function", 3, {"Norm", "mean", "sigma"}}},
    {FitModel::POLY2, {poly2, "Second order polynomial (Bernstein)", 3, {"b0", "b1", "b2"}}},
    {FitModel::POLY3, {poly3, "Third order polynomial (Bernstein)", 4, {"b0", "b1", "b2", "b3"}}},
    {FitModel::POLY5, {poly5, "Fifth order polynomial (Bernstein)", 6, {"b0", "b1", "b2", "b3", "b4", "b5"}}},
    {FitModel::ARGUSROOFIT, {argusRoofit, "Argus Background (same as RooFit)", 4, {"N", "m0", "c", "p"}}},
    {FitModel::ARGUSMODIFIED, {argusModified, "Argus Background (revisited)", 5, {"M", "m_offset", "m_max", "p", "q"}}},
    {FitModel::CHAT1, {chat1, "Modified exponential function", 4, {"A", "Mth", "B", "C"}}},
    {FitModel::CHAT2, {chat2, "Threshold function", 4, {"A", "Mth", "B", "C"}}},
    {FitModel::CHAT3, {chat3, "Polynomial function with threshold", 5, {"A", "Mth", "B", "C", "D"}}}};

  // Function to get model information
  const modelStructure& GetModel(FitModel model) {
    auto it = models.find(model);
    if (it == models.end()) {
      std::cerr << "Error in FitFunctions::GetModel: Requested fit model does not exist.\n";
      exit(1);
    }
    return it->second;
  }
  
} // namespace FitFunctions

  
    
