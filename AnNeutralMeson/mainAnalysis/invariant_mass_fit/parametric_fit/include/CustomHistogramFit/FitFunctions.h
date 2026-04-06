// FitFunctions.h
#ifndef FITFUNCTIONS_H
#define FITFUNCTIONS_H
#include <vector>
#include <string>
#include <map>

namespace FitFunctions {
  
  // Model function typedef
  typedef double (*modelFunction)(const double*, const double*);

  // Model structure
  struct modelStructure {
    modelFunction fun;
    std::string name;
    int nParams;
    std::vector<std::string> paramNames;
  };
  
  // Enum for all available fit models
  enum class FitModel { GAUSS, POLY2, POLY3, POLY5, ARGUSROOFIT, ARGUSMODIFIED, CHAT1, CHAT2, CHAT3 };

  // Get model information
  const modelStructure& GetModel(FitModel model);

  // Model functions prototypes
  double gauss(const double* x, const double* par);
  double poly2(const double* x, const double* par);
  double poly3(const double* x, const double* par);
  double poly5(const double* x, const double* par);
  double argusRoofit(const double* x, const double* par);
  double argusModified(const double* x, const double* par);
  double chat1(const double* x, const double* par);
  double chat2(const double* x, const double* par);
  double chat3(const double* x, const double* par);

  // Map linking model label to their respective information
  extern const std::map<FitModel, modelStructure> models;
};

#endif
