// ConstrainedFitter.h
#ifndef CONSTRAINEDFITTER_H
#define CONSTRAINEDFITTER_H

#include <iostream>
#include <vector>
#include "TF1.h"
#include "TH1.h"
#include "TMinuit.h"

#include "CustomHistogramFit/FitFunctions.h"

class ConstrainedFitter {
public:
  ConstrainedFitter(TH1* hist, std::vector<FitFunctions::FitModel> models);
  ~ConstrainedFitter();
  void SetInitialParams(const std::vector<double>& params);
  void SetParamBounds(const std::vector<std::pair<double, double>>& bounds);
  void SetFitRange(const std::vector<std::pair<double, double>>& fitRange) { fitRange_ = fitRange; }
  void PerformFit();
  void PrintResults();
  
  TF1* GetFitFunction() { return fitFunction_; }
  
 private:
  TH1* hist_;
  TF1* fitFunction_;
  std::string modelName_;
  int nParams_;
  std::vector<double> initialParams_;
  std::vector<std::pair<double, double>> paramBounds_;
  std::vector<std::string> paramNames_;
  double paramStep_;

  // (Potentially discontinuous) fitting range
  std::vector<std::pair<double,double>> fitRange_;

  // Vector of model functions
  std::vector<FitFunctions::modelFunction> singleModelFuns_;
  std::vector<int> singleNParams_;
  
  // Combined model function
  std::function<double(const double*, const double*)> CombinedModelFun_;
  
  void InitializeFunction();

  // Combine all input models in a single sum function
  double CombineModelFunctions(const double *x, const double *par);

  // Should be non-static to access hist_, fitFunction_, nParams_ ...
  double Chi2Function(const double *par);
};

#endif
