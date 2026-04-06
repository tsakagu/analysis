#include "CustomHistogramFit/ConstrainedFitter.h"

#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

#include <algorithm>

ConstrainedFitter::ConstrainedFitter(TH1* hist, std::vector<FitFunctions::FitModel> models)
  : hist_(hist) {

  modelName_ = "";
  singleModelFuns_.resize(models.size());
  singleNParams_.resize(models.size());
  nParams_ = 0;
  int iparGlobal = 0;
  for (int im = 0; im < (int) models.size(); im++) {
    FitFunctions::FitModel model = models[im];
    FitFunctions::modelStructure modelInfo = FitFunctions::GetModel(model);
    modelName_ += (im > 0 ? " and " : "") + modelInfo.name;
    singleModelFuns_[im] = modelInfo.fun;
    singleNParams_[im] = modelInfo.nParams;
    nParams_ += modelInfo.nParams;
    int ipar = 0;
    //std::cout << "singleNParams_[" << im << "] = " << modelInfo.nParams << std::endl;
    for (; ipar < modelInfo.nParams; ipar++, iparGlobal++) {
      /*std::cout << "paramNames_ = [";
      for (int i = 0; i < (int) paramNames_.size() - 1; i++) {
        std::cout << paramNames_[i] << ", ";
      }
      std::cout << paramNames_[(int)paramNames_.size() - 1] << "]\n";*/
      if (std::find(paramNames_.begin(), paramNames_.end(), modelInfo.paramNames[ipar]) == paramNames_.end()) {
        paramNames_.push_back(modelInfo.paramNames[ipar]);
      } else {
        paramNames_.push_back(modelInfo.paramNames[ipar] + "_2");
      }
      //std::cout << "iparGlobal = " << iparGlobal << ", parameter: '" << paramNames_.back() << "'\n";
    }
    //std::cout << "model '" << modelInfo.name << "' " << modelInfo.nParams << " parameters\n";
  }
  //std::cout << "word 0\n";
  CombinedModelFun_ = [this] (const double *x, const double *par) { return this->CombineModelFunctions(x, par); };
  //std::cout << "word 1\n";
  //std::cout << "nParams_ = " << nParams_ << std::endl;
  initialParams_.resize(nParams_, 0.0);
  paramBounds_.resize(nParams_, {0.0, 0.0});
  paramStep_ = 0.001;
  std::cout << "combined model '" << modelName_ << "' " << nParams_ << " parameters\n";
    
  fitFunction_ = new TF1("fitFunction", CombinedModelFun_, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), nParams_);
  for (int ipar = 0; ipar < nParams_; ipar++) fitFunction_->SetParName(ipar, paramNames_[ipar].c_str());
}

ConstrainedFitter::~ConstrainedFitter()
{
  delete fitFunction_;
}

void ConstrainedFitter::SetInitialParams(const std::vector<double>& params) {
  if (params.size() != static_cast<size_t>(nParams_)) {
    std::cerr << "Error in ConstrainedFitter::SetInitialParams: Incorrect number of parameters provided! Expected "
              << nParams_ << " parameters for the model '" << modelName_ << "'. "
              << params.size() << " given instead." << std::endl;
    return;
  }
  initialParams_ = params;

  for (int ipar = 0; ipar < nParams_; ipar++) {
    fitFunction_->SetParameter(ipar, initialParams_[ipar]);
  }
}

void ConstrainedFitter::SetParamBounds(const std::vector<std::pair<double, double>>& bounds) {
  if (bounds.size() != static_cast<size_t>(nParams_)) {
    std::cerr << "Error in ConstrainedFitter::SetParamBounds: Incorrect number of bounds provided! Expected "
              << nParams_ << " parameters for the model '" << modelName_ << "'. "
              << bounds.size() << " given instead." << std::endl;
    return;
  }
  paramBounds_ = bounds;
}

void ConstrainedFitter::PerformFit() {
  if (!hist_) {
    std::cerr << "Error in ConstrainedFitter::PerformFit: No histogram provided!" << std::endl;
    return;
  }

  ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  // Make a lambda function
  // Chi2Function cannot be made static
  std::function<double(double const*)> const Chi2Lambda = [this] (const double *params) { return this->Chi2Function(params); };

  ROOT::Math::Functor functor(Chi2Lambda, nParams_);
  minimizer->SetFunction(functor);

  // Set:
  // - parameter name
  // - initial parameter
  // - parameter step
  // - parameter bounds
  for (int i = 0; i < nParams_; i++) {
    // std::cout << "minimizer->SetLimitedVariable(" << i << ", \"" << paramNames_[i] << "\", " << initialParams_[i] << ", " << paramStep_ << ", " << paramBounds_[i].first << ", " << paramBounds_[i].second << ");\n";
    minimizer->SetLimitedVariable(i, paramNames_[i], initialParams_[i], paramStep_, paramBounds_[i].first, paramBounds_[i].second);
  }

  bool success = minimizer->Minimize();

  if (!success) {
    std::cout << "Minimization failed!" << std::endl;
  } else {
    std::cout << "Minimization successful!" << std::endl;
  }
  std::cout << "Final Chi2: " << Chi2Function(minimizer->X()) << std::endl;

  const double *optimized_params = minimizer->X();
  const double *param_errors = minimizer->Errors();
  
  for (int i = 0; i < nParams_; i++)
  {
    fitFunction_->SetParameter(i, optimized_params[i]);
    fitFunction_->SetParError(i, param_errors[i]);
  }
}

void ConstrainedFitter::PrintResults() {
  if (!fitFunction_) {
    std::cerr << "Error: Fit function not initialized! Cannot print results" << std::endl;
    return;
  }
  
  std::cout << "Fit Results:" << std::endl;
  for (int i = 0; i < nParams_; ++i) {
    double value, error;
    value = fitFunction_->GetParameter(i);
    error = fitFunction_->GetParError(i);
    std::cout << "Parameter " << i << ": " << value << " +- " << error << std::endl;
  }
}

double ConstrainedFitter::CombineModelFunctions(const double* x, const double* par) {
  double sum = 0;
  int nPars = 0;
  int iparGlobal = 0; // parameter index in the combined model
  FitFunctions::modelFunction modelFun;
  for (int ifun = 0; ifun < (int) singleModelFuns_.size(); ifun++) {
    // Iterate over model functions
    modelFun = singleModelFuns_[ifun];
    nPars = singleNParams_[ifun];
    double *params = new double[nPars];
    int ipar = 0; // parameter index in the single model
    for (; ipar < nPars; ipar++, iparGlobal++) {
      params[ipar] = par[iparGlobal];
    }
    sum += modelFun(x, params);
    delete params;
  }
  return sum;
}

double ConstrainedFitter::Chi2Function(const double* par) {
  if (hist_ == nullptr) {
    std::cerr << "Error: Histogram not found in minimization function!" << std::endl;
    exit(1);
  }
  fitFunction_->SetParameters(par);
  
  double chi2 = 0;
  for (int bin = 2; bin <= hist_->GetNbinsX(); ++bin) {
    double x = hist_->GetBinCenter(bin);
    double dx = hist_->GetBinWidth(bin);
    double y_data_error = hist_->GetBinError(bin);
    double y_data = hist_->GetBinContent(bin);
    double y_fit = fitFunction_->Eval(x);

    // Skip bins outside of the fit range
    int nbWindows = fitRange_.size();
    // If no fit range is given, consider that the fit is performed over the full histogram
    if (nbWindows >= 1) {
      bool inFitRange = false;
      for (int iw = 0; iw < nbWindows; iw++) {
        if ((fitRange_[iw].first < x) && (x < fitRange_[iw].second)) {
          inFitRange = true;
          break;
        }
      }
      if (!inFitRange) continue;
    }

    // Decrease error in certain regions to constrain the fit
    /*if (x > 0.08 && x < 0.10) // Improve fit on left tail
      y_data_error /= 1000;
    if (0.18 < x) // Improve fit on right tail
      y_data_error /= 1000;
      if (y_data > 0.5 * peak_height) // Improve fit on the peak
      y_data_error /= 0.01;*/

    // Skip empty_bins
    if (y_data == 0 || y_data_error == 0) continue;
    

    // Check for reasonable values in fit evaluation
    if (std::isnan(y_fit) || std::isinf(y_fit))
    {
      std::cout << "NaN or Inf detected in fit evaluation at bin " << bin << std::endl;
      std::cout << "Parameters were the following:\n";
      for (int i = 0; i < nParams_; i++)
      {
        std::cout << "parameter " << i << ": " << par[i] << std::endl;
      }
      //exit(1);
    }

    double increment = std::pow((y_fit - y_data) / y_data_error, 2);
    if (std::isnan(increment) || std::isinf(increment))
      chi2 += 1e7;
    else
      chi2 += increment;
  }
  return chi2;
}

