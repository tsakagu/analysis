// HistogramUtils.cpp
#include "CustomHistogramFit/HistogramUtils.h"
#include <iostream>

double HistogramUtils::GetMaxBinContentInRange(TH1* hist, double xMin, double xMax) {
  int binMin = hist->FindBin(xMin);
  int binMax = hist->FindBin(xMax);
  double maxContent = 0;
  for (int i = binMin; i <= binMax; ++i) {
    if (hist->GetBinContent(i) > maxContent) {
      /*std::cout << "x = " << hist->GetBinCenter(i) << ", "
        << "content = " << hist->GetBinContent(i) << std::endl;*/
      maxContent = hist->GetBinContent(i);
    }
  }
  return maxContent;
}

int HistogramUtils::FindFirstOccupiedBin(TH1 *hist, double threshold)
{
  if (!hist)
  {
    std::cerr << "FindFirstNonZeroBin -> Histogram not found.\n";
    return -1;
  }

  int nBins = hist->GetNbinsX();
  for (int ibin = 1; ibin <= nBins; ibin++)
  {
    if (hist->GetBinContent(ibin) > threshold) {
      return ibin;
    }
  }
  std::cerr << "FindFirstNonZeroBin -> Empty histogram.\n";
  return -1;
}
