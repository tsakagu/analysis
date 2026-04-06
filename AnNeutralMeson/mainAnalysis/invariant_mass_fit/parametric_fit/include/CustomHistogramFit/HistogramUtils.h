// HistogramUtils.h
#ifndef HISTOGRAMUTILS_H
#define HISTOGRAMUTILS_H

#include "TH1.h"

class HistogramUtils {
public:
    static double GetMaxBinContentInRange(TH1* hist, double xMin, double xMax);
    static int FindFirstOccupiedBin(TH1* hist, double threshold);
};

#endif
