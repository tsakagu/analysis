#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TTree.h"
#include "TVector2.h"

#include "./src/sigmaEff.h"

TGaxis::SetMaxDigits(3);

static double median1D(std::vector<double> v)
{
    if (v.empty())
        return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    if (n % 2 == 1)
        return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// Geometric median via Weiszfeld's algorithm
static TVector2 geometricMedian2D(const std::vector<TVector2> &pts, double tol = 1e-9, int maxIter = 1000, double eps = 1e-12)
{
    if (pts.empty())
        return TVector2(0, 0);

    // Initialize with coordinate-wise median (good robust start)
    std::vector<double> xs, ys;
    xs.reserve(pts.size());
    ys.reserve(pts.size());
    for (const auto &p : pts)
    {
        xs.push_back(p.X());
        ys.push_back(p.Y());
    }
    TVector2 xk(median1D(xs), median1D(ys));

    for (int it = 0; it < maxIter; ++it)
    {
        double numX = 0.0, numY = 0.0, den = 0.0;

        // If close to a data point, return it
        for (const auto &pi : pts)
        {
            const double dx = xk.X() - pi.X();
            const double dy = xk.Y() - pi.Y();
            const double di = std::sqrt(dx * dx + dy * dy);
            if (di < eps)
                return pi;
        }

        for (const auto &pi : pts)
        {
            const double dx = xk.X() - pi.X();
            const double dy = xk.Y() - pi.Y();
            const double di = std::sqrt(dx * dx + dy * dy);
            const double wi = 1.0 / di; // unweighted Weiszfeld
            numX += wi * pi.X();
            numY += wi * pi.Y();
            den += wi;
        }

        TVector2 xnext(numX / den, numY / den);
        const double step = (xnext - xk).Mod();
        xk = xnext;

        if (step < tol)
            break;
    }

    return xk;
}

static void draw1Dhistogram(TH1 *hist,                                 //
                            const bool logy,                           //
                            const std::string &xtitle,                 //
                            const std::string &ytitle,                 //
                            const std::vector<std::string> &addinfo,   //
                            const std::string &drawoption = "hist e1", //
                            const std::string &filename = "plot1d")
{
    if (!hist)
        return;

    TCanvas *c = new TCanvas("c_1d", "c_1d", 900, 800);
    c->SetLogy(logy);
    c->cd();
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.05);
    // gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.08);

    hist->SetTitle("");
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    hist->GetYaxis()->SetTitle(ytitle.c_str());
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->GetYaxis()->SetRangeUser((logy ? hist->GetMinimum(0) * 0.5 : 0.0), hist->GetMaximum() * (logy ? 10.0 : 1.4));
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(kFullCircle);
    hist->SetMarkerSize(1.0);
    hist->SetMarkerColor(kBlack);
    hist->Draw(drawoption.c_str());

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    // set right- and bottom aligned for now, which is only good for 1 line
    latex.SetTextAlign(kHAlignRight + kVAlignBottom);
    for (size_t i = 0; i < addinfo.size(); ++i)
    {
        latex.DrawLatex(1 - gPad->GetRightMargin(), (1 - gPad->GetTopMargin() + 0.01) - 0.045 * i, addinfo[i].c_str());
    }

    c->SaveAs((filename + ".png").c_str());
    c->SaveAs((filename + ".pdf").c_str());
    delete c;
}

static void draw1Dhistogram_GaussianFit(TH1 *hist,                                 //
                                        const bool logy,                           //
                                        const std::string &xtitle,                 //
                                        const std::string &ytitle,                 //
                                        const std::vector<std::string> &addinfo,   //
                                        const std::string &drawoption = "hist e1", //
                                        const std::string &filename = "plot1d")
{
    if (!hist)
        return;

    TCanvas *c = new TCanvas("c_1d", "c_1d", 900, 800);
    c->SetLogy(logy);
    c->cd();
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.05);
    // gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.08);

    hist->SetTitle("");
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    hist->GetYaxis()->SetTitle(ytitle.c_str());
    hist->GetYaxis()->SetTitleOffset(1.4);
    hist->GetYaxis()->SetRangeUser((logy ? hist->GetMinimum(0) * 0.5 : 0.0), hist->GetMaximum() * (logy ? 10.0 : 1.4));
    hist->SetLineWidth(2);
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(kFullCircle);
    hist->SetMarkerSize(1.0);
    hist->SetMarkerColor(kBlack);
    hist->Draw(drawoption.c_str());

    // Fit with Gaussian
    TF1 *gausFit = new TF1("gausFit", "gaus", -15, 15);
    hist->Fit(gausFit, "RQ");
    gausFit->SetLineColor(kRed);
    gausFit->SetLineWidth(2);
    gausFit->Draw("same");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    // set right- and bottom aligned for now, which is only good for 1 line
    latex.SetTextAlign(kHAlignRight + kVAlignBottom);
    for (size_t i = 0; i < addinfo.size(); ++i)
    {
        latex.DrawLatex(1 - gPad->GetRightMargin(), (1 - gPad->GetTopMargin() + 0.01) - 0.045 * i, addinfo[i].c_str());
    }

    TLegend *leg = new TLegend(gPad->GetLeftMargin() + 0.05,    //
                               1 - gPad->GetTopMargin() - 0.15, //
                               gPad->GetLeftMargin() + 0.2,    //
                               1 - gPad->GetTopMargin() - 0.04  //
    );
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(hist, "Data", "pel");
    leg->AddEntry(gausFit, Form("Gaussian Fit: #mu = %.3g [cm], #sigma = %.3f [cm]", gausFit->GetParameter(1), gausFit->GetParameter(2)), "l");
    leg->Draw();

    c->SaveAs((filename + ".png").c_str());
    c->SaveAs((filename + ".pdf").c_str());
    delete c;
}

void beamspot_2DJointMedian() //! TODO
{
    std::string plotdir = "./figure/figure-SvtxMbdVertex";
    system(("mkdir -p " + plotdir).c_str());

    std::string filename = "/sphenix/user/hjheng/sPHENIXRepo/TrackingAnalysis/Silicon_MBD_Vertexing/Silicon_MBD_Comparisons/VertexCompare_run_82405/files/outputVTX_Acts.root"; //! TODO: make it argument for condor job

    std::vector<TVector2> pts;
    std::vector<double> vtxZ;
    TH2 *hM_vertexXY = new TH2F("hM_vertexXY", "Vertex XY positions;vertexX [cm];vertexY [cm]", 200, -1, 1, 200, -1, 1);
    TH1 *hM_vertexZ = new TH1F("hM_vertexZ", "Vertex Z positions;vertexZ [cm];Counts", 60, -30, 30);
    TH1 *hM_vertexZ_TrgCrossing = new TH1F("hM_vertexZ_trgCrossing", "Vertex Z positions for different trigger crossings;vertexZ [cm];Counts", 60, -30, 30);
    TH1 *hM_vertexZ_nonTrgCrossing = new TH1F("hM_vertexZ_nonTrgCrossing", "Vertex Z positions for non-trigger crossings;vertexZ [cm];Counts", 60, -30, 30);

    TFile *file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        std::cout << "Failed to open file: " << filename << std::endl;
        return;
    }

    TTree *tree = static_cast<TTree *>(file->Get("VTX"));
    if (!tree)
    {
        std::cout << "Failed to get TTree VTX from: " << filename << std::endl;
        return;
    }

    int counter;
    std::vector<float> *mbdVertex = 0;
    uint64_t gl1bco, bcotr;
    std::vector<int> *firedTriggers = 0;
    std::vector<float> *trackerVertexX = 0;
    std::vector<float> *trackerVertexY = 0;
    std::vector<float> *trackerVertexZ = 0;
    std::vector<short> *trackerVertexCrossing = 0;
    tree->SetBranchAddress("counter", &counter);
    tree->SetBranchAddress("mbdVertex", &mbdVertex);
    tree->SetBranchAddress("gl1bco", &gl1bco);
    tree->SetBranchAddress("bcotr", &bcotr);
    tree->SetBranchAddress("firedTriggers", &firedTriggers);
    tree->SetBranchAddress("trackerVertexX", &trackerVertexX);
    tree->SetBranchAddress("trackerVertexY", &trackerVertexY);
    tree->SetBranchAddress("trackerVertexZ", &trackerVertexZ);
    tree->SetBranchAddress("trackerVertexCrossing", &trackerVertexCrossing);

    const Long64_t nentries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nentries; ++entry)
    {
        tree->GetEntry(entry);

        size_t n = std::min(trackerVertexX->size(), trackerVertexY->size());
        n = std::min(n, trackerVertexCrossing->size());

        for (size_t i = 0; i < n; ++i)
        {
            // if (trackerVertexCrossing->at(i) != 0)
            // {
            //     continue;
            // }

            pts.emplace_back(trackerVertexX->at(i), trackerVertexY->at(i));
            vtxZ.push_back(trackerVertexZ->at(i));
            hM_vertexXY->Fill(trackerVertexX->at(i), trackerVertexY->at(i));
            hM_vertexZ->Fill(trackerVertexZ->at(i));
            if (trackerVertexCrossing->at(i) == 0)
            {
                hM_vertexZ_TrgCrossing->Fill(trackerVertexZ->at(i));
            }
            else
            {
                hM_vertexZ_nonTrgCrossing->Fill(trackerVertexZ->at(i));
            }
        }
    }

    if (pts.empty())
    {
        std::cout << "No trackerVertex points found with trackerVertexCrossing==0." << std::endl;
        return;
    }

    // Compute coordinate-wise median, mean, and geometric median
    std::vector<double> xs, ys;
    xs.reserve(pts.size());
    ys.reserve(pts.size());
    double mx = 0, my = 0;
    for (const auto &p : pts)
    {
        xs.push_back(p.X());
        ys.push_back(p.Y());
        mx += p.X();
        my += p.Y();
    }
    mx /= pts.size();
    my /= pts.size();

    TVector2 coordMed(median1D(xs), median1D(ys));
    TVector2 geoMed = geometricMedian2D(pts, 1e-10, 2000, 1e-14);

    std::cout << "Mean:                          (" << mx << ", " << my << ")\n";
    std::cout << "Coordinate-wise median:        (" << coordMed.X() << ", " << coordMed.Y() << ")\n";
    std::cout << "Geometric median (2D joint):   (" << geoMed.X() << ", " << geoMed.Y() << ")\n";

    // calculate the effective sigma (1-sigma interval)
    double sigmaXEff = 0.0, sigmaYEff = 0.0, vtxXmin = 0.0, vtxXmax = 0.0, vtxYmin = 0.0, vtxYmax = 0.0;
    sigmaXEff = sigmaEff(xs, 0.6827, vtxXmin, vtxXmax); // 1-sigma interval
    sigmaYEff = sigmaEff(ys, 0.6827, vtxYmin, vtxYmax); // 1-sigma interval
    std::cout << "Effective sigma X (1-sigma): " << sigmaXEff << " (range: " << vtxXmin << " to " << vtxXmax << ")\n";
    std::cout << "Effective sigma Y (1-sigma): " << sigmaYEff << " (range: " << vtxYmin << " to " << vtxYmax << ")\n";

    // calculate the mean, median, and effective sigmal, standard deviation for Z vertex
    double zMean = std::accumulate(vtxZ.begin(), vtxZ.end(), 0.0) / vtxZ.size();
    double zCoordMed = median1D(vtxZ);
    double zStdDev = std::sqrt(std::accumulate(vtxZ.begin(), vtxZ.end(), 0.0, [zMean](double acc, double z) { return acc + (z - zMean) * (z - zMean); }) / vtxZ.size());
    std::cout << "Z Vertex - Mean: " << zMean << ", Coordinate-wise Median: " << zCoordMed << ", Std Dev: " << zStdDev << std::endl;

    TCanvas *c = new TCanvas("c", "c", 800, 700);
    gPad->SetRightMargin(0.15);
    gPad->SetTopMargin(0.07);
    c->cd();
    hM_vertexXY->GetXaxis()->SetTitle("Vertex X [cm]");
    hM_vertexXY->GetYaxis()->SetTitle("Vertex Y [cm]");
    hM_vertexXY->Draw("COLZ");

    // markers: geometric median
    TMarker *mGM = new TMarker(geoMed.X(), geoMed.Y(), 20); // circle
    mGM->SetMarkerSize(1.0);
    mGM->SetMarkerColor(kRed);
    // mGM->Draw("SAME");

    c->Update();

    TLegend *leg = new TLegend(gPad->GetLeftMargin() + 0.02,   //
                               gPad->GetBottomMargin() + 0.03, //
                               gPad->GetLeftMargin() + 0.2,    //
                               gPad->GetBottomMargin() + 0.13  //
    );
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->SetTextAlign(kHAlignLeft + kVAlignBottom);
    // leg->AddEntry(mMean, "Mean", "p");
    // leg->AddEntry(mCM, "Coord-wise Median", "p");
    leg->AddEntry((TObject *)0, Form("Geometric median (X, Y)=(%.3g, %.3g) cm", geoMed.X(), geoMed.Y()), "");
    leg->AddEntry((TObject *)0, Form("Effective #sigma_{X} (1#sigma): %.3g cm, [%.3g, %.3g] cm", sigmaXEff, vtxXmin, vtxXmax), "");
    leg->AddEntry((TObject *)0, Form("Effective #sigma_{Y} (1#sigma): %.3g cm, [%.3g, %.3g] cm", sigmaYEff, vtxYmin, vtxYmax), "");
    leg->Draw("SAME");

    c->SaveAs(Form("%s/beamspot_2DJointMedian.png", plotdir.c_str()));
    c->SaveAs(Form("%s/beamspot_2DJointMedian.pdf", plotdir.c_str()));

    // draw the vertex Z distribution for trigger crossings and non-trigger crossings
    draw1Dhistogram(hM_vertexZ_TrgCrossing,          //
                    false,                           //
                    "Vertex Z [cm]",                 //
                    "Counts",                        //
                    {"Trigger crossings"},           //
                    "hist e1",                       //
                    plotdir + "/vertexZ_trgCrossing" //
    );
    draw1Dhistogram(hM_vertexZ_nonTrgCrossing,          //
                    false,                              //
                    "Vertex Z [cm]",                    //
                    "Counts",                           //
                    {"Non-trigger crossings"},          //
                    "hist e1",                          //
                    plotdir + "/vertexZ_nonTrgCrossing" //
    );

    draw1Dhistogram_GaussianFit(hM_vertexZ,          //
                                false,          //
                                "Vertex Z [cm]", //
                                "Counts",        //
                                {},   //
                                "hist e1",           //
                                plotdir + "/vertexZ_allCrossing" //
    );

    // output tree
    TFile *outFile = TFile::Open(Form("%s/beamspot_2DJointMedian.root", plotdir.c_str()), "RECREATE"); //! TODO: make it argument for condor job
    TTree *outTree = new TTree("vertex", "Vertex Tree with MBD Z vertex and silicon beamspot and per-crossing vertex");
    uint64_t out_gl1bco, out_bcotr;
    float beamspotX = geoMed.X();
    float beamspotY = geoMed.Y();
    float MbdVertexZ;
    bool trigger12; // true if trigger 12 fired (MBD N&S >=2 & |vtx| < 10 cm)
    float SiliconVertexZ_crossing0;
    std::vector<float> SiliconVertexZ_nonzeroCrossing;
    outTree->Branch("counter", &counter);
    outTree->Branch("gl1bco", &out_gl1bco);
    outTree->Branch("bcotr", &out_bcotr);
    outTree->Branch("beamspotX", &beamspotX);
    outTree->Branch("beamspotY", &beamspotY);
    outTree->Branch("MbdVertexZ", &MbdVertexZ);
    outTree->Branch("trigger12", &trigger12);
    outTree->Branch("SiliconVertexZ_crossing0", &SiliconVertexZ_crossing0);
    outTree->Branch("SiliconVertexZ_nonzeroCrossing", &SiliconVertexZ_nonzeroCrossing);
    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
    {
        tree->GetEntry(entry);

        out_gl1bco = gl1bco;
        out_bcotr = bcotr;

        // check if trigger 12 is fired
        trigger12 = false;
        auto it = std::find(firedTriggers->begin(), firedTriggers->end(), 12);
        if (it != firedTriggers->end())
        {
            trigger12 = true;
        }

        // MBD vertex Z
        MbdVertexZ = (mbdVertex->empty()) ? -9999.0 : mbdVertex->at(0); // all the mbdVertex are the same per trigger frame

        // Silicon vertex Z
        SiliconVertexZ_crossing0 = -9999.0;
        SiliconVertexZ_nonzeroCrossing.clear();

        for (size_t i = 0; i < trackerVertexCrossing->size(); ++i)
        {
            if (trackerVertexCrossing->at(i) == 0)
            {
                SiliconVertexZ_crossing0 = trackerVertexZ->at(i); // assuming Z is vertex Z
            }
            else
            {
                SiliconVertexZ_nonzeroCrossing.push_back(trackerVertexZ->at(i)); // assuming Z is vertex Z
            }
        }

        outTree->Fill();
    }
    outFile->cd();
    outTree->Write();
    outFile->Close();

    // clean up
    delete c;
    delete hM_vertexXY;
    file->Close();
}
