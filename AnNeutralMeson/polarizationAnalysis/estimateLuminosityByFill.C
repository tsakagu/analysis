#include <uspin/SpinDBContent.h>
#include <uspin/SpinDBOutput.h>

R__LOAD_LIBRARY(libuspin.so)

void estimateLuminosityRun(int iB, int iScal, std::stringstream& GL1Pscalers, unsigned long& sumScalers, unsigned long runnumber, 
                            unsigned long& LuminosityUp, unsigned long& LuminosityDown,
                            double& pol, double& polerr, double& crossangle, double& crossangle_std);

void estimateLuminosityByFill(const std::string& runListName = "run_to_fill.txt", const std::string& outputfilename = "luminosity_fill.root")
{
  std::ifstream fin;
  fin.open(runListName);
  if (!(fin.is_open()))
  {
    std::cerr << "Error: File " << runListName << " could not be opened.\n";
  }
  std::string line;

  // Polarizations and luminosities are defined differently for each beam and scaler
  const int nBeams = 2;
  const std::string beams[2] = {"blue", "yellow"};
  const int nScalers = 2;
  const std::string scalers[2] = {"mbdns", "zdcns"};

  // Book histograms
  TFile *outfile = new TFile(outputfilename.c_str(), "RECREATE");
  outfile->cd();

  // Polarization graphs
  TH1F *h_polarization[nBeams][nScalers];
  TH1F *h_relative_luminosity[nBeams][nScalers];
  TGraph *graph_polarization_fill[nBeams][nScalers];
  TGraph *graph_relative_luminosity_fill[nBeams][nScalers];
  TGraph *graph_luminosity_up_fill[nBeams][nScalers];
  TGraph *graph_luminosity_down_fill[nBeams][nScalers];
  TGraph *graph_crossangle_fill[nBeams][nScalers];
  int runnumber;
  int fillnumber;
  double Rmbd;
  double Rzdc;
  TTree *tree_relative_luminosity[nBeams];
  for (int iB = 0; iB < nBeams; iB++)
  {
    std::stringstream treename; treename << "tree_relative_luminosity_" << beams[iB]; 
    tree_relative_luminosity[iB] = new TTree(treename.str().c_str(), "relative luminosity tree");
    tree_relative_luminosity[iB]->Branch("fillnumber", &fillnumber, "fillnumber/I"); 
    tree_relative_luminosity[iB]->Branch("mbd", &Rmbd, "mbd/D");
    tree_relative_luminosity[iB]->Branch("zdc", &Rzdc, "zdc/D");
    for (int iScal = 0; iScal < nScalers; iScal++) 
    {
      std::stringstream hname; hname << "h_polarization_" << beams[iB] << "_" << scalers[iScal];
      std::stringstream htitle; htitle << ";" << beams[iB] << " (" << scalers[iScal] << ") " << " beam polarization [%]; luminosity weighted distribution";
      h_polarization[iB][iScal] = new TH1F(hname.str().c_str(), htitle.str().c_str(), 50, 30, 65);
      hname.str(""); hname << "h_relative_luminosity_" << beams[iB] << "_" << scalers[iScal];
      htitle.str(""); htitle << ";" << beams[iB] << " (" << scalers[iScal] << ") " << " beam relative luminosity; luminosity weighted distribution";
      h_relative_luminosity[iB][iScal] = new TH1F(hname.str().c_str(), htitle.str().c_str(), 50, 0.9, 1.1);
      hname.str(""); hname << "graph_polarization_fill_" << beams[iB] << "_" << scalers[iScal];
      htitle.str(""); htitle << ";fill number;" << beams[iB] << " (" << scalers[iScal] << ") " << " beam polarization [%]";
      graph_polarization_fill[iB][iScal] = new TGraph();
      graph_polarization_fill[iB][iScal]->SetName(hname.str().c_str());
      graph_polarization_fill[iB][iScal]->SetTitle(htitle.str().c_str());
      hname.str(""); hname << "graph_relative_luminosity_fill_" << beams[iB] << "_" << scalers[iScal];
      htitle.str(""); htitle << ";fill number;" << beams[iB] << " (" << scalers[iScal] << ") " << " beam relative luminosity";
      graph_relative_luminosity_fill[iB][iScal] = new TGraph();
      graph_relative_luminosity_fill[iB][iScal]->SetName(hname.str().c_str());
      graph_relative_luminosity_fill[iB][iScal]->SetTitle(htitle.str().c_str());
      hname.str(""); hname << "graph_luminosity_up_fill_" << beams[iB] << "_" << scalers[iScal];
      htitle.str(""); htitle << ";fill number;" << beams[iB] << " (" << scalers[iScal] << ") " << " beam luminosity up";
      graph_luminosity_up_fill[iB][iScal] = new TGraph();
      graph_luminosity_up_fill[iB][iScal]->SetName(hname.str().c_str());
      graph_luminosity_up_fill[iB][iScal]->SetTitle(htitle.str().c_str());
      hname.str(""); hname << "graph_luminosity_down_fill_" << beams[iB] << "_" << scalers[iScal];
      htitle.str(""); htitle << ";fill number;" << beams[iB] << " (" << scalers[iScal] << ") " << " beam luminosity down";
      graph_luminosity_down_fill[iB][iScal] = new TGraph();
      graph_luminosity_down_fill[iB][iScal]->SetName(hname.str().c_str());
      graph_luminosity_down_fill[iB][iScal]->SetTitle(htitle.str().c_str());

      hname.str(""); hname << "graph_crossangle_fill_" << beams[iB] << "_" << scalers[iScal];
      htitle.str(""); htitle << ";fill number;" << beams[iB] << " (" << scalers[iScal] << ") " << " beam luminosity down";
      graph_crossangle_fill[iB][iScal] = new TGraph();
      graph_crossangle_fill[iB][iScal]->SetName(hname.str().c_str());
      graph_crossangle_fill[iB][iScal]->SetTitle(htitle.str().c_str());
    }
  }
  
  unsigned long normalization[nBeams][nScalers] = {0};
  long double average_polarization[nBeams][nScalers] = {0};
  long double total_luminosity_up[nBeams][nScalers] = {0};
  long double total_luminosity_down[nBeams][nScalers] = {0};
  long double average_relative_luminosity[nBeams][nScalers] = {0};
  
  //std::cout << "runnumber beam scaler LuminosityUp LuminosityDown RelativeLuminosity" << std::endl;
  std::string entry_content;
  // Read all fill and run numbers
  std::map<int, std::vector<int>> fill2runs;
  while (std::getline(fin, line)) // Run number loop
  {
    std::istringstream entry(line);
    std::getline(entry, entry_content, ' ');
    fillnumber = std::stoi(entry_content);
    std::getline(entry, entry_content, ' ');
    runnumber = std::stoi(entry_content);
    fill2runs[fillnumber].push_back(runnumber);
  }
  
  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iScal = 0; iScal < nScalers; iScal++)
    {
      int count = -1;
      for (const auto& [fill, runs] : fill2runs)
      {
        count++;
        fillnumber = fill;
        std::cout << "fill number = " << fillnumber << std::endl;
        unsigned long LuminosityUp = 0;
        unsigned long LuminosityDown = 0;
        double pol = 0;
        double polerr = 0;
        double crossangle;
        double crossangle_std;
        std::stringstream GL1Pscalers;
        unsigned long sumScalers = 0;
        for (int irun = 0; irun < (int) runs.size(); irun++)
        {
          runnumber = runs[irun];
          estimateLuminosityRun(iB, iScal, GL1Pscalers, sumScalers,
                                runnumber, LuminosityUp, LuminosityDown,
                                pol, polerr, crossangle, crossangle_std);
        }
        
        unsigned long Luminosity = LuminosityUp + LuminosityDown;
        long double RelativeLuminosity = (long double) LuminosityUp / (long double) LuminosityDown;
        if (iScal == 0) Rmbd = RelativeLuminosity;
        else if (iScal == 1) Rzdc = RelativeLuminosity;
        

        average_polarization[iB][iScal] += Luminosity * pol;
        total_luminosity_up[iB][iScal] += (long double) LuminosityUp;
        total_luminosity_down[iB][iScal] += (long double) LuminosityDown;

        //average_relative_luminosity_blue += (float) blueLuminosity * ((float) blueUpLuminosity / (float) blueDownLuminosity);
        normalization[iB][iScal] += Luminosity;

        //std::cout << runnumber << " " << beams[iB] << " " << scalers[iScal] << " " << Luminosity << std::endl;
        //std::cout << runnumber << " " << beams[iB] << " " << scalers[iScal] << " " << LuminosityUp << " " << LuminosityDown << " " << RelativeLuminosity << std::endl;
        
        // fill histograms
        h_polarization[iB][iScal]->Fill(pol * 100, Luminosity);
        graph_polarization_fill[iB][iScal]->SetPoint(count, fillnumber, pol);
        h_relative_luminosity[iB][iScal]->Fill(RelativeLuminosity, Luminosity);
        graph_relative_luminosity_fill[iB][iScal]->SetPoint(count, fillnumber, RelativeLuminosity);
        graph_luminosity_up_fill[iB][iScal]->SetPoint(count, fillnumber, LuminosityUp);
        graph_luminosity_down_fill[iB][iScal]->SetPoint(count, fillnumber, LuminosityDown);
        graph_crossangle_fill[iB][iScal]->SetPoint(count, fillnumber, crossangle);
        //graph_crossangle_run[iB][iScal]->SetPointError(count, 0, crossangle_std);
      }
      //std::cout << runnumber << " " << beams
      tree_relative_luminosity[iB]->Fill();
    }
  }
  for (int iB = 0; iB < nBeams; iB++)
  {
    for (int iScal = 0; iScal < nScalers; iScal++)
    {
      average_polarization[iB][iScal] /= (long double) normalization[iB][iScal];
      average_relative_luminosity[iB][iScal] = total_luminosity_up[iB][iScal] / total_luminosity_down[iB][iScal];
      
      std::cout << std::fixed << std::setprecision(10)
                << "Average values " << beams[iB] << " " << scalers[iScal] << ":\n"
                << "Pol = " << average_polarization[iB][iScal] << "  "
                << "RL = " << average_relative_luminosity[iB][iScal] << "\n";
      
      graph_polarization_fill[iB][iScal]->Write();
      graph_relative_luminosity_fill[iB][iScal]->Write();
      graph_luminosity_up_fill[iB][iScal]->Write();
      graph_luminosity_down_fill[iB][iScal]->Write();
      graph_crossangle_fill[iB][iScal]->Write();
      fin.close();
    }
  }
  
  outfile->Write();
  outfile->Close();
  delete outfile;

  gSystem->Exit(0);

}

void estimateLuminosityRun(int iB, int iScal, std::stringstream& GL1Pscalers, unsigned long& sumScalers, unsigned long runnumber, 
                              unsigned long& LuminosityUp, unsigned long& LuminosityDown,
                              double& pol, double& polerr, double& crossangle, double& crossangle_std)
{
  SpinDBOutput spin_out("phnxrc");
  SpinDBContent *spin_cont = new SpinDBContentv1();

  spin_out.StoreDBContent(runnumber, runnumber);
  spin_out.GetDBContentStore(spin_cont, runnumber);

  // General information
  //std::cout << "Run number: " << runnumber << std::endl;

  // Spin patterns
  int spinpat[120] = {0};
  if (iB == 0)
  {
    for (int i = 0; i < 120; ++i)
    {
      spinpat[i] = -1 * spin_cont->GetSpinPatternBlue(i);
    }
  }
  else if (iB == 1)
  {
    for (int i = 0; i < 120; ++i)
    {
      spinpat[i] = -1 * spin_cont->GetSpinPatternYellow(i);
    }
  }

  crossangle = spin_cont->GetCrossAngle();
  crossangle_std = spin_cont->GetCrossAngleStd();

  // Beam polarization
  //double bluepol, yellpol, bluepolerr, yellpolerr;
  if (iB == 0) spin_cont->GetPolarizationBlue(0, pol, polerr);
  else if (iB == 1) spin_cont->GetPolarizationYellow(0, pol, polerr);
  pol /= 100; polerr /= 100;

  // GL1p Scalers
  // I don't know which one is the most relevant
  unsigned long scaler[120] = {0};
  //GL1Pscalers << "\"";
  sumScalers = 0;
  if (iScal == 0)
  {
    for (int i = 0; i < 120; ++i)
    {
      scaler[i] = spin_cont->GetScalerMbdNoCut(i); sumScalers += scaler[i];
    }
  } else if (iScal == 1)
  {
    for (int i = 0; i < 120; ++i)
    {
      scaler[i] = spin_cont->GetScalerZdcNoCut(i); sumScalers += scaler[i];
    }
  } else if (iScal == 2)
  {
    std::cout << "This should never happen!\n";
    for (int i = 0; i < 120; ++i)
    {
      scaler[i] = spin_cont->GetScalerMbdVertexCut(i); sumScalers += scaler[i];
    }
  }

  // Estimate luminosity:
  // LuminosityUp = 0;
  // LuminosityDown = 0;
  for (int i = 0; i < 120; i++)
  {
    if (spinpat[i] == 1)
      LuminosityUp += scaler[i];
    else if (spinpat[i] == -1)
      LuminosityDown += scaler[i];
  }
}

