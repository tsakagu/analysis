#include <uspin/SpinDBInput.h>
#include <uspin/SpinDBOutput.h>
#include <uspin/SpinDBContent.h>

R__LOAD_LIBRARY(libuspin.so)

bool hasGoodSpinQA(int runnumber);

void macro_select_run_spinqa(const std::string& runListName = "RunList_golden.txt")
{
  std::ifstream runListFile;
  runListFile.open(runListName);

  std::ofstream runListFileOutput;
  runListFileOutput.open("RunList_spinqa.txt");

  std::string line;
  int runnumber;
  while (std::getline(runListFile, line))
  {
    runnumber = std::stoi(line);
    if (hasGoodSpinQA(runnumber))
    {
     runListFileOutput << runnumber << std::endl;
    }
  }
  runListFileOutput.close();
  runListFile.close();

  gSystem->Exit(0);
}

bool hasGoodSpinQA(int runnumber)
{
  SpinDBOutput spin_out("phnxrc");
  SpinDBContent *spin_cont = new SpinDBContentv1();
  spin_out.StoreDBContent(runnumber, runnumber);
  if (!spin_out.GetDBContentStore(spin_cont, runnumber))
  {
    // spin QA row is not available
    return false;
  }
  
  // If badrunqa = 1, the run is discarded
  // It corresponds to either of these conditions:
  // spin pattern must match one of the known patterns given by Main Control Room
  // polarization values for either beam must be between 0 and 1
  // < 10 bunches exhibit a mismatch between intended and measured beam spin pattern
  int badrunqa = spin_cont->GetBadRunFlag();

  if (badrunqa == 1)
  {
    return false;
  }

  // In addition, the GL1P scalers must not be uniformly equal to 0
  const int nbunches = 120;
  for (int i = 0; i < nbunches; i++)
  {
    if (spin_cont->GetScalerMbdVertexCut(i) > 0)
    {
      return true;
    }
  }
  return false;
}
