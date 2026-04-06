#include <string>

#include <phool/recoConsts.h>
#include <ffamodules/CDBInterface.h>

R__LOAD_LIBRARY(libcdbobjects)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4all.so)

void macro_select_run_hotmap(const std::string& runListName = "RunList_fem.txt")
{
  recoConsts *rc = recoConsts::instance();
  std::ifstream runListFile;
  runListFile.open(runListName);

  std::ofstream runListFileOutput;
  runListFileOutput.open("RunList_hotmap.txt");

  std::string line;
  int runnumber;
  std::string calibdir;
  while (std::getline(runListFile, line))
  {
    runnumber = std::stoi(line);
    rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
    rc->set_uint64Flag("TIMESTAMP", runnumber);

    calibdir = CDBInterface::instance()->getUrl("CEMC_BadTowerMap");

    if (!calibdir.empty())
    {
      runListFileOutput << runnumber << std::endl;
    }
  }
  runListFileOutput.close();
  runListFile.close();

  gSystem->Exit(0);
}
