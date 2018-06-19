#include "cancer_growth.h"

int main (int argc, char* argv[]) {
  bdm::InitializeBiodynamo("cancer_growth");

  const int desiredNumberOfBdmSimulations = 1;
  std::array<BDM_Domain, desiredNumberOfBdmSimulations> simulations_list;

  for (int i=0; i<desiredNumberOfBdmSimulations; i++) {
    // TODO: build simulation depending on feb3 info at time 0
    const int n_cell_types =2, n_biochem =1, n_mstruct =1;
    BDM_Domain rve(i, n_cell_types, n_biochem, n_mstruct);
    // set some initial values for these data
    rve.cells_population[0] = 1; // 1000
    rve.cells_population[1] = 1; // 10
    rve.biochemical_level[0] = 1;
    rve.biochemical_level_gradient[0] = RealVector3d({0.0, 0.0, 0.0});
    rve.microstructural_level[0] = 0.0;
    rve.microstructural_level_gradient[0] = RealVector3d({0.0, 0.0, 0.0});
    for (int face = 0; face < 6; face++) {
      rve.face_incomingCells[face]=0;
    }
    // copy data to array structure
    simulations_list[i] = rve;
  }

  const unsigned int numberOfIterations = 1;
  for (unsigned int j=0; j<numberOfIterations; j++) {
    std::cout << "\n\t-- Iteration " << j << " --" << std::endl;
    for (unsigned int i=0; i<simulations_list.size(); i++) {
      //TODO: if need to initialise at each simulation
      // const char *visualisationFileName = "cancer_growth_";
      // char initialisationName [100];
      // strcpy(initialisationName,visualisationFileName);
      // strcat(initialisationName, std::to_string(j).c_str());
      // strcat(initialisationName, "_");
      // strcat(initialisationName, std::to_string(i).c_str());
      // bdm::InitializeBiodynamo(initialisationName);

      // if first iteration, manage folder for visualisation export
      if (j==0) {
        struct stat sb;
        std::string tempsName = "rve_"+std::to_string(simulations_list[i].id())+"_visualisation";
        const char *folderName = tempsName.c_str();
        // if folder exists, remove previous visualisation files
        if (stat(folderName, &sb) == 0 && S_ISDIR(sb.st_mode)) {
          const char *cmd2 = "rm ./";
          const char *cmd3 = "/*";
          char concatenation2 [100];
          strcpy(concatenation2,cmd2);
          strcat(concatenation2, folderName); strcat(concatenation2, cmd3);
          std::system(concatenation2);
        }
        // create folder for visualisatione export
        else {
          const char *cmd = "mkdir ";
          char concatenation [100];
          strcpy(concatenation,cmd); strcat(concatenation, folderName);
          std::system(concatenation);
        }
      }
      bdm::Initialise(simulations_list[i]);
      bdm::Simulate(simulations_list[i]);
      int total = 0;
      for (int k=0; k<6; k++) {
        total = total + simulations_list[i].escaped_cells[k];
      }
      std::cout << "number of escaped_cells: " << total << std::endl;

      // move visualisation files into separated folder for each simulation
      std::string tempsName = "rve_"+std::to_string(simulations_list[i].id())+"_visualisation";
      const char *folderName = tempsName.c_str();
      const char *cmd = "mv -b cancer_growth.pvsm simulation_info.json MyCell_data_* ";
//      const char *cmd = "mv -b MyCell_data_* ";
      char concatenation [100];
      strcpy(concatenation, cmd); strcat(concatenation, folderName);
      std::system(concatenation);

      bdm::ResetBDM();
    }
    // TODO: run feb3 and update BDM_Domain informations
  }

  return 0;
}
