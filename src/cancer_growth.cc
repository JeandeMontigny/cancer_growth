#include "./cancer_growth.h"

int main (int argc, char* argv[]) {

  const int desiredNumberOfBdmSimulations = 4;
  std::array<BDM_Domain, desiredNumberOfBdmSimulations> simulations_list;

  for (int i=0; i<desiredNumberOfBdmSimulations; i++) {
    // TODO: build simulation depending on feb3 info at time 0
    const int n_cell_types =2, n_biochem =1, n_mstruct =1;
    BDM_Domain rve(i, n_cell_types, n_biochem, n_mstruct);
    // set some initial values for these data
    rve.cells_population[0] = 100;
    rve.cells_population[1] = 10;
    rve.biochemical_level[0] = 1.0;
    rve.biochemical_level_gradient[0] = RealVector3d({0.0, 0.0, 0.0});
    rve.microstructural_level[0] = 0.0;
    rve.microstructural_level_gradient[0] = RealVector3d({0.0, 0.0, 0.0});
    // copy data to array structure
    simulations_list[i] = rve;
  }

  const unsigned int numberOfIterations = 5;
  for (unsigned int j=0; j<numberOfIterations; j++) {
    std::cout << "\n\t-- Iteration " << j << " --" << std::endl;
    for (unsigned int i=0; i<simulations_list.size(); i++) {
      bdm::Initialise(simulations_list[i]);
      bdm::Simulate(simulations_list[i]);
      bdm::ResetBDM();
    }
    // TODO: run feb3 and update BDM_Domain informations
  }

  return 0;
}
