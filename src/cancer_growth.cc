#include "./cancer_growth.h"

int main (int argc, char* argv[]) {
    const int n_cell_types =2, n_biochem =1, n_mstruct =1;
    BDM_Domain rve(n_cell_types, n_biochem, n_mstruct);
    // set some initial values for these data
    rve.cells_population[0] = 100;
    rve.cells_population[1] = 10;
    rve.biochemical_level[0] = 1.0;
    rve.biochemical_level_gradient[0] = RealVector3d({0.0, 0.0, 0.0});
    rve.microstructural_level[0] = 0.0;
    rve.microstructural_level_gradient[0] = RealVector3d({0.0, 0.0, 0.0});
    //
    const int n_cell_types2 =4, n_biochem2 =2, n_mstruct2 =2;
    BDM_Domain rve2(n_cell_types, n_biochem, n_mstruct);
    // set some initial values for these data
    rve2.cells_population[0] = 100;
    rve2.cells_population[1] = 10;
    rve2.biochemical_level[0] = 0.2;
    rve2.biochemical_level_gradient[0] = RealVector3d({0.0, 0.0, 0.0});
    rve2.microstructural_level[0] = 0.0;
    rve2.microstructural_level_gradient[0] = RealVector3d({0.0, 0.0, 0.0});

    // list of sinulations to run
    array<BDM_Domain, 2> simulations_list = {rve, rve2};

    for (int i; i < simulation_list.size(); i++) {
      bdm::Initialise(simulation_list[i]);
      bdm::Simulate(simulation_list[i]);
      bdm::ResetBDM();
    }

    return 0;
}
