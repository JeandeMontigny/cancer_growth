#include "./cancer_growth.h"

int main (int argc, char* argv[]) {
    const int n_cell_types =2, n_biochem =1, n_mstruct =1, n_steps = 1000;
    BDM_Domain rve(n_cell_types, n_biochem, n_mstruct);
    // set some initial values for these data
    rve.cells_population[0] = 100;
    rve.cells_population[1] = 10;
    rve.biochemical_level[0] = 1.0;
    rve.biochemical_level_gradient[0] = RealVector3d({0.0, 0.0, 0.0});
    rve.microstructural_level[0] = 0.0;
    rve.microstructural_level_gradient[0] = RealVector3d({0.0, 0.0, 0.0});
    //
    bdm::Initialise(rve);
    bdm::Simulate(rve);
    //
    return 0;
}
