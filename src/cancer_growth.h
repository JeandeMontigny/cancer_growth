#ifndef CANCER_GROWTH_H_
#define CANCER_GROWTH_H_

#include "biodynamo.h"
#include "../bio_module.h"

namespace bdm {

  // 0. Define my custom cell, which extends Cell by adding an extra
// data member cell_type.
  //oxygen_level should be between 0 and 1
  BDM_SIM_OBJECT(MyCell, Cell) {
    BDM_SIM_OBJECT_HEADER(MyCellExt, 1, can_divide_, cell_colour_, oxygen_level_, hypo_division_); // create the header with our new data member

  public:
    MyCellExt() {}
    MyCellExt(const std::array<double, 3>& position) : Base(position) {} // our creator

    // getter and setter for our new data member
    void SetCanDivide(bool d) { can_divide_[kIdx] = d; }
    bool GetCanDivide() { return can_divide_[kIdx]; }
    bool* GetCanDividePtr() { return can_divide_.data(); }

    void SetCellColour(int cellColour) { cell_colour_[kIdx] = cellColour; }
    int GetCellColour() { return cell_colour_[kIdx]; }
    int* GetCellColourPtr() { return cell_colour_.data(); }

    void SetOxygenLevel(double oxygenLevel) {oxygen_level_[kIdx] = oxygenLevel; }
    double GetOxygenLevel() { return oxygen_level_[kIdx]; }
    double* GetOxygenLevelPtr() { return oxygen_level_.data(); }

    void SetHypoDiv(bool div) { hypo_division_[kIdx] = div; }
    bool GetHypoDiv() { return hypo_division_[kIdx]; }
    bool* GetHypoDivPtr() { return hypo_division_.data(); }

  private:
    // declare new data member and define their type
    // private data can only be accessed by public function and not directly
    vec<bool> can_divide_;
    vec<int> cell_colour_;
    vec<double> oxygen_level_;
    vec<bool> hypo_division_;
  };

  
// 2. Define compile time parameter
  template <typename Backend>
    struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
    using BiologyModules = Variant<GrowthModule>;
    using AtomicTypes = VariadicTypedef<MyCell>; // use my custom Cell
  };
  
  template <typename Function, typename TResourceManager = ResourceManager<>>  
// my cell creator
    static void CellCreator(double min, double max, int num_cells, Function cell_builder) {
    auto rm = TResourceManager::Get();

    // Determine simulation object type which is returned by the cell_builder
    using FunctionReturnType = decltype(cell_builder({0, 0, 0}));

    auto container = rm->template Get<FunctionReturnType>();
    container->reserve(num_cells);

    for (int i = 0; i < num_cells; i++) {
      double x = gTRandom.Uniform(min, max);
      double y = gTRandom.Uniform(min, max);
      double z = gTRandom.Uniform(min, max);
      auto new_simulation_object = cell_builder({x, y, z});
      container->push_back(new_simulation_object);
    }
    container->Commit();
  } // end CellCreator


  inline
  double get3DDistSq (array<double, 3> cell1, array<double, 3> cell2) {
    return pow((cell1[0]-cell2[0]),2)+pow((cell1[1]-cell2[1]),2)+pow((cell1[2]-cell2[2]), 2);
  }
// 3.
  template <typename TResourceManager = ResourceManager<>>
  inline
  int Simulate (int argc, const char** argv) {

    size_t nb_of_cells=100; // 10000 // 1000000
    size_t nb_of_cancerous_cells=10; //10
      
    Param::live_visualization_ = false;
    Param::export_visualization_ = true;
    Param::visualization_export_interval_ = 8;
    Param::visualize_sim_objects_["MyCell"] = std::set<std::string>{"diameter_", "cell_colour_", "oxygen_level_"};
   
    // Create an artificial bounds for the simulation space
    Param::bound_space_ = true;
    Param::min_bound_ = 0;
    Param::max_bound_ = 100.0;//1000
    Param::run_mechanical_interactions_ = true;
    gTRandom.SetSeed(2448); // 5807

    auto Construct_regular =  [](const std::array<double, 3>& position) {
      MyCell cell(position);
      cell.SetDiameter(5.0);
      cell.SetCellColour(0);
      // cell.AddBiologyModule(GrowthModule());
      cell.SetCanDivide(false);
      cell.SetOxygenLevel(1.0);
      // cell.SetHypoDiv(true);
      return cell;
    };

    // cell creation
    CellCreator(Param::min_bound_, Param::max_bound_, nb_of_cells, Construct_regular);
    cout << "regular cells created" << endl;

    auto Construct_cancerous =  [](const std::array<double, 3>& position) {
      MyCell cell(position);
      cell.SetDiameter(6.0);
      cell.SetCellColour(1);
      cell.AddBiologyModule(GrowthModule());
      cell.SetCanDivide(true);
      cell.SetOxygenLevel(1.0);
      cell.SetHypoDiv(true);
      return cell;
    };

    // cell creation
    CellCreator(45.0, 55.0, nb_of_cancerous_cells, Construct_cancerous); // 250, 350
    cout << "cancerous cells created" << endl;

    Scheduler<> scheduler;
    const int max_step = 1500;
    // auto rm = TResourceManager::Get();
    
    for (int i=0; i<max_step; i++) {
      if (i%10==0) {
        cout << "step " << i << " out of " << max_step << endl;
      }      
      scheduler.Simulate(1);
    }

    // auto my_cells=rm->template Get<MyCell>();
    // for (size_t i = 0; i < my_cells->size(); i++) {
    //   auto&& cell = (*my_cells)[i];
    //   if (cell.GetCellColour()==2) {
    //     cell.SetOxygenLevel(0.5);
    //   }
    // }
    
    // for (int i=0; i<max_step; i++) {
    //   if (i%10==0) {
    //     cout << "step " << i << " out of " << max_step << endl;
    //   }      
    //   scheduler.Simulate(1);
    // }

    // rm = TResourceManager::Get();
    // my_cells=rm->template Get<MyCell>();
    // for (size_t i = 0; i < my_cells->size(); i++) {
    //   auto&& cell = (*my_cells)[i];
    //   if (cell.GetCellColour()==2) {
    //     cell.SetOxygenLevel(0.1);
    //   }
    // }
    
    // for (int i=0; i<max_step; i++) {
    //   if (i%10==0) {
    //     cout << "step " << i << " out of " << max_step << endl;
    //   }      
    //   scheduler.Simulate(1);
    // }
    
    // int cancerousCellNb=0;
    // rm = TResourceManager::Get();
    // my_cells=rm->template Get<MyCell>();
    // for (size_t i = 0; i < my_cells->size(); i++) {
    //   auto&& cell = (*my_cells)[i];
    //   if (cell.GetCellColour()==2) {
    //     cancerousCellNb++;
    //   }
    // }
    // cout << "nb of cancerous cells in the modelling: " << cancerousCellNb << endl;
    // cout << "done" << endl;

    return 0;
    
  }

} // namespace bdm

#endif // CANCER_GROWTH_H_
