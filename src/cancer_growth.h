#ifndef CANCER_GROWTH_H_
#define CANCER_GROWTH_H_

#include "biodynamo.h"

namespace bdm {


  inline double getRandom(){
    double random_variable =double(rand())/RAND_MAX;
    return random_variable;
  }

  inline double getRandom_m1_1(){
    double random_variable =(double(rand())/RAND_MAX)*2;
    return random_variable-1;
  }

  // 0. Define my custom cell, which extends Cell by adding an extra
// data member cell_type.
  BDM_SIM_OBJECT(MyCell, Cell) {
    BDM_SIM_OBJECT_HEADER(MyCellExt, 1, can_divide_);

  public:
    MyCellExt() {}
  MyCellExt(const std::array<double, 3>& position) : Base(position) {}

    void SetCanDivide(bool d) { can_divide_[kIdx] = d; }

    bool GetCanDivide() { return can_divide_[kIdx]; }
    
    bool* GetCanDividePtr() { return can_divide_.data(); }
    
  private:
    vec<bool> can_divide_;
  };

  
// 1. Define growth behaviour
  struct GrowthModule {
    template <typename T>
    void Run(T* cell) {
//cell grow until it reach a diam of 8
      if (cell->GetDiameter() < 8) {
        cell->ChangeVolume(1000);

        array<double, 3> cell_movements{getRandom_m1_1()*4, getRandom_m1_1()*4, getRandom_m1_1()*4};
        cell->UpdateMassLocation(cell_movements);
        cell->SetPosition(cell->GetMassLocation());
        //Reset biological movement to 0.
        cell->SetTractorForce({0, 0, 0});
      }
// when diam of 8, it has a chance of dividing 
      else {
//        cout << "can divide: " << cell->GetCanDivide() << endl;
        double aNewRandomDouble=getRandom();
//        cout << "chance of: " << aNewRandomDouble  << endl;
      
        if (aNewRandomDouble <= 0.55 && cell->GetCanDivide()==true ){ //0.55
          Divide(*cell);
//          cell->SetCanDivide(false);
          cout << "\tcell divided" << endl;
        }
// if it doesn't divide, it will never be able to divide
        else {
          if (cell->GetCanDivide()==true) {
            cell->SetCanDivide(false);
            cout << "cell can't divide anymore" << endl;
          }
        }
      }
    }
    
    bool IsCopied(BmEvent event) const { return true; }
    ClassDefNV(GrowthModule, 1);
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


  inline double get3DDistSq(array<double, 3> cell1, array<double, 3> cell2) {
    return pow((cell1[0]-cell2[0]),2)+pow((cell1[1]-cell2[1]),2)+pow((cell1[2]-cell2[2]), 2);
  }
// 3.
//  template <typename TResourceManager = ResourceManager<>>
  inline int Simulate(int argc, const char** argv) {

    size_t nb_of_cells=10000;
    size_t nb_of_cancerous_cells=10;
      
    Param::live_visualization_ = true;
   
    // Create an artificial bounds for the simulation space
    Param::bound_space_ = true;
    Param::min_bound_ = 0;

    Param::max_bound_ = 600;//-20 //120 //520
    Param::run_mechanical_interactions_ = true;
    gTRandom.SetSeed(5807);

// cell creation
    auto Construct_regular =  [](const std::array<double, 3>& position) {
      MyCell cell(position);
      cell.SetDiameter(9);
      cell.AddBiologyModule(GrowthModule());
      cell.SetCanDivide(false);
      return cell;
    };
    CellCreator(Param::min_bound_, Param::max_bound_, nb_of_cells, Construct_regular);
    cout << "regular cells created" << endl;

    auto Construct_cancerous =  [](const std::array<double, 3>& position) {
      MyCell cell(position);
      cell.SetDiameter(7);
      cell.AddBiologyModule(GrowthModule());
      cell.SetCanDivide(true);
      return cell;
    };
    CellCreator(250, 350, nb_of_cancerous_cells, Construct_cancerous);
    cout << "cancerous cells created" << endl;

    Scheduler<> scheduler;
    int max_step=2000;
    
    for (int i=0; i<max_step; i++) {
      if (i%10==0) {
        cout << "step " << i << " out of " << max_step << endl;
      }      
/*for (size_t i = 0; i < cells->size(); i++) {
  auto&& cell = (*cells)[i];
  }*/
      scheduler.Simulate(1);
    }
    cout << "done" << endl;
    return 0;
    
//    Exporter exporter;
// exporter.CreatePVDFile("paraviewExport", 1, 1.0);
//    exporter.ToFile(cells, "Positions.dat");
//    exporter.ToVTUFile(cells, "viewExportVTU", 1);
  }

} // namespace bdm

#endif // CANCER_GROWTH_H_
