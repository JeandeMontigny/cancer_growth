#ifndef CANCER_GROWTH_H_
#define CANCER_GROWTH_H_

#include <fstream>
#include "biodynamo.h"
#include "../bio_module.h"
#include "../bdm_feb3_interface.h"

namespace bdm {



  // 0. Define my custom cell, which extends Cell by adding an extra data member cell_type. 
  BDM_SIM_OBJECT(MyCell, Cell) {
    BDM_SIM_OBJECT_HEADER(MyCellExt, 1, can_divide_, cell_colour_, oxygen_level_, hypo_division_); // create the header with our new data member

  public:
    MyCellExt() {}
    MyCellExt(const std::array<double, 3>& position) : Base(position) {} // our creator

    // getter and setter for our new data member
    void  SetCanDivide(bool d) { can_divide_[kIdx] = d; }
    bool  GetCanDivide() { return can_divide_[kIdx]; }
    bool* GetCanDividePtr() { return can_divide_.data(); }

    void SetCellColour(int cellColour) { cell_colour_[kIdx] = cellColour; }
    int  GetCellColour() { return cell_colour_[kIdx]; }
    int* GetCellColourPtr() { return cell_colour_.data(); }

    void    SetOxygenLevel(double oxygenLevel) {oxygen_level_[kIdx] = oxygenLevel; }
    double  GetOxygenLevel() { return oxygen_level_[kIdx]; }
    double* GetOxygenLevelPtr() { return oxygen_level_.data(); }

    void  SetHypoDiv(bool div) { hypo_division_[kIdx] = div; }
    bool  GetHypoDiv() { return hypo_division_[kIdx]; }
    bool* GetHypoDivPtr() { return hypo_division_.data(); }

  private:
    vec<bool> can_divide_;
    vec<int> cell_colour_;
    vec<double> oxygen_level_; // should be between 0 and 1
    vec<bool> hypo_division_;
  };


  
  // 2. Define the compile-time parameter
  template <typename Backend>
  struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
    using BiologyModules = Variant<GrowthModule>;
    using AtomicTypes = VariadicTypedef<MyCell>; // use my custom Cell
  };
  


  // my cell creator
  template <typename Function, typename TResourceManager = ResourceManager<>>  
  static
  void CellCreator (double min, double max, unsigned int n_cells, Function cell_builder) {
    auto rm = TResourceManager::Get();

    // Determine simulation object type which is returned by the cell_builder
    using FunctionReturnType = decltype(cell_builder({0, 0, 0}));

    auto container = rm->template Get<FunctionReturnType>();
    container->reserve(n_cells);

    for (int i=0; i<n_cells; i++) {
      double x = gTRandom.Uniform(min, max);
      double y = gTRandom.Uniform(min, max);
      double z = gTRandom.Uniform(min, max);
      auto new_simulation_object = cell_builder({x, y, z});
      container->push_back(new_simulation_object);
    }
    container->Commit();
  }



  inline
  double pow2 (double a) { return a*a; }



  inline
  double get3DDistSq (std::array<double, 3> cell1, std::array<double, 3> cell2) {
    return pow2((cell1[0]-cell2[0]))
          +pow2((cell1[1]-cell2[1]))
          +pow2((cell1[2]-cell2[2]));
  }



  // 3. Core initialisation routine
  template <typename TResourceManager = ResourceManager<>>
  inline
  void Initialise () {
    const unsigned int n_host_cells = 100;
    const unsigned int n_cancer_cells = 10;

    Param::live_visualization_ = false;
    Param::export_visualization_ = true;
    Param::visualization_export_interval_ = 20;
    Param::visualize_sim_objects_["MyCell"] = std::set<std::string>{"diameter_", "cell_colour_", "oxygen_level_"};
    Param::bound_space_ = true; // create artificial boundary limits for the 3D simulation space
    Param::min_bound_ = 0;
    Param::max_bound_ = 100.0;
    Param::run_mechanical_interactions_ = true;
    gTRandom.SetSeed(2448);

    auto Construct_Host_Cells =  [](const std::array<double, 3>& position) {
      MyCell cell(position);
      cell.SetDiameter(2.0);
      cell.SetCellColour(0);
      // cell.AddBiologyModule(GrowthModule());
      cell.SetCanDivide(false);
      cell.SetOxygenLevel(1.0);
      // cell.SetHypoDiv(true);
      return cell;
    };

    // cell creation (min boundary, max boundary, # of cells, default initialiser for cells)
    CellCreator(0.01, 99.99, n_host_cells, Construct_Host_Cells);

    auto Construct_Cancer_Cells =  [](const std::array<double, 3>& position) {
      MyCell cell(position);
      cell.SetDiameter(6.0);
      cell.SetCellColour(1);
      cell.AddBiologyModule(GrowthModule());
      cell.SetCanDivide(true);
      cell.SetOxygenLevel(1.0);
      cell.SetHypoDiv(true);
      return cell;
    };

    // cell creation (min boundary, max boundary, # of cells, default initialiser for cells)
    CellCreator(45.00, 55.00, n_cancer_cells, Construct_Cancer_Cells);

    cout << "regular cells created = " << n_host_cells << endl;
    cout << "cancerous cells created = " << n_cancer_cells << endl;
  }



  // 4. Core simulation routine
  template <typename TResourceManager = ResourceManager<>>
  inline
  void Simulate () {

    Scheduler<> scheduler;
    const int max_step = 1000;
    // auto rm = TResourceManager::Get();
    
    // create a PVD file for Paraview to process
    {
      std::ofstream fpvd("cells_data.pvd");
      fpvd.setf(std::ios_base::scientific|std::ios_base::unitbuf);
      //
      fpvd << "<?xml version=\"1.0\"?>" << std::endl;
      fpvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      fpvd << "<Collection>" << std::endl;
      for (int i=0; i<=max_step; i++) {
        const double t = (i+0.0)/max_step;
        if ( 0 == i || max_step == i || i%Param::visualization_export_interval_ == 0 ) {
          fpvd << "<DataSet timestep=\"" << t << "\" group=\"\" part=\"0\" file=\"cells_data_" << i << ".pvtu\"/>" << std::endl;
        }
      }
      fpvd << "</Collection>" << std::endl;
      fpvd << "</VTKFile>" << std::endl;
    }

    // iterate for all time-steps
    for (int i=0; i<=max_step; i++) {
      if (0==i%10) cout << "step " << i << " out of " << max_step << endl;
      //
      scheduler.Simulate(1);
      //
    }

    /*
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
    */
  }



} // namespace bdm

#endif // CANCER_GROWTH_H_
