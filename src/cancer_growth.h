#ifndef CANCER_GROWTH_H_
#define CANCER_GROWTH_H_

#include <fstream>
#include <sys/stat.h>
#include "biodynamo.h"
#include "../bdm_feb3_interface.h"

namespace bdm {

  // 0. Define my custom cell, which extends Cell by adding an extra data member cell_type.
  BDM_SIM_OBJECT(MyCell, Cell) {
    BDM_SIM_OBJECT_HEADER(MyCellExt, 1, can_divide_, cell_colour_, oxygen_level_, hypo_division_, is_cancerous_); // create the header with our new data member

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

    void SetIsCancerous(bool cancerous) {is_cancerous_[kIdx] = cancerous; }
    bool GetIsCancerous() { return is_cancerous_[kIdx]; }
    bool* GetIsCancerousPtr() { return is_cancerous_.data(); }

  private:
    vec<bool> can_divide_;
    vec<int> cell_colour_;
    vec<double> oxygen_level_; // should be between 0 and 1
    vec<bool> hypo_division_;
    vec<bool> is_cancerous_;
  };

  struct HostCellBiologyModule : public BaseBiologyModule {
  public:
    HostCellBiologyModule () : BaseBiologyModule(gAllBmEvents) {}

    template <typename T>
    void Run (T* cell) {
      const double currentOxygenLevel = cell->GetOxygenLevel();
      int growthSpeed;
      array<double, 3> cell_movements;
      double divideProba;

      // normoxia: high division rate but low migration
      if (currentOxygenLevel > 0.7) {
        growthSpeed = 0;
        cell_movements = {gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1)};
        divideProba = 0.9;
        cell->SetHypoDiv(true);
      }
      // hypoxia: low division rate but high migration
      else if (currentOxygenLevel > 0.3) {
        growthSpeed = 0;
        cell_movements = {gTRandom.Uniform(-4, 4), gTRandom.Uniform(-4, 4), gTRandom.Uniform(-4, 4)};
        divideProba = 0.4;
        cell->SetHypoDiv(false);
      }
      // necrosis: no division and no migration, just die!
      else {
        return;
        growthSpeed = 0;
        cell_movements = {0.0, 0.0, 0.0};
        divideProba = 0.0;
        cell->SetHypoDiv(false);
      }

      //        cell->ChangeVolume(0.0);
      //        cell->UpdateMassLocation(cell_movements);
      //        cell->SetPosition(cell->GetMassLocation());
      //        cell->SetTractorForce({0, 0, 0});
    }

    ClassDefNV (HostCellBiologyModule, 1);
  }; // end: HostCellBiologyModule

  struct CancerCellBiologyModule : public BaseBiologyModule {
  public:
    CancerCellBiologyModule () : BaseBiologyModule(gAllBmEvents) {}

    template <typename T>
    void Run (T* cell) {
      const double currentOxygenLevel = cell->GetOxygenLevel();
      int growthSpeed;
      array<double, 3> cell_movements;
      double divideProba;

      // normoxia: high division rate but low migration
      if (currentOxygenLevel > 0.7) {
        growthSpeed = 100;
        cell_movements = {gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1)};
        divideProba = 0.9;
        cell->SetHypoDiv(true);
      }
      // hypoxia: low division rate but high migration
      else if (currentOxygenLevel > 0.3) {
        growthSpeed = 40;
        cell_movements = {gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1)};
        divideProba = 0.15;
        cell->SetHypoDiv(false);
      }
      // necrosis: no division and no migration, just die!
      else {
        return;
        growthSpeed = 0;
        cell_movements = {0.0, 0.0, 0.0};
        divideProba = 0.0;
        cell->SetHypoDiv(false);
      }

      // cell grows until it reaches a diameter of ...
      if (cell->GetDiameter() < 8.0) {
        cell->ChangeVolume(growthSpeed);
        cell->UpdatePosition(cell_movements);
        cell->SetPosition(cell->GetPosition());
        cell->SetTractorForce({0, 0, 0});
      }
      else if (cell->GetDiameter() >= 8.0) {
        cell->ChangeVolume(0.0);
        cell->UpdatePosition(cell_movements);
        cell->SetPosition(cell->GetPosition());
        cell->SetTractorForce({0, 0, 0});
      }

      if (cell->GetCanDivide() && cell->GetDiameter() > 7.0) {
        double aNewRandomDouble = gTRandom.Uniform(0.0, 1.0);
        if (aNewRandomDouble <= divideProba) {
          auto daughter = cell->Divide();
          daughter->SetCellColour(cell->GetCellColour()); // daughter takes the cell_colour_ value of her mother
          daughter->SetCanDivide(true); // daughter will be able to divide
          daughter->SetHypoDiv(cell->GetHypoDiv()); // daughter will be able to divide in hypoxy
          daughter->SetOxygenLevel(cell->GetOxygenLevel()); // daughter takes the oxygen_level_ value of her mother
          daughter->SetIsCancerous(true);
        }
      }
    }// end run

    ClassDefNV (CancerCellBiologyModule, 1);
  }; // end: CancerCellBiologyModule


enum Substances { cancerous_diffusion };

  // 1b. Define secretion behavior:
struct CancerousSecretion : public BaseBiologyModule {
  // Daughter cells inherit this biology module
  CancerousSecretion() : BaseBiologyModule(gAllBmEvents) {}

  template <typename T>
  void Run(T* sim_object) {
    if (sim_object->template IsSoType<MyCell>()) {
      auto&& cell = sim_object->template ReinterpretCast<MyCell>();

      if (!init_) {
        dg_ = GetDiffusionGrid(cancerous_diffusion);
        init_ = true;
      }
      auto& secretion_position = cell->GetPosition();
      dg_->IncreaseConcentrationBy(secretion_position, 0.5);
    }
  } // end run

private:
  bool init_ = false;
  DiffusionGrid* dg_ = nullptr;
  ClassDefNV(CancerousSecretion, 1);
}; // end biologyModule SubstanceSecretion

struct CancerousMovement : public BaseBiologyModule {
  // Daughter cells inherit this biology module
  CancerousMovement() : BaseBiologyModule(gAllBmEvents) {}

  template <typename T>
  void Run(T* sim_object) {
    if (sim_object->template IsSoType<MyCell>()) {
      auto&& cell = sim_object->template ReinterpretCast<MyCell>();

      std::array<double, 3> mov = { 0, 0, 0.8 }; // gTRandom.Uniform(-1, 1)

      cell->UpdatePosition(mov);
      cell->SetPosition(cell->GetPosition());
    }
  }

    private:
      ClassDefNV(CancerousMovement, 1);
  };


struct Lymphocytes : public BaseBiologyModule {
public:
  Lymphocytes () : BaseBiologyModule(gAllBmEvents) {}

  template <typename T>
  void Run (T* cell) {

    if (!init_) {
      dg_ = GetDiffusionGrid(cancerous_diffusion);
      init_ = true;
    }
    auto& position = cell->GetPosition();
    std::array<double, 3> gradient;
    double concentration;

    dg_->GetGradient(position, &gradient);
    concentration = dg_->GetConcentration(position);

    std::cout << "concentration: " << concentration << std::endl;

    cell->UpdatePosition(gradient);
    cell->SetPosition(cell->GetPosition());
  }

private:
  bool init_ = false;
  DiffusionGrid* dg_ = nullptr;
  ClassDefNV (Lymphocytes, 1);
}; // end Lymphocytes


  // 2. Define the compile-time parameter
  template <typename Backend>
  struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
    // using BiologyModules = Variant<HostCellBiologyModule>;
    using BiologyModules = Variant<CancerCellBiologyModule, HostCellBiologyModule, CancerousSecretion, CancerousMovement, Lymphocytes>;
    using AtomicTypes = VariadicTypedef<MyCell>;
  };


  // cell creator
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

  // MyCellCreator
    template <typename Function, typename TResourceManager = ResourceManager<>>
    static void MyCellCreator(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, int num_cells, Function cell_builder) {
      auto rm = TResourceManager::Get();
      // Determine simulation object type which is returned by the cell_builder
      using FunctionReturnType = decltype(cell_builder({0, 0, 0}));
      auto container = rm->template Get<FunctionReturnType>();
      container->reserve(num_cells);

      for (int i = 0; i < num_cells; i++) {
        double x = gTRandom.Uniform(x_min, x_max);
        double y = gTRandom.Uniform(y_min, y_max);
        double z = gTRandom.Uniform(z_min, z_max);
        auto new_simulation_object = cell_builder({x, y, z});
        container->push_back(new_simulation_object);
      }
      container->Commit();
  } // end MyCellCreator


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
  void Initialise (BDM_Domain& rve) {

    std::cout << " -- initialization of simulation " << rve.id() << " -- " << std::endl;

    ModelInitializer::DefineSubstance(cancerous_diffusion, "cancerous_diffusion", 1, 0.5, 10);

    // set-up these simulation parameters
    // Param::live_visualization_ = false;
    // Param::export_visualization_ = false;
    // Param::visualization_export_interval_ = 1;
    // Param::visualize_sim_objects_["MyCell"] = std::set<std::string>{"diameter_", "cell_colour_", "oxygen_level_"};
    Param::bound_space_ = false; // create artificial boundary limits for the 3D simulation space
    Param::min_bound_ = 0;
    Param::max_bound_ = 100.0;
    Param::run_mechanical_interactions_ = true;
    gTRandom.SetSeed(2448);

    auto Construct_Host_Cells =  [](const std::array<double, 3>& position) {
      MyCell cell({40, 40, 10});
//      MyCell cell(position);
      cell.SetDiameter(6.5);
      cell.SetCellColour(0);
      cell.AddBiologyModule(Lymphocytes());
      // cell.AddBiologyModule(HostCellBiologyModule());
      // cell.SetCanDivide(false);
      // // cell.SetHypoDiv(true);
      // cell.SetIsCancerous(false);
      return cell;
    };

    // cell creation (min boundary, max boundary, # of cells, default initialiser for cells)
    auto Construct_Cancer_Cells =  [](const std::array<double, 3>& position) {
      MyCell cell({50, 50, 20});
//      MyCell cell(position);
      cell.SetDiameter(6.0);
      cell.SetCellColour(1);
      cell.AddBiologyModule(CancerousSecretion());
      cell.AddBiologyModule(CancerousMovement());
      // cell.AddBiologyModule(CancerCellBiologyModule());
      // cell.SetCanDivide(true);
      // cell.SetHypoDiv(true);
      // cell.SetIsCancerous(true);
      return cell;
    };

    std::string fileName="simulation_"+std::to_string(rve.id())+".txt";
    // if exported external file exist
    if (access(fileName.c_str(), F_OK) != -1) {
      // read outputFile
      std::ifstream file(fileName);
      std::string line;
      int nb_of_cells=0;
      while (std::getline(file, line)) {
        nb_of_cells++;
      }

      // set up resource manager and reserve cells
      auto rm = TResourceManager::Get();
      auto cells = rm->template Get<MyCell>();
      cells->reserve(nb_of_cells);

      array<double, 9> lineInfo;

      std::ifstream fileBis(fileName);
      while (std::getline(fileBis, line)) {
        size_t pos = 0;
        int i = 0;
        // retreived info from line
        while ((pos=line.find(" ")) != std::string::npos) {
          lineInfo[i] = stod(line.substr(0, pos));
          line.erase(0, pos+1);
          i++;
        }
        lineInfo[8] = stod(line);

        // create cells depending on outputFile
        MyCell cell({lineInfo[1], lineInfo[2], lineInfo[3]});
        cell.SetDiameter(lineInfo[4]);
        cell.SetCanDivide(lineInfo[5]);
        cell.SetOxygenLevel(rve.biochemical_level[0]);
        cell.SetHypoDiv(lineInfo[7]);
        if (lineInfo[8] == 1) { // if cancerous cell
          cell.SetIsCancerous(true);
          cell.AddBiologyModule(CancerCellBiologyModule());
        }
        else {
          cell.SetIsCancerous(false);
          cell.AddBiologyModule(HostCellBiologyModule());
        }

        cells->push_back(cell);
      }
      cells->Commit();

      // create cells comming from neighbour cubes. 6 faces
      // create first x half (face 0, x=0) imcoming cells
      if (rve.face_incomingCells[0]!=0) {
        MyCellCreator(Param::min_bound_, Param::max_bound_/4,        Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1),        Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1), rve.face_incomingCells[0], Construct_Cancer_Cells);
      }
      // create second x half (face 1, x=max) imcoming cells
      if (rve.face_incomingCells[1]!=0) {
        MyCellCreator(Param::max_bound_/4, Param::max_bound_,        Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1),        Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1), rve.face_incomingCells[1], Construct_Cancer_Cells);
      }

      // create first y half (face 2, y=0) imcoming cells
      if (rve.face_incomingCells[2]!=0) {
        MyCellCreator(Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1), Param::min_bound_, Param::max_bound_/4, Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1), rve.face_incomingCells[2], Construct_Cancer_Cells);
      }
      // create second y half (face 3, y=max) imcoming cells
      if (rve.face_incomingCells[3]!=0) {
        MyCellCreator(Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1), Param::max_bound_/4, Param::max_bound_, Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1), rve.face_incomingCells[3], Construct_Cancer_Cells);
      }

      // create first z half (face 4, z=0) imcoming cells
      if (rve.face_incomingCells[4]!=0) {
        MyCellCreator(Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1), Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1), Param::min_bound_, Param::max_bound_/4, rve.face_incomingCells[4], Construct_Cancer_Cells);
      }
      // create second z half (face 5, z=max) imcoming cells
      if (rve.face_incomingCells[5]!=0) {
        MyCellCreator(Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1), Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1), Param::max_bound_/4, Param::max_bound_, rve.face_incomingCells[5], Construct_Cancer_Cells);
      }

    }

    else {
      // cell creation (min boundary, max boundary, # of cells, default initialiser for cells)
      CellCreator(Param::min_bound_+(Param::max_bound_*0.1), Param::max_bound_-(Param::max_bound_*0.1), rve.cells_population[0], Construct_Host_Cells);
      CellCreator(5.00, 15.00, rve.cells_population[1], Construct_Cancer_Cells);
    }

    auto rm = TResourceManager::Get();
    auto all_cells = rm->template Get<MyCell>();
    for (unsigned int l=0; l<rve.n_cell_types(); l++) {
      rve.cells_mass[l] = 0.0;
      rve.cells_population[l] = 0;
    }
    for (unsigned int i=0; i<all_cells->size(); i++) {
      auto&& cell = (*all_cells)[i];
      const unsigned int cell_type = cell.GetCellColour();
      //
      cell.SetOxygenLevel(rve.biochemical_level[0]);
      //
      rve.cells_mass[cell_type] += cell.GetMass();
      rve.cells_population[cell_type] += 1;
    }
  } // end Initialise


  template <typename TResourceManager = ResourceManager<>>
  inline void ResetBDM () {
    Param::Reset();
    TResourceManager::Get()->Clear();
  }


  // 4. Core simulation routine
  template <typename TResourceManager = ResourceManager<>>
  inline
  void Simulate (BDM_Domain& rve) {

    Scheduler<> scheduler;
    const int max_step = 81;

    // iterate for all time-steps
    auto rm = TResourceManager::Get();
    auto all_cells = rm->template Get<MyCell>();

    std::cout << " -- running simulation " << rve.id() << " -- " << std::endl;

    const double low_bound = Param::min_bound_ + 0.01*(Param::max_bound_-Param::min_bound_);
    const double up_bound  = Param::max_bound_ - 0.01*(Param::max_bound_-Param::min_bound_);
    unsigned int current_escaped_cells[] = {0, 0, 0, 0, 0, 0};

    for (int i=0; i<=max_step; i++) {
      //
      scheduler.Simulate(1);
      //
      // all_cells = rm->template Get<MyCell>();
      // for (unsigned int i=0; i<all_cells->size(); i++) {
      //   auto thisCell = (*all_cells)[i];
      //   array<double, 3> thisPosition = thisCell.GetPosition();
      //   if (thisPosition[0]<=low_bound) {
      //     thisCell.RemoveFromSimulation();
      //     current_escaped_cells[0] += 1;
      //   } else if (thisPosition[0]>=up_bound ) {
      //     thisCell.RemoveFromSimulation();
      //     current_escaped_cells[1] += 1;
      //   } else if (thisPosition[1]<=low_bound) {
      //     thisCell.RemoveFromSimulation();
      //     current_escaped_cells[2] += 1;
      //   } else if (thisPosition[1]>=up_bound ) {
      //     thisCell.RemoveFromSimulation();
      //     current_escaped_cells[3] += 1;
      //   } else if (thisPosition[2]<=low_bound) {
      //     thisCell.RemoveFromSimulation();
      //     current_escaped_cells[4] += 1;
      //   } else if (thisPosition[2]>=up_bound ) {
      //     thisCell.RemoveFromSimulation();
      //     current_escaped_cells[5] += 1;
      //   }
      // }

      for (int k=0; k<6; k++) {
        rve.escaped_cells[k] = current_escaped_cells[k];
      }
      for (unsigned int l=0; l<rve.n_cell_types(); l++) {
        rve.cells_mass[l] = 0.0;
        rve.cells_population[l] = 0;
      }
      all_cells = rm->template Get<MyCell>();
      for (unsigned int i=0; i<all_cells->size(); i++) {
        auto&& cell = (*all_cells)[i];
        const unsigned int cell_type = cell.GetCellColour();
        //
        rve.cells_mass[cell_type] += cell.GetMass();
        rve.cells_population[cell_type] += 1;
      }
      //
      if (0==i%10) {
        std::cout << " *** Time-step " << i << " out of " << max_step << std::flush;
        std::cout << "; # of cells:" << std::flush;
        for (unsigned int l=0; l<rve.n_cell_types(); l++) {
          std::cout << " " << rve.cells_population[l] << std::flush;
        }
        std::cout << std::endl;
      }

    } // end iterate for all time-steps

    // export simulation in external file
    ofstream outputFile;
    outputFile.open("simulation_" + std::to_string(rve.id()) + ".txt");

    std::cout << " -- exporting simulation " << rve.id() << " -- " << std::endl;

    all_cells = rm->template Get<MyCell>();
    int numberOfCells = all_cells->size();

    for (int i=0; i < numberOfCells; i++) {
      auto thisCell = (*all_cells)[i];
      array<double, 3> thisPosition = thisCell.GetPosition();
      //
      outputFile << i
                 << " " << thisPosition[0] << " " << thisPosition[1] << " " << thisPosition[2]
                 << " " << thisCell.GetDiameter()
                 << " " << thisCell.GetCanDivide()
                 << " " << thisCell.GetOxygenLevel()
                 << " " << thisCell.GetHypoDiv()
                 << " " << thisCell.GetIsCancerous()
                 << std::endl;
    }

    outputFile.close();

  } // end simulate

} // namespace bdm

#endif // CANCER_GROWTH_H_
