#ifndef CANCER_GROWTH_H_
#define CANCER_GROWTH_H_

#include <fstream>
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

        cell->ChangeVolume(0.0);
        cell->UpdateMassLocation(cell_movements);
        cell->SetPosition(cell->GetMassLocation());
        cell->SetTractorForce({0, 0, 0});
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

        // cell grows until it reaches a diameter of ...
        if (cell->GetDiameter() < 8.0) {
          cell->ChangeVolume(growthSpeed);
          cell->UpdateMassLocation(cell_movements);
          cell->SetPosition(cell->GetMassLocation());
          cell->SetTractorForce({0, 0, 0});
        }
        else if (cell->GetDiameter() >= 8.0) {
          cell->ChangeVolume(0.0);
          cell->UpdateMassLocation(cell_movements);
          cell->SetPosition(cell->GetMassLocation());
          cell->SetTractorForce({0, 0, 0});
        }

        if (cell->GetCanDivide() && cell->GetDiameter() > 7.0) {
          double aNewRandomDouble = gTRandom.Uniform(0.0, 1.0);
          if (aNewRandomDouble <= divideProba) {
            auto&& daughter = Divide(*cell);
            daughter.SetCellColour(cell->GetCellColour()); // daughter takes the cell_colour_ value of her mother
            daughter.SetCanDivide(true); // daughter will be able to divide
            daughter.SetHypoDiv(cell->GetHypoDiv()); // daughter will be able to divide in hypoxy
            daughter.SetOxygenLevel(cell->GetOxygenLevel()); // daughter takes the oxygen_level_ value of her mother
            daughter.SetIsCancerous(true);
          }
        }
      }// end run

      ClassDefNV (CancerCellBiologyModule, 1);
    }; // end: CancerCellBiologyModule

  // 2. Define the compile-time parameter
  template <typename Backend>
  struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
    // using BiologyModules = Variant<HostCellBiologyModule>;
    using BiologyModules = Variant<CancerCellBiologyModule, HostCellBiologyModule>;
    using AtomicTypes = VariadicTypedef<MyCell>;
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
  void Initialise (BDM_Domain& rve) {

    std::cout << " -- initialization of simulation " << rve.simulationNb << " --" << std::endl;

    // set-up these simulation parameters
    Param::live_visualization_ = false;
    Param::export_visualization_ = false;
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
      cell.AddBiologyModule(HostCellBiologyModule());
      cell.SetCanDivide(false);
      // cell.SetHypoDiv(true);
      cell.SetIsCancerous(false);
      return cell;
    };

    // cell creation (min boundary, max boundary, # of cells, default initialiser for cells)
    auto Construct_Cancer_Cells =  [](const std::array<double, 3>& position) {
      MyCell cell(position);
      cell.SetDiameter(6.0);
      cell.SetCellColour(1);
      cell.AddBiologyModule(CancerCellBiologyModule());
      cell.SetCanDivide(true);
      cell.SetHypoDiv(true);
      cell.SetIsCancerous(true);
      return cell;
    };

    std::string fileName="simulation_"+std::to_string(rve.simulationNb)+".txt";
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
      cells->reserve(nb_of_cells); // -1?

      array<double, 9> lineInfo;

      std::ifstream fileBis(fileName);
      while (std::getline(fileBis, line)) {
        // retreived info from line
        size_t pos = 0;
        int i = 0;
        while ((pos=line.find(" ")) != std::string::npos) {
//          std::cout << "parsing item: " << line.substr(0, pos) << std::endl;
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
    }

    else {
      // cell creation (min boundary, max boundary, # of cells, default initialiser for cells)
      CellCreator(0.01, 99.99, rve.cells_population[0], Construct_Host_Cells);
      CellCreator(45.00, 55.00, rve.cells_population[1], Construct_Cancer_Cells);
    }

//    cout << "regular cells created = " << rve.cells_population[0] << endl;
//    cout << "cancerous cells created = " << rve.cells_population[1] << endl;

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
    const int max_step = 250;

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
    auto rm = TResourceManager::Get();

    std::cout << " -- running simulation " << rve.simulationNb << " --" << std::endl;

    for (int i=0; i<=max_step; i++) {
      //
      scheduler.Simulate(1);
      //
      auto all_cells = rm->template Get<MyCell>();
      for (unsigned int l=0; l<rve.n_cell_types(); l++) {
        rve.cells_mass[l] = 0.0;
        rve.cells_population[l] = 0;
      }
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
        std::cout << "; # of cells (" << rve.cells_population[0] << "," << rve.cells_population[1] << ")" << std::flush;
        std::cout << std::endl;
      }

    } // end iterate for all time-steps

    // TODO: create separated folder for each simulation, labelled by the BDM_simulation number

    // export simulation in external file
    ofstream outputFile;
    outputFile.open("simulation_" + std::to_string(rve.simulationNb) + ".txt");

    std::cout << " -- exporting simulation " << rve.simulationNb << " --" << std::endl;

    auto my_cells = rm->template Get<MyCell>();
    int numberOfCells = my_cells->size();

    for (int i; i < numberOfCells; i++) {
      auto thisCell = (*my_cells)[i];
      array<double, 3> thisPosition = thisCell.GetMassLocation();

      outputFile << i << " " << thisPosition[0] << " " << thisPosition[1] << " " << thisPosition[2] << " " << thisCell.GetDiameter() << thisCell.GetCanDivide() << " " << thisCell.GetOxygenLevel() << " " << thisCell.GetHypoDiv() << " " << thisCell.GetIsCancerous() << "\n";
    }

    outputFile.close();

  } // end simulate



} // namespace bdm

#endif // CANCER_GROWTH_H_
