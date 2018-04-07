#ifndef BIO_MODULE_H_
#define BIO_MODULE_H_

namespace bdm {



  struct HostCellBiologyModule : public BaseBiologyModule {
  public:
    HostCellBiologyModule () : BaseBiologyModule(gAllBmEvents) {}
    
    template <typename T>
    void Run (T* cell) {
      int growthSpeed;
      array<double, 3> cell_movements;
      double divideProba;

      if (cell->GetOxygenLevel() > 0.8) { // normoxia
        growthSpeed = 0;
        cell_movements = {gTRandom.Uniform(-0.05, 0.05), gTRandom.Uniform(-0.05, 0.05), gTRandom.Uniform(-0.05, 0.05)};
        divideProba = 0.0001;
        cell->SetHypoDiv(true);
      }
      else if (cell->GetOxygenLevel() > 0.6) { // hypoxia
        growthSpeed = 0;
        cell_movements = {gTRandom.Uniform(-0.0001, 0.0001), gTRandom.Uniform(-0.0001, 0.0001), gTRandom.Uniform(-0.0001, 0.0001)};
        divideProba = 0.0001;
        cell->SetHypoDiv(false);
      }
      else { // apotposis/necrosis
        growthSpeed = 0;
        cell_movements = {0.0, 0.0, 0.0};
        divideProba = 0.000;
        cell->SetHypoDiv(false);
      }
        
      cell->ChangeVolume(growthSpeed);
      cell->SetCanDivide(gTRandom.Uniform(0.0, 1.0) < divideProba ? true : false);
      cell->UpdateMassLocation(cell_movements);
      cell->SetPosition(cell->GetMassLocation());
      cell->SetTractorForce({0, 0, 0});

      if (cell->GetCanDivide()) {
        auto&& daughter = Divide(*cell);
        daughter.SetCellColour(cell->GetCellColour());
        daughter.SetCanDivide(cell->GetCanDivide());
        daughter.SetHypoDiv(cell->GetHypoDiv());
        daughter.SetOxygenLevel(cell->GetOxygenLevel());
        daughter.SetTractorForce(cell->GetTractorForce());
      }

    }
    
    ClassDefNV (HostCellBiologyModule, 1);
  }; // end: HostCellBiologyModule



  struct CancerCellBiologyModule : public BaseBiologyModule {
  public:
    CancerCellBiologyModule () : BaseBiologyModule(gAllBmEvents) {}
    
    template <typename T>
    void Run (T* cell) {
      int growthSpeed;
      array<double, 3> cell_movements;
      double divideProba;

      if (cell->GetOxygenLevel() > 0.7) { // normoxia
        growthSpeed = 300;
        cell_movements = {gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1)};
        divideProba = 0.200;
        cell->SetHypoDiv(true);
      }
      else if (cell->GetOxygenLevel() > 0.3) { // hypoxia
        growthSpeed = 40;
        cell_movements = {gTRandom.Uniform(-4, 4), gTRandom.Uniform(-4, 4), gTRandom.Uniform(-4, 4)};
        divideProba = 0.100;
        cell->SetHypoDiv(false);
      }
      else { // apotposis/necrosis
        growthSpeed = 0;
        cell_movements = {0.0, 0.0, 0.0};
        divideProba = 0.000;
        cell->SetHypoDiv(false);
      }
        
      // cell grows until it reaches a diameter of ...
      if (cell->GetDiameter() > 15.0) {
        growthSpeed = 0;
      }
      if (cell->GetDiameter() < 10.0) {
        divideProba = 0.000;
      }

      cell->ChangeVolume(growthSpeed);
      cell->SetCanDivide(gTRandom.Uniform(0.0, 1.0) < divideProba ? true : false);
      cell->UpdateMassLocation(cell_movements);
      cell->SetPosition(cell->GetMassLocation());
      cell->SetTractorForce({0, 0, 0});

      if (cell->GetCanDivide()) {
        auto&& daughter = Divide(*cell);
        daughter.SetCellColour(cell->GetCellColour());
        daughter.SetCanDivide(cell->GetCanDivide());
        daughter.SetHypoDiv(cell->GetHypoDiv());
        daughter.SetOxygenLevel(cell->GetOxygenLevel());
        daughter.SetTractorForce(cell->GetTractorForce());
      }

    }
    
    ClassDefNV (CancerCellBiologyModule, 1);
  }; // end: CancerCellBiologyModule



} // end namespace bdm

#endif
