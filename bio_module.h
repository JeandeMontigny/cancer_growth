#ifndef BIO_MODULE_H_
#define BIO_MODULE_H_

namespace bdm {



  // 1. Define growth behaviour
  struct GrowthModule : public BaseBiologyModule {


    
    GrowthModule() : BaseBiologyModule(gAllBmEvents) {}
    
    template <typename T>
    void Run(T* cell) {
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
        
      //cell grow until it reach a diam of ...
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
        if (aNewRandomDouble <= divideProba)
        {
          auto&& daughter = Divide(*cell);
          daughter.SetCellColour(cell->GetCellColour()); // daughter takes the cell_colour_ value of her mother
          daughter.SetCanDivide(true); // daughter will be able to divide
          daughter.SetHypoDiv(cell->GetHypoDiv()); // daughter will be able to divide in hypoxy
          daughter.SetOxygenLevel(cell->GetOxygenLevel()); // daughter takes the oxygen_level_ value of her mother
        }
      }

      // when diam of 8, it has a chance of dividing (if it can divide)
      /*
      else if (cell->GetCanDivide()==true && cell->GetHypoDiv()==true) {
        double aNewRandomDouble=gTRandom.Uniform(0, 1);
        if (aNewRandomDouble <= divideProba) { //0.55 //0.65
          auto&& daughter = Divide(*cell);
          daughter.SetCellColour(cell->GetCellColour()); // daughter takes the cell_colour_ value of her mother
          daughter.SetCanDivide(true); // daughter will be able to divide
          daughter.SetHypoDiv(true); // daughter will be able to divide in hypoxy
          daughter.SetOxygenLevel(cell->GetOxygenLevel()); // daughter takes the oxygen_level_ value of her mother
        }
        // if it doesn't divide
        else {
          // it won't divide in hypoxy any more
          if (currentOxygenLevel <= 0.7) {
            cell->SetHypoDiv(false);
          }
          // it will never be able to divide any more
          if (currentOxygenLevel > 0.7) {
            cell->SetCanDivide(false);
          }
        }
      }
      */
    }


    
    ClassDefNV(GrowthModule, 1);


    
  }; // end GrowthModule

} // end namespace bdm

#endif
