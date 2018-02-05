#ifndef BIO_MODULE_H_
#define BIO_MODULE_H_

namespace bdm {

// 1. Define growth behaviour
  struct GrowthModule : public BaseBiologyModule {
    
  GrowthModule() : BaseBiologyModule(gAllBmEvents) {}
    
    template <typename T>
      void Run(T* cell) {

      int growthSpeed;
      array<double, 3> cell_movements;
      double divideProba;
      double currentOxygenLevel=cell->GetOxygenLevel();
      // if normoxy: high division rate but low migration
      if (currentOxygenLevel > 0.7) {
        cell->SetHypoDiv(true);
        growthSpeed=100;
        cell_movements={gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1)}; // create an array of 3 ramdom numbers between -1 and 1
        divideProba=0.8;
      }
      // if hypoxy: low division rate but high migration
      else if (currentOxygenLevel > 0.3) {
        growthSpeed = 40;
        //TODO: not a random migration. depending on oxygen gradient. Can get oxygen gradien from feb3
        cell_movements={gTRandom.Uniform(-4, 4), gTRandom.Uniform(-4, 4), gTRandom.Uniform(-4, 4)}; // create an array of 3 ramdom numbers between -4 and 4
        divideProba=0.4;
      }
      // necrose
      else {
        return;
        growthSpeed=0;
        cell_movements={0, 0, 0};
        divideProba=0;
        cell->SetHypoDiv(false);
      }
        
      //cell grow until it reach a diam of 8
      if (cell->GetDiameter() < 8) {
        cell->ChangeVolume(growthSpeed);
        cell->UpdateMassLocation(cell_movements);
        cell->SetPosition(cell->GetMassLocation());
        //Reset biological movement to 0.
        cell->SetTractorForce({0, 0, 0});
      }
      // when diam of 8, it has a chance of dividing (if it can divide)
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
      
    } // end of Run()
    
//    bool IsCopied(BmEvent event) const { return true; }
    ClassDefNV(GrowthModule, 1);
    
  }; // end GrowthModule

} // end namespace bdm

#endif
