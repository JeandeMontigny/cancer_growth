#ifndef BIO_MODULE_H_
#define BIO_MODULE_H_

namespace bdm {

// 1. Define growth behaviour
  struct GrowthModule : public BaseBiologyModule {
    
  GrowthModule() : BaseBiologyModule(gAllBmEvents) {}
    
    template <typename T>
      void Run(T* cell) {

      // if normoxy: high division rate but low migration
      if (cell->GetOxygenLevel() > 0.7) { 
        //cell grow until it reach a diam of 8
        if (cell->GetDiameter() < 8) {
          cell->ChangeVolume(100);
          
          array<double, 3> cell_movements{gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1), gTRandom.Uniform(-1, 1)}; // create an array of 3 ramdom numbers between -1 and 1
          cell->UpdateMassLocation(cell_movements);
          cell->SetPosition(cell->GetMassLocation());
          //Reset biological movement to 0.
          cell->SetTractorForce({0, 0, 0});
        }
        // when diam of 8, it has a chance of dividing 
        else {
          double aNewRandomDouble=gTRandom.Uniform(0, 1);
      
          if (aNewRandomDouble <= 0.8 && cell->GetCanDivide()==true ) { //0.55 //0.65
            auto&& daughter = Divide(*cell);
            daughter.SetCellColour(cell->GetCellColour()); // daughter takes the cell_colour_ value of her mother
            daughter.SetCanDivide(true); // daughter will be able to divide
            daughter.SetOxygenLevel(cell->GetOxygenLevel()); // daughter takes the oxygen_level_ value of her mother
          }
          // if it doesn't divide, it will never be able to divide
          else {
            if (cell->GetCanDivide()==true) {
              cell->SetCanDivide(false);
            }
          }
        }
      }
      
      // if hypoxy: low division rate but high migration
      else if (cell->GetOxygenLevel() > 0.4) {
        // moves in anycase
//TODO: not a random migration. depending on oxygen gradient
        array<double, 3> cell_movements{gTRandom.Uniform(-4, 4), gTRandom.Uniform(-4, 4), gTRandom.Uniform(-4, 4)}; // create an array of 3 ramdom numbers between -1 and 1
        cell->UpdateMassLocation(cell_movements);
        cell->SetPosition(cell->GetMassLocation());
        //Reset biological movement to 0.
        cell->SetTractorForce({0, 0, 0});

        //cell grow until it reach a diam of 8
        if (cell->GetDiameter() < 8) {
          cell->ChangeVolume(40); //100
        }
        // when diam of 8, it has a chance of dividing 
        else {
          double aNewRandomDouble=gTRandom.Uniform(0, 1);
      
          if (aNewRandomDouble <= 0.2 && cell->GetCanDivide()==true ) { // 0.8 for normoxy
            auto&& daughter = Divide(*cell);
            daughter.SetCellColour(cell->GetCellColour()); // daughter takes the cell_colour_ value of her mother
            daughter.SetCanDivide(true); // daughter will be able to divide
            daughter.SetOxygenLevel(cell->GetOxygenLevel()); // daughter takes the oxygen_level_ value of her mother
            daughter.SetHypoDiv(true); // daughter will be able to divide
          }
          // if it doesn't divide, it will never be able to divide
          else {
            if (cell->GetCanDivide()==true && cell->GetHypoDiv()==true) {
              cell->SetHypoDiv(false);
            }
          }
        }
      }
		
      // if necrotic tissue: don't do anything (proba to die?)
//      else if (gTRandom.Uniform(0, 1) <=0.00001) { // proba to die at every time step
//        Delete(*cell); // Delete() method isn't implemented in this version
//      }
      
    } // end of Run()
    
//    bool IsCopied(BmEvent event) const { return true; }
    ClassDefNV(GrowthModule, 1);
  };

} // end namespace bdm

#endif
