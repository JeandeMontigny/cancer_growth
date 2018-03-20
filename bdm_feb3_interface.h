//=================================================================================================
#ifndef __BDM_FEB3_INTERFACE__
#define __BDM_FEB3_INTERFACE__
//=================================================================================================
#include <vector>
#include <array>
//=================================================================================================
namespace bdm_feb3 {
//=================================================================================================
typedef double                 Real;
typedef std::array<double, 3>  RealVector3d;
//=================================================================================================
class InterfaceData {
public:
    InterfaceData () {}
    InterfaceData (int n_cell_types, int n_biochem) { this->init(n_cell_types, n_biochem); }
    ~InterfaceData () {}
    //
    inline
    void init (int n_cell_types, int n_biochem) {
        this->cells_mass.assign(n_cell_types, 0.0);
        this->cells_population.assign(n_cell_types, 0);
        this->ecm_density = 0.0;
        this->ecm_density_gradient = RealVector3d({0.0, 0.0, 0.0});
        this->biochemical_level.assign(n_biochem, 0.0);
        this->biochemical_level_gradient.assign(n_biochem, RealVector3d({0.0, 0.0, 0.0}));
        this->interstitial_fluid_velocity = RealVector3d({0.0, 0.0, 0.0});
    }
public:
    std::vector<Real> cells_mass;
    std::vector<unsigned int> cells_population;
    Real         ecm_density;
    RealVector3d ecm_density_gradient;
    std::vector<Real>         biochemical_level;
    std::vector<RealVector3d> biochemical_level_gradient;
    RealVector3d interstitial_fluid_velocity;
};
//=================================================================================================
} // end of namespace: bdm_feb3
//=================================================================================================
class BDM_CubicDomain {
public:
    BDM_CubicDomain () : _is_init(false) {}
    BDM_CubicDomain (int n_cell_types, int n_biochem) { this->init(n_cell_types, n_biochem); }
    ~BDM_CubicDomain () {}
    //
    inline
    void init (int n_cell_types, int n_biochem) {
        this->data.init(n_cell_types, n_biochem);
        // this->_cells = bdm::daosoa<bdm::Cell>(ncell);
        this->_is_init = true;
    }
    inline
    bool is_init () const { return this->_is_init; }
public:
    bdm_feb3::InterfaceData data;
private:
    bool _is_init;
    // bdm::daosoa<bdm::Cell> _cells;
};
//=================================================================================================
#endif
//=================================================================================================
