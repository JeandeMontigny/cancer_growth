//=================================================================================================
#ifndef __BDM_FEB3_INTERFACE__
#define __BDM_FEB3_INTERFACE__
//=================================================================================================
#include <vector>
#include <array>
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
        this->_n_cell_types = n_cell_types;
        this->_n_biochem = n_biochem;
        this->cells_mass.assign(this->_n_cell_types, 0.0);
        this->cells_population.assign(this->_n_cell_types, 0);
        this->ecm_density = 0.0;
        this->ecm_density_gradient = RealVector3d({0.0, 0.0, 0.0});
        this->biochemical_level.assign(this->_n_biochem, 0.0);
        this->biochemical_level_gradient.assign(this->_n_biochem, RealVector3d({0.0, 0.0, 0.0}));
        this->interstitial_fluid_velocity = RealVector3d({0.0, 0.0, 0.0});
    }
    inline
    unsigned int n_cell_types () const { return this->_n_cell_types; }
    inline
    unsigned int n_biochem () const { return this->_n_biochem; }
public:
    std::vector<Real> cells_mass;
    std::vector<unsigned int> cells_population;
    Real         ecm_density;
    RealVector3d ecm_density_gradient;
    std::vector<Real>         biochemical_level;
    std::vector<RealVector3d> biochemical_level_gradient;
    RealVector3d interstitial_fluid_velocity;
private:
    unsigned int _n_cell_types, 
                 _n_biochem;
};
//-------------------------------------------------------------------------------------------------
class BDM_Domain : public InterfaceData {
public:
    BDM_Domain () : _is_init(false) {}
    BDM_Domain (int n_cell_types, int n_biochem) { BDM_Domain::init(n_cell_types, n_biochem); }
    ~BDM_Domain () {}
    //
    inline
    void init (int n_cell_types, int n_biochem) {
        this->_is_init = true;
        InterfaceData::init(n_cell_types, n_biochem);
    }
    inline
    bool is_init () const { return this->_is_init; }
private:
    bool _is_init;
};
//=================================================================================================
#endif
//=================================================================================================
