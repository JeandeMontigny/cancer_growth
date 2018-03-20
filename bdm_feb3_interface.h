//=================================================================================================
#ifndef __BDM_FEB3_INTERFACE__
#define __BDM_FEB3_INTERFACE__
//=================================================================================================
#include <vector>
//=================================================================================================
namespace bdm_feb3 {
//=================================================================================================
typedef double  Real;
//=================================================================================================
class Vector3d {
public:
    Vector3d () { _comp[0]=_comp[1]=_comp[2]=0.0; }
    Vector3d (const Real& x, const Real& y, const Real& z) { _comp[0]=x; _comp[1]=y; _comp[2]=z; }
    Vector3d (const Vector3d& v) { _comp[0]=v(0); _comp[1]=v(1); _comp[2]=v(2); }
    ~Vector3d () {}
    //
    inline
    Real& operator() (unsigned int i) { return _comp[i]; }
    inline
    const Real& operator() (unsigned int i) const { return _comp[i]; }
private:
    Real _comp[3];
};
//-------------------------------------------------------------------------------------------------
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
        this->ecm_density_gradient = Vector3d();
        this->biochemical_level.assign(n_biochem, 0.0);
        this->biochemical_level_gradient.assign(n_biochem, Vector3d());
        this->interstitial_fluid_velocity = Vector3d();
    }
public:
    std::vector<Real> cells_mass;
    std::vector<unsigned int> cells_population;
    Real     ecm_density;
    Vector3d ecm_density_gradient;
    std::vector<Real>     biochemical_level;
    std::vector<Vector3d> biochemical_level_gradient;
    Vector3d interstitial_fluid_velocity;
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
