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
    InterfaceData (unsigned int id, unsigned int nct, unsigned int nbc, unsigned int nms) {
        this->init(id, nct, nbc, nms);
    }
    ~InterfaceData () {}
    //
    inline
    void init (unsigned int id, unsigned int nct, unsigned int nbc, unsigned int nms) {
        this->_id = id;
        this->_n_cell_types = nct;
        this->_n_biochem = nbc;
        this->_n_mstruct = nms;
        this->cells_mass.assign(this->_n_cell_types, 0.0);
        this->cells_population.assign(this->_n_cell_types, 0);
        this->biochemical_level.assign(this->_n_biochem, 0.0);
        this->biochemical_level_gradient.assign(this->_n_biochem, RealVector3d({0.0, 0.0, 0.0}));
        this->microstructural_level.assign(this->_n_mstruct, 0.0);
        this->microstructural_level_gradient.assign(this->_n_mstruct, RealVector3d({0.0, 0.0, 0.0}));
        this->interstitial_fluid_velocity = RealVector3d({0.0, 0.0, 0.0});
    }
    inline
    unsigned int id () const { return this->_id; }
    inline
    unsigned int n_cell_types () const { return this->_n_cell_types; }
    inline
    unsigned int n_biochem () const { return this->_n_biochem; }
    inline
    unsigned int n_mstruct () const { return this->_n_mstruct; }
public:
    std::vector<Real> cells_mass;
    std::vector<unsigned int> cells_population;
    std::vector<Real>         biochemical_level;
    std::vector<RealVector3d> biochemical_level_gradient;
    std::vector<Real>         microstructural_level;
    std::vector<RealVector3d> microstructural_level_gradient;
    RealVector3d interstitial_fluid_velocity;
    unsigned int escaped_cells[6];
private:
    unsigned int _id, // simulation index
                 _n_cell_types, // number of cell types or/and phenotypes
                 _n_biochem, _n_mstruct; // number of biochemical and micro-structural variables
};
//-------------------------------------------------------------------------------------------------
class BDM_Domain : public InterfaceData {
public:
    BDM_Domain () : _is_init(false) {}
    BDM_Domain (unsigned int id, unsigned int nct, unsigned int nbc, unsigned int nms) {
        BDM_Domain::init(id, nct, nbc, nms);
    }
    ~BDM_Domain () {}
    //
    inline
    void init (unsigned int id, unsigned int nct, unsigned int nbc, unsigned int nms) {
        this->_is_init = true;
        InterfaceData::init(id, nct, nbc, nms);
    }
    inline
    bool is_init () const { return this->_is_init; }
private:
    bool _is_init;
};
//=================================================================================================
#endif
//=================================================================================================
