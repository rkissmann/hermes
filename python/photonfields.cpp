#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hermes/photonfields/PhotonField.h"
#include "hermes/photonfields/CMB.h"
#include "hermes/photonfields/ISRF.h"

namespace py = pybind11;

namespace hermes { namespace photonfields {

void init(py::module &m) {
    
    py::module subm = m.def_submodule("photonfields");
    subm.doc() = "photon fields package";

    // charged gas density models
    py::class_<PhotonField, std::shared_ptr<PhotonField>>(subm, "PhotonField")
	      .def("getEnergyDensity",
		   (QEnergyDensity (PhotonField::*)(const Vector3QLength &, const QEnergy &) const)
		   &PhotonField::getEnergyDensity)
	      .def("getEnergyDensity",
		   (QEnergyDensity (PhotonField::*)(const Vector3QLength &, int) const)
		   &PhotonField::getEnergyDensity)
	      .def("getEnergyAxis", &PhotonField::getEnergyAxis);
    py::class_<CMB, std::shared_ptr<CMB>, PhotonField>(subm, "CMB")
	      .def(py::init<>());
    py::class_<ISRF, std::shared_ptr<ISRF>, PhotonField>(subm, "ISRF")
	      .def(py::init<>());
}

} // namespace photonfields
} // namespace hermes