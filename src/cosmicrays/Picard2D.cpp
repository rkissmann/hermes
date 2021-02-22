#ifdef HERMES_HAVE_HDF5

#include "hermes/cosmicrays/Picard2D.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "hermes/Common.h"

#define DEFAULT_CR_PATH \
	"/Test_Picard_2D/Test_Picard_2D_tfinal"

namespace hermes { namespace cosmicrays {

Picard2D::Picard2D(const std::string &dirname_, const PID &pid_)
    : CosmicRayDensity(pid_), dirname(dirname_) {
	readFile();
}

Picard2D::Picard2D(const PID &pid_)
    : CosmicRayDensity(pid_), dirname(getDataPath(DEFAULT_CR_PATH)) {
	readFile();
}

Picard2D::Picard2D(const std::vector<PID> &pids_)
    : CosmicRayDensity(pids_), dirname(getDataPath(DEFAULT_CR_PATH)) {
	readFile();
}

Picard2D::Picard2D(const std::string &dirname_, const std::vector<PID> &pids_)
    : CosmicRayDensity(pids_), dirname(dirname_) {
	readFile();
}

void Picard2D::readFile() {
	// file is automatically opened by constructor
	h5file = std::make_unique<Hdf5Reader>(Hdf5Reader(filename));

	// read energies
	readEnergyAxis();

	readSpatialGrid2D();
	readDensity2D();
}

QPDensityPerEnergy Picard2D::getDensityPerEnergy(
    const QEnergy &E_, const Vector3QLength &pos_) const {
	return getDensityPerEnergy(static_cast<int>(energyIndex.at(E_)), pos_);
}

QPDensityPerEnergy Picard2D::getDensityPerEnergy(
    int iE_, const Vector3QLength &pos_) const {
	if (pos_.z < zmin || pos_.z > zmax) return QPDensityPerEnergy(0);

	QLength rho = sqrt(pos_.x * pos_.x + pos_.y * pos_.y);
	if (rho > rmax) return QPDensityPerEnergy(0);

	auto pos = Vector3QLength(rho, pos_.z, 0);
	return (grid[iE_])->interpolate(pos);
}

void Picard2D::readEnergyAxis() {
	QEnergy Energy;

	// Determine number of energies stored in file:
	h5file->ReadGlobalAttribute("Entries", dimE);

	// Read individual energies from file:
	for (int iE = 0; iE < dimE; ++iE) {
		// get name of dataset
		std::string dsetName(getDsetName(iE));

		double valE;
		h5file->ReadAttributeFromDataset(dsetName,"Etot",valE);
		Energy = 1_MeV * valE;
		energyRange.push_back(Energy);
		energyIndex[Energy] = iE;
	}


}


std::string Picard2D::getDsetName(std::size_t iE) {
	// Attribute name, where dataset name is stored
	std::stringstream sstream;
	sstream << "Name_om";
	sstream << std::setfill('0') << std::setw(2) << iE;
	std::string attrName = sstream.str();

	// Get name of dataset
	std::string dsetName;
	h5file->ReadGlobalAttribute(attrName, dsetName);

	return dsetName;
}

void Picard2D::readSpatialGrid2D() {

	std::vector<float> rCen, zCen;
	h5file->ReadGlobalAttribute("rGridCentred", rCen);
	h5file->ReadGlobalAttribute("zGridCentred", zCen);

	rmin = rCen[0] * 1_kpc;
	rmax = rCen[rCen.size()-1] * 1_kpc;
	zmin = zCen[0] * 1_kpc;
	zmax = zCen[rCen.size()-1] * 1_kpc;

	dimr = rCen.size();
	dimz = zCen.size();

	QLength deltar = (rmax - rmin) / (dimr - 1);
	QLength deltaz = (zmax - zmin) / (dimz - 1);

	// Vector3d origin(-1*rmax.getValue(), -1*rmax.getValue(),
	// zmin.getValue());
	Vector3d origin(-1 * static_cast<double>(rmax), static_cast<double>(zmin),
	                0);
	Vector3d spacing(static_cast<double>(deltar), static_cast<double>(deltaz),
	                 0);

	for (int i = 0; i < dimE; ++i) {
		grid.push_back(std::make_unique<ScalarGrid2DQPDensityPerEnergy>(
		    ScalarGrid2DQPDensityPerEnergy(origin, dimr, dimz, spacing)));
	}
}

std::size_t Picard2D::calcArrayIndex2D(std::size_t ir, std::size_t iz) {
	return (iz * dimr + ir);
}

void Picard2D::readDensity2D() {

	auto vecSize = dimr * dimz;
	unsigned long nElements = energyRange.size() * vecSize;
	constexpr double fluxToDensity =
			static_cast<double>(4_pi / (c_light * 1_MeV));

	// Read A & Z from file
	int AVal, ZVal;
	int indexAttr = h5file->FindAttrIndex("A of part");
	h5file->ReadGlobalAttribute(indexAttr, AVal);
	indexAttr = h5file->FindAttrIndex("Z of part");
	h5file->ReadGlobalAttribute(indexAttr, ZVal);

	std::vector<int> dims(2);
	std::vector<float> h5data;


	if (isPIDEnabled(PID(ZVal, AVal))) {
		std::cerr << "hermes: info: reading species with Z = " << ZVal
				<< " A = " << AVal << std::endl;
		// loop over all energies and read spatial data.
		for (std::size_t iE = 0; iE < dimE; ++iE) {
			// get name of dataset
			std::string dsetName(getDsetName(iE));

			// Load spatial distribution at energy E
			h5file->ReadDataset(dsetName, dims, h5data);

			// Store data:
			for (std::size_t ir = 0; ir < dimr; ++ir) {
				for (std::size_t iz = 0; iz < dimz; ++iz) {
					grid[iE]->addValue(ir, iz, fluxToDensity * h5data[ calcArrayIndex2D(ir, iz)]);
				}
			}

		} // energy loop

	}
}


}}  // namespace hermes::cosmicrays

#endif  // HERMES_HAVE_HDF5
