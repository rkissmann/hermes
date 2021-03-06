#ifdef HERMES_HAVE_CFITSIO

#include "hermes/neutralgas/RingModel.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>

#include "hermes/Common.h"

namespace hermes { namespace neutralgas {

RingData::RingData(GasType gas) : type(gas) {
	if (gas == GasType::HI) {
		readDataFile("NHrings_Ts300K.fits.gz");
	}
	if (gas == GasType::H2) {
		readDataFile("WCOrings_COGAL.fits.gz");
	}
}

void RingData::readDataFile(const std::string &filename) {
	ffile = std::make_unique<FITSFile>(
	    FITSFile(getDataPath("GasDensity/Remy18/" + filename)));
	ffile->openFile(FITS::READ);

	n_lon = ffile->readKeyValueAsInt("NAXIS1");
	n_lat = ffile->readKeyValueAsInt("NAXIS2");
	n_rings = ffile->readKeyValueAsInt("NAXIS3");

	min_lon = ffile->readKeyValueAsDouble("CRVAL1");
	delta_lon = ffile->readKeyValueAsDouble("CDELT1");
	min_lat = ffile->readKeyValueAsDouble("CRVAL2");
	delta_lat = ffile->readKeyValueAsDouble("CDELT2");

	std::cerr << "Number of rings: " << n_rings << std::endl;

	int firstElement = 1;
	int nElements = n_lon * n_lat * n_rings;
	dataVector = ffile->readImageAsFloat(firstElement, nElements);
}

GasType RingData::getGasType() const { return type; }

int RingData::getRingNumber() const { return n_rings; }

double RingData::getRawValue(int ring, const QDirection &dir) const {
	QAngle lat = 180_deg - dir[0];
	QAngle lon =
	    180_deg - dir[1];  // becaue the galactic centre of
	                       // the ring model is in the middle of the map

	int pxl_lat = static_cast<int>(round(lat / 180_deg * n_lat));
	int pxl_lon = static_cast<int>(round(lon / 360_deg * n_lon));

	// NAXIS1 x NAXIS2 x NAXIS3 => lon x lat x ring
	return dataVector[(ring * n_lat + pxl_lat) * n_lon + pxl_lon];
}

QColumnDensity RingData::getHIColumnDensityInRing(int ring,
                                                  const QDirection &dir) const {
	// the data is given in cm^-2
	return getRawValue(ring, dir) / 1_cm2;
}

QRingCOIntensity RingData::getCOIntensityInRing(int ring,
                                                const QDirection &dir) const {
	// the data is given in K km s^-2
	if (ring == 10 || ring == 11) return QRingCOIntensity(0);
	return getRawValue(ring, dir) * 1_K * 1_km / 1_s;
}

Ring::Ring(std::size_t index_, std::shared_ptr<RingData> dataPtr_,
           QLength innerR_, QLength outerR_)
    : index(index_),
      dataPtr(std::move(dataPtr_)),
      innerR(innerR_),
      outerR(outerR_) {}

Ring::~Ring() {}

std::size_t Ring::getIndex() const { return index; }

std::pair<QLength, QLength> Ring::getBoundaries() const {
	return std::make_pair(innerR, outerR);
}

bool Ring::isInside(const Vector3QLength &pos) const {
	QLength rho = pos.getRho();
	return (rho > innerR && rho < outerR);
}

QRingX0Unit Ring::X0Function(const QDirection &dir_) const {
	return 1.8e20 / (1_cm2 * 1_K * 1_km) * 1_s;
}

GasType Ring::getGasType() const { return dataPtr->getGasType(); }

QColumnDensity Ring::getHIColumnDensity(const QDirection &dir_) const {
	return dataPtr->getHIColumnDensityInRing(index, dir_);
}

QColumnDensity Ring::getH2ColumnDensity(const QDirection &dir_) const {
	return 2 * X0Function(dir_) * dataPtr->getCOIntensityInRing(index, dir_);
}

QColumnDensity Ring::getColumnDensity(const QDirection &dir_) const {
	if (getGasType() == GasType::HI) return getHIColumnDensity(dir_);
	if (getGasType() == GasType::H2) return getH2ColumnDensity(dir_);
	return QColumnDensity(0);
}

RingModel::RingModel(GasType gas)
    : NeutralGasAbstract(), dataPtr(std::make_shared<RingData>(RingData(gas))) {
	std::fill(enabledRings.begin(), enabledRings.end(),
	          true);  // enable all by default
	fillRingContainer();
}

std::array<bool, 12> RingModel::getEnabledRings() const { return enabledRings; }

void RingModel::disableRingNo(int i) {
	if (i >= enabledRings.size()) {
		throw std::runtime_error(
		    "Provided number is bigger than the total "
		    "number of rings. Aborted.");
	}

	enabledRings[i] = false;
}

void RingModel::enableRingNo(int i) {
	if (i >= enabledRings.size()) {
		throw std::runtime_error(
		    "Provided number is bigger than the total "
		    "number of rings. Aborted.");
	}

	enabledRings[i] = true;
}

void RingModel::setEnabledRings(std::array<bool, 12> list) {
	enabledRings = list;
}

bool RingModel::isRingEnabled(int i) const {
	if (i >= enabledRings.size()) return false;
	return enabledRings[i];
}

GasType RingModel::getGasType() const { return dataPtr->getGasType(); }

void RingModel::fillRingContainer() {
	for (std::size_t i = 0; i < dataPtr->getRingNumber(); ++i) {
		ringContainer.push_back(std::make_shared<Ring>(
		    Ring(i, dataPtr, boundaries[i], boundaries[i + 1])));
	}
}

int RingModel::getRingNumber() const { return dataPtr->getRingNumber(); }

std::vector<std::pair<PID, double>> RingModel::getAbundanceFractions() const {
	return abundanceFractions;
}

std::shared_ptr<Ring> RingModel::operator[](const std::size_t i) const {
	return ringContainer[i];
}

std::size_t RingModel::size() const { return ringContainer.size(); }

RingModel::iterator RingModel::begin() { return ringContainer.begin(); }

RingModel::const_iterator RingModel::begin() const {
	return ringContainer.begin();
}

RingModel::iterator RingModel::end() { return ringContainer.end(); }

RingModel::const_iterator RingModel::end() const { return ringContainer.end(); }

}}  // namespace hermes::neutralgas

#endif  // HERMES_HAVE_CFITSIO
