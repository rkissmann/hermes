#ifdef HERMES_HAVE_HDF5

#ifndef HERMES_PICARD2D_H
#define HERMES_PICARD2D_H

#include <map>
#include <memory>
#include <set>

#include "hermes/FITSWrapper.h"
#include "hermes/cosmicrays/Hdf5Reader.h"
#include "hermes/cosmicrays/CosmicRayDensity.h"

namespace hermes { namespace cosmicrays {
/**
 * \addtogroup CosmicRays
 * @{
 */

class Picard2D : public CosmicRayDensity {
  private:
	std::string dirname;
	std::unique_ptr<Hdf5Reader> h5file;

	void readFile();
	void readEnergyAxis();
	void readSpatialGrid2D();
	void readDensity2D();
	std::string getDsetName(std::size_t iE);
	std::size_t calcArrayIndex2D(std::size_t ir, std::size_t iz);

	QLength rmin, rmax, zmin, zmax;
	int dimE;
	int dimz, dimr;
	std::vector<std::unique_ptr<ScalarGrid2DQPDensityPerEnergy>> grid;

	// TODO: implement as std::unordered_map
	std::map<QEnergy, std::size_t> energyIndex;

  public:
	Picard2D(const PID &pid);
	Picard2D(const std::vector<PID> &pids);
	Picard2D(const std::string &filename, const PID &pid_);
	Picard2D(const std::string &filename, const std::vector<PID> &pids);
	QPDensityPerEnergy getDensityPerEnergy(
	    const QEnergy &E_, const Vector3QLength &pos_) const override;
	QPDensityPerEnergy getDensityPerEnergy(int iE_,
	                                       const Vector3QLength &pos_) const;
};

/** @}*/
}}  // namespace hermes::cosmicrays

#endif  // HERMES_DRAGON2D_H

#endif  // HERMES_HAVE_HDF5
