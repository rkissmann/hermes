#ifdef HERMES_HAVE_CFITSIO

#ifndef HERMES_DRAGONCRDENSITY_H
#define HERMES_DRAGONCRDENSITY_H

#include "hermes/cosmicrays/CosmicRayDensity.h"
#include "hermes/FITSWrapper.h"

#include <memory>
#include <map>
#include <set>

namespace hermes { namespace cosmicrays {

class Dragon2DCRDensity: public CosmicRayDensity {
private:
  	std::string filename;
	std::unique_ptr<FITSFile> ffile;

	void readFile();
	void readEnergyAxis();
	void readSpatialGrid2D();
	void readDensity2D();
	std::size_t calcArrayIndex2D(
		std::size_t iE, std::size_t ir, std::size_t iz);

	QLength rmin, rmax, zmin, zmax;
	int dimE;
	int dimz, dimr;
	std::vector<std::unique_ptr<ScalarGrid2DQPDensityPerEnergy> > grid;
	
	//TODO: implement as std::unordered_map
	std::map<QEnergy, std::size_t> energyIndex;
public:
	Dragon2DCRDensity(const PID& pid);
	Dragon2DCRDensity(const std::vector<PID> &pids);
	Dragon2DCRDensity(const std::string &filename,
			const PID &pid_);
	Dragon2DCRDensity(const std::string &filename,
			const std::vector<PID> &pids);
	QPDensityPerEnergy getDensityPerEnergy(const QEnergy& E_,
			const Vector3QLength& pos_) const override;
	QPDensityPerEnergy getDensityPerEnergy(int iE_,
			const Vector3QLength& pos_) const;
};

class Dragon3DCRDensity: public CosmicRayDensity {
private:
  	std::string filename;
	std::unique_ptr<FITSFile> ffile;

	void readFile();
	void readEnergyAxis();
	void readSpatialGrid3D();
	void readDensity3D();
	std::size_t calcArrayIndex2D(
		std::size_t iE, std::size_t ir, std::size_t iz);

	QLength rmin, rmax, zmin, zmax;
	QLength xmin, xmax, ymin, ymax;
	int dimE;
	int dimx, dimy, dimz, dimr;
	std::vector<std::unique_ptr<ScalarGridQPDensityPerEnergy> > grid;
	
	//TODO: implement as std::unordered_map
	std::map<QEnergy, std::size_t> energyIndex;
public:
	Dragon3DCRDensity();
	Dragon3DCRDensity(const std::string &filename_,
			const PID &pid_);
	Dragon3DCRDensity(const std::string &filename_,
			const std::vector<PID> &pids_);
	QPDensityPerEnergy getDensityPerEnergy(const QEnergy& E_,
			const Vector3QLength& pos_) const override;
	QPDensityPerEnergy getDensityPerEnergy(int iE_, const Vector3QLength& pos_) const;
};

} // namespace cosmicrays
} // namespace hermes

#endif // HERMES_DRAGONCRDENSITY_H

#endif // HERMES_HAVE_CFITSIO