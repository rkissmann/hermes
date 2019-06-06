#ifdef HERMES_HAVE_CFITSIO

#ifndef HERMES_DRAGONCRDENSITY_H
#define HERMES_DRAGONCRDENSITY_H

#include "hermes/cosmicRayDensity/CosmicRayDensity.h"
#include "hermes/FITSWrapper.h"
#include <cassert>
#include <memory>
#include <map>
#include <set>

namespace hermes {

class PID {
protected:
	int Z;
	int A;
	int id;
public:
	PID(int Z_, int A_) : Z(Z_), A(A_) {
		assert(A_ >= 0);
		id = Z * 100 + A;
	}

	int getID() const {
		return id;
	}

	virtual ~PID() {
	}

	bool isNucleus() {
		return (A > 0);
	}

	int atomicNr() const {
		return Z;
	}

	int massNr() const {
		return A;
	}
};

static const PID Electron(-1, 0);
static const PID Positron(1, 0);
static const PID Proton(1, 1);
static const PID Helium(2, 4);

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

	void enablePID(const PID &pid_);
	void disablePID(const PID &pid_);
	bool isPIDEnabled(const PID &pid_) const;
  
	QLength rmin, rmax, zmin, zmax;
	int dimE;
	int dimz, dimr;
	std::vector<std::unique_ptr<ScalarGrid2DQPDensityPerEnergy> > grid;
	std::set<int> listOfPIDs;
	
	//TODO: implement as std::unordered_map
	std::map<QEnergy, std::size_t> energyIndex;
public:
	Dragon2DCRDensity();
	Dragon2DCRDensity(const std::string &filename_,
			const PID &pid_);
	Dragon2DCRDensity(const std::string &filename_,
			const std::vector<PID> &pids_);
	QPDensityPerEnergy getDensityPerEnergy(const QEnergy& E_, const Vector3QLength& pos_) const;
	QPDensityPerEnergy getDensityPerEnergy(int iE_, const Vector3QLength& pos_) const;
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

	void enablePID(const PID &pid_);
	void disablePID(const PID &pid_);
	bool isPIDEnabled(const PID &pid_) const;
  
	QLength rmin, rmax, zmin, zmax;
	QLength xmin, xmax, ymin, ymax;
	int dimE;
	int dimx, dimy, dimz, dimr;
	std::vector<std::unique_ptr<ScalarGridQPDensityPerEnergy> > grid;
	std::set<int> listOfPIDs;
	
	//TODO: implement as std::unordered_map
	std::map<QEnergy, std::size_t> energyIndex;
public:
	Dragon3DCRDensity();
	Dragon3DCRDensity(const std::string &filename_,
			const PID &pid_);
	Dragon3DCRDensity(const std::string &filename_,
			const std::vector<PID> &pids_);
	QPDensityPerEnergy getDensityPerEnergy(const QEnergy& E_, const Vector3QLength& pos_) const;
	QPDensityPerEnergy getDensityPerEnergy(int iE_, const Vector3QLength& pos_) const;
};


} // namespace hermes

#endif // HERMES_DRAGONCRDENSITY_H

#endif // HERMES_HAVE_CFITSIO
