#ifndef HERMES_SYNCHROINTEGRATOR_H
#define HERMES_SYNCHROINTEGRATOR_H

#include "hermes/Units.h"
#include "hermes/skymaps/SynchroSkymap.h"
#include "hermes/integrators/Integrator.h"
#include "hermes/magneticField/MagneticField.h"
#include "hermes/gasDensity/GasDensity.h"

#include <memory>
#include <array>

namespace hermes {

class SynchroIntegrator: public Integrator<SynchroSkymap> {
private:
	std::shared_ptr<MagneticField> mfield;
	std::shared_ptr<GasDensity> gdensity;
public:
	SynchroIntegrator(const std::shared_ptr<MagneticField> mfield, const std::shared_ptr<GasDensity> gdensity);
	~SynchroIntegrator();

	void set_output();
	SynchroSkymap::tPixel integral(QDirection iterdir);
};

} // namespace hermes

#endif // HERMES_SYNCHROINTEGRATOR_H
