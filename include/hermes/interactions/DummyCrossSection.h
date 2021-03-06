#ifndef HERMES_DUMMYCROSSSECTION_H
#define HERMES_DUMMYCROSSSECTION_H

#include <memory>

#include "hermes/CacheTools.h"
#include "hermes/interactions/DiffCrossSection.h"

namespace hermes { namespace interactions {
/**
 * \addtogroup Interactions
 * @{
 */

class DummyCrossSection : public DifferentialCrossSection {
  private:
	QDiffCrossSection cs;

  public:
	DummyCrossSection(const QDiffCrossSection &cs_)
	    : DifferentialCrossSection(), cs(cs_){};

	QDiffCrossSection getDiffCrossSection(const QEnergy &E_photon,
	                                      const QEnergy &E_gamma) const {
		return QDiffCrossSection(cs);
	}

	QDiffCrossSection getDiffCrossSection(const QEnergy &E_electron,
	                                      const QEnergy &E_photon,
	                                      const QEnergy &E_gamma) const {
		return QDiffCrossSection(cs);
	};
};

/** @}*/
}}  // namespace hermes::interactions

#endif  // HERMES_DUMMYCROSSSECTION_H
