#ifndef HERMES_CHARGEDGASDENSITY_H
#define HERMES_CHARGEDGASDENSITY_H

#include "hermes/Grid.h"
#include "hermes/Units.h"

namespace hermes { namespace chargedgas {
/**
 * \addtogroup ChargedGas
 * @{
 */

class ChargedGasDensity {
  private:
	QTemperature gasTemp;

  public:
	ChargedGasDensity() : gasTemp(1e4_K) {}
	ChargedGasDensity(QTemperature T) : gasTemp(T) {}
	virtual ~ChargedGasDensity() {}
	virtual QPDensity getDensity(const Vector3QLength &pos) const = 0;

	inline void setTemperature(QTemperature T) { gasTemp = T; }
	inline QTemperature getTemperature() const { return gasTemp; }
};

/** @}*/
}}  // namespace hermes::chargedgas

#endif  // HERMES_CHARGEDGASDENSITY_H
