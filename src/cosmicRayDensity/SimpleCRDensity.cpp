#include "hermes/cosmicRayDensity/SimpleCRDensity.h"
#include <iostream>

namespace hermes {

SimpleCRDensity::SimpleCRDensity() : minE(1_GeV), maxE(1e4_GeV), steps(10) {
	makeEnergyRange();
}

SimpleCRDensity::SimpleCRDensity(QEnergy minE_, QEnergy maxE_, int steps_) {
	minE = minE_;
	maxE = maxE_;
	steps = steps_;
	makeEnergyRange();
}

void SimpleCRDensity::makeEnergyRange() {
	QEnergy energy = minE;
	float energyRatio =
		exp(1./static_cast<double>(steps-1) * log(maxE/minE));

	for (int i = 0; i < steps; ++i) {
		energyRange.push_back(energy);
		energy = energy * energyRatio;
	}
}

QPDensityPerEnergy SimpleCRDensity::getDensityPerEnergy(const QEnergy& E_, const Vector3QLength& pos_) const {
	constexpr int alpha = 3;
	auto Phi0 = 0.1 / (1_GeV * 1_cm*1_cm * 1_s * c_light) * 4_pi; 
	auto E0 = 1_GeV;
	auto E_cutoff = 5_TeV;

	QNumber profile = exp(-1.*fabs(pos_.getZ())/1_kpc);
	QPDensityPerEnergy spectrum = Phi0 * pow<-1*alpha>(E_ / E0) * exp(-E_ / E_cutoff);

	//return profile * spectrum;
	return QPDensityPerEnergy(1);
}

} // namespace hermes