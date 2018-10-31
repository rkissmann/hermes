#include "hermes/example.h"

#include <iostream>
#include <memory>

namespace hermes {

class TestGasDensity: public GasDensity {
public:
        TestGasDensity() { }
        QPDensity getDensity(const Vector3QLength &pos) const {
		return 1.0 / 1_cm3; 
	}
};

void playground() {

	auto B = Vector3QMField(1_muG, 0_muG, 0_muG);
	auto ptr_ufield = std::make_shared<UniformMagneticField>(UniformMagneticField(B));

	int nside = 32;	
	auto ptr_skymap = std::make_shared<RMSkymap>(RMSkymap(nside));
	auto ptr_JF12 = std::make_shared<JF12Field>(JF12Field());
	auto ptr_Gas = std::make_shared<HII_Cordes91>(HII_Cordes91());
	//auto ptr_Gas = std::make_shared<TestGasDensity>(TestGasDensity());
	auto ptr_output = std::make_shared<FITSOutput>(FITSOutput("!example.fits.gz"));

	RMIntegrator RM = RMIntegrator(ptr_JF12, ptr_Gas);
	RM.set_skymap(ptr_skymap);
	RM.compute();

	ptr_output->save(
		ptr_skymap
	);

	
}

} // namespace hermes

int main(void){

	hermes::playground();

	return 0;
}

