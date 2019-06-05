#ifndef HERMES_RADIOSKYMAPRANGE_H
#define HERMES_RADIOSKYMAPRANGE_H

#include <hermes/skymaps/RadioSkymap.h>

namespace hermes {

class RadioSkymapRange {
private:
	typedef std::vector<RadioSkymap> tSkymapsContainer;
	tSkymapsContainer skymaps;
	std::vector<QFrequency> freqs;
	QFrequency minFreq, maxFreq;
	std::size_t nside;
	int freqSteps;
	void initFrequencyRange();
public:
	RadioSkymapRange(std::size_t nside_, QFrequency minFreq_, QFrequency maxFreq_, int freqSteps_);
	~RadioSkymapRange();

	void setIntegrator(std::shared_ptr<IntegratorTemplate<QTemperature> > integrator_);
	void setMask(std::shared_ptr<SkymapMask> mask_);
	void compute();
       
	/** output **/
	void save(std::shared_ptr<Output> output) const;
 
	/** iterator goodies */
        typedef typename tSkymapsContainer::iterator iterator;
        typedef typename tSkymapsContainer::const_iterator const_iterator;
        iterator begin();
        const_iterator begin() const;
        iterator end();
        const_iterator end() const;
};

} // namespace hermes

#endif // HERMES_RADIOSKYMAPRANGE_H