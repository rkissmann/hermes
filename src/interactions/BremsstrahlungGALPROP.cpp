#include "hermes/interactions/BremsstrahlungGALPROP.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include <functional>

#include "hermes/Common.h"

// Following units from Koch and Motz, 1959
#define mc_units (m_electron * c_light)
#define mc2_units (m_electron * c_squared)

#define LIMIT 1000
#define EPSINT 1e-5
#define KEYINT 15

namespace hermes { namespace interactions {

BremsstrahlungGALPROP::BremsstrahlungGALPROP()
    : BremsstrahlungAbstract(),
      cachingEnabled(true),
      cache({std::make_unique<CacheStorageCrossSection>(),
             std::make_unique<CacheStorageCrossSection>(),
             std::make_unique<CacheStorageCrossSection>()}) {
	// Initialize caching for all targets
	for (auto &t : allTargets) {
		cache[static_cast<int>(t)]->setFunction(
		    [t, this](QEnergy T_electron, QEnergy E_gamma) {
			    return this->getDiffCrossSectionForTargetDirectly(t, T_electron,
			                                                      E_gamma);
		    });
	}
}

void BremsstrahlungGALPROP::enableCaching() { cachingEnabled = true; };

void BremsstrahlungGALPROP::disableCaching() { cachingEnabled = false; };

QDiffCrossSection BremsstrahlungGALPROP::getDiffCrossSectionForTarget(
    Target t, const QEnergy &T_electron, const QEnergy &E_gamma) const {
	if (cachingEnabled)
		return cache[static_cast<int>(t)]->getValue(T_electron, E_gamma);
	return getDiffCrossSectionForTargetDirectly(t, T_electron, E_gamma);
}

QNumber BremsstrahlungGALPROP::ElwertFactor(const QNumber &beta_i,
                                            const QNumber &beta_f,
                                            int Z) const {
	return beta_i * (1. - exp(-2.0 * pi * Z * alpha_fine / beta_i)) / beta_f *
	       (1. - exp(-2.0 * pi * Z * alpha_fine / beta_f));
}

QNumber BremsstrahlungGALPROP::xiFunc(const QNumber &T_electron_i,
                                      const QNumber &k, int Z, int N) const {
	constexpr QNumber b = 0.07_MeV / mc2_units;
	constexpr QNumber c = 0.33_MeV / mc2_units;
	return 1. + N / pow<2>(Z) * (1. - exp((b - T_electron_i) / (9. * b))) *
	                (1. - 0.3 * exp(-k / c));
}

QNumber BremsstrahlungGALPROP::Phi_u(const QNumber &gamma_i,
                                     const QNumber &gamma_f,
                                     const QNumber &k) const {
	return 4. * (log(2. * gamma_i * gamma_f / k) - 0.5);
}

inline double BremsstrahlungGALPROP::R_1(double q, double Z) const {
	double F_1 =
	    1. / pow<2>(1. + pow<2>(q) /
	                         pow<2>(2. * static_cast<double>(alpha_fine) * Z));
	return 1 - F_1;
}

inline double BremsstrahlungGALPROP::R_2(double q, double Z) const {
	double F_2 =
	    1. /
	    pow<2>(1. + pow<2>(q) / pow<2>(2. * static_cast<double>(alpha_fine) *
	                                   (Z - 5. / 16.)));
	return 2. * (1. - F_2) - (1. - pow<2>(F_2)) / Z;
}

QNumber BremsstrahlungGALPROP::I_Phi_1(const QNumber &delta_, int Z,
                                       int N) const {
	double result, error;
	double delta = static_cast<double>(delta_);

	auto R_N = [this, delta, Z, N](double q) {
		return ((N == 1) ? R_1(q, Z) : R_2(q, Z)) / pow<3>(q) *
		       pow<2>(q - delta);
	};

	gsl_function_pp<decltype(R_N)> Fp(R_N);
	gsl_function *F = static_cast<gsl_function *>(&Fp);

	gsl_integration_workspace *w = gsl_integration_workspace_alloc(LIMIT);
	gsl_integration_qags(F, delta, 1, 0, EPSINT, LIMIT, w, &result, &error);
	gsl_integration_workspace_free(w);

	return result;
}

QNumber BremsstrahlungGALPROP::I_Phi_2(const QNumber &delta_, int Z,
                                       int N) const {
	double result, error;
	double delta = static_cast<double>(delta_);

	auto R_N = [this, delta, Z, N](double q) {
		return ((N == 1) ? R_1(q, Z) : R_2(q, Z)) / pow<4>(q) *
		       (pow<3>(q) - 6. * pow<2>(delta) * q * std::log(q / delta) +
		        3. * pow<2>(delta) * q - 4. * pow<3>(delta));
	};

	gsl_function_pp<decltype(R_N)> Fp(R_N);
	gsl_function *F = static_cast<gsl_function *>(&Fp);

	gsl_integration_workspace *w = gsl_integration_workspace_alloc(LIMIT);
	gsl_integration_qags(F, delta, 1, 0, EPSINT, LIMIT, w, &result, &error);
	gsl_integration_workspace_free(w);

	return result;
}

QNumber BremsstrahlungGALPROP::Phi_1(const QNumber &gamma_i,
                                     const QNumber &gamma_f, const QNumber &k,
                                     const QNumber &delta, int Z, int N) const {
	QNumber I = I_Phi_1(delta, Z, N);
	return pow<2>(Z - N) * Phi_u(gamma_i, gamma_f, k) +
	       8.0 * Z * (1.0_num - QNumber((N - 1.) / Z) + I);
}

QNumber BremsstrahlungGALPROP::Phi_2(const QNumber &gamma_i,
                                     const QNumber &gamma_f, const QNumber &k,
                                     const QNumber &delta, int Z, int N) const {
	QNumber I = I_Phi_2(delta, Z, N);
	return pow<2>(Z - N) * Phi_u(gamma_i, gamma_f, k) +
	       8.0 * Z * (5.0_num / 6. * (1. - (N - 1.) / Z) + I);
}

QArea BremsstrahlungGALPROP::dsdk_LowEnergy(const QNumber &p_i,
                                            const QNumber &p_f, QNumber k,
                                            int Z) const {
	return 16. * pow<2>(Z * r_electron) * alpha_fine / (3. * k * pow<2>(p_i)) *
	       log((p_i + p_f) / (p_i - p_f));
}

QArea BremsstrahlungGALPROP::dsdk_IntermediateEnergy(
    const QNumber &gamma_i, const QNumber &gamma_f, const QNumber &p_i,
    const QNumber &p_f, const QNumber &k, int Z) const {
	QNumber L = 2. * log((gamma_i * gamma_f + p_i * p_f - 1_num) / k);
	QNumber epsilon_i = log((gamma_i + p_i) / (gamma_i - p_i));
	QNumber epsilon_f = log((gamma_f + p_f) / (gamma_f - p_f));

	QNumber ininner_factor =
	    epsilon_i * (gamma_i * gamma_f + pow<2>(p_i)) / pow<3>(p_i);
	ininner_factor -=
	    epsilon_f * (gamma_i * gamma_f + pow<2>(p_f)) / pow<3>(p_f);
	ininner_factor += 2. * k * gamma_i * gamma_f / pow<2>(p_i * p_f);

	QNumber inner_factor = 8. / 3. * gamma_i * gamma_f / p_i / p_f;
	inner_factor += pow<2>(k) *
	                (pow<2>(gamma_i * gamma_f) + pow<2>(p_i * p_f)) /
	                pow<3>(p_i * p_f);
	inner_factor += k / 2. / p_i / p_f * ininner_factor;

	QNumber factor = 4. / 3.;
	factor -= 2. * gamma_i * gamma_f * (pow<2>(p_f) + pow<2>(p_i)) /
	          pow<2>(p_f * p_i);
	factor +=
	    epsilon_i * gamma_f / pow<3>(p_i) + epsilon_f * gamma_i / pow<3>(p_f);
	factor -= epsilon_i * epsilon_f / p_i / p_f;
	factor += L * inner_factor;

	return pow<2>(Z * r_electron) * alpha_fine * (p_f / p_i) / k * factor;
}

QArea BremsstrahlungGALPROP::dsdk_HighEnergy(const QNumber &gamma_i,
                                             const QNumber &gamma_f,
                                             const QNumber &k, int Z,
                                             int N) const {
	QNumber delta = k / 2. / gamma_i / gamma_f;

	QNumber phi_1, phi_2;
	if (N == 0) {
		phi_1 = pow<2>(Z) * Phi_u(gamma_i, gamma_f, k);
		phi_2 = pow<2>(Z) * Phi_u(gamma_i, gamma_f, k);
	} else {
		phi_1 = Phi_1(gamma_i, gamma_f, k, delta, Z, N);
		phi_2 = Phi_2(gamma_i, gamma_f, k, delta, Z, N);
	}

	QNumber factor = (1_num + pow<2>(gamma_f / gamma_i)) * phi_1 -
	                 2. / 3. * gamma_f / gamma_i * phi_2;
	return pow<2>(r_electron) * alpha_fine / k * factor;
}

QDiffCrossSection BremsstrahlungGALPROP::getDiffCrossSectionForTargetDirectly(
    Target t, const QEnergy &T_electron, const QEnergy &E_gamma) const {
	int Z, N;
	if (t == Target::HII) {
		Z = 1;
		N = 0;
	}
	if (t == Target::HI) {
		Z = 1;
		N = 1;
	}
	if (t == Target::He) {
		Z = 2;
		N = 2;
	}

	QNumber k = E_gamma / mc2_units;
	QNumber T_electron_i = T_electron / mc2_units;
	QNumber T_electron_f = T_electron_i - k;

	if (T_electron_f <= 0_num) return QDiffCrossSection(0);
	if (T_electron_i < 0.01_MeV / mc2_units) return QDiffCrossSection(0);

	QNumber E_electron_i = T_electron_i + 1_num;
	QNumber E_electron_f = T_electron_f + 1_num;

	QNumber p_i = sqrt(T_electron_i * (T_electron_i + 2_num));
	QNumber p_f = sqrt(T_electron_f * (T_electron_f + 2_num));

	QNumber gamma_i = E_electron_i;
	QNumber gamma_f = E_electron_f;

	QNumber beta_i = sqrt(1_num - 1. / (gamma_i * gamma_i));
	QNumber beta_f = sqrt(1_num - 1. / (gamma_f * gamma_f));

	if (T_electron_i < 0.07_MeV / mc2_units) {
		return ElwertFactor(beta_i, beta_f, Z) *
		       dsdk_LowEnergy(p_i, p_f, k, Z) / mc2_units;
	}

	if (T_electron_i < 2.0_MeV / mc2_units) {
		return ElwertFactor(beta_i, beta_f, Z) * xiFunc(T_electron_i, k, Z, N) *
		       dsdk_IntermediateEnergy(gamma_i, gamma_f, p_i, p_f, k, Z) /
		       mc2_units;
	}

	return dsdk_HighEnergy(gamma_i, gamma_f, k, Z, N) / mc2_units;
}

}}  // namespace hermes::interactions
