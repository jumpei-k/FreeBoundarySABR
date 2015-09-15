#pragma once
#define _USE_MATH_DEFINES
#include <boost/math/special_functions/asinh.hpp>
#include <math.h>
#include "IHeatKernel.h"
#include "IIntegrant.h"

namespace ktlib {
	class A2Integrant : public IIntegrant {
	public:
		A2Integrant(const IHeatKernel& heatKernel,
			double initialVolatility,
			double volVol,
			double strike,
			double maturity,
			double eta,
			double qBar,
			double b);
		double operator()(double integrant) const;
	private:
		double indicatorFunction(bool condition) const;
		double calculateSPsi(double psi) const;
		const IHeatKernel& _heatKernel;
		const double _initialVolatility;
		const double _volVol;
		const double _strike;
		const double _maturity;
		const double _eta;
		const double _qBar;
		const double _b;
	};

	A2Integrant::A2Integrant(
		const IHeatKernel& heatKernel,
		double initialVolatility,
		double volVol,
		double strike,
		double maturity,
		double eta,
		double qBar,
		double b) const
		:_heatKernel(heatKernel), _initialVolatility(initialVolatility), _volVol(volVol),
		_strike(strike), _maturity(maturity), _eta(eta), _qBar(qBar), _b(b)
	{
	}

	double A2Integrant::operator() (double integrand) const
	{
		const double stochasticVolTime = _maturity * _volVol * _volVol;
		const double sPsi = calculateSPsi(integrand);
		const double numerator = std::sinh(integrand) *
			(std::cosh(_eta * integrand) * indicatorFunction(_strike >= 0.0) +
			std::sinh(_eta * integrand) * indicatorFunction(_strike < 0.0))
			* _heatKernel(stochasticVolTime, sPsi);
		const double denominator = (_b + std::cosh(integrand)) * std::cosh(sPsi);
		return numerator / denominator;
	}
	double A2Integrant::calculateSPsi(double psi) const
	{
		return asinh(_volVol / _initialVolatility *
			std::sqrt(2.0 * _qBar * (_b + std::cosh(psi))));
	}
	double A2Integrant::indicatorFunction(bool condition) const
	{
		return condition ? 1.0 : 0.0;
	}

}