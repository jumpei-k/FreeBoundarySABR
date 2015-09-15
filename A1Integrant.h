#pragma once
#include <boost/math/special_functions/asinh.hpp>
#include <math.h>
#include "IHeatKernel.h"
#include "IIntegrant.h"

namespace ktlib {
	class A1Integrant : public IIntegrant {
	public:
		A1Integrant(const IHeatKernel& heatKernel,
			double initialVolatility,
			double volVol,
			double maturity,
			double eta,
			double qBar,
			double b);
		double operator() (double integrant) const;
	private:
		double calculateSPhi(double phi) const;
		const IHeatKernel& _heatKernel;
		const double _initialVolatility;
		const double _volVol;
		const double _maturity;
		const double _eta;
		const double _qBar;
		const double _b;
	};
	A1Integrant::A1Integrant(const IHeatKernel& heatKernel,
		double initialVolatility,
		double volVol,
		double maturity,
		double eta,
		double qBar,
		double b)
		:_heatKernel(heatKernel), _initialVolatility(initialVolatility), _volVol(volVol),
		_maturity(maturity), _eta(eta), _qBar(qBar), _b(b)
	{
	}
	double A1Integrant::operator() (
		double integrant) const
	{
		const double stochasticVolTime = _maturity * _volVol * _volVol;
		const double sPhi = calculateSPhi(integrant);
		const double numerator = std::sin(integrant) * std::sin(_eta * integrant)
			* _heatKernel(stochasticVolTime, sPhi);
		const double denominator = (_b - std::cos(integrant)) * std::cosh(sPhi);
		return numerator / denominator;
	}

	double A1Integrant::calculateSPhi(
		double phi) const
	{
		return asinh(_volVol / _initialVolatility *
			std::sqrt(2.0 * _qBar * (_b - std::cos(phi))));
	}
}