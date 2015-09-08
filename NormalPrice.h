#pragma once
#define _USE_MATH_DEFINES
#include <boost/math/distributions/normal.hpp>
#include <math.h>

namespace ktlib {
	class NormalPrice {
	public:
		NormalPrice(
			double initialForward,
			double strike,
			double maturity);
		double operator() (double impliedVolatility) const;
	private:
		double calculateIntegrationRange(double impliedVolatility) const;
		const double _initialForward;
		const double _strike;
		const double _maturity;
	};
	NormalPrice::NormalPrice(
		double initialForward,
		double strike,
		double maturity)
		:
		_initialForward(initialForward),
		_strike(strike),
		_maturity(maturity)
	{
	}
	double NormalPrice::operator() (double impliedVolatility) const
	{
		const double integrationRange = calculateIntegrationRange(impliedVolatility);
		const boost::math::normal normal(0.0, 1.0);
		return (_initialForward - _strike)
			* boost::math::cdf(normal, integrationRange)
			+ impliedVolatility * sqrt(_maturity)
			* exp(-0.5 * pow(integrationRange, 2)) / sqrt(2.0 * M_PI);
	}
	double NormalPrice::calculateIntegrationRange(double impliedVolatility) const
	{
		return (_initialForward - _strike) / (impliedVolatility * sqrt(_maturity));
	}
}