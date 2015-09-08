#pragma once
#include <boost/math/distributions/normal.hpp>
#include <math.h>

namespace ktlib {
	class BlackSholesPrice {
	public:
		BlackSholesPrice(
			double initialForward,
			double strike,
			double maturity);
		double operator() (double impliedVolatility) const;
	private:
		double calculatePositiveD(double impliedVolatility) const;
		double calculateNegativeD(double impliedVolatility) const;
		const double _initialForward;
		const double _strike;
		const double _maturity;
	};
	BlackSholesPrice::BlackSholesPrice(
		double initialForward,
		double strike,
		double maturity)
		:
		_initialForward(initialForward),
		_strike(strike),
		_maturity(maturity)
	{
	}
	double BlackSholesPrice::operator() (double impliedVolatility) const
	{
		const double positiveD = calculatePositiveD(impliedVolatility);
		const double negativeD = calculateNegativeD(impliedVolatility);
		const boost::math::normal normal(0.0, 1.0);
		return _initialForward * boost::math::cdf(normal, positiveD)
			- _strike * boost::math::cdf(normal, negativeD);
	}
	double BlackSholesPrice::calculatePositiveD(double impliedVolatility) const
	{
		return (log(_initialForward / _strike)
				+ 0.5 * pow(impliedVolatility, 2) * _maturity)
			/ (impliedVolatility * sqrt(_maturity));
	}
	double BlackSholesPrice::calculateNegativeD(double impliedVolatility) const
	{
		return (log(_initialForward / _strike)
				- 0.5 * pow(impliedVolatility, 2) * _maturity)
			/ (impliedVolatility * sqrt(_maturity));
	}

}