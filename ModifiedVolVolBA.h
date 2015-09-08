#pragma once
#include <assert.h>
#include <math.h>

namespace ktlib {
	class ModifiedVolVol {
	public:
		ModifiedVolVol(
			double volVol,
			double initialVolatility,
			double skew,
			double correlation);
		double operator() (double initialForward) const;
	private:
		const double _volVol;
		const double _initialVolatility;
		const double _skew;
		const double _correlation;
	};

	ModifiedVolVol::ModifiedVolVol(
		double volVol,
		double initialVolatility,
		double skew,
		double correlation)
		:
		_volVol(volVol),
		_initialVolatility(initialVolatility),
		_skew(skew),
		_correlation(correlation)
	{
	}
	double ModifiedVolVol::operator() (double initialForward) const
	{
		const double modification = -1.5 * (pow(_volVol * _correlation, 2)
			+ _initialVolatility * _volVol * _correlation *
			(1.0 - _skew) * pow(abs(initialForward), _skew - 1.0));
		return sqrt(pow(_volVol, 2) + modification);
	}
}