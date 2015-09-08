#pragma once
#define _USE_MATH_DEFINES
#include <assert.h>
#include <math.h>
#include "ModifiedVolVol.h"

namespace ktlib {
	class ModifiedInitialVolatility {
	public:
		ModifiedInitialVolatility(
			double volVol,
			double initialVolatility,
			double skew,
			double correlation);
		double operator() (double initialForward,
			double strike, double maturity) const;
	private:
		bool isAtTheMoney(double initialForward,
			double strike, double maturity) const;
		double caluculateExceptionAtTheMoney(double initialForward,
			double strike, double maturity) const;
		double caluculateQ(
			double underlying) const;
		double caluculateDq(
			double initialForward,
			double strike) const;
		double caluculateFloaredStrike(
			double initialForward,
			double strike) const;
		double caluculateMinimumInitialVolatility(
			double initialForward,
			double strike) const;
		double caluculateCapitalPhi(
			double initialForward,
			double strike) const;
		double caluculateLeadingVolatility(
			double initialForward,
			double strike) const;
		double caluculateMinimumLeadingVolatility(
			double leadingVolatility,
			double initialForward,
			double strike) const;
		double caluculateAnglePhi(
			double dq,
			double minimumInitialVolatility) const;
		double caluculateL(
			double strike,
			double minimumInitialVolatility) const;
		double caluculateU(
			double dq,
			double minimumInitialVolatility) const;
		double caluculateAngleI(
			double l, double u) const;
		double caluculateMinimumB(
			double angleI, double anglePhi) const;
		double caluculatePerterbationRatio(
			double leadingVolatility,
			double initialForward,
			double strike) const;
		const double _volVol;
		const double _initialVolatility;
		const double _skew;
		const double _correlation;
		const ModifiedVolVol _modifiedVolVol;
	};
	ModifiedInitialVolatility::ModifiedInitialVolatility(
		double volVol,
		double initialVolatility,
		double skew,
		double correlation)
		:
		_volVol(volVol),
		_initialVolatility(initialVolatility),
		_skew(skew),
		_correlation(correlation),
		_modifiedVolVol(
			volVol,
			initialVolatility,
			skew,
			correlation)
	{
	}
	double ModifiedInitialVolatility::operator() (
		double initialForward, double strike, double maturity) const
	{
		if (isAtTheMoney(initialForward, strike, maturity)) {
			return caluculateExceptionAtTheMoney(
				initialForward, strike, maturity);
		}
		const double leadingVolatility = caluculateLeadingVolatility(
			initialForward, strike);
		const double perturbationRatio = caluculatePerterbationRatio(
			leadingVolatility, initialForward, strike);
		return leadingVolatility * (1.0 + perturbationRatio * maturity);
	}
	bool ModifiedInitialVolatility::isAtTheMoney(double initialForward,
		double strike, double maturity) const
	{
		// 0.01 is determined by numerical experiment.
		return abs(initialForward - strike) < 0.01 * sqrt(maturity);
	}
	double ModifiedInitialVolatility::caluculateExceptionAtTheMoney(
		double initialForward, double strike, double maturity) const
	{
		/// koko kara kaku
		const double perturbationRatioVolVolTerm
			= (1.0 - pow(_modifiedVolVol(initialForward) / _volVol, 2)
					- 1.5 * pow(_correlation, 2))
				* pow(_volVol, 2) / 12.0;
		const double perturbationRatioInitialVolTerm = 0.25 * _skew * _correlation
			* _volVol * (abs(initialForward), _skew - 1.0) * _initialVolatility;
		const double perturbationRatio
			= perturbationRatioVolVolTerm + perturbationRatioInitialVolTerm;
		return _initialVolatility * (1.0 + perturbationRatio * maturity);
	}
	double ModifiedInitialVolatility::caluculateQ(double underlying) const
	{
		return pow(abs(underlying), 1.0 - _skew) / (1.0 - _skew);
	}
	double ModifiedInitialVolatility::caluculateDq(
		double initialForward, double strike) const
	{
		const double floardStrike
			= caluculateFloaredStrike(initialForward, strike);
		return caluculateQ(floardStrike) - caluculateQ(initialForward);
	}
	double ModifiedInitialVolatility::caluculateFloaredStrike(
		double initialForward, double strike) const
	{
		return strike > 0.1 * initialForward ? strike : 0.1 * initialForward;
	}
	double ModifiedInitialVolatility::caluculateMinimumInitialVolatility(
		double initialForward, double strike) const
	{
		const double dq = caluculateDq(initialForward, strike);
		const double modification = pow(_volVol * dq, 2)
			+ 2.0 * _correlation * _volVol * dq * _initialVolatility;
		return sqrt(pow(_initialVolatility, 2) + modification);
	}
	double ModifiedInitialVolatility::caluculateCapitalPhi(
		double initialForward, double strike) const
	{
		const double volVolRatio = _modifiedVolVol(initialForward) / _volVol;
		const double dq = caluculateDq(initialForward, strike);
		const double minimum
			= caluculateMinimumInitialVolatility(initialForward, strike);
		const double numerator
			= minimum + _correlation * _initialVolatility + _volVol * dq;
		const double denominator
			= (1.0 + _correlation) * _initialVolatility;
		return pow(numerator / denominator, volVolRatio);
	}
	double ModifiedInitialVolatility::caluculateLeadingVolatility(
		double initialForward, double strike) const
	{
		const double capitalPhi = caluculateCapitalPhi(initialForward, strike);
		const double dq = caluculateDq(initialForward, strike);
		const double numerator
			= 2.0 * capitalPhi * dq * _modifiedVolVol(initialForward);
		const double denominator
			= pow(capitalPhi, 2) - 1.0;
		return numerator / denominator;
	}
	double ModifiedInitialVolatility::caluculateMinimumLeadingVolatility(
		double leadingVolatility,
		double initialForward,
		double strike) const
	{
		const double dq = caluculateDq(initialForward, strike);
		const double modification = _modifiedVolVol(initialForward) * dq;
		return sqrt(pow(leadingVolatility, 2) + pow(modification, 2));
	}
	double ModifiedInitialVolatility::caluculateAnglePhi(
		double dq, double minimumInitialVolatility) const
	{
		const double value = -(dq * _volVol + _initialVolatility * _correlation)
			/ minimumInitialVolatility;
		return acos(value);
	}
	double ModifiedInitialVolatility::caluculateL(
		double strike,
		double minimumInitialVolatility) const
	{
		const double q = caluculateQ(strike);
		return minimumInitialVolatility
			/ (q * _volVol * sqrt(1.0 - pow(_correlation, 2)));
	}
	double ModifiedInitialVolatility::caluculateU(
		double dq,
		double minimumInitialVolatility) const
	{
		const double numerator = dq * _volVol * _correlation
			+ _initialVolatility - minimumInitialVolatility;
		const double denominator = dq * _volVol *  sqrt(1.0 - pow(_correlation, 2));
		return numerator / denominator;
	}
	double ModifiedInitialVolatility::caluculateAngleI(
		double l, double u) const
	{
		const double lConjugate = sqrt(abs(1.0 - pow(l, 2)));
		//Iangel diverges at L = 1
		assert(abs(l - 1.0) > 1.0e-6);
		return l < 1.0
			? 2.0 / lConjugate *
				(atan((u + l) / lConjugate) - atan(l / lConjugate))
			: log((u * (l + lConjugate) + 1.0) / (u * (l - lConjugate) + 1.0))
				/ lConjugate;
	}
	double ModifiedInitialVolatility::caluculateMinimumB(
		double angleI, double anglePhi) const
	{
		const double preFactor = -0.5 * (_skew * _correlation)
			/ ((1.0 - _skew) * sqrt(1.0 - pow(_correlation, 2)));
		const double angleTerm = M_PI - anglePhi - acos(_correlation) - angleI;
		return preFactor * angleTerm;
	}
	double ModifiedInitialVolatility::caluculatePerterbationRatio(
		double leadingVolatility,
		double initialForward,
		double strike) const
	{
		const double dq = caluculateDq(initialForward, strike);
		const double minimumInitialVolatility
			= caluculateMinimumInitialVolatility(initialForward, strike);
		const double minimumLeading = caluculateMinimumLeadingVolatility(
			leadingVolatility, initialForward, strike);
		const double u = caluculateU(dq, minimumInitialVolatility);
		const double l = caluculateL(strike, minimumInitialVolatility);
		const double angleI = caluculateAngleI(l, u);
		const double anglePhi = caluculateAnglePhi(dq, minimumInitialVolatility);
		const double minimumB = caluculateMinimumB(angleI, anglePhi);
		const double volRatio = (_initialVolatility * minimumInitialVolatility)
			/ (leadingVolatility * minimumLeading);
		const double modifiedVolVol = _modifiedVolVol(initialForward);
		const double r = dq * modifiedVolVol / leadingVolatility;
		const double rConjugate = sqrt(1.0 + pow(r, 2));
		const double numerator = pow(modifiedVolVol, 2) * rConjugate
			* (0.5 * log(volRatio) - minimumB);
		const double denominator = r * log(rConjugate + r);
		return numerator / denominator;
	}

}