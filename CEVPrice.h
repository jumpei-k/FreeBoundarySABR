#pragma once
#define _USE_MATH_DEFINES
#include <boost/math/special_functions/asinh.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <assert.h>
#include <math.h>
#include "IHeatKernel.h"
#include "A1Integrant.h"
#include "A2Integrant.h"
#include "TrapezoidalIntegration.h"

namespace ktlib {
	class CevPrice {
	public:
		CevPrice(const IHeatKernel& heatKernel, double requirement);
		double operator() (
			double initialForward,
			double strike,
			double maturity,
			double volVol,
			double initialVolatility,
			double skew,
			double correlation) const;
	private:
		double calculateInterinsicValue(double initialForward, double strike) const;
		double calculateTimeValue(
			double initialForward,
			double strike,
			double maturity,
			double volVol,
			double initialVolatility,
			double skew,
			double correlation) const;
		double indicatorFunction(bool condition) const;
		double calculateIntegrationRange(
			double qBar,
			double maturity) const;
		double calculateQBar(
			double initialForward,
			double strike,
			double skew) const;
		double calculateEta(double skew) const;
		double calculateB(
			double initialForward,
			double strike,
			double skew) const;
		const IHeatKernel& _heatKernel;
		const double _requirement;
	};

	CevPrice::CevPrice(
		const IHeatKernel& heatKernel, double requirement)
		:_heatKernel(heatKernel), _requirement(requirement)
	{
		assert(_requirement > 0.0);
	}

	double CevPrice::operator() (
		double initialForward,
		double strike,
		double maturity,
		double volVol,
		double initialVolatility,
		double skew,
		double correlation) const
	{
		const double interinsicValue
			= calculateInterinsicValue(initialForward, strike);
		const double timeValue
			= calculateTimeValue(initialForward, strike, maturity,
			volVol, initialVolatility, skew, correlation);
		return interinsicValue + timeValue;
	}

	double CevPrice::calculateInterinsicValue(
		double initialForward, double strike) const
	{
		return (initialForward - strike) > 0.0 ? (initialForward - strike) : 0.0;
	}

	double CevPrice::calculateTimeValue(
		double initialForward,
		double strike,
		double maturity,
		double volVol,
		double initialVolatility,
		double skew,
		double correlation) const
	{
		//Definition
		const double qBar = calculateQBar(initialForward, strike, skew);
		const double b = calculateB(initialForward, strike, skew);
		const double eta = calculateEta(skew);
		const std::size_t maxPertition = 1000000;

		//Integration
		const TrapezoidalIntegration integrator(_requirement, maxPertition);
		const boost::shared_ptr<A1Integrant> a1Integrant
			= boost::make_shared<A1Integrant>(_heatKernel, initialVolatility, volVol,
				maturity, eta, qBar, b);
		const boost::shared_ptr<A2Integrant> a2Integrant
			= boost::make_shared<A2Integrant>(_heatKernel, initialVolatility, volVol,strike,
				maturity, eta, qBar, b);
		const double cutOff = calculateIntegrationRange(qBar, maturity);
		const double a1 = integrator(*a1Integrant, 0.0, M_PI);
		const double a2 = integrator(*a2Integrant, 0.0, cutOff);
		const double preFactor = std::sqrt(std::abs(strike * initialForward)) / M_PI;
		return preFactor *
			(indicatorFunction(strike >= 0.0) * a1 + std::sin(eta * M_PI) * a2);
	}

	double CevPrice::indicatorFunction(bool condition) const
	{
		return condition ? 1.0 : 0.0;
	}

	double CevPrice::calculateIntegrationRange(
		double qBar,
		double maturity) const
	{
		assert(qBar * _requirement < maturity);
		const double range = std::log(-maturity / qBar
			* std::log(qBar / maturity * _requirement));
		assert(range > 0.0);
		return range;
	}
	double CevPrice::calculateQBar(
		double initialForward,
		double strike,
		double skew) const
	{
		return pow(abs(strike * initialForward), 1.0 - skew)
			/ pow(1.0 - skew, 2);
	}
	double CevPrice::calculateEta(double skew) const
	{
		return 0.5 / (1.0 - skew);
	}
	double CevPrice::calculateB(
		double initialForward,
		double strike,
		double skew) const
	{
		//divergence at initiailForward = 0 and strike = 0
		assert(abs(initialForward) > 1.0e-10);
		assert(abs(strike) > 1.0e-10);
		//
		const double absoluteForward = std::abs(initialForward);
		const double absoluteStrike = std::abs(strike);
		const double numerator = std::pow(absoluteForward, 2.0 * (1.0 - skew))
			+ std::pow(absoluteStrike, 2.0 * (1.0 - skew));
		const double denominator
			= 2.0 * std::pow((absoluteForward * absoluteStrike), 1.0 - skew);
		return numerator / denominator;
	}
}