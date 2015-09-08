#pragma once
#define _USE_MATH_DEFINES
#include <boost/math/special_functions/asinh.hpp>
#include <assert.h>
#include <math.h>
#include "IHeatKernel.h"

namespace ktlib {
	class ExactCevPrice {
	public:
		ExactCevPrice(const IHeatKernel& heatKernel, double requirement);
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
		double calculateA1Integrand(
			double integrant,
			double initialVol,
			double volVol,
			double maturity,
			double eta,
			double qBar,
			double b) const;
		double calculateA2Integrand(
			double integrant,
			double initialVol,
			double volVol,
			double strike,
			double maturity,
			double eta,
			double qBar,
			double b) const;
		double calculateA1Integral(
			double initialVol,
			double volVol,
			double maturity,
			double eta,
			double qBar,
			double b) const;
		double calculateA2Integral(
			double initialVol,
			double volVol,
			double strike,
			double maturity,
			double eta,
			double qBar,
			double b) const;
		double calculateIntegrationRange(
			double qBar,
			double maturity) const;
		double calculateQBar(
			double initialForward,
			double strike,
			double skew) const;
		double calculateEta(double skew) const;
		double calculateSPhi(
			double phi,
			double initialVol,
			double volVol,
			double qBar,
			double b) const;
		double calculateSPsi(
			double psi,
			double initialVol,
			double volVol,
			double qBar,
			double b) const;
		double calculateB(
			double initialForward,
			double strike,
			double skew) const;
		std::size_t calculateNumberOfIntegralPertition(
			double startPoint,
			double endtPoint) const;
		const IHeatKernel& _heatKernel;
		const double _requirement;
	};

	ExactCevPrice::ExactCevPrice(
		const IHeatKernel& heatKernel, double requirement)
		:_heatKernel(heatKernel), _requirement(requirement)
	{
		assert(_requirement > 0.0);
	}

	double ExactCevPrice::operator() (
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

	double ExactCevPrice::calculateInterinsicValue(
		double initialForward, double strike) const
	{
		return (initialForward - strike) > 0.0 ? (initialForward - strike) : 0.0;
	}

	double ExactCevPrice::calculateTimeValue(
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
		const double a1 = calculateA1Integral(
			initialVolatility, volVol, maturity, eta, qBar, b);
		const double a2 = calculateA2Integral(
			initialVolatility, volVol, strike, maturity, eta, qBar, b);
		const double preFactor = std::sqrt(std::abs(strike * initialForward)) / M_PI;
		return preFactor *
			(indicatorFunction(strike >= 0.0) * a1 + std::sin(eta * M_PI) * a2);
	}

	double ExactCevPrice::indicatorFunction(bool condition) const
	{
		return condition ? 1.0 : 0.0;
	}

	double ExactCevPrice::calculateA1Integrand(
		double integrant,
		double initialVol,
		double volVol,
		double maturity,
		double eta,
		double qBar,
		double b) const
	{
		const double stochasticVolTime = maturity * volVol * volVol;
		const double sPhi = calculateSPhi(integrant, initialVol, volVol, qBar, b);
		const double numerator = std::sin(integrant) * std::sin(eta * integrant)
			* _heatKernel(stochasticVolTime, sPhi);
		const double denominator = (b - std::cos(integrant)) * std::cosh(sPhi);
		return numerator / denominator;
	}
	double ExactCevPrice::calculateA2Integrand(
		double integrant,
		double initialVol,
		double volVol,
		double strike,
		double maturity,
		double eta,
		double qBar,
		double b) const
	{
		const double stochasticVolTime = maturity * volVol * volVol;
		const double sPsi = calculateSPsi(integrant, initialVol, volVol, qBar, b);
		const double numerator = std::sinh(integrant) *
			(std::cosh(eta * integrant) * indicatorFunction(strike >= 0.0) +
				std::sinh(eta * integrant) * indicatorFunction(strike < 0.0))
			* _heatKernel(stochasticVolTime, sPsi);
		const double denominator = (b + std::cosh(integrant)) * std::cosh(sPsi);
		return numerator / denominator;
	}
	double ExactCevPrice::calculateA1Integral(
		double initialVol,
		double volVol,
		double maturity,
		double eta,
		double qBar,
		double b) const
	{
		const double start = 0.0;
		const double end = M_PI;
		const std::size_t numberOfPertition
			= static_cast<std::size_t>((end - start) / _requirement);
		assert(numberOfPertition > 3);
		const double dx = (end - start) / static_cast<double>(numberOfPertition);
		double total = 0.0;
		// use trapezoidal rule with homogeneous pertition
		////evaluate start point
		double valueAtIntegrationPoint = calculateA1Integrand(start + dx,
			initialVol, volVol, maturity, eta, qBar, b);
		total = total + 0.5 * valueAtIntegrationPoint * dx;
		////integrate
		for (std::size_t index = 2; index < (numberOfPertition - 1); ++index) {
			valueAtIntegrationPoint = calculateA1Integrand(start + (dx * index),
				initialVol, volVol, maturity, eta, qBar, b);
			total = total + valueAtIntegrationPoint * dx;
		}
		////evaluate start point
		valueAtIntegrationPoint = calculateA1Integrand(end - dx,
			initialVol, volVol, maturity, eta, qBar, b);
		total = total + 0.5 * valueAtIntegrationPoint * dx;
		return total;
	}
	double ExactCevPrice::calculateA2Integral(
		double initialVol,
		double volVol,
		double strike,
		double maturity,
		double eta,
		double qBar,
		double b) const
	{
		const double start = 0.0;
		const double safetyFactor = 1.5;
		const double end
			= safetyFactor * calculateIntegrationRange(qBar, maturity);
		const std::size_t numberOfPertition
			= static_cast<std::size_t>((end - start) / _requirement);
		assert(numberOfPertition > 3);
		const double dx = (end - start) / static_cast<double>(numberOfPertition);
		double total = 0.0;
		// use trapezoidal rule with homogeneous pertition
		////evaluate start point
		double valueAtIntegrationPoint = calculateA2Integrand(start + dx,
			initialVol, volVol, strike, maturity, eta, qBar, b);
		total = total + 0.5 * valueAtIntegrationPoint * dx;
		////integrate
		for (std::size_t index = 2; index < (numberOfPertition - 1); ++index) {
			valueAtIntegrationPoint = calculateA2Integrand(start + (dx * index),
				initialVol, volVol, strike, maturity, eta, qBar, b);
			total = total + valueAtIntegrationPoint * dx;
		}
		////evaluate start point
		valueAtIntegrationPoint = calculateA2Integrand(end - dx,
			initialVol, volVol, strike, maturity, eta, qBar, b);
		total = total + 0.5 * valueAtIntegrationPoint * dx;
		return total;
	}
	double ExactCevPrice::calculateIntegrationRange(
		double qBar,
		double maturity) const
	{
		assert(qBar * _requirement < maturity );
		const double range = std::log(-maturity / qBar
			* std::log(qBar / maturity * _requirement));
		assert(range > 0.0);
		return range;
	}
	double ExactCevPrice::calculateQBar(
		double initialForward,
		double strike,
		double skew) const
	{
		return pow(abs(strike * initialForward), 1.0 - skew)
			/ pow(1.0 - skew, 2);
	}
	double ExactCevPrice::calculateEta(double skew) const
	{
		return 0.5 / (1.0 - skew);
	}
	double ExactCevPrice::calculateSPhi(
		double phi,
		double initialVol,
		double volVol,
		double qBar,
		double b) const
	{
		return asinh(volVol / initialVol *
			std::sqrt(2.0 * qBar * (b - std::cos(phi))));
	}
	double ExactCevPrice::calculateSPsi(
		double psi,
		double initialVol,
		double volVol,
		double qBar,
		double b) const
	{
		return asinh(volVol / initialVol *
			std::sqrt(2.0 * qBar * (b + std::cosh(psi))));
	}
	double ExactCevPrice::calculateB(
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
	std::size_t ExactCevPrice::calculateNumberOfIntegralPertition(
		double startPoint,
		double endtPoint) const
	{
		return static_cast<std::size_t>((endtPoint - startPoint) / _requirement);
	}
}