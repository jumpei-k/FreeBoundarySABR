#pragma once
#define _USE_MATH_DEFINES
#include "IHeatKernel.h"
#include <math.h>
#include <iostream>

namespace ktlib {
	class NonCorrelatedSabrIntegralKernel : public IHeatKernel {
	public:
		NonCorrelatedSabrIntegralKernel(double requirement);
		double operator() (double time, double space) const;
	private:
		double calculateIntegrant(
			double integrand, double time, double space) const;
		double calculateIntegrationRange(
			double time) const;
		double executeIntegration(double time, double space) const;
		const double _requirement;
	};

	NonCorrelatedSabrIntegralKernel::NonCorrelatedSabrIntegralKernel(double requirement)
		:_requirement(requirement)
	{
	}

	double NonCorrelatedSabrIntegralKernel::operator()(
		double time, double space) const
	{
		const double integration = executeIntegration(time, space);
		const double preFactor = 2.0 * sqrt(2.0) * exp(-0.125 * time)
			/ (time * sqrt(2.0 * M_PI * time));
		return preFactor * integration;
	}
	double NonCorrelatedSabrIntegralKernel::calculateIntegrant(
		double integrand, double time, double space) const
	{
		return integrand * exp(-0.5 * integrand / time)
			* sqrt(cosh(integrand) - cosh(space));
	}
	double NonCorrelatedSabrIntegralKernel::calculateIntegrationRange(
		double time) const
	{
		/// for now, adopt 3 sigma
		return 3.0 * sqrt(time);
	}
	double NonCorrelatedSabrIntegralKernel::executeIntegration(
		double time, double space) const
	{
		const double end = calculateIntegrationRange(time);
		if (space > end) {
			return 0.0;
		}
		const std::size_t numberOfPertition = (end - space) / _requirement;
		double integration = 0.0;
		integration = integration + 0.5 * calculateIntegrant(space + _requirement);
		for (std::size_t index = 1; index < numberOfPertition - 1; ++index){
			integration = integration + calculateIntegrant(
				space + _requirement * index, time, space);
		}
	}

}