#pragma once
#include <assert.h>
#include "IHeatKernel.h"
#include "ExactCevPrice.h"
#include "ModifiedInitialVolatility.h"
#include "ModifiedVolVol.h"

namespace ktlib {
	class CorrelationMappingSabrPrice {
	public:
		CorrelationMappingSabrPrice(
			const IHeatKernel& heatKernel, double requirement);
		double operator() (
			double initialForward,
			double strike,
			double maturity,
			double volVol,
			double initialVolatility,
			double skew,
			double correlation);
	private:
		const ExactCevPrice _exactCev;
	};

	CorrelationMappingSabrPrice::CorrelationMappingSabrPrice(
		const IHeatKernel& heatKernel, double requirement)
		:_exactCev(heatKernel, requirement)
	{
		assert(requirement > 0.0);
	}

	double CorrelationMappingSabrPrice::operator() (
		double initialForward,
		double strike,
		double maturity,
		double volVol,
		double initialVolatility,
		double skew,
		double correlation)
	{
		ModifiedInitialVolatility modifiedInitialVol(
			volVol, initialVolatility, skew, correlation);
		const double modInitial = modifiedInitialVol(
			initialForward, strike, maturity);
		ModifiedVolVol modifiedVolVol(
			volVol, initialVolatility, skew, correlation);
		const double modVol = modifiedVolVol(initialForward);
		return _exactCev(initialForward, strike, maturity,
			modifiedVolVol(initialForward),
			modifiedInitialVol(initialForward, strike, maturity),
			skew, correlation);
	}
}
