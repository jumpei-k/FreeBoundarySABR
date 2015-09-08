#pragma once
#include "IHeatKernel.h"
#include <math.h>

namespace ktlib {
	class NonCorrelatedSabrSimpleKernel : public IHeatKernel {
	public:
		NonCorrelatedSabrSimpleKernel();
		double operator() (double time, double space) const;
	private:
		double calculateR(double time, double space) const;
		double calculateDr(double time) const;
		double calculateG(double space) const;
	};

	NonCorrelatedSabrSimpleKernel::NonCorrelatedSabrSimpleKernel()
	{
	}

	double NonCorrelatedSabrSimpleKernel::operator()(
		double time, double space) const
	{
		const double r = calculateR(time, space);
		const double dr = calculateDr(time);
		const double preFactor = sqrt(sinh(space) / space)
			* exp(-0.5 * pow(space, 2) / time - 0.125 * time);
		return preFactor * (r + dr);
	}
	double NonCorrelatedSabrSimpleKernel::calculateR(
		double time, double space) const
	{
		const double g = calculateG(space);
		// order by time
		const double firstOrder = (3.0 * g) / (8.0 * pow(space, 2));
		const double secondOrder
			= -5.0 * (-8.0 * pow(space, 2) + 3.0 * pow(g, 2) + 24.0 * g)
			/ (128.0 * pow(space, 4));
		const double thirdOrder = 35.0 *
			(-40.0 * pow(space, 2) + 3.0 * pow(g, 3) + 24.0 * pow(g, 2) + 120.0 * g)
			/ (1024.0 * pow(space, 6));
		return 1.0 + firstOrder * time + secondOrder * pow(time, 2)
			+ thirdOrder * pow(time, 3);
	}
	double NonCorrelatedSabrSimpleKernel::calculateDr(double time) const
	{
		return exp(0.125 * time)
			- (3072.0 + 384.0 * time + 24.0 * pow(time, 2) + pow(time, 3))
			/ 3072.0;
	}
	double NonCorrelatedSabrSimpleKernel::calculateG(double space) const
	{
		return space / tanh(space) - 1.0;
	}
}