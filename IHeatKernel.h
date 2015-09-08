#pragma once

namespace ktlib {
	class IHeatKernel {
	public:
		virtual double operator() (double time, double space) const = 0;
	};
}