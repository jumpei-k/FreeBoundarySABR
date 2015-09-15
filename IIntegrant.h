#pragma once
class IIntegrant {
	public:
		virtual double operator() (double integrand) const = 0;
};