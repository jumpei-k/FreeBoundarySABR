#pragma once
#include "IIntegrant.h"
class IntegrantMock : public IIntegrant {
public:
	IntegrantMock();
	double operator() (double integrand) const;
};

double IntegrantMock::operator() (double integrand) const
{
	return 3.0 * integrand * integrand;
}

IntegrantMock::IntegrantMock()
{
}