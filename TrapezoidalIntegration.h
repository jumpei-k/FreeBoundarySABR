#pragma once
#include "IIntegrant.h"
#include <iostream>

namespace ktlib {
	class TrapezoidalIntegration {
	public:
		TrapezoidalIntegration(double accuracy, std::size_t maxPertition);
		double operator()(const IIntegrant& integrant, double start, double end) const;
	private:
		double calculateOddsum(const IIntegrant& integrant, double start, double end,
			std::size_t numberOfPertition) const;
		bool isSurtisfeiedAccuracyCondition(double integration, double oddSum) const;
		bool isSurtisfeiedPertitionCondition(std::size_t numberOfPertition) const;
		const double _accuracy;
		const std::size_t _maxPertition;
	};
}
ktlib::TrapezoidalIntegration::TrapezoidalIntegration(double accuracy, std::size_t maxPertition)
	:_accuracy(accuracy), _maxPertition(maxPertition)
{
}

double ktlib::TrapezoidalIntegration::operator()(
	const IIntegrant& integrant, double start, double end) const
{
	std::size_t numberOfPertition = 2;
	double integration = 0.0;
	double oddSum = calculateOddsum(integrant, start, end, numberOfPertition);
	while (true)
	{
		if (isSurtisfeiedAccuracyCondition(integration, oddSum)
			|| isSurtisfeiedPertitionCondition(numberOfPertition)) {
			return 0.5 *integration + oddSum;
		}
		integration = 0.5 * integration + oddSum;
		numberOfPertition *= 2;
		oddSum = calculateOddsum(integrant, start, end, numberOfPertition);
	} 
}

double ktlib::TrapezoidalIntegration::calculateOddsum(const IIntegrant& integrant,
	double start, double end, std::size_t numberOfPertition) const
{
	const double dx = (end - start) / static_cast<double>(numberOfPertition);
	double sum = 0.0;
	for (std::size_t index = 0; index < numberOfPertition / 2; ++index) {
		sum += dx * integrant.operator()(start + static_cast<double>(2 * index + 1) * dx);
	}
	return sum;
}

bool ktlib::TrapezoidalIntegration::isSurtisfeiedAccuracyCondition(
	double integration, double oddSum) const
{
	return abs((0.5 * integration + oddSum) * _accuracy) > abs(oddSum);
}

bool ktlib::TrapezoidalIntegration::isSurtisfeiedPertitionCondition(
	std::size_t numberOfPertition) const
{
	return numberOfPertition >= _maxPertition;
}