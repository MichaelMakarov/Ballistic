#pragma once
#include <vector>

namespace ball
{
	namespace math
	{
		class LegendrePolynomial
		{
		protected:
			size_t _degree;
			std::vector<double> _values;

		private:
			void CreatePolynom(
				std::vector<double>& values, 
				const size_t degree);
			
			bool _normalized;

		public:
			explicit LegendrePolynomial(
				const size_t degree = 0,
				const bool normalized = false);
			~LegendrePolynomial() {}

			size_t Degree() const;

			bool Normalized() const;

			virtual double operator () (const double x) const;
			double operator [] (const size_t i) const;
		};

		class LegendreFunction : public LegendrePolynomial
		{
		private:
			size_t _derivation;
			bool _normalized;

		public:
			LegendreFunction(
				const size_t degree = 0,
				const size_t derivation = 0,
				const bool normalized = false);
			~LegendreFunction() {}

			size_t Derivation() const;

			bool Normalized() const;

			double operator () (const double x) const override;
		};
	}
}