#pragma once
#include <vector>

namespace ball
{
	namespace math
	{
		class LegendrePolynom
		{
		protected:
			size_t _degree;
			std::vector<double> _values;

		private:
			void CreatePolynom(
				std::vector<double>& values, 
				const size_t degree,
				const size_t size);

			double _norm;

		public:
			explicit LegendrePolynom(const size_t degree = 0);
			~LegendrePolynom() {}

			size_t Degree() const;

			double Norm() const;

			virtual double operator () (const double x) const;
			double operator [] (const size_t i) const;
		};

		class LegendreFunction : public LegendrePolynom
		{
		private:
			size_t _derivation;
			double _norm;

		public:
			LegendreFunction(
				const size_t degree = 0,
				const size_t derivation = 0);
			~LegendreFunction() {}

			size_t Derivation() const;

			double Norm() const;

			double operator () (const double x) const override;
		};
	}
}