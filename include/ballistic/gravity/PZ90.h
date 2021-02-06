#include "Gravity.h"

namespace ball
{
	namespace space
	{
		class PZ90 : public IGravity
		{
		public:
			size_t Count() const override { return 36; }
			double Mu() const override { return 398600.4418e9; }
			double R() const override { return 6378136; }
			double Ec2() const override { return 0.0066943662; }
			double W() const override { return 7.292115e-5; }
			double Fl() const override { return 1.0 / 298.2564151; }
			const std::vector<std::pair<double, double>>& Harmonics() const { return _harmonics; }

		private:
			const static std::vector<std::pair<double, double>> _harmonics;
		};
	}
}