#include "GeoPotential.h"
#include "AstroValues.h"
#include <cmath>
#include <fstream>
#include <iomanip>

namespace ball
{
	namespace tasks
	{
		GeoPotential::GeoPotential(
			const std::unique_ptr<IGravity> pGravity,
			const size_t count)
		{
			_count = count > pGravity->Count() ? pGravity->Count() : count;
			_eR = pGravity->R();
			_eMu = pGravity->Mu();
			size_t dim = 0;
			for (size_t i = 0; i <= _count; ++i) dim += i + 1;
			auto harmonics = pGravity->Harmonics();
			_functions.resize(dim);
			_harmonics.resize(dim);
			dim = 0;
			for (size_t i = 0; i <= _count; ++i)
				for (size_t j = 0; j <= i; ++j)
				{
					_harmonics[dim] = harmonics[dim];
					_functions[dim++] = math::LegendreFunction(i, j, true);
				}
		}

		GeoPotential& GeoPotential::operator = (GeoPotential&& gp) noexcept
		{
			_count = gp._count;
			_eR = gp._eR;
			_eMu = gp._eMu;
			_harmonics = std::move(gp._harmonics);
			_functions = std::move(gp._functions);
			return *this;
		}

		double GeoPotential::operator () (const geometry::RBL& coordinates) const
		{
			double result = 0.0,
				sinphi = std::sin(coordinates.B),
				R_r = _eR / coordinates.R,
				l = coordinates.L,
				mult = 1,
				p;
			size_t k = 0;
			std::vector<std::pair<double, double>> ls(_count + 1);
			for (size_t i = 0; i <= _count; ++i) 
				ls[i] = { std::cos(i * l), std::sin(i * l) };
			for (size_t n = 0; n <= _count; ++n)
			{
				for (size_t m = 0; m <= n; ++m)
				{
					p = _harmonics[k].first * ls[m].first + _harmonics[k].second * ls[m].second;
					result += mult * p * _functions[k](sinphi);
					k++;
				}
				mult *= R_r;
			}
			return _eMu / coordinates.R * result;
		}

		geometry::XYZ GeoPotential::Acceleration(const geometry::XYZ& xyzCoord) const
		{
			using namespace geometry;

			const auto rblCoord{ GCS_OrthoToSpher(xyzCoord) };

			// usefull constants related to current position

			const double sinphi{ std::sin(rblCoord.B) };
			const double cosphi{ std::cos(rblCoord.B) };
			const double tgphi{ sinphi / cosphi };
			const double R_r{ _eR / rblCoord.R };
			const double lambda{ rblCoord.L };
			double mult{ R_r * R_r };

			// the temporary constants

			const double r_2{ rblCoord.R * rblCoord.R };
			const double xy_2{ xyzCoord.X * xyzCoord.X + xyzCoord.Y * xyzCoord.Y };
			const double xy{ std::sqrt(xy_2) };
			const double rxy{ xy * r_2 };

			// the derivatives

			const auto xyzdR{ xyzCoord / rblCoord.R };
			const auto xyzdPhi{
				XYZ(
					-xyzCoord.X * xyzCoord.Z / xy / rblCoord.R,
					-xyzCoord.Y * xyzCoord.Z / xy / rblCoord.R,
					xy / rblCoord.R)
			};
			const auto xyzdLambda{ XYZ(-xyzCoord.Y / xy, xyzCoord.X / xy, 0) };
			/*const auto xyzdPhi{
				XYZ(
					-xyzCoord.X * xyzCoord.Z / rxy,
					-xyzCoord.Y * xyzCoord.Z / rxy,
					std::sqrt(xy_2) / r_2)
			};
			const auto xyzdLambda{ XYZ(-xyzCoord.Y / xy_2, xyzCoord.X / xy_2, 0) };*/
			RBL rbldU_n, rbldU_sum;
			double dPoly;

			// the temporary values

			size_t k{ 3 };
			double poly;
			double kcs, ksc;
			auto cs{ std::vector<std::pair<double, double>>(_count + 1) };
			auto delta = [](const size_t m) { return m == 0 ? 0.5 : 1.0; };

			for (size_t i = 0; i <= _count; ++i)
				cs[i] = { std::cos(i * lambda), std::sin(i * lambda) };
			for (size_t n = 2; n <= _count; ++n)
			{
				rbldU_n.B = rbldU_n.L = rbldU_n.R = 0.0;
				for (size_t m = 0; m <= n; ++m)
				{
					// current Legendre function
					poly = _functions[k](sinphi);
					// a derivative of the current Legendre function
					dPoly = -poly * m * tgphi + (m == n ? 0 : _functions[k + 1](sinphi) * 
								std::sqrt((n - m) * (n + m + 1) * delta(m)));
					// Cnm * cos(m * L) + Snm * sin(m * L)
					kcs = _harmonics[k].first * cs[m].first + _harmonics[k].second * cs[m].second;
					// Snm * cos(m * L) - Cnm * sin(m * L)
					ksc = _harmonics[k].second * cs[m].first - _harmonics[k].first * cs[m].second;
					rbldU_n.R -= poly * kcs;
					rbldU_n.B += dPoly * kcs;
					rbldU_n.L += poly * ksc * m;
					k++;
				}
				rbldU_n.R *= n + 1;
				rbldU_sum += mult * rbldU_n;
				mult *= R_r;
			}
			rbldU_sum.L /= cosphi;
			rbldU_sum *= _eMu / r_2;

			auto centrPot{ XYZ(-_eMu / r_2 * xyzdR.X, -_eMu / r_2 * xyzdR.Y, -_eMu / r_2 * xyzdR.Z) };
			auto harmPot = XYZ(
				rbldU_sum.R * xyzdR.X + rbldU_sum.B * xyzdPhi.X + rbldU_sum.L * xyzdLambda.X,
				rbldU_sum.R * xyzdR.Y + rbldU_sum.B * xyzdPhi.Y + rbldU_sum.L * xyzdLambda.Y,
				rbldU_sum.R * xyzdR.Z + rbldU_sum.B * xyzdPhi.Z
			);
			return centrPot + harmPot;
		}
	}
}