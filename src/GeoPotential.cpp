#include "GeoPotential.h"
#include "Conversions.h"

namespace ball
{
	namespace space
	{
		GeoPotential::GeoPotential(
			const std::shared_ptr<IEarth> pGravity,
			const size_t count)
		{
			_count = count > pGravity->count() ? pGravity->count() : count;
			_eR = pGravity->R();
			_eMu = pGravity->Mu();
			size_t dim = 0;
			for (size_t i = 0; i <= _count; ++i) dim += i + 1;
			_harmonics.resize(dim);
			std::memcpy(_harmonics.data(), pGravity->harmonics().data(), sizeof(std::pair<double, double>) * dim);
		}

		GeoPotential::GeoPotential(const GeoPotential& gp)
		{
			_count = gp._count;
			_eR = gp._eR;
			_eMu = gp._eMu;
			_harmonics = gp._harmonics;
		}

		GeoPotential::GeoPotential(GeoPotential&& gp) noexcept
		{
			_count = gp._count;
			_eR = gp._eR;
			_eMu = gp._eMu;
			_harmonics = std::move(gp._harmonics);
			gp._eR = gp._eMu = gp._count = 0;
		}

		GeoPotential& GeoPotential::operator = (const GeoPotential& gp)
		{
			_count = gp._count;
			_eR = gp._eR;
			_eMu = gp._eMu;
			_harmonics = gp._harmonics;
			return *this;
		}

		GeoPotential& GeoPotential::operator = (GeoPotential&& gp) noexcept
		{
			_count = gp._count;
			_eR = gp._eR;
			_eMu = gp._eMu;
			_harmonics = std::move(gp._harmonics);
			gp._eR = gp._eMu = gp._count = 0;
			return *this;
		}

		double GeoPotential::operator () (const geometry::RBL& coordinates) const
		{
			double result{ 0 };
			const double sinphi{ std::sin(coordinates.B) };
			const double cosphi{ std::cos(coordinates.B) };
			const double R_r{ _eR / coordinates.R };
			const double coslambda{ std::cos(coordinates.L) };
			const double sinlambda{ std::sin(coordinates.L) };
			double mult{ 1 };
			double b;
			size_t k{ 3 };
			thread_local static std::vector<std::pair<double, double>> cs(_count + 1);
			cs[0] = { 1, 0 };
			for (size_t i = 1; i <= _count; ++i)
				cs[i] = {
					cs[i - 1].first * coslambda - cs[i - 1].second * sinlambda,
					cs[i - 1].second * coslambda + cs[i - 1].first * sinlambda
			};
			// legendre functions
			thread_local static auto pnm{ std::vector<double>(_harmonics.size()) };
			pnm[0] = 1;
			pnm[1] = sinphi * std::sqrt(3);
			pnm[2] = cosphi * std::sqrt(3);
			for (size_t n = 2; n <= _count; ++n)
			{
				for (size_t m = 0; m < n; ++m)
				{
					pnm[k] = std::sqrt(2 * n - 1) * sinphi * pnm[k - n] - std::sqrt((n - 1 - m) * (n - 1 + m) / (2.0 * n - 3)) * pnm[(k + 1) - n - n];
					pnm[k++] *= std::sqrt((2.0 * n + 1) / ((n - m) * (n + m)));
				}
				pnm[k++] = std::sqrt(1 + 0.5 / n) * cosphi * pnm[k - n - 1];
			}
			k = 0;
			for (size_t n = 0; n <= _count; ++n)
			{
				for (size_t m = 0; m <= n; ++m)
				{
					b = _harmonics[k].first * cs[m].first + _harmonics[k].second * cs[m].second;
					result += mult * b * pnm[k];
					k++;
				}
				mult *= R_r;
			}
			return _eMu / coordinates.R * result;
		}

		geometry::XYZ GeoPotential::acceleration(const geometry::XYZ& xyzCoord) const
		{
			using namespace geometry;

			// usefull constants related to current position

			const double r{ xyzCoord.length() };
			const double xy{ std::sqrt(xyzCoord.X * xyzCoord.X + xyzCoord.Y * xyzCoord.Y) };
			const double sinphi{ xyzCoord.Z / r };
			const double cosphi{ xy / r };
			const double tgphi{ sinphi / cosphi };
			const double coslambda{ xyzCoord.X / xy };
			const double sinlambda{ xyzCoord.Y / xy };
			const double mu_r_2{ _eMu / r /r };
			const double R_r{ _eR / r };
			double mult{ R_r * R_r };
			const double zxyr{ xyzCoord.Z / xy / r };

			// the derivatives

			const auto xyzdR{ xyzCoord / r };
			const auto xyzdPhi{	XYZ(-xyzCoord.X * zxyr,	-xyzCoord.Y * zxyr,	xy / r) };
			const auto xyzdLambda{ XYZ(-xyzCoord.Y / xy, xyzCoord.X / xy, 0) };

			// the temporary values

			size_t k{ 3 };
			double poly, dpoly;
			double kcs, ksc;
			RBL rbldU_n, rbldU_sum;

			auto delta = [](const size_t m) { return m == 0 ? 0.5 : 1.0; };

			// cosines and sines
			thread_local static auto cs{ std::vector<std::pair<double, double>>(_count + 1) };
			cs[0] = { 1, 0 };
			for (size_t i = 1; i <= _count; ++i)
				cs[i] = {
					cs[i - 1].first * coslambda - cs[i - 1].second * sinlambda,
					cs[i - 1].second * coslambda + cs[i - 1].first * sinlambda 
				};
			// legendre functions
			thread_local static auto pnm{ std::vector<double>(_harmonics.size()) };
			pnm[0] = 1;
			pnm[1] = sinphi * std::sqrt(3);
			pnm[2] = cosphi * std::sqrt(3);
			for (size_t n = 2; n <= _count; ++n)
			{
				for (size_t m = 0; m < n; ++m)
				{
					pnm[k] = std::sqrt(2 * n - 1) * sinphi * pnm[k - n] - std::sqrt((n - 1 - m) * (n - 1 + m) / (2.0 * n - 3)) * pnm[(k + 1) - n - n];
					pnm[k++] *= std::sqrt((2.0 * n + 1) / ((n - m) * (n + m)));
				}
				pnm[k++] = std::sqrt(1 + 0.5 / n) * cosphi * pnm[k - n - 1];
			}
			// calculating the accelerations
			k = 3;
			for (size_t n = 2; n <= _count; ++n)
			{
				rbldU_n.B = rbldU_n.L = rbldU_n.R = 0.0;
				for (size_t m = 0; m <= n; ++m)
				{
					// current Legendre function
					poly = pnm[k];
					// a derivative of the current Legendre function
					dpoly = -poly * m * tgphi + (m == n ? 0 : pnm[k + 1] * std::sqrt((n - m) * (n + m + 1) * delta(m)));
					// Cnm * cos(m * L) + Snm * sin(m * L)
					kcs = _harmonics[k].first * cs[m].first + _harmonics[k].second * cs[m].second;
					// Snm * cos(m * L) - Cnm * sin(m * L)
					ksc = _harmonics[k].second * cs[m].first - _harmonics[k].first * cs[m].second;
					rbldU_n.R -= poly * kcs;
					rbldU_n.B += dpoly * kcs;
					rbldU_n.L += poly * ksc * m;
					k++;
				}
				rbldU_sum.R += (n + 1) * mult * rbldU_n.R;
				rbldU_sum.B += mult * rbldU_n.B;
				rbldU_sum.L += mult * rbldU_n.L;
				mult *= R_r;
			}
			rbldU_sum.L /= cosphi;
			rbldU_sum.R *= mu_r_2;
			rbldU_sum.B *= mu_r_2;
			rbldU_sum.L *= mu_r_2;

			auto centr_pot{ XYZ(-mu_r_2 * xyzdR.X, -mu_r_2 * xyzdR.Y, -mu_r_2 * xyzdR.Z) };
			auto harm_pot = XYZ(
				rbldU_sum.R * xyzdR.X + rbldU_sum.B * xyzdPhi.X + rbldU_sum.L * xyzdLambda.X,
				rbldU_sum.R * xyzdR.Y + rbldU_sum.B * xyzdPhi.Y + rbldU_sum.L * xyzdLambda.Y,
				rbldU_sum.R * xyzdR.Z + rbldU_sum.B * xyzdPhi.Z
			);
			return centr_pot + harm_pot;
		}
	}
}