#include <iostream>
#include <iomanip>
#include <fstream>
#include "general/Mathematics.h"
#include "Integrators.h"
#include "RungeKuttaIntegrator.h"
#include "EverhartIntegrator.h"
#include "AdamsIntegrator.h"
#include "EGM96.h"
#include "general/Mathematics.h"
#include "general/Times.h"

using namespace ball;
using namespace general::math;
using namespace general::time;

void save(const std::string& filename, const std::vector<Vec2>& arr, const std::vector<double>& time)
{
	auto fout{ std::ofstream(filename) };
	for (size_t i = 0; i < arr.size(); ++i) {
		fout << time[i] << " " << arr[i] << std::endl;
	}
	fout.close();
}

void pendulum_integration()
{
	const double dt{ 0.05 };
	constexpr double a{ deg_to_rad(3) };
	const double l{ 0.5 }, g{ 9.81 }, w{ std::sqrt(g / l) };
	const size_t n{ 1000 };
	auto time{ std::vector<double>(n) };
	time[0] = 0;
	for (size_t i = 1; i < n; ++i) time[i] = time[i - 1] + dt;
	double t;

	auto func = [w](const Vec2& x, const double& t) -> Vec2 { return Vec2(x.y(), x.x() * (-w * w)); };
	auto invoker{ make_classfunc<Vec2, double, decltype(func)>(&func) };

	std::cout << "RK\n";
	auto rkangles{ std::vector<Vec2>(n) };
	rkangles[0] = Vec2(a, 0);
	auto rkint{ RKIntegrator<Vec2, double>() };
	//rkint.func = func;
	for (size_t i = 1; i < n; ++i) rkint.integrate(rkangles[i - 1], time[i - 1], dt, rkangles[i], t, invoker);
	save("RK results.txt", rkangles, time);

	std::cout << "Everhart\n";
	auto evangles{ std::vector<Vec2>(n) };
	evangles[0] = Vec2(a, 0);
	auto evint{ EverhartIntegrator<Vec2, double, 4>() };
	//evint.func = func;
	for (size_t i = 1; i < n; ++i) evint.integrate(evangles[i - 1], time[i - 1], dt, evangles[i], t, invoker);
	save("Ev results.txt", evangles, time);

	std::cout << "Analitical\n";
	auto anangles{ std::vector<Vec2>(n) };
	anangles[0] = Vec2(a, 0);
	for (size_t i = 0; i < n; ++i) {
		anangles[i].x() = a * std::cos(w * time[i]);
		anangles[i].y() = -a * w * std::sin(w * time[i]);
	}
	save("analitical.txt", anangles, time);
}

template<Arithmetic X, ball::Time T, class Int, class F>
std::vector<std::pair<X, T>> integrate(
	const X& x0, const T& t0, const double dt, const size_t n,
	const SinglestepIntegrator<Int, X, T>& integrator, 
	const Func<F, X, T>& func)
{
	if (dt <= 0) throw std::invalid_argument("Invalid dt <= 0!");
	auto data = std::vector<std::pair<X, T>>(n);
	data[0] = std::make_pair(x0, t0);
	for (size_t i = 1; i < data.size(); ++i) {
		integrator.integrate(data[i - 1].first, data[i - 1].second, dt, data[i].first, data[i].second, func);
	}
	return data;
}
template<Arithmetic X, ball::Time T, class F>
void integrate_by_adams(
	std::vector<X>& xlist, std::vector<T>& tlist, const double dt,
	const AdamsIntegrator<X, T>& integrator,
	const  Func<F, X, T>& func)
{
	if (dt <= 0) throw std::invalid_argument("Invalid dt <= 0!");
	auto ptr_x = xlist.cbegin();
	auto ptr_t = tlist.cbegin();
	for (size_t i = integrator.degree(); i < xlist.size(); ++i) {
		integrator.integrate(ptr_x._Ptr, ptr_t._Ptr, dt, xlist[i], tlist[i], func);
		++ptr_x; ++ptr_t;
	}
}


void orbit_integration()
{
	constexpr const double r0{ 7e6 }, mu = EGM96::Mu();
	const double v0 = std::sqrt(mu / r0);
	Vec3 vel = Vec3(1 / std::sqrt(2), 0, 1 / std::sqrt(2)) * v0;
	Vec3 pos = Vec3(0, r0, 0);
	PV x0 = PV(pos, vel);
	auto func = [mu, r0](const PV& vec, const double& t) {
		return PV(vec.vel, normalize(vec.pos) * (-mu / r0 / r0));
	};
	auto invoker{ make_classfunc<PV, double>(&func) };
	const double tk{ 86400 }, t0{ 0 };
	const double dt{ 30 };
	size_t n = 1 + static_cast<size_t>((tk - t0) / dt);
	auto calc_traj = [mu, r0, vel, pos](const PV& x0, const double t0, const double dt, const size_t n) {
		auto data = std::vector<std::pair<Vec3, double>>(n);
		const double w{ vel.length() / r0 };
		double angle{ 0 };
		Vec3 axis = normalize(cross(pos, vel));
		for (size_t i = 0; i < n; ++i) {
			data[i].first = rotate_vector(pos, axis, angle);
			data[i].second = t0 + dt * i;
			angle += w * dt;
		}
		return data;
	};

	Stopwatch sw;
	sw.start();
	auto traj0 = calc_traj(x0, t0, dt, n);
	sw.finish();
	std::cout << sw.duration() << std::endl;
	sw.start();
	auto traj1 = integrate<PV, double, RKIntegrator<PV, double>, ClassFunc<PV, double, decltype(func)>>(x0, 0, dt, n, RKIntegrator<PV, double>(), invoker);
	sw.finish();
	std::cout << sw.duration() << std::endl;
	sw.start();
	auto traj2 = integrate<PV, double, EverhartIntegrator<PV, double, 7>, ClassFunc<PV, double, decltype(func)>>(x0, 0, dt, n, EverhartIntegrator<PV, double, 7>(), invoker);
	sw.finish();
	std::cout << sw.duration() << std::endl;
	AdamsIntegrator<PV, double> adams;
	auto traj3 = std::vector<PV>(traj1.size());
	auto tlist = std::vector<double>(traj3.size());
	for (size_t i = 0; i < adams.degree(); ++i) {
		traj3[i] = traj1[i].first;
		tlist[i] = traj1[i].second;
	}
	sw.start();
	integrate_by_adams(traj3, tlist, dt, adams, invoker);
	sw.finish();
	std::cout << sw.duration() << std::endl;
	

	auto fout = std::ofstream("orbit integration.txt");
	fout << std::setprecision(16);
	for (size_t i = 0; i < n; ++i) {
		fout << traj0[i].second << " " << traj0[i].first << " " << 
			traj1[i].first.pos << " " << 
			traj2[i].first.pos << " " <<
			traj3[i].pos << std::endl;
	}
	fout.close();
}

int main()
{
	orbit_integration();

	return 0;
}