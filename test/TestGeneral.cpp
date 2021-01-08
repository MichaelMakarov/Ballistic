#include <iostream>
#include "Times.h"
#include "GeneralConstants.h"
#include "Mathematics.h"
#include "Structures.h"
#include "Geometry.h"
#include "PZ90.h"
#include <fstream>

using namespace ball::geometry;
using namespace ball::time;

void TestDateTime();
void TestLegendrePolynom();
void TestMatrix();
void TestTle();

int main()
{
	TestDateTime();
	TestLegendrePolynom();
	TestMatrix();
	TestTle();
	return 0;
}

void TestTle()
{
	std::cout << "\n...Test Tle...\n";

	auto fread = std::ifstream("TleFile.tle");
	TLE tle1;
	if (!LoadTle(fread, tle1))
		std::cout << "Failed to load tle content from file!\n";
	fread.close();
	auto op = tle1.ToOscul(std::make_unique<ball::space::PZ90>()->Mu());
	std::cout << op << std::endl;
}

void TestDateTime()
{
	std::cout << "\n...Time tests...\n";

	auto dt0{ DateTime() };
	auto dt1{ DateTime(2000, 1, 3, 23, 23, 23, 423) };
	auto dt2{ DateTime(2000, 1, 3, 0, 0, 1, 0) };
	auto dt3{ DateTime(2013, 12, 23, 1, 56, 34, 100) };
	auto dt4{ DateTime(2013, 12, 23, 11, 5, 7, 0) };
	auto dt5{ DateTime(2013, 12, 23, 20, 23, 19, 0) };

	std::cout << "dt0: " << dt0 << std::endl;
	std::cout << "dt1: " << dt1 << std::endl;
	std::cout << "dt2: " << dt2 << std::endl;
	std::cout << "dt3: " << dt3 << std::endl;
	std::cout << "dt4: " << dt4 << std::endl;
	std::cout << "dt5: " << dt5 << std::endl;

	auto jd0{ JD(JD2000) };
	auto jd1{ JD(dt1) };
	auto jd2{ JD(dt2) };
	auto jd3{ JD(dt3) };
	auto jd4{ JD(dt4) };
	auto jd5{ JD(dt5) };
	auto jd6{ jd5 + 4 / HOURS_PER_DAY };

	std::cout << "JD0 = " << jd0 << std::endl;
	std::cout << "JD1 = " << jd1 << std::endl;
	std::cout << "JD2 = " << jd2 << std::endl;
	std::cout << "JD3 = " << jd3 << std::endl;
	std::cout << "JD4 = " << jd4 << std::endl;
	std::cout << "JD5 = " << jd5 << std::endl;
	std::cout << "JD6 = " << jd6 << std::endl;

	std::cout << "JD0 to DateTime: " << jd0.ToDateTime() << std::endl;
	std::cout << "JD1 to DateTime: " << jd1.ToDateTime() << std::endl;
	std::cout << "JD2 to DateTime: " << jd2.ToDateTime() << std::endl;
	std::cout << "JD3 to DateTime: " << jd3.ToDateTime() << std::endl;
	std::cout << "JD4 to DateTime: " << jd4.ToDateTime() << std::endl;
	std::cout << "JD5 to DateTime: " << jd5.ToDateTime() << std::endl;
	std::cout << "JD6 to DateTime: " << jd6.ToDateTime() << std::endl;

}

void TestLegendrePolynom()
{
	using namespace ball::math;

	std::cout << "\n...Legendre polynoms tests...\n";

	auto p1{ LegendrePolynomial(3) };
	auto p2{ LegendrePolynomial(10) };
	auto p3{ LegendrePolynomial(15) };
	
	double v1 = -0.5, v2 = 0.5;

	std::cout << "P3 degree: " << p1.Degree() << "; P10 degree: " << p2.Degree() << std::endl;

	std::cout << "P3(" << v1 << ") = " << p1(v1) << std::endl;
	std::cout << "P10(" << v1 << ") = " << p2(v1) << std::endl;
	std::cout << "P15(" << v1 << ") = " << p3(v1) << std::endl;
	std::cout << "P3(" << v2 << ") = " << p1(v2) << std::endl;
	std::cout << "P10(" << v2 << ") = " << p2(v2) << std::endl;
	std::cout << "P15(" << v2 << ") = " << p3(v2) << std::endl;

	std::cout << "...Legendre functions tests...\n";

	auto f1{ LegendreFunction(3, 1) };
	auto f2{ LegendreFunction(10, 5) };
	auto f3{ LegendreFunction(15, 0) };
	auto f4{ LegendreFunction(3, 2) };
	auto f5{ LegendreFunction(4, 4) };

	std::cout << "P3_1(" << v1 << ") = " << f1(v1) << std::endl;
	std::cout << "P10_5(" << v1 << ") = " << f2(v1) << std::endl;
	std::cout << "P15_0(" << v1 << ") = " << f3(v1) << std::endl;
	std::cout << "P3_2(" << v1 << ") = " << f4(v1) << std::endl;
	std::cout << "P3_1(" << v2 << ") = " << f1(v2) << std::endl;
	std::cout << "P10_5(" << v2 << ") = " << f2(v2) << std::endl;
	std::cout << "P15_0(" << v2 << ") = " << f3(v2) << std::endl;
	std::cout << "P3_2(" << v2 << ") = " << f4(v2) << std::endl;
	std::cout << "P4_4(" << v2 << ") = " << f5(v2) << std::endl;
}

void TestMatrix()
{
	using namespace ball::math;

	std::cout << "\n...Matrix tests...\n";

	const double a[3]{ 0.212340538, 0.590533136, 0.911412040 };
	auto m1{ Matrix3x3() };
	auto m2{ Matrix3x3(1, 2, 3, 1, 2, 4, 3, 2, 1) };
	auto m3{ Matrix3x3::Inv(m2) };
	auto m4{ Matrix3x3::Eye() };
	auto m5{ Matrix3x3(
		a[0], a[0] * a[0], a[0] * a[0] * a[0],
		a[1], a[1] * a[1], a[1] * a[1] * a[1],
		a[2], a[2] * a[2], a[2] * a[2] * a[2]) };
	auto m6{ Matrix3x3::Inv(m5) };
	
	std::cout << "M1: " << m1 << std::endl;
	std::cout << "M2: " << m2 << std::endl;
	std::cout << "M4: " << m4 << std::endl;
	std::cout << "M3 = M2^-1: " << m3 << std::endl;
	std::cout << "M2 * M3: " << m2 * m3 << std::endl;
	std::cout << "M2 + M4: " << m2 + m4 << std::endl;
	std::cout << "M2 - M4: " << m2 - m4 << std::endl;
	std::cout << "M2 * 3.5: " << m2 * 3.5 << std::endl;
	std::cout << "M5: " << m5 << std::endl;
	std::cout << "M6: " << m6 << std::endl;
	std::cout << "M5 * M6: " << m5 * m6 << std::endl;
}
