#include <iostream>
#include "DateTime.h"
#include "Constants.h"
#include "LegendrePolynom.h"
#include "MathFunctions.h"
#include "Matrix3x3.h"
#include <fstream>

using namespace ball;

void TestDateTime();
void TestLegendrePolynom();
void TestFactorial();
void TestMatrix();

int main()
{
	std::cout << "...Testing factorials...\n";
	TestFactorial();
	std::cout << "...Testing DateTime...\n";
	TestDateTime();
	std::cout << "\n...Testing Legendre polynoms...\n";
	TestLegendrePolynom();
	std::cout << "\n...Testing matrices...\n";
	TestMatrix();
	return 0;
}

void TestDateTime()
{
	ball::time::DateTime dt1,
		dt2(2013, 12, 23, 10, 56, 34, 100),
		dt3(2013, 8, 31, 9, 5, 7, 0),
		dt4(2064, 5, 7, 10, 14, 19, 0);
	std::cout << "DateTime1: " << dt1 << std::endl;
	std::cout << "DateTime2: " << dt2 << std::endl;
	std::cout << "DateTime3: " << dt3 << std::endl;
	std::cout << "DateTime4: " << dt4 << std::endl;
	std::cout << "DateTime1 > DateTime2: " << (dt3 > dt2) << std::endl;
	std::cout << "DateTime1 < DateTime2: " << (dt3 < dt2) << std::endl;
	ball::time::JD jd2(dt2), jd3(dt3), jd4 = ball::time::JD2000,
		jd5(dt4);
	std::cout << "JD2 = " << jd2 << std::endl;
	std::cout << "JD3 = " << jd3 << std::endl;
	std::cout << "JD4 = " << jd4 << std::endl;
	std::cout << "JD5 = " << jd5 << std::endl;
	std::cout << "JD2 to DateTime: " << jd2.ToDateTime() << std::endl;
	std::cout << "JD3 to DateTime: " << jd3.ToDateTime() << std::endl;
	std::cout << "JD4 to DateTime: " << jd4.ToDateTime() << std::endl;
	std::cout << "JD5 to DateTime: " << jd5.ToDateTime() << std::endl;
}

void TestLegendrePolynom()
{
	std::cout << "Legendre polynoms\n";
	math::LegendrePolynom p1(3), 
		p2(10), 
		p3(15);
	double v1 = -0.5, v2 = 0.5;
	std::cout << "P3 degree: " << p1.Degree() << "; P10 degree: " << p2.Degree() << std::endl;
	std::cout << "P3(" << v1 << ") = " << p1(v1) << std::endl;
	std::cout << "P10(" << v1 << ") = " << p2(v1) << std::endl;
	std::cout << "P15(" << v1 << ") = " << p3(v1) << std::endl;
	std::cout << "P3(" << v2 << ") = " << p1(v2) << std::endl;
	std::cout << "P10(" << v2 << ") = " << p2(v2) << std::endl;
	std::cout << "P15(" << v2 << ") = " << p3(v2) << std::endl;

	std::cout << "Legendre functions\n";
	math::LegendreFunction f1(3, 1),
		f2(10, 5),
		f3(15, 0),
		f4(3, 2);
	std::cout << "P3_1(" << v1 << ") = " << f1(v1) << std::endl;
	std::cout << "P10_5(" << v1 << ") = " << f2(v1) << std::endl;
	std::cout << "P15_0(" << v1 << ") = " << f3(v1) << std::endl;
	std::cout << "P3_2(" << v1 << ") = " << f4(v1) << std::endl;
	std::cout << "P3_1(" << v2 << ") = " << f1(v2) << std::endl;
	std::cout << "P10_5(" << v2 << ") = " << f2(v2) << std::endl;
	std::cout << "P15_0(" << v2 << ") = " << f3(v2) << std::endl;
	std::cout << "P3_2(" << v2 << ") = " << f4(v2) << std::endl;
}

void TestFactorial()
{
	std::ofstream fout("factorials.txt");
	for (size_t i = 21; i <= 30; ++i)
		fout << i << " : " << math::Factorial(i) << std::endl;
	fout.close();
}

void TestMatrix()
{
	const double a[3]{ 0.212340538, 0.590533136, 0.911412040 };
	math::Matrix3x3 m1, m2(1, 2, 3, 1, 2, 4, 3, 2, 1),
		m3 = math::Matrix3x3::Inv(m2),
		m4 = math::Matrix3x3::Eye(),
		m5 = math::Matrix3x3(
			a[0], a[0] * a[0], a[0] * a[0] * a[0],
			a[1], a[1] * a[1], a[1] * a[1] * a[1],
			a[2], a[2] * a[2], a[2] * a[2] * a[2]),
		m6 = math::Matrix3x3::Inv(m5);
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
