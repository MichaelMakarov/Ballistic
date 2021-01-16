#include "Times.h"
#include "GeneralConstants.h"
#include <cmath>

namespace ball
{
	namespace time
	{
		bool Date::Check(
			const size_t year,
			const unsigned short month,
			const unsigned short day)
		{
			const unsigned short daysOrder1[12]
			{
				31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
			};
			const unsigned short daysOrder2[12]
			{
				31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
			};
			const bool leap = year % 400 == 0 ? true :
				(year % 100 == 0 ? false :
					(year % 4 == 0 ? true : false));
			if (day > 0)
			{
				if (leap) return day <= daysOrder2[month - 1];
				else return day <= daysOrder1[month - 1];
			}
			return false;
		}

		Date::Date(
			const size_t year,
			const unsigned short month,
			const unsigned short day)
		{
			if (month < 1 || month > 12)
				throw std::invalid_argument("Invalid month!");
			if (!Check(year, month, day))
				throw std::invalid_argument("Invalid day!");
			_year = year;
			_month = month;
			_day = day;
		}

		bool operator > (const Date& f, const Date& s)
		{
			if (f._year > s._year) return true;
			else if (f._year == s._year)
			{
				if (f._month > s._month) return true;
				else if (f._month == s._month)
				{
					return f._day > s._day;
				}
			}
			return false;
		}
		bool operator < (const Date& f, const Date& s)
		{
			if (f._year < s._year) return true;
			else if (f._year == s._year)
			{
				if (f._month < s._month) return true;
				else if (f._month == s._month)
				{
					return f._day < s._day;
				}
			}
			return false;
		}
		bool operator == (const Date& f, const Date& s)
		{
			return f._year == s._year && f._month == s._month && f._day == s._day;
		}

		std::ostream& operator << (std::ostream& o, const Date& d)
		{
			o << d._year << '/' << d._month << '/' << d._day;
			return o;
		}


		Time::Time(
			const unsigned short hour,
			const unsigned short minute,
			const unsigned short second,
			const unsigned short millisecond)
		{
			if (hour > 23)
				throw std::invalid_argument("Invalid hour > 23!");
			if (minute > 59)
				throw std::invalid_argument("Invalid minute > 59!");
			if (second > 59)
				throw std::invalid_argument("Invalid second > 59!");
			if (millisecond > 999)
				throw std::invalid_argument("Invalid millisecond > 999!");
			_hour = hour;
			_minute = minute;
			_second = second;
			_millisec = millisecond;
		}

		std::ostream& operator << (std::ostream& o, const Time& t)
		{
			o << t._hour << ':' << t._minute << ':' << t._second << '.' << t._millisec;
			return o;
		}

		bool operator > (const Time& f, const Time& s)
		{
			if (f._hour > s._hour) return true;
			else if (f._hour == s._hour)
			{
				if (f._minute > s._minute) return true;
				else if (f._minute == s._minute)
				{
					if (f._second > s._second) return true;
					else if (f._second == s._second)
					{
						return f._millisec > s._millisec;
					}
				}
			}
			return false;
		}
		bool operator < (const Time& f, const Time& s)
		{
			if (f._hour < s._hour) return true;
			else if (f._hour == s._hour)
			{
				if (f._minute < s._minute) return true;
				else if (f._minute == s._minute)
				{
					if (f._second < s._second) return true;
					else if (f._second == s._second)
					{
						return f._millisec < s._millisec;
					}
				}
			}
			return false;
		}
		bool operator == (const Time& f, const Time& s)
		{
			return f._hour == s._hour && f._minute == s._minute && f._second == s._second && f._millisec == s._millisec;
		}


		std::ostream& operator << (std::ostream& o, const DateTime& d)
		{
			o << d._date << ' ' << d._time;
			return o;
		}
		bool operator > (const DateTime& f, const DateTime& s)
		{
			if (f._date > s._date) return true;
			if (f._date == s._date)
				return f._time > s._time;
			return false;
		}
		bool operator < (const DateTime& f, const DateTime& s)
		{
			if (f._date < s._date) return true;
			else if (f._date == s._date)
				return f._time < s._time;
			return false;
		}
		bool operator == (const DateTime& f, const DateTime& s)
		{
			return f._date == s._date && f._time == s._time;
		}

		JD::JD(const double jd)
		{
			double jdn;
			_time = std::modf(jd, &jdn);
			_day = static_cast<size_t>(jdn);
		}
		JD::JD(const DateTime& dt)
		{
			Date d = dt.get_date();
			Time t = dt.get_time();
			unsigned int a = (14 - d.get_month()) / 12,
				y = 4800 + d.get_year() - a,
				m = d.get_month() + 12 * a - 3;
			_day = d.get_day() + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
			_time = (t.get_hour() + (t.get_minute() + (t.get_second() + t.get_millisecond() * 1e-3) / 60) / 60) / 24;
		}

		DateTime JD::to_datetime() const
		{
			size_t
				a = _day + 32044,
				b = (4 * a + 3) / 146097,
				c = a - (146097 * b) / 4,
				d = (4 * c + 3) / 1461,
				e = c - (1461 * d) / 4,
				m = (5 * e + 2) / 153,
				day = e - (153 * m + 2) / 5 + 1,
				month = m + 3 - 12 * (m / 10),
				year = 100 * b + d - 4800 + (m / 10),
				t = static_cast<size_t>(_time * SEC_PER_DAY * 1e3);
			auto millisec = static_cast<unsigned short>(t % 1000);
			t /= 1000;
			auto second = static_cast<unsigned short>(t % 60);
			t /= 60;
			auto minute = static_cast<unsigned short>(t % 60);
			t /= 60;
			auto hour = static_cast<unsigned short>(t % 24);
			return DateTime(year, month, day, hour, minute, second, millisec);
		}

		double JD::to_double() const
		{
			return _time + _day;
		}

		JD& JD::operator = (const double jd)
		{
			double jdn;
			_time = std::modf(jd, &jdn);
			_day = static_cast<size_t>(jdn);
			return *this;
		}

		void JD::add_days(const int n)
		{
			_day += n;
		}
		void JD::add_hours(const int n)
		{
			double day;
			_time = std::modf(_time + n / HOURS_PER_DAY, &day);
			_day += static_cast<size_t>(day);
		}
		void JD::add_minutes(const int n)
		{
			double day;
			_time = std::modf(_time + n / MIN_PER_DAY, &day);
			_day += static_cast<size_t>(day);
		}
		void JD::add_seconds(const int n)
		{
			double day;
			_time = std::modf(_time + n / SEC_PER_DAY, &day);
			_day += static_cast<size_t>(day);
		}

		JD operator + (const JD& jd, const double dt)
		{
			auto result{ jd };
			double day;
			result._time = std::modf(jd._time + dt, &day);
			result._day += static_cast<size_t>(day);
			return result;
		}
		JD operator - (const JD& jd, const double dt)
		{
			auto result{ jd };
			double day;
			result._time = std::modf(jd._time - dt, &day);
			result._day += static_cast<size_t>(day);
			return result;
		}

		double operator - (const JD& f, const JD& s)
		{
			return (f._time - s._time) + (f._day - s._day);
		}

		bool operator < (const JD& f, const JD& s)
		{
			return f._day < s._day && f._time < s._time;
		}
		bool operator > (const JD& f, const JD& s)
		{
			return f._day > s._day && f._time > s._time;;
		}
		bool operator == (const JD& f, const JD& s)
		{
			return f._day == s._day && std::abs(f._time - s._time) < (1e-3 / SEC_PER_DAY);
		}
		bool operator <= (const JD& f, const JD& s)
		{
			return f._day < s._day || (f._day == s._day && f._time <= s._time);
		}
		bool operator >= (const JD& f, const JD& s)
		{
			return f._day > s._day || (f._day == s._day && f._time >= s._time);
		}

		std::ostream& operator << (std::ostream& o, const JD& jd)
		{
			o << jd._day << std::ends << jd._time;
			return o;
		}
	}
}