#include <stdexcept>
#include "DateTime.h"
#include "Constants.h"

namespace ball
{
    namespace time
    {
        const unsigned int Date::DaysCountByMonthesOrd[12] =
        {
            31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
        };
        const unsigned int Date::DaysCountByMonthesLeap[12] =
        {
            31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
        };
        const unsigned int Date::DaysSumByMonthesOrd[12] =
        {
            31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365
        };
        const unsigned int Date::DaysSumByMonthesLeap[12] =
        {
            31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366
        };

        Date::Date(
            const unsigned int year,
            const unsigned int month,
            const unsigned int day)
        {
            if (!CheckMonth(month))
                throw std::runtime_error("Invalid month!");
            if (!CheckDay(day, month, year))
                throw std::runtime_error("Invalid day!");
            _year = year;
            _month = month;
            _day = day;
        }

        bool Date::IsLeapYear(
            const unsigned int year)
        {
            return year % 400 ?
                true : (year % 4 ?
                    (year % 100 ?
                        false : true) : false);
        }

        bool Date::CheckMonth(
            const unsigned int month)
        {
            if (month <= 12 && month >= 1) return true;
            return false;
        }

        bool Date::CheckDay(
            const unsigned int day,
            const unsigned int month,
            const unsigned int year)
        {
            if (day >= 1 && day <= 31)
            {
                if (IsLeapYear(year))
                {
                    if (day <= DaysCountByMonthesLeap[month - 1])
                        return true;
                }
                else {
                    if (day <= DaysCountByMonthesOrd[month - 1])
                        return true;
                }
            }
            return false;
        }

        unsigned int Date::Year() const { return _year; }
        unsigned int Date::Month() const { return _month; }
        unsigned int Date::Day() const { return _day; }

        std::ostream& operator << (std::ostream& o, const Date& d)
        {
            o << d._year << '/' << d._month << '/' << d._day;
            return o;
        }

        bool operator > (const Date& f, const Date& s)
        {
            if (f._year > s._year) return true;
            if (f._year == s._year)
            {
                if (f._month > s._month) return true;
                if (f._month == s._month)
                    if (f._day > s._day) return true;
            }
            return false;
        }
        bool operator < (const Date& f, const Date& s)
        {
            if (f._year < s._year) return true;
            if (f._year == s._year)
            {
                if (f._month < s._month) return true;
                if (f._month == s._month)
                    if (f._day < s._day) return true;
            }
            return false;
        }
        bool operator == (const Date& f, const Date& s)
        {
            return f._year == s._year &&
                f._month == s._month &&
                f._day == s._day;
        }

        Time::Time(
            const unsigned int hour,
            const unsigned int minute,
            const unsigned int second,
            const unsigned int millisec)
        {
            if (hour > MaxHour) throw std::runtime_error("Invalid hour!");
            if (minute > MaxMinute) throw std::runtime_error("Invalid minute!");
            if (second > MaxSecond) throw std::runtime_error("Invalid second!");
            if (millisec > MaxMillisec) throw std::runtime_error("Invalid mili seconds!");
            _hour = hour;
            _minute = minute;
            _second = second;
            _millisec = millisec;
        }

        unsigned int Time::Hour() const { return _hour; }
        unsigned int Time::Minute() const { return _minute; }
        unsigned int Time::Second() const { return _second; }
        unsigned int Time::Millisecond() const { return _millisec; }

        std::ostream& operator << (std::ostream& o, const Time& t)
        {
            o << t._hour << ':' << t._minute << ':' << t._second << '.' << t._millisec;
            return o;
        }
        bool operator > (const Time& f, const Time& s)
        {
            if (f._hour > s._hour) return true;
            if (f._hour == s._hour)
            {
                if (f._minute > s._minute) return true;
                if (f._minute == s._minute)
                {
                    if (f._second > s._second) return true;
                    if (f._second == s._second)
                        if (f._millisec > s._millisec)
                            return true;
                }
            }
            return false;
        }
        bool operator < (const Time& f, const Time& s)
        {
            if (f._hour < s._hour) return true;
            if (f._hour == s._hour)
            {
                if (f._minute < s._minute) return true;
                if (f._minute == s._minute)
                {
                    if (f._second < s._second) return true;
                    if (f._second == s._second)
                        if (f._millisec < s._millisec)
                            return true;
                }
            }
            return false;
        }
        bool operator == (const Time& f, const Time& s)
        {
            return f._hour == s._hour &&
                f._minute == s._minute &&
                f._second == s._second &&
                f._millisec == s._millisec;
        }

        DateTime::DateTime(
            const unsigned int year,
            const unsigned int month,
            const unsigned int day,
            const unsigned int hour,
            const unsigned int minute,
            const unsigned int second,
            const unsigned int millisec)
        {
            _date = time::Date(year, month, day);
            _time = time::Time(hour, minute, second, millisec);
        }
        DateTime::DateTime(
            const time::Date& date,
            const time::Time& time)
        {
            _date = date;
            _time = time;
        }

        const time::Date& DateTime::Date() const { return _date; }
        const time::Time& DateTime::Time() const { return _time; }

        std::ostream& operator << (std::ostream& o, const DateTime& d)
        {
            o << d._date << ' ' << d._time;
            return o;
        }
        bool operator > (const DateTime& f, const DateTime& s)
        {
            if (f._date > s._date) return true;
            if (f._date == s._date)
                if (f._time > s._time) return true;
            return false;
        }
        bool operator < (const DateTime& f, const DateTime& s)
        {
            if (f._date < s._date) return true;
            if (f._date == s._date)
                if (f._time < s._time) return true;
            return false;
        }
        bool operator == (const DateTime& f, const DateTime& s)
        {
            return f._date == s._date && f._time == s._time;
        }

        JD::JD(const double jd)
        {
            if (jd < 0.0) throw std::runtime_error("Invalid JD!");
            _jd = jd;
        }
        JD::JD(const DateTime& dt)
        {
            Date d = dt.Date();
            Time t = dt.Time();
            unsigned int a = (14 - d.Month()) / 12,
                y = d.Year() + 4800 - a,
                m = d.Month() + 12 * a - 3;
            _jd = d.Day() + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
            _jd += ((static_cast<double>(t.Hour()) - 12) +
                (static_cast<double>(t.Minute()) +
                    (static_cast<double>(t.Second()) +
                        static_cast<double>(t.Millisecond()) * 1e-3) / 60) / 60) / 24;
        }

        double JD::Get() const { return _jd; }
        unsigned int JD::JDN() const { return static_cast<unsigned int>(_jd); }

        DateTime JD::ToDateTime() const
        {
            unsigned int jdn = JDN(),
                a = jdn + 32044,
                b = (4 * a + 3) / 146097,
                c = a - (146097 * b) / 4,
                d = (4 * c + 3) / 1461,
                e = c - (1461 * d) / 4,
                m = (5 * e + 2) / 153,
                day = e - (153 * m + 2) / 5 + 1,
                month = m + 3 - 12 * (m / 10),
                year = 100 * b + d - 4800 + (m / 10);
            double p = _jd - jdn;
            unsigned int hour = p * HoursPerDay,
                minute = p * MinPerDay - hour * MinPerHour,
                second = p * SecPerDay - hour * SecPerHour - minute * SecPerMin,
                millisec = static_cast<unsigned int>(p * SecPerDay * 1000) % 1000;
            if (hour >= 12) {
                day += 1;
                hour -= 12;
            }
            return DateTime(year, month, day, hour, minute, second, millisec);
        }

        void JD::AddDays(const int n) { _jd += n; }
        void JD::AddHours(const int n)
        {
            _jd += static_cast<double>(n) / HoursPerDay;
        }
        void JD::AddMinutes(const int n)
        {
            _jd += static_cast<double>(n) / MinPerDay;
        }
        void JD::AddSeconds(const int n)
        {
            _jd += static_cast<double>(n) / SecPerDay;
        }

        std::ostream& operator << (std::ostream& o, const JD& jd)
        {
            o << jd._jd;
            return o;
        }
    }
}