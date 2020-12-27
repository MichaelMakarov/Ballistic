#pragma once
#include <ostream>

namespace ball
{
    namespace time
    {
        struct Date
        {
        private:

            static const unsigned int DaysCountByMonthesOrd[12];
            static const unsigned int DaysCountByMonthesLeap[12];
            static const unsigned int DaysSumByMonthesOrd[12];
            static const unsigned int DaysSumByMonthesLeap[12];

            unsigned int _day, _month, _year;

            bool CheckMonth(
                const unsigned int month);
            bool CheckDay(
                const unsigned int day,
                const unsigned int month,
                const unsigned int year);
            bool IsLeapYear(
                const unsigned int year);

        public:
            Date() : _year(0), _month(0), _day(0) {}
            Date(
                const unsigned int year,
                const unsigned int month,
                const unsigned int day);
            ~Date() {}

            unsigned int Day() const;
            unsigned int Month() const;
            unsigned int Year() const;

            friend std::ostream& operator << (std::ostream& o, const Date& d);

            friend bool operator > (const Date& f, const Date& s);
            friend bool operator < (const Date& f, const Date& s);
            friend bool operator ==(const Date& f, const Date& s);
        };

        struct Time
        {
        private:
            static const unsigned int MaxHour = 23;
            static const unsigned int MaxMinute = 59;
            static const unsigned int MaxSecond = 59;
            static const unsigned int MaxMillisec = 999;

            unsigned int _hour, _minute, _second, _millisec;

        public:
            Time() : _hour(0), _minute(0), _second(0), _millisec(0) {}
            Time(
                const unsigned int hour,
                const unsigned int minute,
                const unsigned int second,
                const unsigned int millisec = 0);
            ~Time() {}

            unsigned int Hour() const;
            unsigned int Minute() const;
            unsigned int Second() const;
            unsigned int Millisecond() const;

            friend std::ostream& operator << (std::ostream& o, const Time& d);

            friend bool operator > (const Time& f, const Time& s);
            friend bool operator < (const Time& f, const Time& s);
            friend bool operator ==(const Time& f, const Time& s);
        };

        struct DateTime
        {
        private:
            Date _date;
            Time _time;

        public:
            DateTime() {}
            DateTime(const Date& date, const Time& time);
            DateTime(
                const unsigned int year,
                const unsigned int month,
                const unsigned int day,
                const unsigned int hour,
                const unsigned int minute,
                const unsigned int second,
                const unsigned int millisec);
            ~DateTime() {}

            const Date& GetDate() const { return _date; }
            const Time& GetTime() const { return _time; }

            friend std::ostream& operator << (std::ostream& o, const DateTime& d);

            friend bool operator > (const DateTime& f, const DateTime& s);
            friend bool operator < (const DateTime& f, const DateTime& s);
            friend bool operator ==(const DateTime& f, const DateTime& s);
        };

        struct JD
        {
        private:
            size_t _jd;
            double _jt;

        public:
            JD() {}
            explicit JD(const DateTime& dt);
            JD(const size_t jdn, const double time);
            JD(const double jd);
            ~JD() {}

            JD& operator = (const double jd);
            JD& operator = (const JD& jd);

            double Get() const;
            size_t JDN() const;
            DateTime ToDateTime() const;

            void AddDays(const int n);
            void AddHours(const int n);
            void AddMinutes(const int n);
            void AddSeconds(const int n);

            operator double() const { return _jt; }

            friend std::ostream& operator << (std::ostream& o, const JD& jd);
        };
    }
}