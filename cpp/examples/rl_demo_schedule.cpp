// rl_demo_schedule EFFECTIVE TERMINATION FREQ CALENDAR BDC STUB EOM PAYLAG
//   e.g. rl_demo_schedule 2026-09-08 2036-09-08 1Y USD-SIFMA MF ShortFront 1 2
// Prints one CSV row per period: index,unadjStart,unadjEnd,start,end,payment,regular,act360,thirty360
#include "rl/core/daycount.hpp"
#include "rl/core/error.hpp"
#include "rl/core/schedule.hpp"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    using namespace rl::core;
    try {
        if (argc != 9) {
            std::cerr << "usage: rl_demo_schedule EFFECTIVE TERMINATION FREQ CALENDAR BDC STUB EOM PAYLAG\n";
            return 2;
        }
        ScheduleSpec spec;
        spec.effective = Date::fromIso(argv[1]);
        spec.termination = Date::fromIso(argv[2]);
        spec.frequency = Tenor::parse(argv[3]);
        spec.calendar = Calendar::get(argv[4]);
        spec.convention = bdcFromString(argv[5]);
        spec.terminationConvention = spec.convention;
        spec.stub = stubRuleFromString(argv[6]);
        spec.endOfMonth = std::atoi(argv[7]) != 0;
        spec.paymentLag = std::atoi(argv[8]);

        Schedule s = Schedule::generate(spec);
        std::printf("index,unadjStart,unadjEnd,start,end,payment,regular,act360,thirty360\n");
        for (std::size_t i = 0; i < s.size(); ++i) {
            const Period& p = s[i];
            std::printf("%zu,%s,%s,%s,%s,%s,%d,%.10f,%.10f\n", i, p.unadjustedStart.toIso().c_str(), p.unadjustedEnd.toIso().c_str(),
                        p.start.toIso().c_str(), p.end.toIso().c_str(), p.payment.toIso().c_str(), p.regular ? 1 : 0,
                        yearFraction(DayCount::Act360, p.start, p.end), yearFraction(DayCount::Thirty360US, p.start, p.end));
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
