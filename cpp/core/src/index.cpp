#include "rl/core/index.hpp"

#include "rl/core/error.hpp"

namespace rl::core {

std::string_view toString(IndexKind k) { return k == IndexKind::Overnight ? "Overnight" : "Term"; }

std::string_view toString(Compounding c) {
    switch (c) {
    case Compounding::Simple: return "Simple";
    case Compounding::CompoundedDaily: return "CompoundedDaily";
    case Compounding::AveragedDaily: return "AveragedDaily";
    }
    return "?";
}

namespace {
Index overnight(std::string name, Currency ccy, DayCount dc, std::string cal, Compounding comp, int payDelay) {
    Index ix;
    ix.name = std::move(name);
    ix.currency = ccy;
    ix.kind = IndexKind::Overnight;
    ix.tenor = Tenor::days(1);
    ix.dayCount = dc;
    ix.fixingCalendar = std::move(cal);
    ix.fixingLag = 0;
    ix.compounding = comp;
    ix.paymentDelayDays = payDelay;
    return ix;
}

Index term(std::string name, Currency ccy, Tenor tenor, DayCount dc, std::string cal, int fixingLag) {
    Index ix;
    ix.name = std::move(name);
    ix.currency = ccy;
    ix.kind = IndexKind::Term;
    ix.tenor = tenor;
    ix.dayCount = dc;
    ix.fixingCalendar = std::move(cal);
    ix.fixingLag = fixingLag;
    ix.compounding = Compounding::Simple;
    return ix;
}
} // namespace

// Conventions: ISDA 2021 Definitions floating-rate matrix and the CME/LCH clearing specs for the
// standard OIS in each currency. Payment delay is the swap convention that travels with the index.
Index Index::sofr() { return overnight("USD-SOFR", Currency::USD, DayCount::Act360, "USD-SIFMA", Compounding::CompoundedDaily, 2); }
Index Index::fedFunds() { return overnight("USD-FEDFUNDS", Currency::USD, DayCount::Act360, "USD-FEDWIRE", Compounding::AveragedDaily, 2); }
Index Index::estr() { return overnight("EUR-ESTR", Currency::EUR, DayCount::Act360, "EUR-TARGET2", Compounding::CompoundedDaily, 1); }
Index Index::sonia() { return overnight("GBP-SONIA", Currency::GBP, DayCount::Act365Fixed, "GBP", Compounding::CompoundedDaily, 0); }
Index Index::tona() { return overnight("JPY-TONA", Currency::JPY, DayCount::Act365Fixed, "JPY", Compounding::CompoundedDaily, 2); }
Index Index::euribor3m() { return term("EUR-EURIBOR-3M", Currency::EUR, Tenor::months(3), DayCount::Act360, "EUR-TARGET2", 2); }
Index Index::euribor6m() { return term("EUR-EURIBOR-6M", Currency::EUR, Tenor::months(6), DayCount::Act360, "EUR-TARGET2", 2); }
Index Index::termSofr3m() { return term("USD-TERMSOFR-3M", Currency::USD, Tenor::months(3), DayCount::Act360, "USD-SIFMA", 2); }

Index Index::get(std::string_view name) {
    if (name == "USD-SOFR") return sofr();
    if (name == "USD-FEDFUNDS") return fedFunds();
    if (name == "EUR-ESTR") return estr();
    if (name == "GBP-SONIA") return sonia();
    if (name == "JPY-TONA") return tona();
    if (name == "EUR-EURIBOR-3M") return euribor3m();
    if (name == "EUR-EURIBOR-6M") return euribor6m();
    if (name == "USD-TERMSOFR-3M") return termSofr3m();
    throw InputError("unknown index: " + std::string(name));
}

} // namespace rl::core
