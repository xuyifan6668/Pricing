#include "rl/bootstrap/instrument_set.hpp"

#include "rl/core/error.hpp"
#include "rl/core/index.hpp"
#include "rl/instruments/builder.hpp"
#include "rl/instruments/cashflow.hpp"
#include "rl/instruments/conventions.hpp"

#include <algorithm>

namespace rl::bootstrap {

using market::Quote;
using market::QuoteType;

namespace {

core::Date lastAccrualEnd(const instruments::Leg& leg) { return leg.lastAccrualEnd(); }

BootstrapInstrument buildOne(const Quote& q, core::Date asOf, const InstrumentSetOptions& opt) {
    const std::string curve = q.curve.empty() ? q.index : q.curve; // basis swaps pin the spread-index curve unless told otherwise
    switch (q.type) {
    case QuoteType::OisSwap: {
        auto b = instruments::OisSwapBuilder(instruments::Conventions::oisSwap(q.index), asOf).id(q.id).fixedRate(q.asRate()).notional(opt.notional);
        if (q.tenor) b.tenor(*q.tenor);
        else b.effective(*q.start).maturity(*q.end);
        instruments::Swap s = b.build();
        const core::Date maturity = lastAccrualEnd(s.fixedLeg());
        return BootstrapInstrument{q.id, q.type, curve, maturity, q.asRate(), q, BootstrapTrade(std::move(s))};
    }
    case QuoteType::BasisSwap: {
        auto b = instruments::BasisSwapBuilder(instruments::Conventions::basisSwap(q.index, q.referenceIndex), asOf).id(q.id).spread(q.asRate()).notional(opt.notional);
        if (q.tenor) b.tenor(*q.tenor);
        else b.effective(*q.start).maturity(*q.end);
        instruments::BasisSwap s = b.build();
        const core::Date maturity = lastAccrualEnd(s.legs()[0]);
        return BootstrapInstrument{q.id, q.type, curve, maturity, q.asRate(), q, BootstrapTrade(std::move(s))};
    }
    case QuoteType::Future: {
        instruments::Future f = instruments::FutureBuilder(instruments::Conventions::future(q.index))
                                    .id(q.id).start(*q.start).end(*q.end).price(q.asPrice()).convexityAdjustmentBp(q.convexityAdjustmentBp).build();
        const core::Date maturity = f.end();
        return BootstrapInstrument{q.id, q.type, curve, maturity, q.asPrice(), q, BootstrapTrade(std::move(f))};
    }
    case QuoteType::Deposit: {
        auto b = instruments::DepositBuilder(instruments::Conventions::deposit(q.index), asOf).id(q.id).rate(q.asRate()).notional(opt.notional);
        if (q.tenor) b.tenor(*q.tenor);
        else b.start(*q.start).end(*q.end);
        instruments::Deposit d = b.build();
        const core::Date maturity = d.end();
        return BootstrapInstrument{q.id, q.type, curve, maturity, q.asRate(), q, BootstrapTrade(std::move(d))};
    }
    case QuoteType::Fra: {
        instruments::Fra f = instruments::FraBuilder(core::Index::get(q.index)).id(q.id).start(*q.start).end(*q.end).rate(q.asRate()).notional(opt.notional).build();
        const core::Date maturity = f.end();
        return BootstrapInstrument{q.id, q.type, curve, maturity, q.asRate(), q, BootstrapTrade(std::move(f))};
    }
    case QuoteType::Turn: {
        TurnOfYear t{q.id, q.index, *q.start, *q.end, q.asRate()};
        const core::Date maturity = t.end;
        return BootstrapInstrument{q.id, q.type, curve, maturity, t.jump, q, BootstrapTrade(std::move(t))};
    }
    }
    throw InputError("quote " + q.id + ": unknown quote type");
}

} // namespace

InstrumentSet buildInstrumentSet(const market::QuoteSet& quotes, const InstrumentSetOptions& options) {
    InstrumentSet out;
    std::vector<std::string> order;
    for (const Quote& q : quotes) {
        if (q.type == QuoteType::Turn && !options.includeTurns) continue;
        out.push_back(withContext("quote " + q.id, [&] { return buildOne(q, quotes.asOf(), options); }));
        if (std::find(order.begin(), order.end(), out.back().curveName) == order.end()) order.push_back(out.back().curveName);
    }
    auto rank = [&](const std::string& c) { return static_cast<std::size_t>(std::find(order.begin(), order.end(), c) - order.begin()); };
    std::stable_sort(out.begin(), out.end(), [&](const BootstrapInstrument& a, const BootstrapInstrument& b) {
        if (a.curveName != b.curveName) return rank(a.curveName) < rank(b.curveName);
        return a.maturity < b.maturity;
    });
    for (std::size_t i = 1; i < out.size(); ++i) {
        if (out[i].curveName == out[i - 1].curveName && out[i].maturity == out[i - 1].maturity)
            throw InputError("instruments " + out[i - 1].id + " and " + out[i].id + " pin the same pillar " + out[i].maturity.toIso() + " on " + out[i].curveName +
                             " (singular Jacobian)");
        RL_REQUIRE(out[i].maturity > quotes.asOf(), "instrument " + out[i].id + " matures on or before as-of");
    }
    if (!out.empty()) RL_REQUIRE(out.front().maturity > quotes.asOf(), "instrument " + out.front().id + " matures on or before as-of");
    return out;
}

std::vector<std::string> curveNames(const InstrumentSet& instruments) {
    std::vector<std::string> out;
    for (const auto& bi : instruments)
        if (std::find(out.begin(), out.end(), bi.curveName) == out.end()) out.push_back(bi.curveName);
    return out;
}

std::vector<const BootstrapInstrument*> instrumentsOn(const InstrumentSet& instruments, std::string_view curveName) {
    std::vector<const BootstrapInstrument*> out;
    for (const auto& bi : instruments)
        if (bi.curveName == curveName) out.push_back(&bi);
    std::stable_sort(out.begin(), out.end(), [](const BootstrapInstrument* a, const BootstrapInstrument* b) { return a->maturity < b->maturity; });
    return out;
}

std::string_view tradeKind(const BootstrapTrade& t) noexcept {
    switch (t.index()) {
    case 0: return "deposit";
    case 1: return "fra";
    case 2: return "future";
    case 3: return "swap";
    case 4: return "basis_swap";
    case 5: return "turn";
    default: return "?";
    }
}

} // namespace rl::bootstrap
