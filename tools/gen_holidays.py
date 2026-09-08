#!/usr/bin/env python3
"""Generate the holiday files under data/holidays/ from published rules.

Run:  uv run tools/gen_holidays.py          (regenerates all files, byte-identical if rules unchanged)
Test: uv run pytest tools/

Output format (one file per calendar): one ISO date per line, optional "# comment" after it, blank
lines and full-line comments ignored. The C++ Calendar loader (cpp/core/src/calendar.cpp) only ever
reads these files; it never computes a holiday itself. That is deliberate: production systems load
holiday lists from a reference-data source, they do not encode "third Monday of February" in C++.

Calendars:
  us_sifma    US bond market (SIFMA recommended full closes). SOFR fixings and USD Treasury settlement.
  us_fedwire  Federal Reserve / Fedwire. USD cash settlement; Fed Funds fixings.
  uk          England & Wales bank holidays. SONIA, GBP settlement.
  target2     TARGET2 (Eurosystem). ESTR, EUR settlement.
  tokyo       Japanese national holidays plus bank holidays (Dec 31, Jan 2, Jan 3). TONA, JPY settlement.

Every departure from the exact published rule is tagged SIMPLIFICATION below.
"""
from __future__ import annotations

import datetime as dt
from pathlib import Path

FIRST_YEAR = 2000
LAST_YEAR = 2075

D = dt.date
MON, TUE, WED, THU, FRI, SAT, SUN = range(7)


# ----------------------------------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------------------------------
def easter(year: int) -> D:
    """Gregorian Easter Sunday. Meeus/Jones/Butcher algorithm (Meeus, Astronomical Algorithms, ch. 8)."""
    a = year % 19
    b, c = divmod(year, 100)
    d, e = divmod(b, 4)
    f = (b + 8) // 25
    g = (b - f + 1) // 3
    h = (19 * a + b - d - g + 15) % 30
    i, k = divmod(c, 4)
    l = (32 + 2 * e + 2 * i - h - k) % 7
    m = (a + 11 * h + 22 * l) // 451
    month, day = divmod(h + l - 7 * m + 114, 31)
    return D(year, month, day + 1)


def nth_weekday(year: int, month: int, weekday: int, n: int) -> D:
    """n-th (1-based) given weekday of the month; n = -1 means the last one."""
    if n > 0:
        first = D(year, month, 1)
        delta = (weekday - first.weekday()) % 7
        return first + dt.timedelta(days=delta + 7 * (n - 1))
    last = D(year + (month == 12), month % 12 + 1, 1) - dt.timedelta(days=1)
    delta = (last.weekday() - weekday) % 7
    return last - dt.timedelta(days=delta)


def observed_sun_to_mon(d: D) -> D:
    return d + dt.timedelta(days=1) if d.weekday() == SUN else d


def observed_nearest_weekday(d: D) -> D:
    if d.weekday() == SAT:
        return d - dt.timedelta(days=1)
    if d.weekday() == SUN:
        return d + dt.timedelta(days=1)
    return d


# ----------------------------------------------------------------------------------------------------
# United States
# ----------------------------------------------------------------------------------------------------
def _us_fixed_holidays(year: int) -> list[tuple[str, D]]:
    hols = [
        ("New Year's Day", D(year, 1, 1)),
        ("Martin Luther King Jr. Day", nth_weekday(year, 1, MON, 3)),
        ("Presidents' Day", nth_weekday(year, 2, MON, 3)),
        ("Memorial Day", nth_weekday(year, 5, MON, -1)),
        ("Independence Day", D(year, 7, 4)),
        ("Labor Day", nth_weekday(year, 9, MON, 1)),
        ("Columbus Day", nth_weekday(year, 10, MON, 2)),
        ("Veterans Day", D(year, 11, 11)),
        ("Thanksgiving Day", nth_weekday(year, 11, THU, 4)),
        ("Christmas Day", D(year, 12, 25)),
    ]
    if year >= 2022:  # Federal holiday from 17 Jun 2021; Fed and SIFMA first observed it in 2022.
        hols.append(("Juneteenth", D(year, 6, 19)))
    return hols


def us_fedwire(year: int) -> list[tuple[D, str]]:
    """Federal Reserve Bank holidays. Rule: Sunday holidays are observed Monday; Saturday holidays are
    NOT observed (Fedwire is open the preceding Friday). Source: federalreserve.gov/aboutthefed/k8.htm."""
    out = []
    for name, d in _us_fixed_holidays(year):
        if d.weekday() == SAT:
            continue
        if d.weekday() == SUN:
            out.append((d + dt.timedelta(days=1), f"{name} (observed)"))
        else:
            out.append((d, name))
    return out


_SIFMA_ONE_OFFS = {
    D(2001, 9, 11): "September 11 attacks",
    D(2001, 9, 12): "September 11 attacks",
    D(2001, 9, 13): "September 11 attacks",
    D(2001, 9, 14): "September 11 attacks",
    D(2004, 6, 11): "National day of mourning (Reagan)",
    D(2012, 10, 30): "Hurricane Sandy",
    D(2018, 12, 5): "National day of mourning (G.H.W. Bush)",
    D(2025, 1, 9): "National day of mourning (Carter)",
}


def us_sifma(year: int) -> list[tuple[D, str]]:
    """SIFMA recommended full-close days for the US fixed-income market.

    Rules: Fed holidays; Sunday -> Monday; Saturday -> preceding Friday (SIFMA differs from Fedwire
    here: e.g. Fri 2026-07-03 and Fri 2021-12-24 are SIFMA closes with Fedwire open); Good Friday.

    SIMPLIFICATION: (1) New Year's Day on a Saturday is treated as not observed (matches SIFMA's actual
    2021-12-31 recommendation, which was a full open day). (2) Good Friday is always a full close.
    SIFMA in fact recommends only an early close (12:00) when Good Friday coincides with the monthly
    employment report (2010, 2012, 2015, 2021, 2023); production keeps a separate early-close list.
    (3) One-off closes are hard-coded above and only known through 2025.
    """
    out = []
    for name, d in _us_fixed_holidays(year):
        if d.weekday() == SAT:
            if name == "New Year's Day":
                continue
            out.append((d - dt.timedelta(days=1), f"{name} (observed)"))
        elif d.weekday() == SUN:
            out.append((d + dt.timedelta(days=1), f"{name} (observed)"))
        else:
            out.append((d, name))
    out.append((easter(year) - dt.timedelta(days=2), "Good Friday"))
    out += [(d, n) for d, n in _SIFMA_ONE_OFFS.items() if d.year == year]
    return out


# ----------------------------------------------------------------------------------------------------
# United Kingdom (England & Wales)
# ----------------------------------------------------------------------------------------------------
_UK_ONE_OFFS = {
    D(2002, 6, 3): "Spring bank holiday (moved for Golden Jubilee)",
    D(2002, 6, 4): "Golden Jubilee",
    D(2011, 4, 29): "Royal wedding",
    D(2012, 6, 4): "Spring bank holiday (moved for Diamond Jubilee)",
    D(2012, 6, 5): "Diamond Jubilee",
    D(2020, 5, 8): "Early May bank holiday (moved to VE Day 75)",
    D(2022, 6, 2): "Spring bank holiday (moved for Platinum Jubilee)",
    D(2022, 6, 3): "Platinum Jubilee",
    D(2022, 9, 19): "State funeral of Queen Elizabeth II",
    D(2023, 5, 8): "Coronation of King Charles III",
}
# Regular holidays that were *moved* in those years (so the regular date is not a holiday).
_UK_SUPPRESSED = {D(2002, 5, 27), D(2012, 5, 28), D(2020, 5, 4), D(2022, 5, 30)}


def uk(year: int) -> list[tuple[D, str]]:
    """England & Wales bank holidays. Source: gov.uk/bank-holidays and the Banking and Financial
    Dealings Act 1971. Substitute-day rule: a weekend holiday moves to the next weekday not already a
    holiday, so Christmas/Boxing Day resolve jointly.

    SIMPLIFICATION: one-offs known through 2023 only; future royal or jubilee days are unknowable.
    """
    out = []
    out.append((observed_sun_to_mon(D(year, 1, 1)) if D(year, 1, 1).weekday() != SAT
                else D(year, 1, 3), "New Year's Day" + (" (substitute)" if D(year, 1, 1).weekday() >= SAT else "")))
    e = easter(year)
    out.append((e - dt.timedelta(days=2), "Good Friday"))
    out.append((e + dt.timedelta(days=1), "Easter Monday"))
    for d, name in [(nth_weekday(year, 5, MON, 1), "Early May bank holiday"),
                    (nth_weekday(year, 5, MON, -1), "Spring bank holiday"),
                    (nth_weekday(year, 8, MON, -1), "Summer bank holiday")]:
        if d not in _UK_SUPPRESSED:
            out.append((d, name))
    xmas, boxing = D(year, 12, 25), D(year, 12, 26)
    if xmas.weekday() == SAT:        # Sat 25, Sun 26 -> Mon 27, Tue 28
        out += [(D(year, 12, 27), "Christmas Day (substitute)"), (D(year, 12, 28), "Boxing Day (substitute)")]
    elif xmas.weekday() == SUN:      # Sun 25, Mon 26 -> Mon 26 Boxing, Tue 27 Christmas subst
        out += [(D(year, 12, 26), "Boxing Day"), (D(year, 12, 27), "Christmas Day (substitute)")]
    elif xmas.weekday() == FRI:      # Fri 25, Sat 26 -> Mon 28
        out += [(xmas, "Christmas Day"), (D(year, 12, 28), "Boxing Day (substitute)")]
    else:
        out += [(xmas, "Christmas Day"), (boxing, "Boxing Day")]
    out += [(d, n) for d, n in _UK_ONE_OFFS.items() if d.year == year]
    return out


# ----------------------------------------------------------------------------------------------------
# TARGET2
# ----------------------------------------------------------------------------------------------------
def target2(year: int) -> list[tuple[D, str]]:
    """TARGET2 closing days (ECB): 1 Jan, Good Friday, Easter Monday, 1 May, 25 Dec, 26 Dec.
    No weekend substitution. SIMPLIFICATION: the same rule is applied to 2000-2001, when TARGET (the
    predecessor) had a slightly different list."""
    e = easter(year)
    return [
        (D(year, 1, 1), "New Year's Day"),
        (e - dt.timedelta(days=2), "Good Friday"),
        (e + dt.timedelta(days=1), "Easter Monday"),
        (D(year, 5, 1), "Labour Day"),
        (D(year, 12, 25), "Christmas Day"),
        (D(year, 12, 26), "Christmas Holiday"),
    ]


# ----------------------------------------------------------------------------------------------------
# Tokyo
# ----------------------------------------------------------------------------------------------------
def _equinox_day(year: int, vernal: bool) -> int:
    """Vernal/autumnal equinox day (JST), valid 2000-2099. Standard approximation used by the Japanese
    National Astronomical Observatory tables. SIMPLIFICATION: the official date is announced each
    February for the following year from astronomical computation; this formula matches it for all
    years 2000-2099 in published tables but is not the legal source."""
    base = 20.8431 if vernal else 23.2488
    return int(base + 0.242194 * (year - 1980)) - int((year - 1980) / 4)


def tokyo(year: int) -> list[tuple[D, str]]:
    """Japanese national holidays (Act on National Holidays) plus bank holidays, as observed by the
    Tokyo money market. Rules current as of 2026:
      substitute holiday: a national holiday on Sunday moves to the next day that is not a holiday;
      citizens' holiday: a weekday between two national holidays is itself a holiday.
    Bank holidays: 31 Dec, 2 Jan, 3 Jan (not national holidays; markets closed).
    Olympic-year moves (2020, 2021) and 2019 enthronement days are hard-coded."""
    nat: dict[D, str] = {}

    def add(d: D, name: str) -> None:
        nat[d] = name

    add(D(year, 1, 1), "New Year's Day")
    add(nth_weekday(year, 1, MON, 2), "Coming of Age Day")
    add(D(year, 2, 11), "National Foundation Day")
    if year >= 2020:
        add(D(year, 2, 23), "Emperor's Birthday")
    elif year <= 2018:
        add(D(year, 12, 23), "Emperor's Birthday")
    add(D(year, 3, _equinox_day(year, True)), "Vernal Equinox Day")
    add(D(year, 4, 29), "Showa Day" if year >= 2007 else "Greenery Day")
    add(D(year, 5, 3), "Constitution Memorial Day")
    add(D(year, 5, 4), "Greenery Day" if year >= 2007 else "Citizens' Holiday")
    add(D(year, 5, 5), "Children's Day")
    if year == 2020:
        add(D(2020, 7, 23), "Marine Day (Olympics)")
        add(D(2020, 7, 24), "Sports Day (Olympics)")
        add(D(2020, 8, 10), "Mountain Day (Olympics)")
    elif year == 2021:
        add(D(2021, 7, 22), "Marine Day (Olympics)")
        add(D(2021, 7, 23), "Sports Day (Olympics)")
        add(D(2021, 8, 8), "Mountain Day (Olympics)")
    else:
        add(nth_weekday(year, 7, MON, 3) if year >= 2003 else D(year, 7, 20), "Marine Day")
        if year >= 2016:
            add(D(year, 8, 11), "Mountain Day")
        add(nth_weekday(year, 10, MON, 2), "Sports Day" if year >= 2020 else "Health and Sports Day")
    add(nth_weekday(year, 9, MON, 3) if year >= 2003 else D(year, 9, 15), "Respect for the Aged Day")
    add(D(year, 9, _equinox_day(year, False)), "Autumnal Equinox Day")
    add(D(year, 11, 3), "Culture Day")
    add(D(year, 11, 23), "Labour Thanksgiving Day")
    if year == 2019:
        add(D(2019, 5, 1), "Enthronement Day")
        add(D(2019, 10, 22), "Enthronement Ceremony")

    out: dict[D, str] = dict(nat)
    # Citizens' holiday: weekday sandwiched between two national holidays.
    for d in list(nat):
        mid = d + dt.timedelta(days=1)
        if mid not in nat and (mid + dt.timedelta(days=1)) in nat and mid.weekday() != SUN:
            out[mid] = "Citizens' Holiday"
    # Substitute holidays (rule since 2007; before that only the immediate Monday).
    for d in sorted(nat):
        if d.weekday() == SUN:
            s = d + dt.timedelta(days=1)
            while s in out:
                s += dt.timedelta(days=1)
            out[s] = f"{nat[d]} (substitute)"
    for d, name in [(D(year, 1, 2), "Bank holiday"), (D(year, 1, 3), "Bank holiday"),
                    (D(year, 12, 31), "Bank holiday")]:
        out.setdefault(d, name)
    return sorted(out.items())


# ----------------------------------------------------------------------------------------------------
CALENDARS = {
    "us_sifma": (us_sifma, "US bond market (SIFMA). SOFR fixing calendar."),
    "us_fedwire": (us_fedwire, "Federal Reserve / Fedwire. USD settlement, Fed Funds fixing calendar."),
    "uk": (uk, "England & Wales bank holidays. SONIA / GBP settlement."),
    "target2": (target2, "TARGET2. ESTR / EUR settlement."),
    "tokyo": (tokyo, "Tokyo: Japanese national holidays + bank holidays. TONA / JPY settlement."),
}


def generate(name: str) -> str:
    fn, desc = CALENDARS[name]
    lines = [f"# {name}: {desc}",
             f"# Generated by tools/gen_holidays.py for {FIRST_YEAR}-{LAST_YEAR}. Do not edit by hand.",
             "# Format: YYYY-MM-DD [# comment]. Weekends are never listed; the loader treats Sat/Sun as non-business.",
             ""]
    seen = set()
    for year in range(FIRST_YEAR, LAST_YEAR + 1):
        for d, comment in sorted(fn(year)):
            if d in seen:
                continue
            seen.add(d)
            lines.append(f"{d.isoformat()}  # {comment}")
    return "\n".join(lines) + "\n"


def main() -> None:
    out_dir = Path(__file__).resolve().parent.parent / "data" / "holidays"
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in CALENDARS:
        (out_dir / f"{name}.txt").write_text(generate(name), encoding="utf-8", newline="\n")
        print(f"wrote {out_dir / (name + '.txt')}")


if __name__ == "__main__":
    main()
