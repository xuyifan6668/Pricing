"""Pin the holiday generator to independently known dates. Run: uv run pytest tools/"""
import datetime as dt

import gen_holidays as g

D = dt.date


def dates(fn, year):
    return {d for d, _ in fn(year)}


def test_easter_known_values():
    # Published Easter Sundays.
    assert g.easter(2024) == D(2024, 3, 31)
    assert g.easter(2025) == D(2025, 4, 20)
    assert g.easter(2026) == D(2026, 4, 5)
    assert g.easter(2000) == D(2000, 4, 23)
    assert g.easter(2038) == D(2038, 4, 25)  # latest possible Easter this century


def test_sifma_vs_fedwire_saturday_rule_2026():
    # 4 July 2026 is a Saturday: SIFMA closes Friday 3 July, Fedwire does not.
    assert D(2026, 7, 3) in dates(g.us_sifma, 2026)
    assert D(2026, 7, 3) not in dates(g.us_fedwire, 2026)
    assert D(2026, 7, 4) not in dates(g.us_fedwire, 2026)


def test_sifma_christmas_2021():
    # Christmas 2021 was a Saturday: SIFMA closed Friday 24 Dec; Fedwire open.
    assert D(2021, 12, 24) in dates(g.us_sifma, 2021)
    assert D(2021, 12, 24) not in dates(g.us_fedwire, 2021)
    # New Year 2022 fell on Saturday: neither closed Friday 31 Dec 2021.
    assert D(2021, 12, 31) not in dates(g.us_sifma, 2021)


def test_us_2026_fixed_dates():
    s = dates(g.us_sifma, 2026)
    for d in [D(2026, 1, 1), D(2026, 1, 19), D(2026, 2, 16), D(2026, 4, 3), D(2026, 5, 25),
              D(2026, 6, 19), D(2026, 9, 7), D(2026, 10, 12), D(2026, 11, 11), D(2026, 11, 26),
              D(2026, 12, 25)]:
        assert d in s, d
    f = dates(g.us_fedwire, 2026)
    assert D(2026, 4, 3) not in f  # Good Friday is not a Fed holiday
    assert D(2026, 6, 19) in f


def test_juneteenth_starts_2022():
    assert D(2021, 6, 18) not in dates(g.us_fedwire, 2021)
    assert D(2022, 6, 20) in dates(g.us_fedwire, 2022)  # 19 June 2022 is a Sunday


def test_uk_2026_and_one_offs():
    u = dates(g.uk, 2026)
    assert u == {D(2026, 1, 1), D(2026, 4, 3), D(2026, 4, 6), D(2026, 5, 4), D(2026, 5, 25),
                 D(2026, 8, 31), D(2026, 12, 25), D(2026, 12, 28)}
    assert D(2023, 5, 8) in dates(g.uk, 2023)
    assert D(2022, 9, 19) in dates(g.uk, 2022)
    assert D(2022, 5, 30) not in dates(g.uk, 2022)  # Spring holiday moved to 2 June
    assert D(2022, 6, 2) in dates(g.uk, 2022) and D(2022, 6, 3) in dates(g.uk, 2022)


def test_uk_christmas_substitutes():
    # 2021: 25 Sat, 26 Sun -> 27, 28
    assert {D(2021, 12, 27), D(2021, 12, 28)} <= dates(g.uk, 2021)
    # 2022: 25 Sun -> Mon 26 Boxing, Tue 27 Christmas
    assert {D(2022, 12, 26), D(2022, 12, 27)} <= dates(g.uk, 2022)
    # 2020: 25 Fri, 26 Sat -> Mon 28
    assert {D(2020, 12, 25), D(2020, 12, 28)} <= dates(g.uk, 2020)


def test_target2_2026():
    assert dates(g.target2, 2026) == {D(2026, 1, 1), D(2026, 4, 3), D(2026, 4, 6), D(2026, 5, 1),
                                      D(2026, 12, 25), D(2026, 12, 26)}


def test_tokyo_2026():
    t = dates(g.tokyo, 2026)
    expected = {D(2026, 1, 1), D(2026, 1, 2), D(2026, 1, 3), D(2026, 1, 12), D(2026, 2, 11), D(2026, 2, 23),
                D(2026, 3, 20), D(2026, 4, 29), D(2026, 5, 3), D(2026, 5, 4), D(2026, 5, 5), D(2026, 5, 6),
                D(2026, 7, 20), D(2026, 8, 11), D(2026, 9, 21), D(2026, 9, 22), D(2026, 9, 23),
                D(2026, 10, 12), D(2026, 11, 3), D(2026, 11, 23), D(2026, 12, 31)}
    assert t == expected


def test_tokyo_equinox_formula():
    assert g._equinox_day(2024, True) == 20 and g._equinox_day(2024, False) == 22
    assert g._equinox_day(2025, True) == 20 and g._equinox_day(2025, False) == 23
    assert g._equinox_day(2027, True) == 21 and g._equinox_day(2027, False) == 23


def test_tokyo_olympic_years():
    assert D(2020, 7, 23) in dates(g.tokyo, 2020) and D(2020, 7, 24) in dates(g.tokyo, 2020)
    assert D(2021, 8, 9) in dates(g.tokyo, 2021)  # Mountain Day 8 Aug 2021 (Sun) -> Mon 9 Aug


def test_files_are_current():
    """The committed files must equal what the generator produces."""
    from pathlib import Path
    root = Path(__file__).resolve().parent.parent / "data" / "holidays"
    for name in g.CALENDARS:
        assert (root / f"{name}.txt").read_text(encoding="utf-8") == g.generate(name), name
