#include "rl/market/curve_interpolation.hpp"

#include "rl/core/error.hpp"

#include <algorithm>
#include <cctype>
#include <string>

namespace rl::market {

std::string_view toString(CurveInterpolation k) {
    switch (k) {
    case CurveInterpolation::LinearLogDf: return "linear_logdf";
    case CurveInterpolation::MonotoneConvex: return "monotone_convex";
    }
    return "?";
}

CurveInterpolation curveInterpolationFromString(std::string_view s) {
    std::string u;
    for (char c : s) {
        if (c == '-' || c == ' ') c = '_';
        u.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    if (u == "linear_logdf" || u == "linear" || u == "loglinear" || u == "log_linear" || u == "linear_log_df")
        return CurveInterpolation::LinearLogDf;
    if (u == "monotone_convex" || u == "monotoneconvex" || u == "mc" || u == "hagan_west")
        return CurveInterpolation::MonotoneConvex;
    throw InputError("unknown curve interpolation: '" + std::string(s) + "'");
}

} // namespace rl::market
