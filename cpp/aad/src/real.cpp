#include "rl/aad/real.hpp"

#include <ostream>

namespace rl::aad {

std::ostream& operator<<(std::ostream& os, const Real& r) {
    os << r.value();
    if (r.isActive()) os << "[#" << r.node() << "]";
    return os;
}

} // namespace rl::aad
