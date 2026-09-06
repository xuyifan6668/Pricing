#include "rl/core/error.hpp"

#include <sstream>

namespace rl {

Error::Error(std::string message) : std::runtime_error(message), message_(std::move(message)) {}

void Error::addContext(std::string ctx) {
    context_.push_back(std::move(ctx));
    what_.clear();
}

const char* Error::what() const noexcept {
    if (what_.empty()) {
        try {
            std::ostringstream os;
            os << message_;
            if (!context_.empty()) {
                os << " [context: ";
                // Stored innermost-first; print outermost-first so it reads like a call stack.
                for (auto it = context_.rbegin(); it != context_.rend(); ++it) {
                    if (it != context_.rbegin()) os << " > ";
                    os << *it;
                }
                os << "]";
            }
            what_ = os.str();
        } catch (...) {
            return message_.c_str();
        }
    }
    return what_.c_str();
}

namespace detail {
void throwInput(const char* expr, const char* file, int line, const std::string& msg) {
    std::ostringstream os;
    os << msg << " (failed: " << expr << " at " << file << ":" << line << ")";
    throw InputError(os.str());
}
} // namespace detail

} // namespace rl
