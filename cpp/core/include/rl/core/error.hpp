#pragma once
// Error hierarchy per docs/DESIGN.md §2.7 / D6.
//
//   rl::Error                 base; carries a message and a context chain ("bootstrap > USD-SOFR > 2y")
//     rl::InputError          bad arguments, unsorted nodes, unknown convention string, ...
//     rl::NumericalError      non-convergence, NaN, bracket failure
//     rl::CalibrationError    carries the residual vector of the failed fit
//     rl::NotImplemented      for a declared feature that has no implementation
//
// Inner loops are noexcept and validate before the loop. Nothing here is thread-shared.

#include <exception>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace rl {

class Error : public std::runtime_error {
public:
    explicit Error(std::string message);

    /// The message without the context chain.
    const std::string& message() const noexcept { return message_; }
    /// Innermost context first, e.g. {"node 2y", "USD-SOFR", "bootstrap"}.
    const std::vector<std::string>& context() const noexcept { return context_; }
    /// Append one level of context (called while unwinding through withContext()).
    void addContext(std::string ctx);
    /// "message [context: bootstrap > USD-SOFR > node 2y]"
    const char* what() const noexcept override;

private:
    std::string message_;
    std::vector<std::string> context_;
    mutable std::string what_;
};

class InputError : public Error {
public:
    using Error::Error;
};

class NumericalError : public Error {
public:
    using Error::Error;
};

class CalibrationError : public Error {
public:
    CalibrationError(std::string message, std::vector<double> residuals)
        : Error(std::move(message)), residuals_(std::move(residuals)) {}
    const std::vector<double>& residuals() const noexcept { return residuals_; }

private:
    std::vector<double> residuals_;
};

class NotImplemented : public Error {
public:
    explicit NotImplemented(std::string_view function)
        : Error("not implemented: " + std::string(function)) {}
};

/// Run f(); if an rl::Error escapes, push `ctx` onto its context chain and rethrow the same object.
template <class F>
auto withContext(std::string_view ctx, F&& f) -> decltype(f()) {
    try {
        return f();
    } catch (Error& e) {
        e.addContext(std::string(ctx));
        throw;
    }
}

namespace detail {
[[noreturn]] void throwInput(const char* expr, const char* file, int line, const std::string& msg);
}

} // namespace rl

/// Precondition check that throws rl::InputError with file:line. Use for argument validation
/// *outside* hot loops only; the loop itself is noexcept.
#define RL_REQUIRE(cond, msg)                                                                          \
    do {                                                                                               \
        if (!(cond)) ::rl::detail::throwInput(#cond, __FILE__, __LINE__, (msg));                       \
    } while (false)

/// Body of a declared but unimplemented function.
#define RL_NOT_IMPLEMENTED(name) throw ::rl::NotImplemented(name)
