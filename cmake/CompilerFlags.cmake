# Compiler flags.
#
# The one rule that matters: NO FAST-MATH ANYWHERE. Reassociation makes results depend on inlining
# decisions, so Windows and Linux CI disagree, the AAD tangent-vs-adjoint identity fails at 1e-12,
# and -ffast-math's no-NaN assumption silently deletes isnan guards. /fp:precise and
# -ffp-contract=off keep every double operation IEEE-754 as written.

add_library(rl_flags INTERFACE)
add_library(rl::flags ALIAS rl_flags)

if(MSVC)
  target_compile_options(rl_flags INTERFACE
    /W4 /WX /permissive- /Zc:__cplusplus /Zc:preprocessor /EHsc /utf-8
    /fp:precise            # and never /fp:fast; /fp:contract is deliberately NOT set
    /bigobj                # explicit instantiation TUs get large
    /wd4127                # conditional expression is constant (if constexpr patterns in tests)
    $<$<CONFIG:Release,RelWithDebInfo>:/O2>
  )
  target_compile_definitions(rl_flags INTERFACE
    _CRT_SECURE_NO_WARNINGS NOMINMAX WIN32_LEAN_AND_MEAN _USE_MATH_DEFINES)
  if(RL_SANITIZE)
    # MSVC has ASan only; UBSan/TSan live in the Linux clang presets.
    target_compile_options(rl_flags INTERFACE /fsanitize=address)
    target_link_options(rl_flags INTERFACE /INCREMENTAL:NO)
  endif()
else()
  target_compile_options(rl_flags INTERFACE
    -Wall -Wextra -Wpedantic -Werror -Wshadow -Wno-unused-parameter
    -fno-fast-math -ffp-contract=off
    $<$<CONFIG:Release,RelWithDebInfo>:-O2>
  )
  if(RL_SANITIZE)
    target_compile_options(rl_flags INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(rl_flags INTERFACE -fsanitize=address,undefined)
  endif()
endif()
