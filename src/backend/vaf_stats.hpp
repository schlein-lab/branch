#pragma once

// BRANCH v0.2 — VAF statistics helpers.
//
// Wilson score interval for a binomial proportion. Chosen over the
// normal-approximation (Wald) interval because Wald degenerates at the
// [0, 1] boundaries — VAF distributions live exactly there for
// homozygous calls, so Wald would report useless 0-width CIs.
//
// Reference: Wilson, E.B. (1927). "Probable inference, the law of
// succession, and statistical inference." JASA 22(158), 209-212.

#include <cmath>
#include <cstdint>

#include "backend/backend_vtable.hpp"

namespace branch::backend {

// Wilson 95% confidence interval for binomial proportion p = k / n.
//
//     centre = (p + z²/2n) / (1 + z²/n)
//     half   = z * sqrt( p(1-p)/n + z²/4n² ) / (1 + z²/n)
//     lo     = centre - half
//     hi     = centre + half
//
// Edge case: total == 0 returns the uninformative prior
// {point = 0.5, ci_low = 0.0, ci_high = 1.0}. The caller should treat
// this as "no data" rather than "VAF = 0.5".
inline VAFEstimate wilson_ci(std::uint32_t successes,
                             std::uint32_t total) noexcept {
    if (total == 0) {
        return VAFEstimate{.point = 0.5f, .ci_low = 0.0f, .ci_high = 1.0f};
    }

    // z = 1.96 for a two-sided 95% CI (normal-approximation quantile).
    constexpr double z = 1.96;
    constexpr double z2 = z * z;

    const double n = static_cast<double>(total);
    const double p = static_cast<double>(successes) / n;

    const double denom = 1.0 + z2 / n;
    const double centre = (p + z2 / (2.0 * n)) / denom;
    const double half =
        z * std::sqrt(p * (1.0 - p) / n + z2 / (4.0 * n * n)) / denom;

    double lo = centre - half;
    double hi = centre + half;
    // Numerical clamp: on unusual inputs (k=0 or k=n) the half-width
    // can push the interval infinitesimally past [0, 1]. Clip so the
    // returned CI is always a valid probability range.
    if (lo < 0.0) lo = 0.0;
    if (hi > 1.0) hi = 1.0;

    return VAFEstimate{
        .point = static_cast<float>(p),
        .ci_low = static_cast<float>(lo),
        .ci_high = static_cast<float>(hi),
    };
}

}  // namespace branch::backend
