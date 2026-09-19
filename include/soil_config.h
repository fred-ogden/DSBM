#ifndef SOIL_CONFIG_H
#define SOIL_CONFIG_H

// Compile-time discretization count.
#ifndef NDISC
#error "NDISC must be defined at compile time (number of soil discs)."
#endif

// Fictitious lower bound on water content (m3/m3).
#ifndef THETA_MIN
#define THETA_MIN 1.0e-03
#endif

// Tiny epsilon to avoid divide-by-zero.
#ifndef SOIL_EPS
#define SOIL_EPS 1.0e-12
#endif

// Inline helper macro
#ifndef SOIL_INLINE
#define SOIL_INLINE static inline
#endif

#endif // SOIL_CONFIG_H
