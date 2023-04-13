# Python Codes for Avocado #

## ntt ##
- `ntt_tablegen`: Generate zeta table for ntt
- `primitive_root`: Find primitive root by exhaustive search

---

## pack ##
- `pack_codegen`: Generate codes doing packing & unpacking variables to `uint8_t`
  - How to use: Change `s`, `e`, `ename` in the top of the code
    - `s`: start of variable range (ex) -2
    - `e`: end of variable range (ex) 2
    - `ename` (if range contains negative value): variable max value (ex) "ETA"

---

## plot_distribution ##
- `plot_gaussian`: Plot data ("gaussian.txt") and gaussian distribution with $\mu = 0, \sigma = \sqrt{1 / 2}$
- `plot_hyperball_radius`: Plot data ("hyperball_radius.txt") and radius distribution
- `plot_uniformdbl`: Plot data ("uniformdbl.txt") and uniform distribution using 32bit precision with python
- `plot_z_radius`: Plot data ("z_radius.txt") and radius distribution

---

## ziggurat ##
All files are about gaussian distribution with $\mu = 0, \sigma = 1$
- `binsearch_x1`: Find ziggurat start position `x_1` by binary search
- `tablegen`: Generate table according to `x_1`. This generates variables `x`, `y`, `k`
  - `x`: x coordinate multiplied by $2^{-31}$
  - `y`: y coordinate
  - `k`: ratio of `x`s multiplied by $2^{31}$
- `ziggurat`: Compute gaussian distribution with ziggurat algorithm and compare it with gaussian distribution expected by plotting

---

## gen ##
- `genJ`: Generate public nonzero `j` with packed unsigned 32bit number

---
