# Verification scripts: E1-E2 interference in X-ray resonant scattering from powders

Two short Python scripts that check, independently, the closed-form expressions for the E1-E2 (dipole-quadrupole) interference term of the X-ray resonant scattering cross section from a powder sample. The same expressions appear in:

- the article *Circular dichroism in the X-ray resonant scattering from powders* (D. Sier, M. O. J. Y. Hunault, A. Juhin, C. Brouder and V. Vercamer);
- Appendix F, *Theoretical developments with LFMT*, section *Isotropic RIXS: the powder formula*, of V. Vercamer's PhD thesis, *Spectroscopic and Structural Properties of Iron in Silicate Glasses*.

The scripts do not reuse either derivation: they rebuild each quantity from its definition (Wigner symbols, recoupling of spherical tensors) and compare it with the closed form.

## Requirements

Python 3.9 or later, with `sympy` and `numpy`:

```sh
pip install sympy numpy
python check9j.py
python checkgamma.py
```

or, without installing anything globally:

```sh
uv run --with sympy --with numpy python check9j.py
uv run --with sympy --with numpy python checkgamma.py
```

Each script runs in a few seconds and takes no argument.

## check9j.py: the coefficients sigma(b,u), from 9j to 6j symbols (exact)

Formula checked: the reduction of the general coefficient, a sum over g of products of three Wigner 9j symbols with the reduced matrix elements h(l) and their phase, to the closed form with two 6j symbols,

sigma(b,u) = sum over g = 1, 2 of (i sqrt(5)/2) sqrt((2b+1)(2g+1)(2u+1)) {1 1 g; 2 1 b} {2 b 1; u 1 1} [S(g)* + (-1)^(b+u) S(g)],

for the five couples (b,u) = (1,0), (1,1), (1,2), (2,1), (2,2).

How: `sympy.physics.wigner` evaluates every 9j and 6j symbol exactly. The structure factors stay symbolic, S(g) = R_g + i I_g for g = 1, 2. The script also checks that the g = 0 and g = 3 terms vanish, and prints a `g ... nonzero` line if one does not.

Expected output: one line per couple, the general expression on the left (`full`, printed in floating point) and the closed form on the right (`closed-form`, exact). They must be equal:

```
b=1 u=0  full: 0.166666666666667*sqrt(5)*I1 - 0.166666666666667*sqrt(15)*I2   closed-form: sqrt(5)*I1/6 - sqrt(15)*I2/6
b=1 u=1  full: 0.0833333333333333*sqrt(15)*I*R1 - 0.25*sqrt(5)*I*R2   closed-form: sqrt(15)*I*R1/12 - sqrt(5)*I*R2/4
b=1 u=2  full: 0.0833333333333333*I1 - 0.0833333333333333*sqrt(3)*I2   closed-form: I1/12 - sqrt(3)*I2/12
b=2 u=1  full: 0.75*I1 + 0.25*sqrt(3)*I2   closed-form: 3*I1/4 + sqrt(3)*I2/4
b=2 u=2  full: 0.25*sqrt(3)*I*R1 + 0.25*I*R2   closed-form: sqrt(3)*I*R1/4 + I*R2/4
```

This is the check that fixes the signs of sigma(2,1) and of sigma(1,0) + sigma(1,2).

## checkgamma.py: angular factors and polarisation averages (numerical)

Formulas checked:

1. The angular factors gamma(b,u), built by recoupling the spherical components of the incident polarisation e, its conjugate e* and the wave vector k with those of the scattered polarisation e_s, against their explicit forms:
   - gamma(1,0) = P_c^s (k.k_s) / (3 sqrt(2)), gamma(1,1) = 0, gamma(1,2) = -P_c^s (k.k_s) / (3 sqrt(10));
   - gamma(2,1) = P_c / sqrt(10) (1/3 - |k.e_s|^2);
   - gamma(2,2) = -i / (2 sqrt(30)) times the sum of four mixed products (e.e_s*)(e*×e_s).k + ...;
   where P_c and P_c^s are the circular polarisation rates of the incident and scattered beams.
2. The averages over the scattered polarisation, taken as the mean over two orthogonal linear polarisations perpendicular to k_s: <gamma(1,0)> = 0, <gamma(2,1)> = -P_c (1 - 3 (k.k_s)^2) / (6 sqrt(10)), and <gamma(2,2)> = i / (2 sqrt(30)) [(e.k_s)(e*×k_s).k + (e*.k_s)(e×k_s).k].
3. The same bracket in the laboratory frame (k along z, incident polarisation set by chi and psi, k_s by the angles alpha and beta): it equals cos(2 chi) sin^2(alpha) sin(2(beta - psi)), and e × e* = i sin(2 chi) z.

How: random unit vectors k and k_s and random complex polarisations perpendicular to them (fixed seed, so every run gives the same numbers), Clebsch-Gordan coefficients from sympy, everything else in numpy.

Expected output: three trials; on each line the recoupled value comes first and the explicit formula second. They agree to about 1e-16, the leftover real or imaginary parts of order 1e-18 being rounding. The last two lines check the laboratory-frame formula and e × e*.

## Reading a mismatch

A difference on one line points to the matching closed form (a sign, a factor sqrt(2) or sqrt(10), a missing term). With the expressions as they stand in the article and in Appendix F, every line agrees.

## License

MIT, see [LICENSE](LICENSE).
