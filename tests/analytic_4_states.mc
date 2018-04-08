/* Maxima batch file used to show analytic results for various pH
values. This is a particular solution to the diamond configuration
using values found in the configuration file 'graph.csv'. In order to
use:

$ maxima
> batch ("analytic_4_states.mc");
> float(ev(sol, pH=5));

This will give the results for [g1,g2,g3,g4]. Note that the output of
the python application gives its results as [g1,g2,g4,g3] due to the
definition of the 'state.csv' file.

`eqns` was computed manually, and can be found in the derivations
documentation. 
*/

energy(pKa) := log(10) * (pH - pKa);
variance(pKa_var) := log(10)^2 * pKa_var;

eqns : [s12 * ((g2 - g1) - d12) + s14 * ((g4 - g1) - d14) = 0, -s12 * ((g2 - g1) - d12) + s23 * ((g3 - g2) - d23) = 0, -s23 * ((g3 - g2) - d23) - s43 * ((g3 - g4) - d43) = 0, -s14 * ((g4 - g1) - d14) + s43 * ((g3 - g4) - d43) = 0];

sol : solve(eqns,[g1,g2,g3,g4]);

sol : ev(sol, d12=energy(6), d23=energy(7) , d14=energy(5.7) , d43=energy(6));
sol : ev(sol, s12=1/variance(0.1), s23=1/variance(0.5), s14=1/variance(5), s43=1/variance(0.1));
sol : ev(sol, %r1 = 0);