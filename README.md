# ElementalBalancing

[![Build Status](https://github.com/dfabianus/ElementalBalancing.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dfabianus/ElementalBalancing.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/dfabianus/ElementalBalancing.jl.svg?branch=main)](https://travis-ci.com/dfabianus/ElementalBalancing.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/dfabianus/ElementalBalancing.jl?svg=true)](https://ci.appveyor.com/project/dfabianus/ElementalBalancing-jl)
[![Coverage](https://codecov.io/gh/dfabianus/ElementalBalancing.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/dfabianus/ElementalBalancing.jl)
[![Coverage](https://coveralls.io/repos/github/dfabianus/ElementalBalancing.jl/badge.svg?branch=main)](https://coveralls.io/github/dfabianus/ElementalBalancing.jl?branch=main)

![](fig/application.svg)

This package implements the elemental balancing concept from [1] in order to calculate unmeasured reaction rates based on other measured rates.

### The algorithm
An excerpt from the source code, where the calculation functions are described.
```julia
# Equations from the book chapter: Gregory N., Material balances and data consistency [1]:
# https://edisciplinas.usp.br/pluginfile.php/5381165/mod_resource/content/1/Chap4.pdf
solve_rc(Ec, Em, rm) = -pinv(Ec) * Em * rm              # eq. 4.14
redundancy(Ec, Em) = Em-Ec*pinv(Ec)*Em                  # eq. 4.17
eps(Rred, rm) = Rred * rm                               # eq. 4.20
P(Rred, F) = Rred * F *Rred'                            # eq. 4.24
delta(Rred, rm, F, P) = (F*Rred'*inv(P) * Rred)* rm     # eq. 4.26
rm_best(rm, delta) = rm-delta                           # eq. 4.27
h(eps, P) = eps' * inv(P) * eps                         # eq. 4.29

# Transfering the equations to the EBALProblem type
solve_rc(prob::EBALProblem, rm)     = solve_rc(Ec(prob), Em(prob), rm)
redundancy(prob::EBALProblem)       = redundancy(Ec(prob), Em(prob))
eps(prob::EBALProblem, rm)          = eps(reduced(redundancy(prob)), rm)
P(prob::EBALProblem, F)             = P(reduced(redundancy(prob)), F)
delta(prob::EBALProblem, rm, F)     = delta(reduced(redundancy(prob)), rm, F, P(prob, F))
rm_best(prob::EBALProblem, rm, F)   = rm_best(rm, delta(prob, rm, F))
h(prob::EBALProblem, rm, F)         = h(eps(prob, rm), P(prob, F))
```

### Application example
First, we define the parameters of the elemental balance $E$, $M$ and $F$ and put them into an `EBALProblem` object.
```julia
# define the parameters of the balancing problem
E = [1 1 1 0;               # elemental matrix E
    4.113 4 0 -4]
i_known = [2,3,4]           # indices of known rates
i_unknown = [1]             # indices of unknown rates
M = [26.5, 30, 44, 32]      # molecular weights of the species
F = diagm([1, 1, 3])        # variance covariance matrix

prob = EBALProblem(E, M, i_known, i_unknown)
```
Then we load the test data, extract the measured reaction rates, in this case $Q_S$, $Q_{CO2} and $Q_{O2}$ and `solve` the problem for each measurement point in `r`.
```julia
raw_data = CSV.read(path_to_data, DataFrame) 
r = [[-QS, QCO2, QO2] for (QS, QCO2, QO2) in zip(raw_data.Q_S, raw_data.Q_CO2, raw_data.Q_O2)]
sol = map(x -> solve(prob, x, F), r)
```
Extracting the results and plotting the resulting estimated rates against the raw rate measurements gives:
```julia
rm_best = [measurement.rm_best_mol for measurement in sol]
rc_best = [measurement.rc_best_mol for measurement in sol]
plot(raw_data.time, [ri[1] for ri in r], label="Q_S", linewidth = 2,
    size=(500,350),
    legendfontsize = 9,
    titlelocation = :left,
    bottom_margin=10Plots.px,
    left_margin=10Plots.px,
    tickfontsize = 10,
    xlabelfontsize = 10,
    ylabelfontsize = 10,
    grid = false,
    framestyle = :box,
    legend=:bottomleft,
    ylabel = "molar rates",
    xlabel = "Time in hours",
    title = "Elemental balancing rate estimation")
plot!(raw_data.time, [ri[1] for ri in rm_best], label="Q_S_reconciled", c=:black, linestyle = :dash)
plot!(raw_data.time, [ri[2] for ri in r], label="Q_CO2", linewidth = 2)
plot!(raw_data.time, [ri[2] for ri in rm_best], label="Q_CO2_reconciled", c=:black, linestyle = :dash)
plot!(raw_data.time, [ri[3] for ri in r], label="Q_O2", linewidth = 2)
plot!(raw_data.time, [ri[3] for ri in rm_best], label="Q_O2_reconciled", c=:black, linestyle = :dash)
plot!(raw_data.time, [ri[1] for ri in rc_best], label="Q_X_reconciled", c=:grey, linewidth = 2)
```

![plot_1](fig/soft_sensing_test_1.svg)

**Further reading**:
[1] https://edisciplinas.usp.br/pluginfile.php/5381165/mod_resource/content/1/Chap4.pdf 
