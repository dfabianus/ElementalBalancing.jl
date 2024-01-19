using ElementalBalancing
using Test
using LinearAlgebra
using CSV
using DataFrames
using Plots

path_to_data = joinpath(@__DIR__, "..", "data", "online_preprocessed.csv")
raw_data = CSV.read(path_to_data, DataFrame) 

# define the parameters of the balancing problem
E = [1 1 1 0;               # elemental matrix E
    4.113 4 0 -4]
i_known = [2,3,4]           # indices of known rates
i_unknown = [1]             # indices of unknown rates (automatically inferred)
M = [26.5, 30, 44, 32]      # molecular weights of the species
F = diagm([1, 1, 3])        # variance covariance matrix

prob = EBALProblem(E, M, i_known)

r = [[-QS, QCO2, QO2] for (QS, QCO2, QO2) in zip(raw_data.Q_S, raw_data.Q_CO2, raw_data.Q_O2)]
sol = map(x -> solve(prob, x, F), r)
rm_best = [measurement.rm_best_mol for measurement in sol]
rc_best = [measurement.rc_best_mol for measurement in sol]

@test sum([ri[1] for ri in rm_best]) < 4749
@test sum([ri[2] for ri in rm_best]) > 2858
@test sum([ri[3] for ri in rm_best]) < -2804
@test sum([ri[1] for ri in rc_best]) > 1891

### PLOTS ###
#plot(raw_data.time, [(ri./ri')[1,3] for ri in rm_best], label="Y_CO", linewidth = 2)

# p = plot(raw_data.time, [ri[1] for ri in r], label="Q_S", linewidth = 2,
#     size=(500,350),
#     legendfontsize = 9,
#     titlelocation = :left,
#     bottom_margin=10Plots.px,
#     left_margin=10Plots.px,
#     tickfontsize = 10,
#     xlabelfontsize = 10,
#     ylabelfontsize = 10,
#     grid = false,
#     framestyle = :box,
#     legend=:bottomleft,
#     ylabel = "molar rates",
#     xlabel = "Time in hours",
#     title = "Elemental balancing rate estimation")
# plot!(raw_data.time, [ri[1] for ri in rm_best], label="Q_S_reconciled", c=:black, linestyle = :dash)
# plot!(raw_data.time, [ri[2] for ri in r], label="Q_CO2", linewidth = 2)
# plot!(raw_data.time, [ri[2] for ri in rm_best], label="Q_CO2_reconciled", c=:black, linestyle = :dash)
# plot!(raw_data.time, [ri[3] for ri in r], label="Q_O2", linewidth = 2)
# plot!(raw_data.time, [ri[3] for ri in rm_best], label="Q_O2_reconciled", c=:black, linestyle = :dash)
# plot!(raw_data.time, [ri[1] for ri in rc_best], label="Q_X_reconciled", c=:grey, linewidth = 2)
# savefig(p, joinpath(@__DIR__, "..", "fig", "soft_sensing_test_1.svg"))

# p2 = plot(raw_data.time, -[ri[1] for ri in rc_best]./[ri[1] for ri in rm_best], 
#     label="Y_XS", linewidth = 2, ylim = (0, 1),
#     size=(500,350),
#     legendfontsize = 9,
#     titlelocation = :left,
#     bottom_margin=10Plots.px,
#     left_margin=10Plots.px,
#     tickfontsize = 10,
#     xlabelfontsize = 10,
#     ylabelfontsize = 10,
#     grid = false,
#     framestyle = :box,
#     legend=:bottomleft,
#     ylabel = "yield [mol X/mol S]",
#     xlabel = "Time in hours",
#     title = "Elemental balancing rate estimation")
# savefig(p2, joinpath(@__DIR__, "..", "fig", "soft_sensing_test_1_yield.svg"))