# Toxic oxygen Model, Quasi steady state oxygen
# 2/16/21

using Plots, DifferentialEquations

function cf_ode!(yp, y, p, t)
    r, d, η, k, λ, β, q, b = p
    c, f = y

    # climax has non constant growth
    # yp[1] = β*λ*c/(b+η*c)*(1 - (c + f)/k) - d*c
    # yp[2] = r*f*(1 - (f + c)/k) - d*f - q*λ*f/(b+η*c)

    w = λ / (b + η * c)
    yp[1] = β * w * c * (1 - (c + f) / k) - d * c
    yp[2] = r * f * (1 - (f + c) / k) - d * f - q * f * w

end


# ode parameters
β = 40.5 * rand()
r = 6.2 * rand()
d = 4.58 * rand()
b = 0.28 * rand()
k = 1.0e10 * rand()
η = 9e-8 * rand()
λ = 0.05 * rand()
q = 13.994 * rand()

# b = 0.999*q*λ/(r-d)

p = [r, d, η, k, λ, β, q, b]

# y0 = [0.00000008 * k * 0.000000001, 0.2 * k]

y0 = [0.00000008 * k * 0.000000001*rand(), 0.2 * k]
# rn = rand()
# if rn > 0.1
#     y0[1] = 0
# end

tspan = (0.0, 2000.0)

prob = ODEProblem(cf_ode!, y0, tspan, p)
sol = solve(prob)

c = sol[1, :]
f = sol[2, :]

# p1 = plot!(f,c, xlabel = "f", ylabel = "c", legend=:false)
# p1 = scatter!([f[end]], [c[end]])
# display(p1)






# p1 = plot()
# for i = 1:50
#     global p; global tspan
#     # p = [r*rand(),d*rand(),μ*rand(),η*rand(),k,λ*rand(),β*rand(),q*rand(),b*rand()]
#
#     fstar = (k/r)*(r - d - q*λ/b)
#     c0 = rand()
#     f0 = fstar*rand()
#     y0 = [c0, f0, 0.22]
#
#     prob = ODEProblem(cf_ode!,y0,tspan,p)
#     sol = solve(prob)
#
#     c = sol[1,:]
#     f = sol[2,:]
#
#     p1 = plot!(f,c, xlabel = "f", ylabel = "c", legend=:false)
#     p1 = scatter!([f[end]], [c[end]])
#     # display(p1)
#
# end
# display(p1)



# these are exclusion steady states
cstar = (k * β * λ - b * d * k) / (d * k * η + β * λ)
fstar = (k / r) * (r - d - q * λ / b)

# coexistence steady states
cx1 =
    (
        -2 * b * d * r * η +
        d * β * η * λ +
        sqrt(d) * sqrt(β) * sqrt(4 * q * r + d * β) * η * λ
    ) / (2 * d * r * η^2)
cx2 =
    (
        -2 * b * d * r * η + d * β * η * λ -
        sqrt(d) * sqrt(β) * sqrt(4 * q * r + d * β) * η * λ
    ) / (2 * d * r * η^2)


p = plot(sol.t[c.>0], c[c.>0], label = "Climax", lw = 2, yaxis = :log)
p = plot!(
    sol.t[f.>0],
    f[f.>0],
    label = "Attack",
    lw = 2,
    legend = :right,
    xlabel = "Time (days)",
    ylabel = "Population",
    yaxis = :log,
)
# p = plot(sol.t[c.>0],c[c.>0], label = "Climax", lw = 2)
# p = plot!(sol.t[f.>0],f[f.>0], label = "Attack", lw = 2, legend=:topright, xlabel = "Time (days)", ylabel = "Population Density")
