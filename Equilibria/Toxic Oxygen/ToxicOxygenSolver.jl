# Toxic oxygen Model
# 2/12/21

using Plots, DifferentialEquations

function cf_ode!(yp,y,p,t)
    Ec,Ac,nc,rf,d,ϵ,μ,η,k,λ,β,q = p
    c,f,w = y

    t_treat = 28
    ϵ = 0
    if t >= t_treat
        ϵ = p[6]
    end

    # climax has non constant growth
    yp[1] = β*w*c*(1 - (c + f)/k) - d*c
    yp[2] = rf*f*(1 - (f + c)/k) - d*f - ϵ*f - q*f*w
    yp[3] = λ - μ*w - η*c*w

end

# ode parameters
Ec = 15.98736; Ac = 0.0139; nc = 1.0
rf = 1.0036; d = 0.6016; μ = 0.07273;
k = 1.0e8
η = 4.1176e-6
λ = 0.22
ϵ = 0.0
β = 50.0
q = 20.0 

p = [Ec,Ac,nc,rf,d,ϵ,μ,η,k,λ,β,q]

y0 = [0.00009*k, 0.02*k, 0.22] # climax only equilibrium

tspan = (0.0,250.0)

prob = ODEProblem(cf_ode!,y0,tspan,p)
sol = solve(prob)

c = sol[1,:]
f = sol[2,:]
w = sol[3,:]


p = plot(sol.t[c.>0],c[c.>0], label = "C ODE", lw = 2, yaxis=:log)
p = plot!(sol.t[f.>0],f[f.>0], label = "F ODE", lw = 2, legend=:right, xlabel = "Time (days)", ylabel = "Population Density",yaxis=:log)
# p = plot(sol.t[c.>0],c[c.>0], label = "C ODE", lw = 2)
# p = plot!(sol.t[f.>0],f[f.>0], label = "F ODE", lw = 2, legend=:right, xlabel = "Time (days)", ylabel = "Population Density")


display(p)
