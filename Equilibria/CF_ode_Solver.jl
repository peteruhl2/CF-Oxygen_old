### Simple solver for the CF ode model
### 1/27/21

using Plots, DifferentialEquations

function cf_ode!(yp,y,p,t)
    Ec,Ac,nc,rf,d,ϵ,μ,η,k = p
    c,f,w = y

    α = 1
    β = 1
    λ = 0.22

    t_treat = 28
    ϵ = 0
    if t >= t_treat
        ϵ = p[6]
    end

    # climax has non constant growth
    yp[1] = (Ec*w^nc/(Ac^nc + w^nc))*c*(1 - c - α*f) - d*c
    yp[2] = rf*f*(1 - f - β*c) - d*f - ϵ*f
    yp[3] = λ - μ*w - η*k*c*w

end

# ode parameters
Ec = 16.0049; Ac = 0.0139; nc = 1.0
rf = 14.0436; d = 0.7016; ϵ = 0.0; μ = 1.0273;
k = 10000
η = 0.8176/k
λ = 0.22

p = [Ec,Ac,nc,rf,d,ϵ,μ,η,k]

# y0 = [0.8283,0.0165,0.1388] # best fitting ode parameters
y0 = [0.8789,0.0165,0.1388]
tspan = (0.0,10000.0)
prob = ODEProblem(cf_ode!,y0,tspan,p)
sol = solve(prob)

c = sol[1,:]
f = sol[2,:]

### stability conditions
C_invade = Ec*λ - rf*λ  - Ac*rf*μ

B = Ac*d*k*η + Ec*λ
D = rf*(d*λ + Ac*λ*μ - Ec*λ)
F_invade = B*(rf-d) + D


p = plot(sol.t[c.>0],c[c.>0], label = "C ODE", lw = 2, yaxis=:log)
p = plot!(sol.t[f.>0],f[f.>0], label = "F ODE", lw = 2, legend=:left, xlabel = "Time (days)", ylabel = "Population Density",yaxis=:log)
display(p)

println("C_invade is ", C_invade," F_invade is ", F_invade)
