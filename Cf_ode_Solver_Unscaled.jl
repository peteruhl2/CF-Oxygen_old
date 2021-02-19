### Simple solver for the CF ode model, unscaled model
### 2/10/21

using Plots, DifferentialEquations

function cf_ode!(yp,y,p,t)
    Ec,Ac,nc,rf,d,ϵ,μ,η,k,λ = p
    c,f,w = y

    α = 1
    β = 1
    # λ = 0.22

    t_treat = 28
    ϵ = 0
    if t >= t_treat
        ϵ = p[6]
    end

    # climax has non constant growth
    yp[1] = (Ec*w^nc/(Ac^nc + w^nc))*c*(1 - (c + α*f)/k) - d*c
    yp[2] = rf*f*(1 - (f + β*c)/k) - d*f - ϵ*f
    yp[3] = λ - μ*w - η*c*w

end

# ode parameters
Ec = 15.98736; Ac = 0.0139; nc = 1.0
rf = 0.6036; d = 0.6016; μ = 1.27273;
k = 1.0e8
η = 8.1176e-5
λ = 0.25
ϵ = 0.0

p = [Ec,Ac,nc,rf,d,ϵ,μ,η,k,λ]

### different initial condtions
# y0 = [0.8283,0.0165,0.1388] # best fitting ode parameters
# y0 = [0.0000008789,(rf-d)/rf,λ/μ] # attack only equilibrium

c_only = (Ec*λ - d*λ -Ac*d*μ)/(Ac*d*k*η + Ec*λ)
w_c = (Ac*d*k*η + Ec*λ)/(Ec*k*η + Ec*μ - d*k*η)
y0 = [0.8*k, 0.2*k, 0.22] # climax only equilibrium

tspan = (0.0,1000000.0)

prob = ODEProblem(cf_ode!,y0,tspan,p)
sol = solve(prob)

c = sol[1,:]
f = sol[2,:]
w = sol[3,:]

### stability conditions
x1 = Ec - rf*(1 + Ac*(μ/λ)) # this one is good I think, c invasion condition
x2 = rf*Ac*k*η + rf*λ + rf*Ac*μ - Ac*d*k*η - Ec*λ # this might be right for the f condition

### stability another way
y1 = Ec/(1 + Ac*(μ/λ))
y2 = (Ac*d*k*η + Ec*λ)/(Ac*k*η + λ + Ac*μ)


p = plot(sol.t[c.>0],c[c.>0], label = "C ODE", lw = 2, yaxis=:log)
p = plot!(sol.t[f.>0],f[f.>0], label = "F ODE", lw = 2, legend=:right, xlabel = "Time (days)", ylabel = "Population Density",yaxis=:log)
# p = plot(sol.t[c.>0],c[c.>0], label = "C ODE", lw = 2)
# p = plot!(sol.t[f.>0],f[f.>0], label = "F ODE", lw = 2, legend=:right, xlabel = "Time (days)", ylabel = "Population Density")


display(p)

# println("C_invade is ", C_invade," F_invade is ", F_invade)
println(x1)
println(x2)
println()
println(y1 > rf)
println(rf > y2)
println()
println(y1 - y2)
println()
