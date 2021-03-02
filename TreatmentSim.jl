### Solver for messing with different treatment times
### 2/8/21

using Plots, DifferentialEquations


function Treat(f, t, t_start, t_end)
    if t > t_start && t < t_end
        return 0.8
    else
        return 0.0
    end
end


function cf_ode!(yp,y,p,t)
    Ec,Ac,nc,rf,d,ϵ,μ,η,k,q = p
    c,f,w = y

    α = 1
    β = 1
    λ = 0.0094

    t_start = 25
    t_end = t_start + 5
    ϵ = Treat(f, t, t_start, t_end)

    # climax has non constant growth
    yp[1] = (Ec*w^nc/(Ac^nc + w^nc))*c*(1 - (f + c)) - d*c
    yp[2] = rf*f*(1 - (f + c)) - d*f - ϵ*f - q*w*f
    yp[3] = λ - μ*w - η*k*c*w
    # println(Treat(t,f))
end

# ode parameters
Ec = 11.9435; Ac = 0.0007; nc = 1.0
rf = 26.2009; d = 0.2249; ϵ = 0.0; μ = 1.0;
q = 0.2117
k = 1e8
η = 3.3e-16
λ = 0.0094

p = [Ec,Ac,nc,rf,d,ϵ,μ,η,k,q]

### different initial condtions
y0 = [0.8913,0.0028,0.0115] # best fitting ode parameters
# y0 = [0.0000008789,(rf-d)/rf,λ/μ] # attack only equilibrium

c_only = (Ec*λ - d*λ -Ac*d*μ)/(Ac*d*k*η + Ec*λ)
w_c = (Ac*d*k*η + Ec*λ)/(Ec*k*η + Ec*μ - d*k*η)
# y0 = [c_only, 0.2, w_c] # climax only equilibrium

tspan = (0.0,40.0)

prob = ODEProblem(cf_ode!,y0,tspan,p)
sol = solve(prob)

c = sol[1,:]
f = sol[2,:]

### stability conditions
x1 = Ec - rf*(1 + Ac*(μ/λ)) # this one is good I think, c invasion condition
x2 = rf*Ac*k*η + rf*λ + rf*Ac*μ - Ac*d*k*η - Ec*λ # this might be right for the f condition

### stability another way
y1 = rf*((Ac/λ)*(k*η + μ) + 1) - (Ac/λ)*d*k*η
y2 = rf*(1 + Ac*(μ/λ))


# p = plot(sol.t[c.>0],c[c.>0], label = "C ODE", lw = 2, yaxis=:log)
# p = plot!(sol.t[f.>0],f[f.>0], label = "F ODE", lw = 2, legend=:right, xlabel = "Time (days)", ylabel = "Population Density",yaxis=:log)
p = plot!(sol.t[c.>0],c[c.>0], label = "C ODE", lw = 2)
p = plot!(sol.t[f.>0],f[f.>0], label = "F ODE", lw = 2, legend=:left, xlabel = "Time (days)", ylabel = "Population Density")


display(p)

# println("C_invade is ", C_invade," F_invade is ", F_invade)
println(x1)
println(x2)
println()
println(y1 > Ec)
println(Ec > y2)
println()
