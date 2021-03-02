# Toxic oxygen Model, phase plane
# 2/18/21

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
Ec = 42.98736*rand()
β = 30.0*rand()
Ac = 0.0139
nc = 1.0
rf = 18.0036*rand()
d = 0.8016*rand()
μ = 5.0273*rand()
k = 1.0e9
η = 6.76e-9*rand()
λ = 0.22
ϵ = 0.0
q = 0.5*rand()

p = [Ec,Ac,nc,rf,d,ϵ,μ,η,k,λ,β,q]

y0 = [0.9*k, 0.02*k, 0.22] # climax only equilibrium

tspan = (0.0,250.0)

prob = ODEProblem(cf_ode!,y0,tspan,p)
sol = solve(prob)

c = sol[1,:]
f = sol[2,:]
w = sol[3,:]

p1 = plot()
for i = 1:50
    global p; global tspan
    # p = [r*rand(),d*rand(),μ*rand(),η*rand(),k,λ*rand(),β*rand(),q*rand(),b*rand()]

    c0 = k*rand()
    f0 = k*rand()
    y0 = [c0, f0, 0.22]

    prob = ODEProblem(cf_ode!,y0,tspan,p)
    sol = solve(prob)

    c = sol[1,:]
    f = sol[2,:]
    w = sol[3,:]

    # p1 = plot!(f,c, xlabel = "f", ylabel = "c", legend=:false)
    p1 = plot!(sol,vars=(1,2,3), legend=:false)
    # p1 = scatter!([f[end]], [c[end]], [w[end]])
    display(p1)

end
display(p1)
