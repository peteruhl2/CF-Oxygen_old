# Full CF ode solver for treatment simulations
# 4/23/21

cd("C:\\Users\\peter\\OneDrive\\Documents\\GitHub\\CF-Oxygen\\Simulations")

using Plots, DifferentialEquations

function cf_ode!(yp,y,p,t)
    # display(t)
    β,b,n,r,k,dbs,dn,α,ϵ,q,μ,η,x0 = p
    c,f,x = y

    λ = μ*x0
    dbs = BrSpec(t,p)

    ### total death rates
    dc = dn + dbs
    df = dn + α*dbs

    ## antibiotics targeting attack
    ϵ = 0
    if t >= t_c
        ϵ = p[9];
    end

    rc = (β*x^n/(b^n + x^n))
    rf = (r + β*(1 - x^n/(b^n + x^n)))
    # display(rc)

    # climax has non constant growth
    yp[1] = rc*c*(1 - (c+f)/k) - dc*c
    yp[2] = rf*f*(1 - (c+f)/k) - df*f - ϵ*f - q*f*x
    yp[3] = λ - μ*x - η*c*x

end

function BrSpec(t,p)
    global t_b
    if t < t_b
        dbs = p[6];
    else
        dbs = 0;
    end
end

# treatment times
global t_b
global t_c

t_b = 19 + 10
t_c = 33*Inf
N0 = 6.7e8

# parameters
r = 0.0046
β = 16.6388
b = 13.4256
n = 2.6626
# n = 2

dn = 0.6045
dbs = 6.7686
α = 0.8976

ϵ = 1.2124
μ = 200*23*60*24

k = 10^10
η = 3.1611e-4
q = 3.2747e-5

frac = 0.8659
x0 = 14.6287
λ = μ*x0

p = [β,b,n,r,k,dbs,dn,α,ϵ,q,μ,η,x0]

c0 = frac*N0
f0 = (1 - frac)*N0


tspan = (0.0,40.0)
y0 = [Complex(c0), Complex(f0), Complex(x0)]

prob = ODEProblem(cf_ode!,y0,tspan,p)
sol = solve(prob)

c = real.(sol[1,:])
f = real.(sol[2,:])

Ct = real.(c./(c .+ f))
Ft = real.(f./(c .+ f))

# p = plot(sol.t[c.>0],c[c.>0], label = "Climax", lw = 2, yaxis=:log)
# p = plot!(sol.t[f.>0],f[f.>0], label = "Attack", lw = 2, legend=:right, xlabel = "Time (days)", ylabel = "Population Density",yaxis=:log)
# p = plot(sol.t[c.>0],c[c.>0], label = "Climax", lw = 2)
# p = plot!(sol.t[f.>0],f[f.>0], label = "Attack", lw = 2, legend=:topright, xlabel = "Time (days)", ylabel = "Population Density")

p1 = plot(sol.t, Ct, label = "C", lw = 2)
p1 = plot!(sol.t, Ft, label = "F", lw = 2)
display(p1)
