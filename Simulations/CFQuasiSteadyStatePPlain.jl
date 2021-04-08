# Full CF quasi steady state phase plane
# 4/1/21

cd("C:\\Users\\peter\\OneDrive\\Documents\\GitHub\\CF-Oxygen\\Simulations")

using Plots, DifferentialEquations

function cf_ode!(yp,y,p,t)
    β,b,n,r,k,d,ϵ,q,μ,η,x0 = p
    c,f = y

    λ = μ*x0
    d = BrSpec(t,p)

    # climax has non constant growth
    x = λ/(μ + η*c)
    yp[1] = β*x*c*(1 - (c+f)/k) - d*c
    yp[2] = r*f*(1 - (c+f)/k) - d*f - ϵ*f - q*f*x

end

function BrSpec(t,p)
    global t_b
    if t < t_b
        d = p[6];
    else
        d = 0;
    end
end

global t_b
# t_b = 19 + 7*10
t_b = Inf
N0 = 6.7e8

# parameters
r = 24.6262
β = 18.7624
b = 12.4
n = 1.0

d = 0.1030

ϵ = 1.1688*0
μ = 200*23*60*24

k = 10^10
η = 4.1834e-4
q = 4.9954e-5

frac = 0.5

x0 = 17.74
λ = μ*x0

p = [β,b,n,r,k,d,ϵ,q,μ,η,x0]

c0 = frac*N0
f0 = (1 - frac)*N0

tspan = (0.0,10000.0)
y0 = [c0, f0]


p1 = plot()
for i = 1:50
    global p; global tspan; global N0
    # p = [r*rand(),d*rand(),μ*rand(),η*rand(),k,λ*rand(),β*rand(),q*rand(),b*rand()]

    # initials conditions
    c0 = rand()*N0*frac
    f0 = rand()*N0*(1 - frac)

    # # maybe take c0 = 0
    # rn = rand()
    # if rn < 0.1
    #     c0 = 0
    # end
    #
    # # maybe take f0 = 0
    # rn = rand()
    # if rn < 0.1
    #     f0 = 0
    # end

    y0 = [c0, f0]

    prob = ODEProblem(cf_ode!,y0,tspan,p)
    sol = solve(prob)

    c = sol[1,:]
    f = sol[2,:]

    p1 = plot!(f,c, xlabel = "f", ylabel = "c", legend=:false)
    p1 = scatter!([f[end]], [c[end]])
    # display(p1)

end
display(p1)
