# Toxic oxygen Model, Quasi steady state oxygen
# 2/16/21

using Plots, DifferentialEquations

function cf_ode!(yp,y,p,t)
    r,d,η,k,λ,β,q,b = p
    c,f = y

    # climax has non constant growth
    # yp[1] = β*λ*c/(b+η*c)*(1 - (c + f)/k) - d*c
    # yp[2] = r*f*(1 - (f + c)/k) - d*f - q*λ*f/(b+η*c)

    w = λ/(b + η*c)
    yp[1] = β*w*c*(1 - (c + f)/k) - d*c
    yp[2] = r*f*(1 - (f + c)/k) - d*f - q*f*w

end


# ode parameters
β = 30.5*rand()
r = 30.9*rand()
d = 0.58*rand()
b = 0.28*rand()
k = 1.0e8
η = 9e-8*rand()
λ = 0.05
q = 5.794*rand()

# b = 0.999*q*λ/(r-d)

p = [r,d,η,k,λ,β,q,b]

y0 = [rand()*k, rand()*k, 0.22*rand()]

tspan = (0.0,25000.0)

prob = ODEProblem(cf_ode!,y0,tspan,p)
sol = solve(prob)

c = sol[1,:]
f = sol[2,:]

p1 = plot!(f,c, xlabel = "f", ylabel = "c", legend=:false)
p1 = scatter!([f[end]], [c[end]])
display(p1)






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

    p1 = plot!(f,c, xlabel = "f", ylabel = "c", legend=:false)
    p1 = scatter!([f[end]], [c[end]])
    display(p1)

end
display(p1)



# these are exclusion steady states
cstar = (k*β*λ - b*d*k)/(d*k*η + β*λ)
fstar = (k/r)*(r - d - q*λ/b)

# coexistence steady states
cx1 = (-2*b*d*r*η + d*β*η*λ + sqrt(d)*sqrt(β)*sqrt(4*q*r + d*β)*η*λ)/(2*d*r*η^2)
cx2 = (-2*b*d*r*η + d*β*η*λ - sqrt(d)*sqrt(β)*sqrt(4*q*r + d*β)*η*λ)/(2*d*r*η^2)

top1 = - (-2b*r + d*k*η - 2k*r*η + (sqrt(d)*k*sqrt(4q*r + d*β)*η)/sqrt(β) + β*λ + sqrt(β)*sqrt(4q*r+d*β)/sqrt(d))
fx1 = top1/2r*η

top2 = 2b*r - d*k*η + 2k*r*η + (sqrt(d)*k*sqrt(4q*r+d*β)*η)/sqrt(β) - β*λ + (sqrt(β)*sqrt(4q*r+d*β*λ))/sqrt(d)
fx2 = top2/2r*η
