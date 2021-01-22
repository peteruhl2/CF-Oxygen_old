# read CF averages and 2.5th and 97.5th percentile data from
# 'AB_avgs.csv'
# 12/17/2020

using CSV, Plots, DifferentialEquations

cd("C:\\Users\\peter\\Onedrive\\Desktop\\cyst fib\\OxygenModels")

# data =========================================================================
tdata = [0; 14; 19; 26; 28; 31; 33; 35; 38]
cdata = [0.8683; 0.6897; 0.6986; 0.3270; 0.1840; 0.4025; 0.3302; 0.8594; 0.9077]
fdata = [0.1317; 0.3103; 0.3014; 0.6730; 0.8160; 0.5975; 0.6698; 0.1406; 0.0923]
# ==============================================================================

# stuff = CSV.read("Simulations\\AB_avgs.csv", Matrix, skipto = 2)
stuff = CSV.read("Simulations\\AB_avgs.csv")
stuff = convert(Matrix{Float64}, stuff)
stuff = stuff

c_25 = stuff[:,1]
c_avg = stuff[:,2]
c_975 = stuff[:,3]

f_25 = stuff[:,4]
f_avg = stuff[:,5]
f_975 = stuff[:,6]


# Solve ODE ===================================================================#
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
Ec = 11.7549; Ac = 0.0139; nc = 1.4408
rf = 14.3436; d = 0.7016; ϵ = 0.5428; μ = 1.4273;
k = 10000
η = 0.8176/k

p = [Ec,Ac,nc,rf,d,ϵ,μ,η,k]

# y0 = [0.8283,0.0165,0.1388] # best fitting ode parameters
y0 = [0.8789,0.0165,0.1388]
tspan = (0.0,40.0)
prob = ODEProblem(cf_ode!,y0,tspan,p)
sol = solve(prob)

#==============================================================================#

# Plots =======================================================================#

p = plot(c_avg./k, color = :blue, label = "C ABM Mean",lw = 2)
p = plot!(c_25./k, color = :blue, label = "C ABM 95th percentile",lw = 0.5, linestyle = :dashdot)
p = plot!(c_975./k, color = :blue, label = "",lw = 0.5, linestyle = :dashdot)
p = plot!(f_avg./k, color = :red, label = "F ABM Mean",lw = 2)
p = plot!(f_25./k, color = :red, label = "F ABM 95th percentile", lw = 0.5, linestyle = :dashdot)
p = plot!(f_975./k, color = :red, label = "",lw = 0.5, linestyle = :dashdot)
p = plot!(24*sol.t,sol[1,:], label = "C ODE", lw = 2)
p = plot!(24*sol.t,sol[2,:], label = "F ODE", lw = 2, legend=:left, xlabel = "Time (days)", ylabel = "Population Density")
p = plot!(xticks = (0:120:960, ["0","5","10","15","20","25","30","35","40"]))
p = scatter!(tdata*24, cdata, label = "C Data", markershape = :hexagon, color = :blue)
p = scatter!(tdata*24, fdata, label = "F Data", markershape = :hexagon, color = :red, legendfontsize=7)

display(p)

vline!([24*28], label = "Start of treatment")
