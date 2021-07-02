# read CF averages and 2.5th and 97.5th percentile data from
# 'AB_avgs.csv'
# updated 6/2020

using CSV, Plots, DifferentialEquations

cd("C:\\Users\\peter\\OneDrive\\Documents\\GitHub\\CF-Oxygen\\Non-spatial Model")

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
#=============================================================================#

function get_rc(ox)
    β = 16.64/36; b = 13.4256; n = 2.6626; rcmin = 0.0;
    return rc = β*ox.^n./(b^n .+ ox.^n) + rcmin
end

function get_rf(ox)
    β = 16.64/36; b = 13.4256; n = 2.6626; rcmin = 0.0;
    return rc = β*(1 - ox.^n./(b^n .+ ox.^n)) + rcmin
end

function get_dc(t)
    t_b = 19*36; dn = 0.6045/36; dbs = 6.7686/36

    if t < t_b
        return dn + dbs
    else
        return dn
    end
end

function get_df(t,w)
    t_b = 19*36; t_c = 33*36; γ = 0.8976; dn = 0.6045/36
    dbs = 6.7686/36; q = 3.27e-5/36; ϵ = 1.21/36

    # base rate
    df = dn + q*w

    # broad spectrum rate
    if t < t_b
        df += γ*dbs
    end

    # clindamycin rate
    if t > t_c
        df += ϵ
    end

    return df
end

# broad spectrum antibiotic function
function BrSpec(t,p)
    global t_b
    if t < t_b
        return p[5];
    else
        return 0.;
    end
end

#==============================================================================#
# Solve ODE ===================================================================#
function cf_ode(yp,y,p,t)
    # Ec,Ac,nc,rf,d,ϵ,μ,η,k = p
    β,b,n,dn,dbs,γ,ϵ,μ,k,η,λ = p
    c,f,w = y

    # get broad spectrum antibiotic amount
    dbs = BrSpec(t,p)

    t_c = 33
    ϵ = 0
    if t >= t_c
        ϵ = p[7]
    end

    # set up death rates
    dc = dn + dbs;
    df = dn + γ*dbs;

    # climax has non constant growth
    yp[1] = (β*w^n/(b^n + w^n))*c*(1 - c - f) - dc*c
    yp[2] = (β*(1 - w^n/(b^n + w^n)))*f*(1 - f - c) - df*f - ϵ*f - q*w*f
    yp[3] = λ - μ*w - η*k*c*w

    # println(q*w*f)
    # println(get_rc(x)," ", get_rf(x))
end

# =============================================================================#
# treatment times and initial ode population
t_b = 19
t_c = 33
N0 = 6.7e8

# ode parameters
# Ec = 11.7549; Ac = 0.0139; nc = 1.4408
# rf = 14.3436; d = 0.7016; ϵ = 0.5428; μ = 1.4273;
# k = 10000
# η = 0.8176/k

# p = [Ec,Ac,nc,rf,d,ϵ,μ,η,k]

β = 16.64; b = 13.4256; n = 2.6626;
dn = 0.6045; dbs = 6.7686; γ = 0.8976;
ϵ = 1.21; μ = 200*23*60*24;
q = 3.27e-5
# k = 10^10
# η = 3.16e-4

# adjusted k and η, did this so that k*η is constant
k = 10000
η = 3.16e2
w0 = 14.63
λ = μ*w0

p = [β,b,n,dn,dbs,γ,ϵ,μ,k,η,λ]

frac = 0.8654
# c0 = frac*N0/k
# f0 = (1 - frac)*N0/k

# hard coded initial conditions
# did this for the adjusted η
c0 = 0.0579818
f0 = 0.0090182
y0 = [Complex(c0), Complex(f0), Complex(w0)]
tspan = (0.0,40.0)
prob = ODEProblem(cf_ode,y0,tspan,p)
sol = solve(prob)

# get time series
C_temp = real.(sol[1,:])
F_temp = real.(sol[2,:])
Ct = C_temp./(C_temp + F_temp)
Ft = F_temp./(C_temp + F_temp)

# get relative abundances from ABM numbers
C_25 = c_25./(c_25 + f_25)
C_avg = c_avg./(c_avg + f_avg)
C_975 = c_975./(c_975 + f_975)

F_25 = f_25./(c_25 + f_25)
F_avg = f_avg./(c_avg + f_avg)
F_975 = f_975./(c_975 + f_975)



### Plots =====================================================================#
### adjust ode time
t_ode = sol.t*36

p = plot(C_avg, color = :blue, label = "C ABM Mean",lw = 2)
p = plot!(C_25, color = :blue, label = "C ABM 95th percentile",lw = 0.5, linestyle = :dashdot)
p = plot!(C_975, color = :blue, label = "",lw = 0.5, linestyle = :dashdot)
p = plot!(F_avg, color = :red, label = "F ABM Mean",lw = 2)
p = plot!(F_25, color = :red, label = "F ABM 95th percentile", lw = 0.5, linestyle = :dashdot)
p = plot!(F_975, color = :red, label = "",lw = 0.5, linestyle = :dashdot)
# p = plot!(36*sol.t,Ct, label = "C ODE", lw = 2)
# p = plot!(36*sol.t,Ft, label = "F ODE", lw = 2, legend=:left, xlabel = "Time (days)", ylabel = "Population Density")
p = plot!(xticks = (0:180:1440, ["0","5","10","15","20","25","30","35","40"]))
# p = scatter!(tdata*36, cdata, label = "C Data", markershape = :hexagon, color = :blue)
# p = scatter!(tdata*36, fdata, label = "F Data", markershape = :hexagon, color = :red, legendfontsize=7)
p = plot!(legend=:left)

display(p)


p = plot(c_avg, color = :blue, label = "C ABM Mean",lw = 2)
p = plot!(c_25, color = :blue, label = "C ABM 95th percentile",lw = 0.5, linestyle = :dashdot)
p = plot!(c_975, color = :blue, label = "",lw = 0.5, linestyle = :dashdot)
p = plot!(f_avg, color = :red, label = "F ABM Mean",lw = 2)
p = plot!(f_25, color = :red, label = "F ABM 95th percentile", lw = 0.5, linestyle = :dashdot)
p = plot!(f_975, color = :red, label = "",lw = 0.5, linestyle = :dashdot)
p = plot!(xticks = (0:180:1440, ["0","5","10","15","20","25","30","35","40"]))
p = plot!(legend=:left)
display(p)
