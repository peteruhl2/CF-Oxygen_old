# Read simulation averages from files and plot
# 1/13/21

using CSV, Plots, Statistics

# get average time data ========================================================
cd("C:\\Users\\peter\\OneDrive\\Documents\\GitHub\\CF-Oxygen\\Spatial Model\\Simulations")
ExtTimes = CSV.read("AvgTimes.csv"; header=true)
times = names(ExtTimes)
# ExtTimes = convert(Matrix{Float64}, AvgTimes)


# Make dictionary of average times to extinction
MTimes = Dict()

for i in 1:length(times)
    MTimes[times[i]] = mean(ExtTimes[!, i])
end


# get simulation data ==========================================================
path = "C:\\Users\\peter\\OneDrive\\Documents\\GitHub\\CF-Oxygen\\Spatial Model\\Simulations"
G = ["0.0005","0.005","0.9","0.35","0.65","0.75","0.95","1.1","1.2","0.05","0.5","0.8"]
g = G[12]
# sim_dir = "\\" * "Sim g = 0.0005, k = 1600"
sim_dir = "\\" * "Sim g = " * g * ", k = 1600"

cd(path*sim_dir)

# stuff = CSV.read("AB_avgs.csv", Matrix, skipto = 2)
stuff = CSV.read("AB_avgs.csv")
stuff = convert(Matrix{Float64}, stuff)

c_25 = stuff[:,1]
c_avg = stuff[:,2]
c_975 = stuff[:,3]

f_25 = stuff[:,4]
f_avg = stuff[:,5]
f_975 = stuff[:,6]


# plots here ===================================================================
# mean series
p = plot!(c_avg, color = :blue, label = "C ABM Mean",lw = 2)
p = plot!(f_avg, color = :red, label = "F ABM Mean",lw = 2, title = "g = $g")
p = plot!(xlabel = "Time (hours)", ylabel = "Population Density")

# # 95th percentiles
# p = plot!(c_25, color = :blue, label = "C ABM 95th percentile",lw = 0.5, linestyle = :dashdot)
# p = plot!(c_975, color = :blue, label = "",lw = 0.5, linestyle = :dashdot)
# p = plot!(f_25, color = :red, label = "F ABM 95th percentile", lw = 0.5, linestyle = :dashdot)
# p = plot!(f_975, color = :red, label = "",lw = 0.5, linestyle = :dashdot, legend=:right)
p = vline!([MTimes[g]], label = "Mean extinction time")
p = plot!(legend = false)
display(p)
