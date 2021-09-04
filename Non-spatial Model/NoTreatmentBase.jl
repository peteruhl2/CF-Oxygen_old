### spatially homogeneous ABM with no treatment
### no ode stuff in this version
### 6/23/2021

# get to home directory
cd(@__DIR__)
# cd("C:\\Users\\peter\\Onedrive\\Desktop\\cyst fib\\OxygenModels")

using Plots, DifferentialEquations

#=============================================================================#

function get_rc(ox,p)
    β,b,n,rcmin,dn,q = p
    return rc = β*ox.^n./(b^n .+ ox.^n) + rcmin
end

function get_rf(ox,p)
    β,b,n,rcmin,dn,q = p
    return rf = β*(1 - ox.^n./(b^n .+ ox.^n))
end

function get_dc(t,p)
    β,b,n,rcmin,dn,q = p
    return dn
end

function get_df(t,w,p)
    β,b,n,rcmin,dn,q = p

    return dn + q*w
end

#==============================================================================#
# Start ABM debugging here ====  ==============================================#
#  ============================================================================#

# Agent based parameters ====================================================#
tmax = 36*50

β = 16.64/36
b = 13.4256
n = 2.6626
rcmin = 0.0

dn = 0.6045/36
q = 3.27e-5/36

p_fixed = [β,b,n,rcmin,dn,q]

# amount of oxygen
w = 14.63

# size of domain
k = 10000
n = convert(Int64,floor(sqrt(k)))
D = zeros(n,n,2)
D[:,:,2] .= w

### Oxygen ode stuff ========================================================#
fw(w,p,t) = p[1] - p[2]*w - p[3]*p[4]*w

# ode parameters
λ = 9.690912e7/36; μ = 6624000/36; η = (3.16e2/36); C_tot = 0.;
p = [λ,μ,η,C_tot]

# results arrays (temp)
C = []
F = []
P = [] # total population
ox = []

### initial amounts of c and f
c1 = ceil(0.0590*length(D[:,:,1]))
f1 = ceil(0.008*length(D[:,:,1]))
### hard coded initial conditions from ode
# c1 = 580.0
# f1 = 91.0

#populate for f1
f0 = 0
while (f0 < f1)
    xdim = rand(1:n)
    ydim = rand(1:n)
    if D[xdim,ydim,1] == 0
        D[xdim,ydim,1] = 2
        global f0 += 1
    end
end

#populate for c1
c0 = 0
while (c0 < c1)
    xdim = rand(1:n)
    ydim = rand(1:n)
    if D[xdim,ydim,1] == 0
        D[xdim,ydim,1] = 1
        global c0 += 1
    end
end

# first time results
pop = c1 + f1
append!(C,c1)
append!(F,f1)
append!(P,pop)
append!(ox,w)

t = 0
# time loop
@time begin
while (true)
    global w; global pop;
    samp = 0;

    global t += 1
    if t > tmax break end
    println("$t out of ",tmax)

    ### set treatment death rate here
    dc = get_dc(t,p_fixed)
    df = get_df(t,w,p_fixed)

    ### this will be stuff for the oxygen ode
    # count c's for ode
    C_tot = 0
    for j = 1:length(D[:,:,1])
        @isdefined(j)
        # count c's
        if D[j] == 1 C_tot += 1 end
    end

    # solve ode here and update w
    p = [λ,μ,η,float(C_tot)]
    tspan = (0.0,1.0)
    u0 = w
    prob = ODEProblem(fw,u0,tspan,p)
    sol = solve(prob)
    w = sol.u[end]

    # agent based loop
    while samp < pop
        # take samples
        idim = rand(1:n)
        jdim = rand(1:n)

        # if sampled spot is c
        if D[idim,jdim,1] == 1
            samp += 1 #increment sample
            rn = rand() # get random number
            rc = get_rc(w,p_fixed) # get oxygen dependent growth rate

            # division event for c
            if rn < rc
                inew = rand(1:n)
                jnew = rand(1:n)

                # check if new spot is empty
                if D[inew,jnew,1] == 0
                    D[inew,jnew,1] = 1
                end

            # death event for c
            elseif rn < (rc + dc)
                D[idim,jdim,1] = 0
            end
        end # end if sampled spot is c

        # if sampled spot is f
        if D[idim,jdim,1] == 2
            samp += 1 # increment samp
            rn = rand() # get random number
            rf = get_rf(w,p_fixed) # get oxygen dependent growth rate

            # division event for f
            if rn < rf
                # find new spot to populate
                inew = rand(1:n)
                jnew = rand(1:n)

                # check if new spot is empty
                if D[inew,jnew,1] == 0
                    D[inew,jnew,1] = 2
                end

                # death event for f
            elseif rn < (rf + df)
                D[idim,jdim,1] = 0
            end
        end # end if sampled spot is f
    end #end agent based loop

    # count c's and f's
    c = 0; f = 0;
    for j = 1:length(D[:,:,1])
        @isdefined(j)
        # count c's
        if D[j] == 1 c += 1 end
        # count f's
        if D[j] == 2 f += 1 end
    end

    # ### this does the movie ===================================================#
    # BS = "Broad-spectrum antibiotics"
    # Clin = "Antibiotics targeting attack"
    # title = "No treatment"
    # t_b = 19*36; t_c = 33*36;
    
    # # set title
    # if t < t_b
    #     title = BS
    # end
    # if t > t_c
    #     title = Clin
    # end
    # if t%1 == 0
    #     p = heatmap(D[:,:,1],title = title*" $t",legend=true,clims=(0,2))
    #     display(p)
    # end
    # ### end of movie stuff ====================================================#


    # update populations
    pop = c + f
    append!(C,c)
    append!(F,f)
    append!(P,pop)
    append!(ox,w)

    # end if one goes extint
    if c == 0 || f == 0 break end
    if f > c break end

end # end timer
end # end time loop

# # End ABM debugging ===========================================================#
# # =============================================================================#

# absolute
p1 = plot((C)[C.>0],label = "C ABM", lw = 2)
p1 = plot!((F)[F.>0],label = "F ABM", lw = 2)
p2 = plot(ox)
p = plot(p1,p2,layout = (2,1),legend=false, xlabel = "t (days)")
p = plot!(xticks = (0:180:1440, ["0","5","10","15","20","25","30","35","40"]))
display(p)

# relative
p1 = plot((C./P)[C.>0],label = "C ABM", lw = 2)
p1 = plot!((F./P)[F.>0],label = "F ABM", lw = 2)
# p2 = plot(ox)
# p = plot(p1,p2,layout = (2,1),legend=false, xlabel = "t (hours)")
p = plot!(xticks = (0:180:1800, ["0","5","10","15","20","25","30","35","40","45","50"]))
p = plot!(legend=:right)
p = plot!(xlabel = "Time (days)")
p = plot!(ylabel = "Relative Abundance")
p = plot!(title = "Spatially homogeneous model")
p = plot!(ylims = (0,1))
display(p)

# find crossing point
tol = 10
switchpt = findall(abs.(C.-F).<20)

vline!([switchpt],linecolor=:black, label = "Population switch")
plot!(legend=:bottomright)
