# Non-spatial agent based model with oxygen
# the oxygen ode acts on the whole grid with number of C cells
# 12/14/2020

# get to home directory
cd("C:\\Users\\peter\\Onedrive\\Desktop\\cyst fib\\OxygenModels")

using Plots, DifferentialEquations
using Statistics

#=============================================================================#

function get_rc(ox)
    Ec = 11.7549/24; Ac = 0.0139; nc = 1.4408; rcmin = 0.0;
    return rc = Ec*ox.^nc./(Ac^nc .+ ox.^nc) + rcmin
end

#==============================================================================#
# Solve ODE ===================================================================#
function cf_ode(yp,y,p,t)
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

    # println(get_rc(x)," ", get_rf(x))
end

# ode parameters
Ec = 11.7549; Ac = 0.0139; nc = 1.4408
rf = 14.3436; d = 0.7016; ϵ = 0.5428; μ = 1.4273;
k = 10000
η = 0.8176/k

p = [Ec,Ac,nc,rf,d,ϵ,μ,η,k]

y0 = [0.8283,0.0165,0.1388]
tspan = (0.0,40.0)
prob = ODEProblem(cf_ode,y0,tspan,p)
sol = solve(prob)

#==============================================================================#

# Agent based parameters ======================================================#
tmax = 24*40

# amount of oxygen
w = 0.1388

# size of domain
n = convert(Int64,floor(sqrt(k)))
D = zeros(n,n,2)
D[:,:,2] .= w

### Oxygen ode stuff ==========================================================#
fw(w,p,t) = p[1] - p[2]*w - p[3]*p[4]*w

# ode parameters
λ = 0.22/24; μ = 1.4273/24; η = (0.8176/24)/n^2; C_tot = 0.;
p = [λ,μ,η,C_tot]

# treatment parameter
ϵ = 0.5428/24
t_treat = 24*28 # starts at 28 days


# attack growth rate and common death rate
dc = 0.7016/24
df = 0.7016/24
rf = 14.3436/24


# results arrays (temp)
C = []
F = []
P = [] # total population
ox = []

# initial amounts of c and f
c1 = ceil(0.8283*length(D[:,:,1]))
f1 = ceil(0.0165*length(D[:,:,1]))

# OLD =========================================================================#
# #populate for f1
# for i = 1:f1
#     xdim = rand(1:n)
#     ydim = rand(1:n)
#     if D[xdim,ydim,1] == 0
#         D[xdim,ydim,1] = 2
#     end
# end

# # populate for c1
# for i = 1:c1
#     xdim = rand(1:n)
#     ydim = rand(1:n)
#     if D[xdim,ydim,1] == 0
#         D[xdim,ydim,1] = 1
#     end
# end
# END OLD =====================================================================#

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

# # get actual initial c and f
# c1 = 0; f1 = 0
# for i=1:length(D[:,:,1])
#     #count number of c
#     if D[i] == 1 global c1 += 1 end
#     #count number of f
#     if D[i] == 2 global f1 += 1 end
# end

# first time results
pop = c1 + f1
append!(C,c1)
append!(F,f1)
append!(P,pop)

t = 0
# time loop
@time begin
while (true)
    global w; global pop;
    samp = 0;

    global t += 1
    if t > tmax break end
    println("$t out of ",tmax)

    # set treatment death rate here
    if t < t_treat
        df = 0.7037/24
    else
        df = 0.7037/24 + ϵ
    end


    # this will be stuff for the oxygen ode
    u0 = w

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
            rc = get_rc(w)

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
        end #end if sampled spot is c

        # if sampled spot is f
        if D[idim,jdim,1] == 2
            samp += 1 # increment samp
            rn = rand(); # get random number
            # rf = get_rf(w)

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

    # # this does the movie
    # if t%1 == 0
    #     p = heatmap(D[:,:,1],title = "Cells",legend=true,clims=(0,2))
    #     display(p)
    # end

    # update populations
    pop = c + f
    append!(C,c)
    append!(F,f)
    append!(P,pop)
    append!(ox,w)

    # end if one goes extint
    if c == 0 || f == 0 break end

end # end timer
end # end time loop

p1 = plot((C./n^2)[C.>0],label = "C ABM", lw = 2)
p1 = plot!((F./n^2)[F.>0],label = "F ABM", lw = 2)
p1 = plot!(24*sol.t,sol[1,:], label = "C ODE")
p1 = plot!(24*sol.t,sol[2,:], label = "F ODE", legend=:left)
p2 = plot(ox)
p2 = plot!(24*sol.t,sol[3,:])
p = plot(p1,p2,layout = (2,1),legend=false)
display(p)

println(C[1], " ", F[1])
