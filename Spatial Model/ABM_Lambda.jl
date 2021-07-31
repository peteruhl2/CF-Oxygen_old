# Updated spatial ABM to get an ansemble and count
# the time it takes for the populations to switch
# July 29, 2021

# get to home directory
cd(@__DIR__)

# # for making a gif for VII meeting
# cd("C:\\Users\\peter\\OneDrive\\Desktop\\Sims")

using Plots, Statistics, DifferentialEquations

# Rate functions ===============================================================
# ===============================================================================

function get_rc(ox,p)
    β,b,n_param,rcmin,dn,q = p
    return rc = β*ox.^n_param./(b^n_param .+ ox.^n_param) + rcmin
end

function get_rf(ox,p)
    β,b,n_param,rcmin,dn,q = p
    return rf = β*(1 - ox.^n_param./(b^n_param .+ ox.^n_param))
end

function get_dc(t,p)
    β,b,n_param,rcmin,dn,q = p
    return dn
end

function get_df(t,w,p)
    β,b,n_param,rcmin,dn,q = p

    return dn + q*w
end

# ABM functions ================================================================
# ==============================================================================
function birth_event!(D,i,j,n,bug)
    # randomly sample neighbors and populate if empty
    rn = rand()
    # println(rn)
    # rn = 0.1

    if rn < 1/8 #1
        if (i-1) >= 1 && D[i-1,j,1] == 0
            D[i-1,j,1] = bug
        end
    elseif rn < 2/8 #2
        if (i-1) >= 1 && (j+1) <= n && D[i-1,j+1,1] == 0
            D[i-1,j+1,1] = bug
        end
    elseif rn < 3/8 #3
        if (j+1) <= n && D[i,j+1,1] == 0
            D[i,j+1,1] = bug
        end
    elseif rn < 4/8 #4
        if (i+1) <= n && (j+1) <= n && D[i+1,j+1,1] == 0
            D[i+1,j+1,1] = bug
        end
    elseif rn < 5/8 #5
        if (i+1) <= n && D[i+1,j,1] == 0
            D[i+1,j,1] = bug
        end
    elseif rn < 6/8 #6
        if (i+1) <= n && (j-1) >= 1 && D[i+1,j-1,1] == 0
            D[i+1,j-1,1] = bug
        end
    elseif rn < 7/8 # 7
        if (j-1) >=1 && D[i,j-1,1] == 0
            D[i,j-1,1] = bug
        end

    else #8
        if (i-1) >= 1 && (j-1) >= 1 && D[i-1,j-1,1] == 0
            D[i-1,j-1,1] = bug
        end
    end
end

function get_neigh(D,i,j,n)
    # add up number of neighbors and amount of oxygen they have
    X = 0
    neigh = 0
    for ii = i-1:i+1
        for jj = j-1:j+1
            if (i == ii && j == jj) continue end
            if (ii>0) && (jj>0) && (ii<=n) && (jj<=n)
                X += D[ii,jj,2]
                neigh += 1
            end
        end
    end
    return [X,neigh]
end

#=============================================================================#


# Main funciton ===============================================================#
function ABM_switch(frac, id)
    # size of domain ==========================================================#
    k = 40.0^2
    n = convert(Int64,floor(sqrt(k)))
    D = zeros(n,n,2)

    # parameters
    β = 16.64/36
    b = 13.64
    n_param = 2.66
    rcmin = 0.0
    dn = 0.6045/36
    q = 3.27e-5/36
    λ = 9.7e7/36*frac
    μ = 6.62e6/36
    g = 0.0 # diffustion rate
    Cyn = 0. # C in the spot?
    neigh = 0. # number of neighbors
    X = 0. # total amount of oxygen in neighbors

    # best fitting value of k*η
    keta = 3.16e6
    # η = keta/k
    η = 3.16e6/36

    p_ = [β,b,n_param,rcmin,dn,q,λ,μ,η,g,Cyn,neigh,X]


    # initial oxygen
    w = 14.6
    D[:,:,2] .= w

    # time stuff
    t = 0
    tmax = 50*36*Inf

    ### Oxygen ode ================================================================#
    # fx(x,X,p,t) = λ - μ*x - η*Cyn*x - g*neigh*x + g*X
    # fx(x,p,t) = p[7] - p[8]*x - p[9]*p[11]*x - p[10]*p[12]*x + p[10]*p[13]

    p_ox = [λ,μ,η,Cyn,g,neigh,X]
    fx(x,p,t) = p[1] - p[2]*x - p[3]*p[4]*x - p[5]*p[6]*x + p[5]*p[7]

    ### ===========================================================================#

    # results arrays (temp)
    C = []
    F = []
    P = [] # total population
    ox = []

    # initial amounts of c and f
    c1 = ceil(2*0.0590*length(D[:,:,1]))
    f1 = ceil(0.008*length(D[:,:,1]))

    #populate for f1
    f0 = 0
    while (f0 < f1)
        xdim = rand(1:n)
        ydim = rand(1:n)
        if D[xdim,ydim,1] == 0
            D[xdim,ydim,1] = 2
            f0 += 1
        end
    end

    #populate for c1
    c0 = 0
    while (c0 < c1)
        xdim = rand(1:n)
        ydim = rand(1:n)
        if D[xdim,ydim,1] == 0
            D[xdim,ydim,1] = 1
            c0 += 1
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
    # @time begin
    while (true)
        samp = 0; # global pop; global t; global tmax; global p_ox
        t += 1
        if t > tmax break end

        # check the time now and then
        if t%100 == 0
            println("Thread$id, $t out of ",tmax)
        end

        for i = 1:n # loop over cell array
            for j = 1:n
                # global Cyn = D[i,j,1]
                if D[i,j,1] == 1
                    Cyn = 1
                else Cyn = 0
                end

                p_ox
                stuff = get_neigh(D,i,j,n)
                X = stuff[1]
                neigh = stuff[2]
                # X,neigh = get_neigh(D,i,j,n)

                p_ox[4] = Cyn
                p_ox[6] = neigh
                p_ox[7] = X

                tspan = (0.0,1.0)
                u0 = D[i,j,2]
                prob = ODEProblem(fx,u0,tspan,p_ox)
                sol = solve(prob)
                D[i,j,2] = sol.u[end]

                # display(p_ox)

                # xij = D[i,j,2]
                # k1 = fx(xij,X)
                # k2 = fx(xij + 0.5*step*k1,X)
                # k3 = fx(xij + 0.5*step*k2,X)
                # k4 = fx(xij + step*k3,X)
                # D[i,j,2] = xij + (step/6)*(k1 + 2*k2 + 2*k3 + k4)

                # println("ode is at: ", t2*step)
            end
        end # loop over cell array



        # agent based loop
        while samp < pop
            # take samples
            idim = rand(1:n)
            jdim = rand(1:n)

            # get death rates
            # dc = get_dc(t)
            # df = get_df(t,w)

            # if sampled spot is c
            if D[idim,jdim,1] == 1
                samp += 1 #increment sample
                rn = rand() # get random number
                rc = get_rc(D[idim,jdim,2],p_)
                dc = get_dc(t,p_)
                # println("rc = ",rc)

                # division event for c
                if rn < rc
                    birth_event!(D, idim, jdim, n, 1)

                # death event for c
                elseif rn < (rc + dc)
                    D[idim,jdim,1] = 0
                end
            end #end if sampled spot is c

            # if sampled spot is f
            if D[idim,jdim,1] == 2
                samp += 1 # increment samp
                rn = rand(); # get random number
                rf = get_rf(D[idim,jdim,2],p_) # get oxygen dependent growth rate
                df = get_df(t,D[idim,jdim,2],p_) # need to get the death rate too

                # division event for f
                if rn < rf
                    birth_event!(D, idim, jdim, n, 2)

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

        # this does the movie
        if t%1 == 0
            p1 = heatmap(D[:,:,1],title = "Cells",legend=true,clims=(0,2))
            p2 = heatmap(D[:,:,2],title = "Oxygen",legend=true, clims=(0.05, maximum(D[:,:,2])))
            # p2 = heatmap(D[:,:,2],title = "Oxygen",legend=true, clims=(5.0,15))
            p = plot(p1,p2,layout = (1,2),legend=true)
            display(p)
        end

        # # just to save figures
        # # for a gif 1/13/21 and 6/18/21
        # if (t-1) % 1 == 0
        #     if t < 10
        #         println("00000$t")
        #         savefig("patch00000$t.png")
        #     elseif t < 100
        #         println("0000$t")
        #         savefig("patch0000$t.png")
        #     elseif t < 1000
        #         println("000$t")
        #         savefig("patch000$t.png")
        #     elseif t < 10000
        #         println("00$t")
        #         savefig("patch00$t.png")
        #     end
        #     # println("$t")
        #     # savefig("patch$t.png")
        # end

        # update populations
        # pop = c + f
        # append!(C,c)
        # append!(F,f)
        # append!(P,pop)
        # append!(ox,mean(D[:,:,2]))

        # end if one goes extint
        if c == 0 || f == 0 break end
        if f > c
            return t
        end

    # end # end timer
    end # end time loop

end # end ABM_switch

# # absolute
# p1 = plot((C)[C.>0],label = "C ABM", lw = 2)
# p1 = plot!((F)[F.>0],label = "F ABM", lw = 2)
# p2 = plot(ox)
# p = plot(p1,p2,layout = (2,1),legend=false, xlabel = "t (days)")
# p = plot!(xticks = (0:180:1440, ["0","5","10","15","20","25","30","35","40"]))
# display(p)
#
# # relative
# p1 = plot((C./P)[C.>0],label = "C ABM", lw = 2)
# p1 = plot!((F./P)[F.>0],label = "F ABM", lw = 2)
# # p2 = plot(ox)
# # p = plot(p1,p2,layout = (2,1),legend=false, xlabel = "t (hours)")
# p = plot!(xticks = (0:180:1800, ["0","5","10","15","20","25","30","35","40","45","50"]))
# p = plot!(legend=:right)
# p = plot!(xlabel = "Time (days)")
# p = plot!(ylabel = "Relative Abundance")
# p = plot!(title = "Spatially dependent model")
# p = plot!(ylims = (0,1))
#
# display(p)
#
# # find crossing point
# tol = 10
# switchpt = findall(abs.(C.-F).<10)
# switchpt = switchpt[1]
#
# vline!([switchpt[1]],linecolor=:black, label = "Population switch")

# go to simulations folder
cd("C:\\Users\\peter\\OneDrive\\Documents\\GitHub\\CF-Oxygen\\Spatial Model\\Spatial Simulations")

n_sims = 1
Fracs = 0.5

# keep track of sims run
sims_run = []

io = open("frac = $Fracs.csv", "w")

@time begin
# Threads.@threads for i = 1:n_sims
for i = 1:n_sims

    # add i to list when it runs
    if i ∉ sims_run
        append!(sims_run,i)
    end

    # which simulations its on and how many remaining
    println("Thread number ", Threads.threadid(), " lambda frac = ", Fracs, " Simulation ",i)
    println("Remaining: ", n_sims - length(sims_run))

    # write switch time to file
    t = ABM_switch(Fracs,Threads.threadid())
    write(io, "$t\n")
end
end # end timer

close(io)
