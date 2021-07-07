# Updated spatial ABM with toxic oxygen
# June 16, 2021

# get to home directory
cd(@__DIR__)

# # for making a gif for VII meeting
# cd("C:\\Users\\peter\\OneDrive\\Documents\\GitHub\\CF-Oxygen\\Spatial Model\\figs")
cd("C:\\Users\\peter\\OneDrive\\Desktop\\Sims")

using Plots, Statistics

# Rate functions ===============================================================
# ===============================================================================

function get_rc(ox)
    global b
    β = 26.64/36; n = 2.6626; rcmin = 0.0;
    return rc = β*ox.^n./(b^n .+ ox.^n) + rcmin
end

function get_rf(ox)
    global b
    β = 26.64/36; n = 2.6626; rcmin = 0.0;
    return rc = β*(1 - ox.^n./(b^n .+ ox.^n)) + rcmin
end

function get_dc(t)
    t_b = 19*36; dn = 0.6045/36; dbs = 6.7686/36

    # if t < t_b
    #     return dn + dbs
    # else
    #     return dn
    # end

    return dn
end

function get_df(t,w)
    global q
    t_b = 19*36; t_c = 33*36; γ = 0.8976; dn = 0.6045/36
    dbs = 6.7686/36; ϵ = 1.21/36
    # q = 1.2;

    # base rate
    df = dn + q*w

    # # broad spectrum rate
    # if t < t_b
    #     df += γ*dbs
    # end
    #
    # # clindamycin rate
    # if t > t_c
    #     df += ϵ
    # end

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
# oxygen growth parameters
b = 0.03
q = 4.5

# size of domain and initial values
w = 0.0388 # initial oxygen
k = 40^2 # carrying capacity
n = convert(Int64,floor(sqrt(k)))
D = zeros(n,n,2)
D[:,:,2] .= w

# time stuff
t = 0
tmax = 1000
t_treat = 24*28

# ode stuff
λ = 0.42/24; μ = 1.0273/24; Cyn = 0.; neigh = 0; X = 0;
η = 22.0
g = 1.2 # diffustion rate

fx(x,X) = λ - μ*x - η*Cyn*x - g*neigh*x + g*X
step = 0.05 # as large as possible w/o blowing up the ode

# death rates constant
dc = 0.7016/24
df = 0.7016/24

# treatment parameter
ϵ = 0.5428/24

# results arrays (temp)
C = []
F = []
P = [] # total population
ox = []

# initial amounts of c and f
c1 = ceil(0.0683*length(D[:,:,1]))
f1 = ceil(0.05*length(D[:,:,1]))

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

# time loop
@time begin
while (true)
    samp = 0; global pop; global t; global tmax
    t += 1
    if t > tmax break end
    println("$t out of ",tmax)

    global step; t2 = 0
    while (t2*step <= 1) # will be ode loop
        t2 += 1

        for i = 1:n # loop over cell array
            for j = 1:n
                # global Cyn = D[i,j,1]
                if D[i,j,1] == 1
                    global Cyn = 1
                else global Cyn = 0
                end

                stuff = get_neigh(D,i,j,n)
                global X = stuff[1]
                global neigh = stuff[2]
                # X,neigh = get_neigh(D,i,j,n)

                xij = D[i,j,2]
                k1 = fx(xij,X)
                k2 = fx(xij + 0.5*step*k1,X)
                k3 = fx(xij + 0.5*step*k2,X)
                k4 = fx(xij + step*k3,X)
                D[i,j,2] = xij + (step/6)*(k1 + 2*k2 + 2*k3 + k4)

                # println("ode is at: ", t2*step)
            end
        end # loop over cell array

    end # end ode loop


    # agent based loop
    while samp < pop
        # take samples
        idim = rand(1:n)
        jdim = rand(1:n)

        # get death rates
        dc = get_dc(t)
        # df = get_df(t,w)

        # if sampled spot is c
        if D[idim,jdim,1] == 1
            samp += 1 #increment sample
            rn = rand() # get random number
            rc = get_rc(D[idim,jdim,2])

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
            rf = get_rf(D[idim,jdim,2]) # get oxygen dependent growth rate
            df = get_df(t,D[idim,jdim,2]) # need to get the death rate too

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
        # p2 = heatmap(D[:,:,2],title = "Oxygen",legend=true) #,clims=(0.05,maximum(D[:,:,2])))
        p2 = heatmap(D[:,:,2],title = "Oxygen",legend=true, clims=(0.0,0.1))
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
    pop = c + f
    append!(C,c)
    append!(F,f)
    append!(P,pop)
    append!(ox,mean(D[:,:,2]))

    # end if one goes extint
    # if c == 0 || f == 0 break end

end # end timer
end # end time loop


p1 = plot((C./n^2)[C.>0],label = "C ABM", lw = 2)
p1 = plot!((F./n^2)[F.>0],label = "F ABM", lw = 2)
p2 = plot(ox)
p = plot(p1,p2,layout = (2,1),legend=false, xlabel = "t (hours)")
display(p)

# println(C[end]," ",F[end])
