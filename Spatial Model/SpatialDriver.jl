# Driver program for the spatially heterogenous ABM
# 12/22/2020

# get to home directory
cd(@__DIR__)

using Plots, Statistics

#=============================================================================#
function birth_event!(D, i, j, n, bug)
    # randomly sample neighbors and populate if empty
    rn = rand()
    # println(rn)
    # rn = 0.1

    if rn < 1 / 8 #1
        if (i - 1) >= 1 && D[i-1, j, 1] == 0
            D[i-1, j, 1] = bug
        end
    elseif rn < 2 / 8 #2
        if (i - 1) >= 1 && (j + 1) <= n && D[i-1, j+1, 1] == 0
            D[i-1, j+1, 1] = bug
        end
    elseif rn < 3 / 8 #3
        if (j + 1) <= n && D[i, j+1, 1] == 0
            D[i, j+1, 1] = bug
        end
    elseif rn < 4 / 8 #4
        if (i + 1) <= n && (j + 1) <= n && D[i+1, j+1, 1] == 0
            D[i+1, j+1, 1] = bug
        end
    elseif rn < 5 / 8 #5
        if (i + 1) <= n && D[i+1, j, 1] == 0
            D[i+1, j, 1] = bug
        end
    elseif rn < 6 / 8 #6
        if (i + 1) <= n && (j - 1) >= 1 && D[i+1, j-1, 1] == 0
            D[i+1, j-1, 1] = bug
        end
    elseif rn < 7 / 8 # 7
        if (j - 1) >= 1 && D[i, j-1, 1] == 0
            D[i, j-1, 1] = bug
        end

    else #8
        if (i - 1) >= 1 && (j - 1) >= 1 && D[i-1, j-1, 1] == 0
            D[i-1, j-1, 1] = bug
        end
    end
end

function get_rc(ox)
    Ec = 11.7549 / 24
    Ac = 0.0139
    nc = 1.4408
    rcmin = 0.0
    return rc = Ec * ox .^ nc ./ (Ac^nc .+ ox .^ nc) + rcmin
end

function get_neigh(D, i, j, n)
    # add up number of neighbors and amount of oxygen they have
    X = 0
    neigh = 0
    for ii = i-1:i+1
        for jj = j-1:j+1
            if (i == ii && j == jj)
                continue
            end
            if (ii > 0) && (jj > 0) && (ii <= n) && (jj <= n)
                X += D[ii, jj, 2]
                neigh += 1
            end
        end
    end
    return [X, neigh]
end

function sim(num::Int64, g::Float64, k::Int64, movie::Bool = false)
    # make file for results
    io = open("sim_$num.csv", "w")

    # size of domain and initial values
    w = 0.1388 # initial oxygen
    # k = 25^2 # carrying capacity, getting from input now
    n = convert(Int64, floor(sqrt(k)))
    D = zeros(n, n, 2)
    D[:, :, 2] .= w

    # time stuff
    t = 0
    tmax = Inf #24*40
    t_treat = 24 * 28

    # ode stuff
    λ = 0.22 / 24
    μ = 1.4273 / 24
    Cyn = 0.0
    neigh = 0
    X = 0
    # η = (0.8176/24)/n^2 # this might not be per capita for this model
    η = (0.8176 / 24)
    # g = 0.05 # diffustion rate, getting from input now

    fx(x, X) = λ - μ * x - η * Cyn * x - g * neigh * x + g * X
    step = 0.005 # as large as possible w/o blowing up the ode

    # attack growth rate and common death rate
    dc = 0.7016 / 24
    df = 0.7016 / 24
    rf = 14.3436 / 24

    # treatment parameter
    ϵ = 0.5428 / 24

    # results arrays (temp)
    C = []
    F = []
    P = [] # total population
    ox = []

    # initial amounts of c and f
    c1 = ceil(0.8283 * length(D[:, :, 1]))
    f1 = ceil(0.0165 * length(D[:, :, 1]))

    #populate for f1
    f0 = 0
    while (f0 < f1)
        xdim = rand(1:n)
        ydim = rand(1:n)
        if D[xdim, ydim, 1] == 0
            D[xdim, ydim, 1] = 2
            f0 += 1
        end
    end

    #populate for c1
    c0 = 0
    while (c0 < c1)
        xdim = rand(1:n)
        ydim = rand(1:n)
        if D[xdim, ydim, 1] == 0
            D[xdim, ydim, 1] = 1
            c0 += 1
        end
    end

    # first time results
    pop = c1 + f1
    append!(C, c1)
    append!(F, f1)
    append!(P, pop)
    append!(ox, w)

    # time loop
    @time begin
        while (true)
            samp = 0 #global pop; global t; global tmax
            t += 1
            if t > tmax
                break
            end
            println("$t out of ", tmax)

            t2 = 0
            while (t2 * step <= 1) # will be ode loop
                t2 += 1

                for i = 1:n # loop over cell array
                    for j = 1:n
                        # global Cyn = D[i,j,1]
                        if D[i, j, 1] == 1
                            Cyn = 1
                        else
                            Cyn = 0
                        end

                        stuff = get_neigh(D, i, j, n)
                        X = stuff[1]
                        neigh = stuff[2]
                        # X,neigh = get_neigh(D,i,j,n)

                        xij = D[i, j, 2]
                        k1 = fx(xij, X)
                        k2 = fx(xij + 0.5 * step * k1, X)
                        k3 = fx(xij + 0.5 * step * k2, X)
                        k4 = fx(xij + step * k3, X)
                        D[i, j, 2] =
                            xij + (step / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

                        # println("ode is at: ", t2*step)
                    end
                end # loop over cell array

            end # end ode loop


            # agent based loop
            while samp < pop
                # take samples
                idim = rand(1:n)
                jdim = rand(1:n)

                # if sampled spot is c
                if D[idim, jdim, 1] == 1
                    samp += 1 #increment sample
                    rn = rand() # get random number
                    rc = get_rc(D[idim, jdim, 2])

                    # division event for c
                    if rn < rc
                        birth_event!(D, idim, jdim, n, 1)

                        # death event for c
                    elseif rn < (rc + dc)
                        D[idim, jdim, 1] = 0
                    end
                end #end if sampled spot is c

                # if sampled spot is f
                if D[idim, jdim, 1] == 2
                    samp += 1 # increment samp
                    rn = rand() # get random number

                    # division event for f
                    if rn < rf
                        birth_event!(D, idim, jdim, n, 2)

                        # death event for f
                    elseif rn < (rf + df)
                        D[idim, jdim, 1] = 0
                    end
                end # end if sampled spot is f
            end #end agent based loop

            # count c's and f's
            c = 0
            f = 0
            for j = 1:length(D[:, :, 1])
                @isdefined(j)
                # count c's
                if D[j] == 1
                    c += 1
                end
                # count f's
                if D[j] == 2
                    f += 1
                end
            end

            # this does the movie
            if movie == true
                if t % 1 == 0
                    p1 = heatmap(
                        D[:, :, 1],
                        title = "Cells",
                        legend = true,
                        clims = (0, 2),
                    )
                    p2 = heatmap(D[:, :, 2], title = "Oxygen", legend = true) #,clims=(0.05,maximum(D[:,:,2])))
                    p = plot(p1, p2, layout = (1, 2), legend = true)
                    display(p)
                end
            end

            # # just to save figures
            # if (t-1) % 10 == 0
            #     println("$t")
            #     savefig("patch$t.png")
            # end

            # update populations
            pop = c + f
            append!(C, c)
            append!(F, f)
            append!(P, pop)
            append!(ox, mean(D[:, :, 2]))

            new_w = mean(D[:, :, 2])
            # write number to file
            write(io, "$c,$f,$pop,$new_w\n")

            # end if one goes extint
            if c == 0 || f == 0
                break
            end

        end # end timer
    end # end time loop

    # close file
    close(io)
end # sim function
#==============================================================================#

function main()
    # model parameters
    g = 5.95
    k = 40^2
    runs = 1

    # make director for simulations
    dir = "Sim g = $g, k = $k"

    # go to simulation folder
    cd("C:\\Users\\peter\\OneDrive\\Documents\\GitHub\\CF-Oxygen\\Spatial Model\\Simulations")
    mkdir(dir)
    cd(dir)

    for i = 1:runs
        sim(i, g, k, false)
    end

    # go back up to the simulations folder
    cd("C:\\Users\\peter\\OneDrive\\Documents\\GitHub\\CF-Oxygen\\Spatial Model\\Simulations")
end

main()
