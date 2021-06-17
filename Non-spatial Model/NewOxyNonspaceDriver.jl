### Driver for spatially homogeneous agent-based model with most up to date
### parameters as of June 2021
### 6/6/2021

# get to directory for simulations
cd("C:\\Users\\peter\\OneDrive\\Documents\\GitHub\\CF-Oxygen\\Non-spatial Model")

using Plots
using Statistics

### Functions =================================================================#
###============================================================================#

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

### Simulation function =======================================================#
function sim(num::Int64, movie::Bool=false)

    # move to directery and make file
    cd("C:\\Users\\peter\\OneDrive\\Documents\\GitHub\\CF-Oxygen\\Non-spatial Model\\Simulations")
    io = open("sim_$num.csv","w")

    # Agent based parameters ====================================================#
    tmax = 36*40

    # amount of oxygen
    # w = 0.1388
    w = 14.63

    # size of domain
    n = convert(Int64,floor(sqrt(k)))
    D = zeros(n,n,2)
    D[:,:,2] .= w

    ### Oxygen ode stuff ========================================================#
    fw(w,p,t) = p[1] - p[2]*w - p[3]*p[4]*w

    # ode parameters
    # λ = 0.22/24; μ = 1.4273/24; η = (0.8176/24)/n^2; C_tot = 0.;
    λ = 9.690912e7/36; μ = 6624000/36; η = (3.16e2/36); C_tot = 0.;
    p = [λ,μ,η,C_tot]

    # treatment parameter(s)
    ϵ = 1.21/36
    t_b = 19*36
    t_c = 33*36

    # attack growth rate and death rates
    # adjusted for new death terms 6/4/21
    dn = 0.6045/36
    dbs = 6.7686/36

    # results arrays (temp)
    C = []
    F = []
    P = [] # total population
    ox = []

    ### initial amounts of c and f
    # c1 = ceil(0.8283*length(D[:,:,1]))
    # f1 = ceil(0.165*length(D[:,:,1]))
    ### hard coded initial conditions from ode
    c1 = 580.0
    f1 = 91.0

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
    @time begin
    while (true)
        # global w; global pop;
        samp = 0;

        t += 1
        if t > tmax break end
        println("$t out of ",tmax)

        ### set treatment death rate here
        dc = get_dc(t)
        df = get_df(t,w)

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
                rc = get_rc(w) # get oxygen dependent growth rate

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
                rf = get_rf(w) # get oxygen dependent growth rate

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
        #
        # # set title
        # if t < t_b
        #     title = BS
        # end
        # if t > t_c
        #     title = Clin
        # end
        # if t%1 == 0
        #     p = heatmap(D[:,:,1],title = title,legend=true,clims=(0,2))
        #     display(p)
        # end
        # ### end of movie stuff ====================================================#


        # update populations
        pop = c + f
        append!(C,c)
        append!(F,f)
        append!(P,pop)
        append!(ox,w)

        write(io, "$c $f $pop $w\n")

        # # end if one goes extint, turned off to always get the full time range
        # if c == 0 || f == 0 break end

    end # end timer
    end # end time loop

    # close file
    close(io)
end

function main()
    mkdir("Simulations")

    runs = 100

    for i = 1:runs
        sim(i)
    end

end

main()
