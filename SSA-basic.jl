module SSAnonVD

export SSA

"""
propensity: function to output the propensities given the state of the system.
args:
- n: the current state vector.
- pars: the system parameters.
"""
function propensity(n::Vector, pars::Vector, hypars::Vector)
    # Reaction rates
    kA = pars[1]; kB = pars[2]; kC = pars[3];
    hS = hypars[1]; hR = hypars[2]; 
    f_r = zeros(hS+hR+1);

    for i in 1:hS-1
        f_r[i] = kA*n[i];
    end
    f_r[hS] = kB*n[hS]
    f_r[hS+1] = kC*n[hS+1];
    for i in hS+2:hS+hR+1
        f_r[i] = n[i]*hR;
    end
    return f_r::Vector{Float64}
end

"""
SSA: function to perform the SSA (with the delayed degradation of nRNA) given the necessary parameters.
args:
- S_time: the number of individual simulations in the ensemble (set to 1 for single trajectory).
- pars: the parameters for the ensemble simulations.
- hpars: the hyperparams describing the number of 
- tol_time: the total simulation time to run for.
- sp: the storage time period, i.e., if sp = 1.0 then final state vector stored every 1.0s.
returns:
- the state vector at the specified times.
"""
function SSA(S_time::Int, pars::Vector, hypars::Vector{Int}, tol_time::Float64, sp::Float64)

    sp <= tol_time || error("The storage time period must be less than the total simulation time!")

    hS = hypars[1]; hR = hypars[2];

    # M = Number of reactions, N = Number of reactants
    M = hS+hR+1::Int; N=hS+hR+1::Int;

    # Define stoichiometry matrix
    S_mat = zeros(M,N);
    # construct row-by-row
    S_mat[1,1] = -1; S_mat[1,hS] = +1; # only first row differs
    for i in 2:M
        S_mat[i,i] = -1; S_mat[i,i-1] = +1;
    end
    S_mat = transpose(S_mat)

    times = convert(Array{Float64,1},LinRange(tol_time,0.0,floor(Int,tol_time/sp)+1));

    # Define reactants trjatory vector
    n = zeros(N,S_time,length(times));

    # define the means of nuclear and cyto for the IC
    kA = pars[1]; kB = pars[2]; kC = pars[3];
    muN = kA*kB/(kC*(kA+(hS-1)kB)); muC = kA*kB/(kA+(hS-1)kB);
    muCi = muC/hR;

    for sim in 1:S_time
        n_temp = convert(Vector{Int64}, vcat([1],zeros(hS-1),[round(Int,muN)],round(Int,muC) .*ones(hR))); # gene starts U1 state, rest near SS values
        T = 0;
        sim_times = copy(times);

        # define counter m for updating storage.
        m = 1;
        while T < tol_time
            # Step 1: Calculate propensity
            f_r = propensity(n_temp, pars, hypars); # propensity of each reaction.
            lambda = sum(f_r); # total propensity of any reaction.

            # Step 2: Calculate tau and mu using random number genrators
            r1 = rand(2,1);
            tau = (1/lambda)*log(1/r1[1]);
            next_r = findfirst(x -> x>=r1[2]*lambda,cumsum(f_r));

            while T+tau >= sim_times[end]
                n[1:N,sim,m] = n_temp; # m used here.
                pop!(sim_times);
                m += 1;
                if length(sim_times) == 0
                    break
                end
            end

            # update the system time
            T += tau;
            # update the state vector
            prod = S_mat[next_r,1:N];
            for i in 1:N
                n_temp[i] += prod[i]
            end

        end
        # if mod(sim,1000) == 0
        #     println(sim)
        # end
    end
    # let's only return the numbers of nuclear and cyto
    nuc_data = n[hS+1,:,:]
    cyto_data = n[hS+1:end,:,:]
    return nuc_data, sum(cyto_data,dims=(1))[1,:,:]
end

end
# # data is a 2-D array:
# # first dim is constant time, second dim diff. times. (cols of constant time)
# function mean(data)
#     means = zeros(size(data)[2])
#     for m in 1:length(means)
#         sum_all = sum(data[:,m]);
#         avg = sum_all / length(data[:,m]);
#         means[m] = avg;
#     end
#     return means
# end

# # data is a 2-D array:
# # first dim is constant time, second dim diff. times.
# function var(data, means)
#     var = zeros(length(data[1,:]))
#     for v in 1:length(var)
#         m = means[Int(v)]
#         squares = zeros(length(data[:,v]))
#         for s in 1:length(squares)
#             squares[s] = data[s,v] ^ 2;
#         end
#         sum_all_sq = sum(squares);
#         va = ((sum_all_sq) / length(data[:,v])) - m^2;
#         var[v] = va;
#     end
#     return var
# end

# end