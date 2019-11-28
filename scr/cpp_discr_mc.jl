cd("/Users/Frank/Sync/DOCUMENTS/onderzoek/CPP-discrete/")

using RCall
using Distributions
using DataFrames
using DelimitedFiles
using SpecialFunctions
using TimerOutputs
using Random
using RCall

to = TimerOutput()

@rlibrary nilde
@rlibrary ggplot2
include("fundefs.jl")

function gendata(n,setting)
    Δ = ones(n)

    if setting==1
        λ0 = 2.0
        p0 = [1,0,0,1,0,1]/3
        m0 = length(p0)
        obs =zeros(Int64,n)
        for i in 1:n
            n_jumps = rand(Poisson(λ0*Δ[i]))
            obs[i] = sum(wsample(1:m0,p0,n_jumps))
        end
    else
        sp = 1/3 # succes probability
        λ0 = -log(sp)
        m0 = 30
        p0 = -(1-sp).^(1:m0)./((1:m0)*log(sp))
        obs = rand(Geometric(sp),n)
    end

    Δ, vec(map(x->Int.(x),obs)) # ensure observations are integer valued
end

function postmean_and_bgtrunc(obs,maxsupp,maxeval,Δ,IT,BI,a,c)
    n = length(obs)
    # Compute the vector that maps each element in obs to the vector z,
    # which consists of the unique values of z, sorted ascending
    z = sort(unique(obs))
    m = min(maxsupp, maximum(z)) #note that maximum(z)=maximum(obs)

    # compute solutions to Diophantine equation (using R)
    @timeit to "compute sols Diophantine equation"  sols, num_sols = diophantine(z,m)

    # For each element in obs, we need to find its location in z, so
    # this is the vector "ind" computed here: so we have z[ind[k]] = obs[k] for k in 1:n
    ind = zeros(Int64,n)
    for i=1:n  ind[i] = searchsorted(z,obs[i])[1]  end

    if !(sum(abs.([z[ind[i]]-obs[i]<0.0001 for i in 1:n]))==250)  error("something wrong") end

    μ = zeros(Int64,IT,m)
    ν = zeros(IT,m)
    b = ones(IT,m)
    γ = ones(IT)

    sum_acc = 0
    proposaltype = 0.2 #0.1  # prob ot picking from the uniform distribution on all Diophantine solutions

    # initialisation at random
    imp = zeros(Int64,IT,n) # present set of pointers to imputed data  element (it,i) contains the pointers at iteration "it" for increment "i"
    for i in 1:n    imp[1,i] = sample(1:num_sols[ind[i]])  end
    μ[1,:] = imp2mu(imp[1,:],sols,ind,m,n)

    impute_ind = findall(obs.>1.1)
    nonimpute_ind = findall(obs.<1.1) # if the jump is 0 or 1, then there is just one way to impute
    imp[:,nonimpute_ind] .= 1

    sum_Δ = sum(Δ)

    for it in 2:IT
        global it
        global acc
        for k in 1:m
            ν[it,k] = rand(Gamma(a + μ[it-1,k], 1/(1/b[it-1,k] + sum_Δ)))
            b[it,k] = rand(InverseGamma(a+c,γ[it-1]+ν[it,k]))
        end
        γ[it] = rand(Gamma(c*m+1,1 ./ (1+sum(1 ./ b[it,:]))))
        for i in impute_ind
          @timeit to "impute step" imp[it,i], acc  = update_imp2(imp[it-1,i],ν[it,:],Δ[i],sols[ind[i]],num_sols[ind[i]],m,proposaltype)
          sum_acc += acc
        end
        μ[it,:] = imp2mu(imp[it,:],sols,ind,m,n)

        if mod(it,25000)==0    println(it)    end

    end
    ν_postmean = vec(mean(ν[BI:IT,:],dims=1))
    ν_postmedian = vec(median(ν[BI:IT,:],dims=1))

    @timeit to "comp bg trunc" λ_bg, p_bg = bg_truncated_estimator(obs,maximum(z))

    # return estimates of length maxsupp by augmenting with zeros
    (aug(ν_postmean,maxeval), aug(ν_postmedian,maxeval), aug(λ_bg * p_bg,maxeval))

end

function aug(x, maxval) # augment vector x with zeros
    if maxval < length(x)
        return(x)
    else
        return(vcat(x, fill(0.0, maxval-length(x))) )# augment x with zeros such
    end
end

function l1error(resi,ν_true)
    [sum(abs.(resi[1]-ν_true)),sum(abs.(resi[2]-ν_true)),sum(abs.(resi[3]-ν_true))]
end

################################## prior specification for ν--------------------------
a = 0.01
c = 2.0 # settings 1 and 2
#c = 0.01 # for setting 3


################################## Set simulation pars --------------------------
# setting=1 corresponds to BG example with uniform distribution on {1,4,6}
# setting=2 corresponds to Geometric distribution
n = 250 # sample size
setting = 1

maxeval = 50
# in the following, maxsupp is the number that ensures that probabilities will be computed
if setting==1
    maxsupp = 15 # ensure >=6
    ν_true = aug([2/3, 0, 0, 2/3, 0, 2/3], maxeval)#vcat([2/3, 0, 0, 2/3, 0, 2/3], fill(0.0,maxsupp-6))
else
    maxsupp = 20
    sp = 1/3 # succes probability
    λ0 = -log(sp)
    p0 = -(1-sp).^(1:maxeval)./((1:maxeval)*log(sp))
    ν_true = λ0 * p0
end

Random.seed!(123) #  set RNG
IT = Int64(4*10^5)
BI = Int64(2*10^5)

B = 50 # MC sample size
res = Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}}[]
for b in 1:B
        println(b)
        Δ, obs = gendata(n,setting)
        push!(res,postmean_and_bgtrunc(obs, maxsupp, maxeval, Δ, IT, BI, a, c))
end

println("Timing:")
show(to, allocations = false, compact = true)

############################### Compute errors -------------------------------------
err = map(x-> l1error(x,ν_true), res)
err_bayesmean = map(x->x[1],err)
err_bayesmedian = map(x->x[2],err)
err_bg = map(x->x[3],err)

f = open("./mcstudy/l1errors.csv","w")
write(f, "bayesmean, bayesmedian,bg\n")
writedlm(f, hcat(err_bayesmean, err_bayesmedian, err_bg), ',')
close(f)

# visualisation
errdf = DataFrame(bayesmean=err_bayesmean, bayesmedian=err_bayesmedian,bg=err_bg)
@rput errdf
R"""
library(tidyverse)
errdf %>% gather(value=error, key=type, bayesmean, bayesmedian, bg) %>%
ggplot(aes(x=type, y=error)) +
geom_boxplot(fill = "lightsteelblue1") + xlab("") +
 ylab("absolute (L1) error")+ coord_flip() + theme_minimal()
"""
