using RCall
using Distributions
using DataFrames
using DelimitedFiles
using SpecialFunctions
using Random
@rlibrary nilde
@rlibrary ggplot2

Random.seed!(34) #  set RNG
include("fundefs.jl")

#-------------------------- load / generate data
scenarios = ["generated1","generated2",
            "testdata_A1","testdata_A2","testdata_A3",
            "testdata_B1","testdata_B2","testdata_C",
            "horsekicks","plantpopulation"]
data_choice = scenarios[9]
generated = scenarios[1:8]  # set of scenarios that are generated and the true Levy measure is known
println(data_choice)
include("setdata.jl")

#-------------------------- prior specification for ν
a = 0.01;  c = 2.0

#-------------------------- preprocessing of data, computing solutions of Diophantine equation
z = sort(unique(obs))
maxsupp = 25
m = min(maxsupp,maximum(z))
# For each element in obs, we need to find its location in z
# ind is such that z[ind[k]] = obs[k] for k in 1:m
ind = zeros(Int64,n)
for i=1:n  ind[i] = searchsorted(z,obs[i])[1]  end
sols, num_sols = diophantine(z,m)

#-------------------------- set nr of mcmc iterations
IT = Int64(5*10^5)  # nr of iterations
BI = Int(IT/2)      # nr of burnin iterations
proposaltype = 0.2 #0.1  # prob ot picking from the uniform distribution on all Diophantine solutions

#-------------------------- initialisation of mcmc algorithm
μ = zeros(Int64,IT,m)
ν = zeros(IT,m)
b = ones(IT,m)
γ = ones(IT)

sum_acc = 0 # count accepted mh steps

# initialisation of imp is at random
imp = zeros(Int64,IT,n) # present set of pointers to imputed data  element (it,i) contains the pointers at iteration "it" for increment "i"
for i in 1:n    imp[1,i] = sample(1:num_sols[ind[i]])  end
μ[1,:] = imp2mu(imp[1,:],sols,ind,m,n)
impute_ind = findall(obs.>1.1)
nonimpute_ind = findall(obs.<1.1) # if the jump is 0 or 1, then there is just one way to impute
imp[:,nonimpute_ind] .= 1
sum_Δ = sum(Δ)

#-------------------------- mcmc algorithm
for it in 2:IT
  global it
  global sum_acc
  global acc
  for k in 1:m
    ν[it,k] = rand(Gamma(a + μ[it-1,k], 1/(1/b[it-1,k] + sum_Δ)))
    b[it,k] = rand(InverseGamma(a+c,γ[it-1]+ν[it,k]))
  end
  γ[it] = rand(Gamma(c*m+1,1 ./ (1+sum(1 ./ b[it,:]))))
  for i in impute_ind
      imp[it,i], acc  = update_imp2(imp[it-1,i],ν[it,:],Δ[i],sols[ind[i]],num_sols[ind[i]],m,proposaltype)
      sum_acc += acc
  end
  μ[it,:] = imp2mu(imp[it,:],sols,ind,m,n)

  if mod(it,10000)==0    println(it)    end
end

#-------------------------- extract iterates for λ and p
λ = sum(ν,dims=2)
p = zeros(IT,m)
for i in 1:IT p[i,:] = ν[i,:]/λ[i]  end

#-------------------------- compute Buchmann-Grubel truncated estimator
λ_bg, p_bg = bg_truncated_estimator(obs,m)
ν_bg = λ_bg * p_bg

#-------------------------- write results to file
include("process_output.jl")
