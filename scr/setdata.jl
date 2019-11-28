if data_choice=="generated1"
  # Buchmann-Grubel example
  n = 250 # sample size
  Δ = ones(n)   #Δ = rand(Uniform(0,3.0),n)
  λ0 = 2.0
  p0 = [1,1,1]/3
  m0 = length(p0)
  obs =zeros(Int64,n)
  for i in 1:n
     n_jumps = rand(Poisson(λ0*Δ[i]))
     obs[i] = sum(wsample(1:m0,p0,n_jumps))
   end
elseif data_choice=="generated2"
  # geometric distr
  sp = 1/3 # succes probability
  λ0 = -log(sp)
  m0 = 30
  p0 = -(1-sp).^(1:m0)./((1:m0)*log(sp))
  obs = rand(Geometric(sp),n)
elseif data_choice=="testdata_A1"
  obs =readdlm("testdata_A1.csv")
  n = length(obs)
  Δ = ones(n)
  λ0 = 2.0
  p0 = [1,0,0,1,0,1]/3
  m0 = length(p0)
elseif data_choice=="testdata_A2"
  obs =readdlm("testdata_A2.csv")
  n = length(obs)
  Δ = ones(n)
  λ0 = 2.0
  p0 = [1,0,0,1,0,1]/3
  m0 = length(p0)
elseif data_choice=="testdata_A3"
  dat =readdlm("testdata_A3.csv")
  obs = dat[:,1]
  Δ = dat[:,2]
  n = length(obs)
  λ0 = 2.0
  p0 = [1,0,0,1,0,1]/3
  m0 = length(p0)
elseif data_choice=="testdata_B1"
  obs =readdlm("testdata_B1.csv")
  n = length(obs)
  Δ = ones(n)
  sp = 1/3 # succes probability
  λ0 = -log(sp)
  m0 = 30
  p0 = -(1-sp).^(1:m0)./((1:m0)*log(sp))
elseif data_choice=="testdata_B2"
  obs =readdlm("testdata_B2.csv")
  n = length(obs)
  Δ = ones(n)
  sp = 1/6 # succes probability
  λ0 = -log(sp)
  m0 = 30
  p0 = -(1-sp).^(1:m0)./((1:m0)*log(sp))
elseif data_choice=="testdata_C"
  obs =readdlm("testdata_C.csv")
  n = length(obs)
  Δ = ones(n)
  λ0 = 5
  p0 = [0.5, 0, 0.25, 0,  0.25]
  m0 = length(p0)
elseif data_choice=="horsekicks"
  obs = [fill(0,109); fill(1,65); fill(2,22);fill(3,3);fill(4,1)]
  n = length(obs)
  Δ = ones(1,n)
  m0 = maximum(obs)
else
  obs =readdlm("plantpopulation.csv")[2:end,1]
  n = length(obs)
  Δ = ones(1,n)
  m0 = maximum(obs)
end

# ensure observations are integer valued
obs = vec(map(x->Int.(x),obs))
