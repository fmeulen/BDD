############# Print some output to screen
ν_postmean = mean(ν[BI:IT,:],dims=1)'
λ_postmean = mean(λ[BI:IT],dims=2)[1]
p_postmean = mean(p[BI:IT,:],dims=1)'

# compute for each ν_k 2.5% and 97.5% quantiles
α = 0.05
ν_95low = mapslices(v-> quantile(v,  α/2), ν[BI:IT,:], dims=1)'
ν_95up = mapslices(v-> quantile(v,  1-α/2), ν[BI:IT,:], dims=1)'

acc_rate = sum_acc/(length(impute_ind)*(IT-1))
println("Acceptance percentage of imputed segments equals: ",
            round(100*acc_rate;digits=2),"%")

if data_choice in generated
  if m>=m0
    p00 = [p0; zeros(m-m0)]
  else
    p00 = p0
  end
  ν0 = λ0 * p00
  if m0>=m
    ν_postmean = [ν_postmean; zeros(m0-m)]
    ν_95low = [ν_95low; zeros(m0-m)]
    ν_95up = [ν_95up; zeros(m0-m)]
    ν_bg= [ν_bg; zeros(m0-m)]
    p_postmean = [p_postmean; zeros(m0-m)]
    p_bg = [p_bg; zeros(m0-m)]
  end
  error_postmean = sum(abs.(ν_postmean - ν0 ))
  error_bg = sum(abs.( ν_bg - ν0 ))
end

if data_choice in generated
  λ_out = [λ_postmean  λ_bg λ0]
  νp_out = [ν_postmean ν_bg ν0 ν_95low ν_95up p_postmean  p_bg  p00]
  [error_postmean  error_bg]
else
  λ_out = [λ_postmean  λ_bg]
  νp_out = [ν_postmean ν_bg  ν_95low ν_95up   p_postmean   p_bg]
end

############# Write output to csv files (ensure a relative directory with name "out" exists)

# write info about the simulation
facc = open("./out/info.txt","w")
write(facc, "Data choice: ",data_choice,"\n")
write(facc, "Number of iterations: ",IT,"\n")
if data_choice in generated
    write(facc,"lambda0 = ", string(λ0),"\n")
    write(facc,"p0 = ", string(p0),"\n\n")
    write(facc,"nu0 = ", string(ν0),"\n\n")
    write(facc,"Error for Buchmann-Grubel= ", string(error_bg),"\n")
    write(facc,"Error for posterior mean of nu= ", string(error_postmean),"\n\n")
end
write(facc,"Prob ot picking from the uniform distribution on all Diophantine solutions= ",string(proposaltype),"\n")
write(facc, "Average acceptance probability equals: ",string(round(mean(acc_rate);digits=3)),"\n\n")
write(facc,"n is number of observations, support of jumps size distribution is {1,...,m}","\n")
write(facc, "[ n,  m] = ",string([n, m]),"\n\n")
#write(facc, "elapsed time ",string(elapsed_time), "\n\n")
write(facc, "---- Prior specification ----","\n")
write(facc, "a= ",string(a),"\n")
write(facc, "c= ",string(c,"\n"))
write(facc, "average value of gamma= ",string(round(mean(γ[BI:IT]);digits=2),"\n\n"))
write(facc, "[lambda_postmean lambda_bg lambda0]=", string(λ_out),"\n")
close(facc)

# write observations to csv file
writedlm("./out/observations.csv",obs,',')   # observations

# write mcmc iterates to csv file
f = open("./out/iterates.csv","w")
indices = 1:m
h0 = "iterate, lambda,gamma,"
h1 = prod("nu$i,"  for i in indices)
h2 = prod("b$i,"  for i in indices)
h3 = prod("p$i"*(i == indices[end] ? "\n" : ", ")  for i in indices)
head = h0 * h1 * h2 * h3
write(f, head)
iterates = [γ λ ν b p]
subsamp = 2:50:IT
writedlm(f,[subsamp iterates[subsamp,:]],',')
close(f)

# write estimates for ν and p to filemode
f = open("./out/nu_p.csv","w")
if data_choice in generated
    head = "nu_postmean,nu_bg,nu0,  nu_95low, nu_95up,  p_postmean,   p_bg, p0"
else
    head = "nu_postmean,nu_bg, nu_95low, nu_95up,  p_postmean, p_bg"
end
write(f,head,"\n")
writedlm(f,νp_out,',')
close(f)

########## plotting
if data_choice in generated
    nu_p = DataFrame(nu_postmean=νp_out[:,1],
                nu_bg=νp_out[:,2],nu0=νp_out[:,3],
                nu_95low=νp_out[:,4], nu_95up=νp_out[:,5],
                p_postmean=νp_out[:,6], p_bg=νp_out[:,7],
                p0=νp_out[:,8])
    @rput nu_p
    R"""
        #pdf('estimates.pdf',width=7, height=3)
        library(tidyverse)
        nu_p %>% mutate(nu=factor(1:nrow(nu_p)))%>% filter(as.numeric(nu)<11) %>%
        ggplot(aes(x=nu, y=nu_postmean)) +
        geom_point(aes(x=nu,y=nu0),colour="gray71",shape=16, size=4)+
        geom_point(shape=8, size=2,colour="black")+
        geom_errorbar(width=.2, aes(ymin=nu_95low, ymax=nu_95up),colour="black")+
        geom_point(aes(x=nu,y=nu_bg),colour="red2",size=2,shape=17) +
        xlab("")+ylab("")
        #+  theme(plot.title = element_text(hjust = 1))+ggtitle("a")
        #dev.off()
    """
end
println("error postmean: ",   error_postmean)
println("error bg: ", error_bg)
