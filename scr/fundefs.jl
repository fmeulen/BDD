function diophantine(z,m)
    lz = size(z)[1]
    @rput z
    @rput lz
    @rput m

    R"""
    library(nilde)  

    num_sols <-rep(0,lz)
    sols <- list(lz)

    # check if the first element of z equals zero, in that case manually add the trivial solution
    if (z[1]==0)
    {
      sols[[1]] <- matrix(0,nrow=m,ncol=1)
      num_sols[[1]] <- 1
      start_i <- 2
    }  else {
      start_i <- 1
    }

    for (i in start_i:lz)
    {
      U <- nlde(1:m,z[i])
      sols[[i]] <- U$solutions
      num_sols[i] <- U$p.n
    }
    """

    # import variables in Julia
    @rget sols
    @rget num_sols

    # convert to integer valued arays
    sols = map(x->Int.(x),sols)
    num_sols = map(x->Int.(x),num_sols)
    sols, num_sols
end


function imp2mu(impvec,sols,ind,m,n)
  # in: impvec is a row of imp, and hence a vector of integers
  # out: compute mu_1,...,mu_m for impvec
  μ =zeros(Int64,m)
  for i in 1:n
    col_i = impvec[i]  # which column to select
    μ += sols[ind[i]][:,col_i] # select imputed solution
  end
  μ
end



function update_imp2(impval,ν,Δi,sols_i,num_sols_i,m,proposaltype)
  # sols_i is the matrix of Diophantine solutions; num_sols_i the number of these solutions
  # propose new index
  if rand()> proposaltype
        if impval==1
            impval_circ = sample([2,num_sols_i])
        elseif impval==num_sols_i
            impval_circ = sample([1,num_sols_i-1])
        else
            impval_circ = sample(impval .+ [-1,1])
        end
   else
      impval_circ = sample(1:num_sols_i)  # randomly choose a solution
      end
      μi = sols_i[:,impval]
      μi_circ = sols_i[:,impval_circ]
      A  = 0
      for j in 1:m
        A += (μi_circ[j] - μi[j]) * log(Δi * ν[j]) + lfactorial(μi[j]) - lfactorial(μi_circ[j])
      end
      if log(rand()) <  A
        impval = impval_circ
        acc = 1
      else
        acc =0
      end
      impval, acc
end


"""
Compute Buchmann-Grubel truncated plug-in estimator.
"""
function bg_truncated_estimator(obs,m)
  qhat = table(obs)
  λhat = -log(qhat[1])
  #m = maximum(obs)
  p_BG = zeros(m)
  p_BG[1] = max(0, min(qhat[2]/(qhat[1]*λhat),1))
  for k=2:m
    ss = 0
    for j=1:k-1
      ss += j * p_BG[j] * qhat[k+1-j]
    end
    x_nk = (qhat[k+1]/λhat - ss/k)/qhat[1]   # sticking to BG notation
    p_BG[k] = max(0,min(x_nk,1-sum(p_BG[1:k-1])))
  end
  λhat,  p_BG
end

function table(X)     map(x->count(y->x==y,X),0:maximum(X))/length(X)  end
