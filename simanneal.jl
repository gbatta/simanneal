function simanneal(obj_fn,obj_fn_args::Tuple,lb::Vector,ub::Vector,c::Float64,nn_g_max::Int64,nn_a_max::Int64,del::Vector,max_eval::Int64,func_tol::Float64,func_tol_array_size::Int64,noisy::Boolean,stall_iter_lim_sc=5)
  init_cost=obj_fn(obj_fn_args...)[1]
  moment_size=length(obj_fn_args[4])-1
  num_args=length(obj_fn_args)
  tee_accept_zero=init_cost
  tee_accept=tee_accept_zero
  num_params=length(obj_fn_args[1])
  sp=length(obj_fn_args[1])
  tee_gen=Array{Float64}(num_params)
  fill!(tee_gen,1.0)
  tee_gen_zero=Array{Float64}(num_params)
  fill!(tee_gen_zero,1.0)
  kay_i=zeros(Float64,sp)
  kay_a=0.0
  stall_iter_lim=stall_iter_lim_sc*sp
  func_tol_array=zeros(Float64,func_tol_array_size)
  func_tol_ave=Inf
  converge=0
  nn_r=0
  nn_a=0
  nn_g=0
  wwold   = copy(obj_fn_args[1])
  wwbest  = copy(wwold)
  cost_old  =init_cost
  cost_new  =init_cost
  cost_best =init_cost
  bounds    =ub-lb
  nn_eval=0
  s=enumerate(zeros(Float64,sp))
  s_done=Array{Float64}(sp)
  ave_diff=0.0
  srand(500990)
  p_accept=0
  v=0
  tt=0
  if noisy==1
      junk=vcat("vvv","p_accept","nn_eval","nn_g","nn_a","nn_r",fill("wwold",sp),fill("wwnew",sp),fill("wwbest",sp),"cost_old","cost_new","cost_best",fill("s_done",sp),"tee_accept","tee_accept_zero",
      fill("tee_gen",sp),fill("tee_gen_zero",sp),fill("kay_i",sp),"kay_a",fill("func_tol_array",func_tol_array_size),fill("inner_costs",moment_size))
  end

  wwdelta=Array{Float64}(sp)
  while converge==0
    while nn_a<=nn_a_max
      wwnew   = Array{Float64}(num_params)
      for L=1:sp
        tt+=1
        srand(L+L*tt+tt)
        v=rand()
        q=sign(v-0.5)*tee_gen[L]*((1+(1/tee_gen[L]))^(abs(2*v-1))-1)
        test=wwold[L]+q*bounds[L]

        while test<lb[L] || test>ub[L]
          tt+=1
          srand(L*100+L*tt+tt)
          vv=rand()
          q=sign(vv-0.5)*tee_gen[L]*((1+(1/tee_gen[L]))^(abs(2*vv-1))-1)
          test=wwold[L]+q*bounds[L]

        end
        splice!(wwnew,L,test)

      end

#      writecsv("testnew.csv",wwnew)

      nn_g+=1

      # [1] references the value of the cost function; [2] would deliver the Mx1 inner_cost_function array
      COSTS_new=obj_fn(wwnew,obj_fn_args[2:num_args]...)
      cost_new=COSTS_new[1]
      inner_costs=COSTS_new[2]

      p_accept=1/(1+exp((cost_new-cost_old)/tee_accept))
      srand(650000+nn_g)
      vvv=rand()

      nn_eval+=1

      if nn_eval>=max_eval
        converge=1
      end



      if vvv<=p_accept
        for L=1:length(wwnew)
          splice!(wwold,L,wwnew[L])
        end
        cost_old  = cost_new

        if cost_new < cost_best
          cost_best    = cost_new
          for ii=1:length(wwnew)
            splice!(wwbest,ii,wwnew[ii])
          end
        end
        nn_a+=1
      end

      # Check for Annealing
      if nn_g==nn_g_max
        kay_i=kay_i+ones(Float64,sp)

        for L=1:sp
          splice!(tee_gen,L,tee_gen_zero[L]*exp(-c*kay_i[L]^(1/sp)))
        end

        kay_a=kay_a+1
        tee_accept=tee_accept_zero*exp(-c*kay_a^(1/sp))
        nn_g=0
      end
      if noisy==1
          junk=hcat(junk,vcat(vvv,p_accept,nn_eval,nn_g,nn_a,nn_r,wwold,wwnew,wwbest,cost_old,cost_new,cost_best,s_done,tee_accept,tee_accept_zero,tee_gen,tee_gen_zero,kay_i,kay_a,func_tol_array,inner_costs))
          writecsv("simanneal_output.csv",junk)
      end
    end  # end of reannealing loop

    # If converge==1, make sure to break the outer 'while' loop, so we don't go through (time-consuming) reannealing process
    if converge==1
      break
    end

    wwdelta=copy(wwbest)
     # Have to construct matrix equal to number of moments x number of moments, because the wwwdelta_mat
     # function is modified every time it's run
    wwdelta_mat=repmat(wwdelta,1,sp)

    function sensitivities(e,del=del)
        wwdelta_mat[e[1],e[1]]=min(wwbest[e[1]] + del[e[1]],ub[e[1]]) #make sure increment doesn't exceed upper bound
        e=obj_fn(wwdelta_mat[:,e[1]],obj_fn_args[2:num_args]...)[1]
    end

    s_done=map(sensitivities,s)


    s_max=maximum(s_done)
    for i=1:sp
      splice!(kay_i,i,((1/c)*log((tee_gen_zero[i]/tee_gen[i])*(s_max/s_done[i])))^sp)
      splice!(tee_gen,i,(maximum(s_done)/s_done[i])*tee_gen[i])
      splice!(tee_gen_zero,i,tee_gen[i])
    end

    tee_accept_zero =  cost_new     # Reset T_accept_zero to last acepted value, which
                                    # is cost_new if we have entered reannealing
    tee_accept      =  cost_best
    kay_a=((1/c)*log(tee_accept_zero/tee_accept))^sp
    #Reset reannealing cycle
    nn_a=0
    #Increment reannealing counter
    nn_r+=1
    insert!(func_tol_array,1,cost_best) # add new best cost and push out old best cost , for averaging
    pop!(func_tol_array)
    comp=fill(cost_best,size(func_tol_array))
    #Compare latest cost_best to cost_bests from previous 2 reannealing cycles
    func_tol_ave=abs(mean(func_tol_array-comp))

    if (func_tol_ave<=func_tol && nn_r>=stall_iter_lim)
      converge=1
    end
  end
  return (cost_best,wwbest)
end
