pwSEM.prepare.data.set<-function(data,grouping.variables=NULL){
  #This function takes an original data set (data), removes any lines with
  #missing values, and sorts it according to the grouping variables.  If there
  #are no grouping variables, it adds a single grouping variable called "group"

  #This line removes lines having any missing values
  data<-data[complete.cases(data),]
  #If there are no grouping variables, then a column called "group" is created
  #with the same value "g" for all lines.
  if(is.null(grouping.variables)){
    data$group=rep("g",dim(data)[1])
    return(data)
  }
  #This line sorts the data frame according to the grouping variables
  data[do.call(what=order,args=
                 data[,c(grouping.variables,grouping.variables)]),]
}

#Create a new "pwSEM" class
pwSEM.class<-function(x){
  structure(x,class="pwSEM.class")
}
#' @export
pwSEM<-function(sem.functions,dependent.errors=NULL,data,
                use.permutations=F,n.perms=5000,do.smooth=F,
                all.grouping.vars=NULL){
  #sem.functions is a list giving the gamm4 or gam models associated with each
  #variable in the sem, including exogenous variables.
  #dependent.errors is a list giving any free covariances, given
  #in the form of X~~y
  #all.grouping.vars is a character vector giving the names of all
  #variables involved in the sem functions that define groups for
  #random effects.
  #remove lines with missing values, sort by grouping variables,
  #and add a new single grouping variable if there are no groups

  data<-pwSEM.prepare.data.set(data=data,grouping.variables=all.grouping.vars)
  n.data.lines<-dim(data)[1]

  n.functions<-length(sem.functions)
  for(i in 1:n.functions){
    if(!(inherits(sem.functions[i][[1]],"gam") |
         inherits(sem.functions[i][[1]]$mer, "lmerMod") |
         inherits(sem.functions[i][[1]]$mer, "glmerMod") |
         inherits(sem.functions[i][[1]],"gamm")))
      stop("Only gam or gamm4 functions can be used in pwSEM")
  }
  dag<-get.dag.from.sem(sem.functions)
  #This gives the names of the variables that are not latents
  not.latent.vars<-row.names(dag)
  if(is.null(dependent.errors))equivalent.mag<-mag<-dag
  if(!is.null(dependent.errors)){
    #This adds free covariances to the adjacency matrix as "100"
    mag<-add.dependent.errors(DAG=dag,dependent.errors=dependent.errors)
    x2<-MAG.to.DAG.in.pwSEM(mag)
    #This gets the names of the latent variables that have been
    #added in the extended DAG (x) to represent the free covariances
    latents<-extract.latents(dag.with.latents=x2,not.latent.vars=not.latent.vars)
    #This gets the d-separation equivalent DAG for this MAG
    equivalent.mag<-DAG.to.MAG.in.pwSEM(full.DAG=x2,latents=latents,
                                        conditioning.latents=NULL)
  }
  basis.set<-basiSet.MAG(equivalent.mag)
  #  basis.set<-basiSet(equivalent.dag)
  if(!is.null(basis.set)){
    out.dsep<-test.dsep.claims(my.list=sem.functions,my.basis.set=basis.set,
                               data=data,use.permutations=use.permutations,do.smooth=do.smooth,
                               n.perms=n.perms,all.grouping.vars=all.grouping.vars,
                               observed.vars=not.latent.vars)
    C.stat<--2*sum(log(out.dsep$null.probs))
    p.C.stat<-1-pchisq(C.stat,df=2*length(out.dsep$null.probs))
    dsep.null.probs<-out.dsep$null.probs
  }
  else C.stat<-p.C.stat<-dsep.null.probs<-NULL
  x<-list(causal.graph=mag,dsep.equivalent.causal.graph=equivalent.mag,basis.set=basis.set,
          dsep.probs=dsep.null.probs,sem.functions=sem.functions,
          C.statistic=C.stat,prob.C.statistic=p.C.stat,
          AIC=get.AIC(sem.functions),n.data.lines=n.data.lines,
          use.permutations=use.permutations,n.perms=n.perms)
  class(x)<-"pwSEM.class"
  x
}


get.dag.from.sem<-function(sem.functions){
  #sem.functions is a list giving the gamm4 or gam models associated with each
  #variable in the sem, including exogenous variables.
  #This function obtains the DAG from this list
  n.vars<-length(sem.functions)
  n.list<-rep(NA,n.vars)
  fo<-rep(list(y~x),n.vars)
  for(i in 1:n.vars){
    if(inherits(sem.functions[i][[1]],"gam"))
      x<-sem.functions[[i]]$formula
    else
      if(inherits(sem.functions[i][[1]],"gamm"))
        x<-sem.functions[[i]]$gam$formula
      else
        if(inherits(sem.functions[i][[1]]$mer, "lmerMod") |
           inherits(sem.functions[i][[1]]$mer, "glmerMod"))
          x<-sem.functions[[i]]$gam$formula
        else stop("Error in get.dag.from.sem")

        if(x[3]!="1()"){
          n.list[i]<-i
          fo[[i]]<-x
        }
  }
  n.list<-n.list[!is.na(n.list)]
  #  strip.formula removes and smoothing calls like s(x)
  for(i in 1:n.vars)fo[[i]]<-strip.formula(fo[[i]])
  #Remember this function!
  do.call(DAG,args=fo[n.list])
}

basiSet.MAG<-function(cgraph){
  #gives the basis set of a DAG or MAG
  #cgraph has existing edges of 0-->1, 100<-->100 or 10--10
  #resulting mag will have existing edges of 0-->1, 100<-->100
  #or 10--10

  mag<-cgraph
  nod<-rownames(mag)
  dv<-length(nod)
  ind<-NULL
  test<-NULL
  for(r in 1:dv){
    for(s in r:dv){
      #if there is an element>0 in column s then r & s are adjacent
      #in mag
      if((mag[r,s]!=0) | (mag[s,r]!=0) | r==s)
        next
      else{
        test<-1
        ed<-nod[c(r,s)]
        pa.r<-nod[mag[,r]==1]
        pa.s<-nod[mag[,s]==1]
        msep<-union(pa.r,pa.s)
        msep<-setdiff(msep,ed)
        b<-list(c(ed,msep))
        ind<-c(ind,b)
      }
    }
  }
  if(is.null(test))cat("No elements in the basis set","\n")
  return(ind)
}

test.dsep.claims<-function(my.list,my.basis.set,data,use.permutations=F,
                           n.perms=5000,do.smooth=F,all.grouping.vars,
                           observed.vars){
  #dsep.claims are the elements returned in basiSet()
  n.claims<-length(my.basis.set)
  out<-data.frame(probs=rep(NA,n.claims))
  n.lines<-dim(data)[1]
  for(i in 1:n.claims){
    y<-get.residuals(my.list=my.list,dsep=my.basis.set[[i]],data=data,
                     do.smooth=do.smooth,all.grouping.vars=all.grouping.vars,
                     observed.vars=observed.vars)
    if(!use.permutations){
      out$probs[i]<-generalized.covariance(R1=y$residuals1,R2=y$residuals2)$prob
    }
    if(use.permutations){
      out$probs[i]<-perm.generalized.covariance(R1=y$residuals1,
                                                R2=y$residuals2,nperm=n.perms)$permutation.prob
    }
  }
  list(basis.set=my.basis.set,null.probs=out$probs)
}

summary.pwSEM.class<-function(x,structural.equations=F){
  #x is an object produced by pwSEM()
  #structural.equations=T to output each structural equation
  cat("Causal graph:","\n")
  var.names<-row.names(x$causal.graph)
  n.vars<-dim(x$causal.graph)[1]
  var.nums<-1:n.vars
  for(i in 1:(n.vars-1)){
    for(j in (i+1):n.vars){
      if(x$causal.graph[i,j]==1)
        cat(var.names[i]," ->",var.names[j],sep="",fill=T)
      if(x$causal.graph[j,i]==1)
        cat(var.names[j]," ->",var.names[i],sep="",fill=T)
      if(x$causal.graph[i,j]==100)
        cat(var.names[i],"<->",var.names[j],sep="",fill=T)
    }
  }
  if(any(x$causal.graph!=x$dsep.equivalent.causal.graph)){
    cat("d-separation equivalent DAG or MAG",fill=T)
    for(i in 1:(n.vars-1)){
      for(j in (i+1):n.vars){
        if(x$dsep.equivalent.causal.graph[i,j]==1)
          cat(var.names[i]," ->",var.names[j],sep="",fill=T)
        if(x$dsep.equivalent.causal.graph[j,i]==1)
          cat(var.names[j]," ->",var.names[i],sep="",fill=T)
        if(x$dsep.equivalent.causal.graph[i,j]==100)
          cat(var.names[i],"<->",var.names[j],sep="",fill=T)
      }
    }
  }
  n<-length(x$basis.set)
  n2<-1:n
  cat("\n")
  cat("Basis Set","\n")
  for(i in 1:n){
    cat("(",n2[i],") ",x$basis.set[[i]][1],"_||_",
        x$basis.set[[i]][2],"| {",x$basis.set[[i]][-c(1,2)],"}",
        sep=" ",fill=T)
  }
  cat("\n")
  if(x$use.permutations){
    cat("Null probabilities are based on permutation method",fill=T)
    cat("with",x$n.perms,"random permutations",fill=T)
    cat("number of observations in data set:",x$n.data.lines,fill=T,"\n")
  }
  if(!x$use.permutations & x$n.data.lines<100)
    cat("Given small sample size, you should specify use.permutations=T",fill=T)
  cat("Null probabilities of d-separation claims in basis set",fill=T)
  for(i in 1:n){
    cat("(",n2[i],") ",x$dsep.probs[i],
        sep="",fill=T)
  }
  cat("\n")
  cat("C-statistic:",x$C.statistic,", df =",2*n,
      ", null probability:",x$prob.C.statistic,fill=T,"\n")
  cat("AIC statistic:",x$AIC,fill=T,"\n")
  if(structural.equations){
    n.funs<-length(x$sem.functions)
    cat("_______Piecewise Functions__________",fill=T)
    cat("\n")

    for(i in 1:n.funs){
      if(inherits(x$sem.functions[[i]],"gam")){
        for.i<-as.character(x$sem.functions[[i]]$formula)
        cat("(",i,"):",for.i[2],for.i[1],for.i[3],fill=T)
        cat("         Parametric coefficients:",fill=T)
        print(summary(x$sem.functions[[i]])$p.table)
        cat("         Smoother terms:","\n",fill=T)

        print(summary(x$sem.functions[[i]])$s.table)
        cat("___________________","\n")
      }
      else{
        for.i<-as.character(x$sem.functions[[i]]$gam$formula)
        cat("(",i,"):",for.i[2],for.i[1],for.i[3],fill=T)
        cat("         Parametric coefficients:",fill=T)
        print(summary(x$sem.functions[[i]]$gam)$p.table)
        cat("         Smoother terms:","\n",fill=T)
        print(summary(x$sem.functions[[i]]$gam)$s.table)
        cat("___________________","\n")
      }
    }


  }
}

add.dependent.errors<-function(DAG,dependent.errors){
  #This function takes a DAG and adds the free covariances to the
  #adjacency matrix as "100" to form a MAG
  #dependent.errors is a list containing free covariances in the
  #form of formulae: x~~y
  n.free<-length(dependent.errors)
  latent.names<-rep(NA,n.free)
  for(i in 1:n.free)
    latent.names[i]<-paste("L",as.character(i),sep="")
  var.names<-row.names(DAG)
  full.names<-c(var.names,latent.names)
  n.vars<-length(var.names)
  MAG<-DAG
  var.nums<-1:n.vars
  for(i in 1:n.free){
    x<-as.character(dependent.errors[[i]])[2]
    y<-gsub("~","",as.character(dependent.errors[[i]])[3])
    x.index<-var.nums[full.names==x]
    y.index<-var.nums[full.names==y]
    MAG[x.index,y.index]<-MAG[y.index,x.index]<-100
  }
  MAG
}

MAG.to.DAG.in.pwSEM<-function(x){
  index = which(x==100,arr.ind=T)
  n <- nrow(index)/2 # pairs of dependent errors
  new.DAG <- matrix(0,ncol=ncol(x)+n,nrow=nrow(x)+n) # new DAG including latent variables
  new.DAG[1:ncol(x),1:nrow(x)]<-x
  # add names to matrix
  colnames(new.DAG) <- c(colnames(x),paste("L",1:n,sep=""))
  rownames(new.DAG) <- c(colnames(x),paste("L",1:n,sep=""))
  new.DAG[(nrow(x)+1):(nrow(x)+n),] <- 0 # the latents do not have causal parents, so set to zero
  # find pairs of observed variables caused by the same latent
  cat <- apply(index,1,paste,collapse="")
  cat2 <- apply(index[,c(2,1)],1,paste,collapse="")
  match <- match(cat,cat2) # matrix row number that has the second pair of the dependent error
  # replace dependent error by directed path from Latent
  no = 1
  for (i in 1:n){
    new.DAG[ncol(x)+no,index[i,1]] <- 1
    new.DAG[ncol(x)+no,index[match[i],1]] <- 1
    no = no+1
    index <- index[-match[i],]
    if (no > n){stop}
  }
  new.DAG[new.DAG==100] = 0
  return(new.DAG)

}

find.possible.Q.in.pwSEM<-function(nvars, x, y) {
  z <- 1:nvars
  z[x] <- z[y] <- 0
  z[z > 0]
}

test.conditioning.latents.in.pwSEM<-
function(latents,conditioning.latents){
  if(is.null(conditioning.latents))return(0)
  if(setequal(conditioning.latents,intersect(latents,conditioning.latents)))return(0)
  return(1)
}

pairs.without.edge.in.pwSEM<-
function(my.graph) {
  nvars<-dim(my.graph)[2]
  com <- combn(1:nvars, 2)
  ncombs <- dim(com)[2]
  keep <- rep(T, ncombs)
  for (i in 1:ncombs) {
    # if(there is an edge between this pair) then remove from com
    if (my.graph[com[1, i], com[2, i]] != 0 |
        my.graph[com[2, i], com[1, i]]!=0) {
      com[1, i] <- com[2, i] <- 0
      keep[i]<-F
    }
  }
  matrix(com[, keep],ncol=sum(keep))
}

full.conditioning.set.in.pwSEM<-
function(observed.conditioning,
         conditioning.latents){
  union(observed.conditioning,conditioning.latents)
}

dag.name.in.pwSEM<-
function (amat,n)
{
  rownames(amat)[n]
}

orient.MAG<-
function(full.DAG,latent.conditioners,cgraph,observed.vars){
  # RETURNS the oriented cgraph

  #This function implements the orientation rules of Richardson & Spirtes
  #full.DAG is the original DAG and latent.conditioners is the names of
  #any latent conditioning variables in it.
  #cgraph is input as a matrix with directed edges (X-->Y),coded as 0-->1,
  #or undirected edges (X--Y), coded as 1--1.
  #This is the matrix after step 3 of Richardson & Spirtes, before the
  #undirected edges are oriented
  #observed.vars is the names of the observed variables in the DAG
  n.observed<-length(observed.vars)
  for(i in 1:n.observed){
    for(j in 1:n.observed){
      # test if there is an undirected edge between variables
      if(cgraph[i,j]==1 & cgraph[j,i]==1){
        test<-0
        # is i ancestral to j?
        if(is.directed.path.in.pwSEM(use.dag=full.DAG,start.var=observed.vars[i],
                            end.var=observed.vars[j])){
          # there is a directed path from i to j
          cgraph[i,j]<-1
          cgraph[j,i]<-0
          next
        }
        # is j ancestral to i?
        if(is.directed.path.in.pwSEM(use.dag=full.DAG,start.var=observed.vars[j],
                            end.var=observed.vars[i])){
          # there is a directed path from i to j
          cgraph[j,i]<-1
          cgraph[i,j]<-0
          next
        }
        # can the edge be turned into correlated errors between i and j?
        if(is.null(latent.conditioners)){
          cgraph[i,j]<-cgraph[j,i]<-100
          next
        }

        else
          n.lc<-length(latent.conditioners)
        for(k in 1:n.lc){
          if(is.directed.path.in.pwSEM(use.dag=full.DAG,start.var=observed.vars[i],
                              end.var= latent.conditioners[k]) &
             is.directed.path.in.pwSEM(use.dag=full.DAG,start.var=observed.vars[j],

                              end.var=latent.conditioners[k])){
            test<-1
            cgraph[i,j]<-cgraph[j,i]<-10
          }
        }
        if(test==0){cgraph[i,j]<-cgraph[j,i]<-100
        }
      }
    }
  }
  cgraph
}

is.directed.path.in.pwSEM<-
function(use.dag,start.var,end.var){
  # start.var=character name of first variable in use.dag
  # end.var=character name of second variable in use.dag
  # if there is a directed path between start.var and end.var, and these two variables
  # are not d-separated in use.dag without conditioning, then there is a directed path between them
  # returns as TRUE or FALSE
  var.names<-colnames(use.dag)
  #start.node is a single number giving the column of use.dag containing start.var
  start.node<-(1:length(var.names))[colnames(use.dag)==start.var]
  #end.node is a single number giving the column of use.dag containing end.var
  end.node<-(1:length(var.names))[colnames(use.dag)==end.var]
  #findPath is a function in ggm that finds one path between two nodes of a graph
  #it returns a vector of numbers giving the sequence of nodes of this path
  #starting at st and going to en
  test1<-length(findPath(amat=use.dag,st=start.node,en=end.node))>0
  test2<-!dSep(amat=use.dag,first=start.var,second=end.var,cond=NULL)
  #if TRUE, there is a path from start.var to end.var and this path has no colliders
  test1 & test2
  return(test1 & test2)
}

DAG.to.MAG.in.pwSEM<-function (full.DAG, latents = NA,
                               conditioning.latents=NULL)
{
  # full.DAG is a binary (0/1) matrix produced from DAG() function in ggm
  # latents is a character vector giving the names of the latents
  # conditioning.latents is a character vector giving the names
  # of those latents that serve as conditioning variables for
  # sampling (i.e. sampling bias)
  # The final Mixed Acyclic graph is returned.

  ####################

  #Requires the ggm library
  #First, make sure that the conditioning latents, if present, are a subset
  # of the declared latents.  If not, stop.
  if(test.conditioning.latents.in.pwSEM(latents,conditioning.latents)!=0){
    stop("Conditioning latents must be a proper subset of all latents")
  }
  #########################################
  # main function
  #
  full.vars<-row.names(full.DAG)
  full.vars.index<-1:length(full.vars)
  n.observed<-length(full.vars)-length(latents)
  observed.DAG<-full.DAG
  observed.vars<-full.vars
  observed.vars.index<-full.vars.index
  for(i in 1:length(latents)){
    observed.vars[latents[i]==full.vars]<-NA
    observed.vars.index[latents[i]==full.vars]<-NA
    observed.DAG[latents[i]==full.vars,]<-NA
    observed.DAG[,latents[i]==full.vars]<-NA
  }
  latent.vars.index<-match(latents,full.vars)
  total.n.vars<-dim(full.DAG)[2]

  if(sum(is.na(latents))>0){
    return(full.DAG)
  }

  if(sum(is.na(latents))==0){
    n.latents<-length(latents)
    for(i in 1:n.latents){
      ok<-F
      for(j in 1:length(full.vars))if(latents[i]==full.vars[j])ok<-T
      if(!ok)stop("ERROR: latent variable name not in the DAG")
    }
  }
  observed.vars<-observed.vars[!is.na(observed.vars)]
  observed.vars.index<-observed.vars.index[!is.na(observed.vars.index)]
  #
  #STEP 2 of Shipley & Douma
  # construct initial observed DAG by removing latents and conserving directed
  # edges between pairs of observed variables
  #
  if(n.observed<=0)stop(cat("No observed variables","\n"))
  if(n.observed==1)stop(cat("Only one observed variable","\n"))
  if(n.observed==2)stop(cat("Only two observed variables","\n"))

  observed.DAG<-observed.DAG[observed.vars.index,observed.vars.index]

  if(n.observed<=0){
    stop(cat("All variables are latent; there is no equivalent observed DAG","\n"))
  }

  #Finished STEP 2 of Shipley & Douma.
  #Start STEP 3
  pairs.to.test<-pairs.without.edge.in.pwSEM(observed.DAG)
  n.pairs.to.test<-dim(pairs.to.test)[2]
  n.remaining<-length(observed.vars)-2
  # if all observed variables share an edge then return...
  if(n.pairs.to.test<=0){
    stop(cat("Since there are only two observed variables, nothing further will be done","\n"))
  }
  add.edge<-matrix(NA,nrow=2,ncol=n.pairs.to.test)
  # for each pair (i) to test, determine dsep in full graph given only the observed variables
  # plus all conditioning latents.

  kount<-0
  # i cycles over each pair that are not adjacent...
  for(i in 1:n.pairs.to.test){
    is.pair.dsep<-F
    # get those other observed variables in graph except this pair...
    possible.Q<-find.possible.Q.in.pwSEM(n.observed,pairs.to.test[1,i],pairs.to.test[2,i])

    # Do do unconditional dseparation...
    # i.e. conditional order=0
    first.var<-observed.vars.index[pairs.to.test[1,i]]
    second.var<-observed.vars.index[pairs.to.test[2,i]]
    test<-dSep(amat=full.DAG,first=dag.name.in.pwSEM(full.DAG,first.var),
               second=dag.name.in.pwSEM(full.DAG,second.var),
               cond=full.conditioning.set.in.pwSEM(NULL,
                                          conditioning.latents))
    # if first.var is dsep from second.var then there is no edge between them;
    if(test){
      is.pair.dsep<-T
      next
    }
    # if here then there are potential conditional variables to consider
    # so cycle through all possible conditional orders...
    if(sum(is.na(possible.Q)==0)){
      n.possible.Q<-length(possible.Q)
      #now, determine, using observed.vars.index[possible.Q], if the pair are dsep
      # in the full graph
      # j gives the conditional order for a given pair
      for(j in 1:n.possible.Q){

        # Q has column = different combinations and rows=elements in each combination
        dQ<-combn(possible.Q,j)

        if(j==n.possible.Q) dQ<-matrix(possible.Q,nrow=j,ncol=1)

        n.Q<-dim(dQ)[2]

        first.var<-observed.vars.index[pairs.to.test[1,i]]

        #   pairs.to.test[1,i],"dag name=",dag.name.in.pwSEM(full.DAG,first.var),"\n")
        second.var<-observed.vars.index[pairs.to.test[2,i]]

        #    pairs.to.test[2,i],"dag.name=",dag.name.in.pwSEM(full.DAG,second.var),"\n")
        # k cycles through these different combinations
        for(k in 1:n.Q){
          cond.vars<-as.vector(observed.vars.index[dQ[,k]])
          test<-dSep(amat=full.DAG,first=dag.name.in.pwSEM(full.DAG,first.var),second=dag.name.in.pwSEM(full.DAG,second.var),
                     cond=full.conditioning.set.in.pwSEM(dag.name.in.pwSEM(full.DAG,cond.vars),
                                                conditioning.latents)
          )
          # if first.var dsep from second.var then there is no edge...
          if(test){
            is.pair.dsep<-T
            break
          }
        }
      }
    }
    if(!is.pair.dsep){
      kount<-kount+1
      add.edge[1,kount]<-pairs.to.test[1,i]
      add.edge[2,kount]<-pairs.to.test[2,i]
    }
  }

  #Add undirected edges to non-adjacent pairs in
  #the observed graph that are not d-separated in the DAG given
  #any combination of other observed variables PLUS any latent conditioning
  #variables.

  # convert observed DAG to a partially oriented graph
  cgraph<-matrix(0,n.observed,n.observed,dimnames=list(observed.vars,observed.vars))
  for(i in 1:(n.observed-1)){
    for(j in (i+1):n.observed){
      if(observed.DAG[i,j]==1 & observed.DAG[j,i]==0){
        cgraph[i,j]<-1
        cgraph[j,i]<-0
      }
      if(observed.DAG[j,i]==1 & observed.DAG[i,j]==0){
        cgraph[j,i]<-1
        cgraph[i,j]<-0
      }
    }
  }
  for(i in 1:kount){
    cgraph[add.edge[1,i],add.edge[2,i]]<-cgraph[add.edge[2,i],add.edge[1,i]]<-1
  }

  #cgraph now holds the partially oriented inducing graph, with X--Y if there is an inducing
  # path between observed variables X & Y.

  #Now, orient these if there are directed paths from i to j
  cgraph<-orient.MAG(full.DAG=full.DAG,latent.conditioners=conditioning.latents,
                     cgraph=cgraph,observed.vars=observed.vars)
  return(cgraph)
}


extract.latents<-function(dag.with.latents,not.latent.vars){
  #extracts the latent variables from a DAG
  vars<-row.names(dag.with.latents)
  vars[-which(vars %in% not.latent.vars)]
}

strip.formula<-function(fo){
  #This takes a formula (fo),  strips out any smoothing
  #calls, and then returns it as a formula
  x<-as.character(fo)
  dep<-x[2]
  temp.ind<-x[3]
  temp.ind<-gsub("\\)","",gsub("s\\(","",temp.ind))
  temp.ind
  formula(paste(dep,"~",temp.ind))
}

get.AIC<-function(sem.functions){
  #calculates the AIC statistic for the set of sem.functions
  #based on Shipley, B. & Douma, J.C. 2020. Generalized AIC
  #and chi-squared statistics for path models consistent with
  #directed acyclic graphs.  Ecology 101(3):e02960.
  n.functions<-length(sem.functions)
  AIC.vector<-rep(NA,n.functions)
  for(i in 1:n.functions){
    if(inherits(sem.functions[i][[1]],"gam")){
      AIC.vector[i]<-AIC(sem.functions[i][[1]])
    }
    else
      if(inherits(sem.functions[i][[1]]$mer,"lmerMod") |
         inherits(sem.functions[i][[1]],"gamm") |
         inherits(sem.functions[i][[1]]$mer,"glmerMod")){
        AIC.vector[i]<-AIC(sem.functions[i][[1]]$mer)
      }
  }
  sum(AIC.vector)
}

get.residuals<-function(my.list,dsep,data,do.smooth,
                        all.grouping.vars,observed.vars){
  #my.list is the list of gam or gamm functions for the dag
  #dsep is the dsep claim X_||Y|Z for which the residuals of X~Z and Y~Z are
  #returned.
  #data is the data set to be used
  #observed.vars holds the names of all observed variables in MAG
  info<-set.up.info.for.dsep.regressions(my.list,
                                         all.grouping.vars=all.grouping.vars)
  n.vars<-length(info$var.name)
  if(n.vars!=length(observed.vars))
    stop("You have not modelled all variables in the DAG/MAG.
         You must also explicitly model all exogenous variables")
  else
    if(sum(sort(info$var.name)!=sort(observed.vars))>0)
      stop("You have not modelled all variables in the DAG/MAG.
         You must also explicitly model all exogenous variables")

  #This holds the variable names and their number in the dag

  set.var.nums<-data.frame(var.name=info$var.name,var.num=1:n.vars)
  #HERE, WE HAVE TO GET THE VARIABLE NUMBER (v1, v2) FOR EACH OF THE TWO VARIABLES IN THE
  #DSEP CLAIM

  v1<-set.var.nums$var.num[set.var.nums$var.name==dsep[1]]
  v2<-set.var.nums$var.num[set.var.nums$var.name==dsep[2]]
  #fo holds the two formulae X~Z and Y~Z
  fo<-create.formulae.from.basis.set(dsep.claim=dsep,do.smooth=do.smooth)
  #n.levels are the number of nesting levels for X and Y
  n.levels1<-sum(!is.na(info$grouping.structure[,set.var.nums$var.name==dsep[1]]))
  n.levels2<-sum(!is.na(info$grouping.structure[,set.var.nums$var.name==dsep[2]]))
  #if n.levels==0 then use gam, since no nesting structure
  #if n.levels>0 then use gamm for X
  if(n.levels1==0){
    r1<-residuals(gam(formula=fo$formula1,family=info$family[[v1]],data=data))
  }
  if(n.levels2==0){
    r2<-residuals(gam(formula=fo$formula2,family=info$family[[v2]],data=data))
  }
  if(n.levels1>0){
    fit<-gamm4(formula=fo$formula1,random=info$random[[v1]],
               family=info$family[[v1]],data=data)
    r1<-residuals(fit$mer,type="response")
  }
  if(n.levels2>0){
    fit<-gamm4(formula=fo$formula2,random=info$random[[v2]],
               family=info$family[[v2]],
               data=data)
    r2<-residuals(fit$mer,type="response")
  }
  data.frame(residuals1=r1,residuals2=r2)
}

generalized.covariance<-function(R1,R2){
  #generalized covariance measure from Shah, R.D. & Peters,J. 2020.
  #The hardness of conditional independence testing and the generalized
  #covariance measure. The Annals of Statistics 48:1514-1538.
  #
  #R1 & R2 are the residuals of two regressions and n is the sample size
  if(length(R1)!=length(R2))stop("In generalized.covariance, lengths of residual vectors not equal")
  n<-length(R1)
  num<-sqrt(n)*(1/n)*sum(R1*R2)
  a<-sum((R1*R2)^2)/n
  b<-(sum(R1*R2)/n)^2
  denom<-sqrt(a-b)
  T.stat<-num/denom
  list(T.stat=T.stat,prob=2*(1-pnorm(abs(T.stat))))
}

perm.generalized.covariance<-function(R1,R2,nperm=5000,plot.hist=F){
  #R1, R2 are vectors holding the residuals; i.e. E[function]-observed
  #The R1 vector is permuted each time.
  T.stat<-generalized.covariance(R1,R2)$T.stat
  perm.T<-rep(NA,nperm)
  if(length(R1)!=length(R2))stop("In perm.generalized.covariance, residual vectors not same length")
  n<-length(R1)
  for(i in 1:nperm){
    for(j in n){
      sel<-sample.int(n)
      perm.T[i]<-generalized.covariance(R1[sel],R2)$T.stat
    }
  }
  prob<-(sum(abs(as.numeric(perm.T))>=abs(T.stat))+1)/length(as.numeric(perm.T))
  CI.prob<-1.96*sqrt((prob*(1-prob)/nperm))
  if(plot.hist){
    hist(as.numeric(perm.T))
    abline(v=T.stat,col="red")
  }
  list(permutation.prob=prob,lower95.CI=prob-CI.prob,upper95.CI=prob+CI.prob)
}

set.up.info.for.dsep.regressions<-function(fun.list,
                                           all.grouping.vars){
  #fun.list is the list of "n" gam(m) functions describing the links between variables
  #in the DAG
  #all.grouping.vars is a character vector containing the names
  #of all the grouping variables included in fun.list
  n<-length(fun.list)
  #"n" is the number of variables, and gamm4 or gam models, in the DAG
  out<-list(var.name=rep(NA,n),family=rep(NA,n),link=rep(NA,n),
            grouping.structure=matrix(NA,nrow=50,ncol=n),
            random=list())
  for(i in 1:n){
    x<-NA
    if(inherits(fun.list[i][[1]], "gam")){
      x<-extract.variable.info.from.gam(fun.list[i][[1]])
    }
    else
      if(inherits(fun.list[i][[1]]$mer, "lmerMod") |
         inherits(fun.list[i][[1]]$mer, "glmerMod") |
         inherits(fun.list[i][[1]],"gamm")){
        x<-extract.variable.info.from.gamm4(fun.list[i][[1]],
                                            all.grouping.vars=all.grouping.vars)
      }
    out$var.name[i]<-x$var.name[1]
    out$family[i]<-x$family[1]
    out$link[i]<-x$family[2]
    #"n2" is the number of nesting levels for that variable
    n2<-length(x$grouping.structure)
    if(inherits(fun.list[i][[1]]$mer, "lmerMod") |
       inherits(fun.list[i][[1]]$mer, "glmerMod")){
      out$random[[i]]<-x$random.formula
      out$grouping.structure[1:n2,i]<-x$grouping.structure
    }
    if(inherits(fun.list[i][[1]], "gam"))
      out$grouping.structure[1:n2,i]<-NA
  }

  r<-rep(NA,50)
  for(i in 1:50)r[i]<-paste("level",as.character(i))
  dimnames(out$grouping.structure)<-list(r,out$var.name)
  out
}

create.formulae.from.basis.set<-function(dsep.claim,do.smooth=F){
  #This function creates 2 formulae from the 2 elements of a dsep claim
  #dsep.claim is a dsep claim taken from the basiSet() function
  ncond<-length(dsep.claim)-2
  #If ncond>0 then add the first conditioning variable to cond
  if(ncond>0 & do.smooth)cond<-paste("s(",dsep.claim[3],")")
  if(ncond>0 & !do.smooth)cond<-paste(dsep.claim[3])
  #if ncond>1 then add these to the conditioning set
  if(ncond>1){
    for(i in 2:ncond){
      if(do.smooth)cond<-paste(cond,"+","s(",dsep.claim[i+2],")")
      if(!do.smooth)cond<-paste(cond,"+",dsep.claim[i+2])
    }
  }
  if(ncond>0){
    fo1<-as.formula(paste(dsep.claim[1],"~",cond))
    fo2<-as.formula(paste(dsep.claim[2],"~",cond))
  }
  else if(ncond==0){
    fo1<-as.formula(paste(dsep.claim[1],"~",as.character(1)))
    fo2<-as.formula(paste(dsep.claim[2],"~",as.character(1)))
  }
  else stop("Error")
  list(formula1=fo1,formula2=fo2)
}

extract.variable.info.from.gamm4<-function(fo,all.grouping.vars){
  #Given a gamm object,this function extracts the following information:
  #the dependent variable name
  #the distribution and link function
  #the names of the nesting variables
  #the random part of the function
  #fo is the gamm4 formula
  #all.grouping.vars is a character vector giving all of the random
  #grouping variables over all functions in the SEM
  if(!inherits(fo$mer, "lmerMod") &
     !inherits(fo$mer, "glmerMod") &
     !inherits(fo,"gamm"))stop("ERROR. The function passed to
    extract.variable.info.from.gamm() is not an allowed model object")
  #This takes a gamm() object (fo) and extracts the following information:
  #var.name: the name of the dependent variable
  var.name<-as.character(fo$gam$formula[2][[1]])
  #family: distribution & link function
  family<-fo$gam$family

  #grouping.structure: names of the grouping variables
  if(inherits(fo,"gamm")){
    grouping.structure<-names(fo$lme$groups)
  }
  else{
    if(inherits(fo$mer, "lmerMod") |
       inherits(fo$mer, "glmerMod"))
      grouping.structure<-names(summary(fo$mer)$ngrps)
  }
  grouping.structure<-grouping.structure[
    which(grouping.structure %in% all.grouping.vars)]

  # making random formula: the random part of the formula for the dependent variable
  n.groups<-length(grouping.structure)
  if(n.groups>0){
    grouping.structure<-grouping.structure[
      !is.na(match(grouping.structure,table=all.grouping.vars))]
  }
  if(n.groups==1)random.fo<-paste("~(1 |",grouping.structure[1],")")
  if(n.groups>1){
    random.fo<-paste("~(1 |",grouping.structure[1],") +")
    for(i in 2:n.groups){
      if(i<max(n.groups))random.fo<-paste(random.fo,
                                          "(1 | ",grouping.structure[i],") +")
      if(i==max(n.groups))random.fo<-paste(random.fo,
                                           "(1 | ",grouping.structure[i],")")
    }
  }
  list(var.name=var.name,family=family,grouping.structure=grouping.structure,
       random.formula=formula(random.fo))
}

extract.variable.info.from.gam<-function(fo){
  if(!inherits(fo, "gam"))stop("ERROR. The function passed to
    extract.variable.info.from.gam() is not a gam object")
  #This takes a gam() object (fo) and extracts the following information:

  #var.name: the name of the dependent variable
  var.name<-as.character(fo$formula[2])
  #family: distribution & link function
  family<-fo$family
  #The name of the single grouping structure is "no.groups"
  grouping.structure<-NA
  names(grouping.structure)<-"no.groups"
  list(var.name=var.name,family=family,grouping.structure=grouping.structure)
}
