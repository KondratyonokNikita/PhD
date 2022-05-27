Sys.setlocale("LC_TIME", "C")
library(data.table)
options(java.parameters = "-Xmx10000m")

Sys.setlocale("LC_TIME", "C")
options(stringsAsFactors = FALSE)

# a + b omega
# (1, omega) is integral basis
# omega = \sqrt{\theta} if \theta \not\equiv 1(\mod 4)
# omega = \frac{1 + \sqrt{\theta}}{2} if \theta \equiv 1(\mod 4)

sum <- function(x, y){
  res <- list()
  res[["a"]] <- x[["a"]] + y[["a"]]
  res[["b"]] <- x[["b"]] + y[["b"]]
  res[["theta"]] <- x[["theta"]]
  return(res)
}

Prod <- function(x, y){
  res <- list()
  if (x[["theta"]] %% 4 != 1) {
    res[["a"]]<-x[["a"]]*y[["a"]]+x[["b"]]*y[["b"]]*x[["theta"]]
    res[["b"]]<-x[["a"]]*y[["b"]]+x[["b"]]*y[["a"]]
  } else {
    res[["a"]]<-x[["a"]]*y[["a"]]+x[["b"]]*y[["b"]]*(x[["theta"]]-1)/4
    res[["b"]]<-x[["a"]]*y[["b"]]+x[["b"]]*y[["a"]]+x[["b"]]*y[["b"]]
  }
  res[["theta"]]<-x[["theta"]]
  return(res)
}

IsEqual <- function(x, y){
  eps<-10^{-6}
  return(as.logical(abs(x[["a"]]-y[["a"]])<eps & abs(x[["b"]]-y[["b"]])<eps & x[["theta"]]==y[["theta"]]))
}

IsEquiv <- function(x, y, m){ # m is rational integer
  res<-sum(x, Inverse(y))
  res[["a"]]<-round(res[["a"]])
  res[["b"]]<-round(res[["b"]])
  return(as.logical(res[["a"]]%%m==0 & res[["b"]]%%m==0))
}

Norm <- function(x){
  if (x[["theta"]] %% 4 != 1) {
    return(
      x[["a"]]*x[["a"]]+x[["b"]]*x[["b"]]*(-x[["theta"]])
    )
  } else {
    return(
      ((2*x[["a"]]+x[["b"]])*(2*x[["a"]]+x[["b"]])+x[["b"]]*x[["b"]]*(-x[["theta"]]))/4
    )
  }
}

NormIdeal <- function(x){
	return(abs(Norm(x)))
}

GetAdjoint <- function(x){
  if (x[["theta"]] %% 4 != 1) {
  	res<-x
  	res[["b"]]<--res[["b"]]
  	return(res)
  } else {
    res<-x
    res[["a"]]<-res[["a"]]+res[["b"]]
    res[["b"]]<--res[["b"]]
  return(res)
  }
}

Div <- function(x, y){
	res<-Prod(x, GetAdjoint(y))
	res[["a"]]<-res[["a"]]/Norm(y)
	res[["b"]]<-res[["b"]]/Norm(y)
	return(res)
}

Unit <- function(theta){
	res<-list()
	if (theta==2){# done
		res[["a"]]<-1
		res[["b"]]<-1
		res[["theta"]]<-theta
	} else if (theta==3) {# done
	  res[["a"]]<-2
	  res[["b"]]<-1
	  res[["theta"]]<-theta
	} else if (theta==5) {# using omega # done
	  res[["a"]]<-0
	  res[["b"]]<-1
	  res[["theta"]]<-theta
	} else if (theta==6) {# done
	  res[["a"]]<-5
	  res[["b"]]<-2
	  res[["theta"]]<-theta
	} else if (theta==7) {# done
	  res[["a"]]<-8
	  res[["b"]]<-3
	  res[["theta"]]<-theta
	} else if (theta==11) {# done
	  res[["a"]]<-10
	  res[["b"]]<-3
	  res[["theta"]]<-theta
	} else if (theta==13) {# using omega # done
	  res[["a"]]<-1
	  res[["b"]]<-1
	  res[["theta"]]<-theta
	} else if (theta==17) {# using omega # done
	  res[["a"]]<-3
	  res[["b"]]<-2
	  res[["theta"]]<-theta
	} else if (theta==19) {# done
	  res[["a"]]<-170
	  res[["b"]]<-39
	  res[["theta"]]<-theta
	} else if (theta==21) {# using omega # done
	  res[["a"]]<-2
	  res[["b"]]<-1
	  res[["theta"]]<-theta
	} else if (theta==29) {# using omega # done
	  res[["a"]]<-2
	  res[["b"]]<-1
	  res[["theta"]]<-theta
	} else if (theta==33) {# using omega # done
	  res[["a"]]<-19
	  res[["b"]]<-8
	  res[["theta"]]<-theta
	} else if (theta==37) {# using omega # done
	  res[["a"]]<-5
	  res[["b"]]<-2
	  res[["theta"]]<-theta
	} else if (theta==41) {# using omega # done
	  res[["a"]]<-27
	  res[["b"]]<-10
	  res[["theta"]]<-theta
	} else if (theta==57) {# using omega # done
	  res[["a"]]<-131
	  res[["b"]]<-40
	  res[["theta"]]<-theta
	} else if (theta==73) {# using omega
	  res[["a"]]<-943
	  res[["b"]]<-250
	  res[["theta"]]<-theta
	}
	return(res)
}

IntUnit <-function(theta){
	res<-list()
	res[["a"]]<-1
	res[["b"]]<-0
	res[["theta"]]<-theta
	return(res)
}

NegIntUnit <-function(theta){
  res<-list()
  res[["a"]]<--1
  res[["b"]]<-0
  res[["theta"]]<-theta
  return(res)
}

Inverse <- function(x){
	res<-list()
	res[["a"]]<--x[["a"]]
	res[["b"]]<--x[["b"]]
	res[["theta"]]<-x[["theta"]]
	return(res)
}

FractionalPartRational <- function(x, y){ # x is integer elements of ring, y is positive rational integer
	res<-list()
	res[["a"]]<-x[["a"]]-floor(x[["a"]]/y)*y
	res[["b"]]<-x[["b"]]-floor(x[["b"]]/y)*y
	res[["theta"]]<-x[["theta"]]
	return(list(x=res, y=y))
}

Abs <- function(x){
  if (x[["theta"]] %% 4 != 1) {
  	return(abs(x[["a"]]+x[["b"]]*sqrt(x[["theta"]])))
  } else {
    return(abs(x[["a"]]+x[["b"]]*(1+sqrt(x[["theta"]]))/2))
  }
}

GetOrbit <- function(a, b){
	base.unit<-Unit(a[["theta"]])
	norm.b<-Norm(b)
	initial.frac.part<-FractionalPartRational(Prod(a, GetAdjoint(b)), norm.b)
	curr.iteration<-initial.frac.part[["x"]]
	curr.degree<-0
	
	orbit.a<-curr.iteration[["a"]]
	orbit.b<-curr.iteration[["b"]]
	orbit.sign<-1
	orbit.degree<-curr.degree
	
	while (!IsEquiv(x=initial.frac.part[["x"]],
	                y=FractionalPartRational(Prod(curr.iteration, base.unit), norm.b)[["x"]],
	                m=norm.b)){
		curr.iteration<-FractionalPartRational(Prod(curr.iteration, base.unit), norm.b)[["x"]]
		curr.degree<-curr.degree + 1
		orbit.a<-c(orbit.a, curr.iteration[["a"]])
		orbit.b<-c(orbit.b, curr.iteration[["b"]])
		orbit.sign<-c(orbit.sign, 1)
		orbit.degree<-c(orbit.degree, curr.degree)
		
	}
	
	base.neg_unit<-Prod(Unit(a[["theta"]]), NegIntUnit(a[["theta"]]))
	initial.frac.part<-FractionalPartRational(Prod(a, GetAdjoint(b)), norm.b)
	curr.iteration<-initial.frac.part[["x"]]
	curr.degree<-0
	
	orbit.a<-c(orbit.a, curr.iteration[["a"]])
	orbit.b<-c(orbit.b, curr.iteration[["b"]])
	orbit.sign<-c(orbit.sign, -1)
	orbit.degree<-c(orbit.degree, curr.degree)
	
	while (!IsEquiv(x=initial.frac.part[["x"]],
	                y=FractionalPartRational(Prod(curr.iteration, base.neg_unit), norm.b)[["x"]],
	                m=norm.b)){
	  curr.iteration<-FractionalPartRational(Prod(curr.iteration, base.neg_unit), norm.b)[["x"]]
	  curr.degree<-curr.degree + 1
	  orbit.a<-c(orbit.a, curr.iteration[["a"]])
	  orbit.b<-c(orbit.b, curr.iteration[["b"]])
	  orbit.sign<-c(orbit.sign, -1)
	  orbit.degree<-c(orbit.degree, curr.degree)
	  
	}
	orbit.a<-orbit.a/norm.b
	orbit.b<-orbit.b/norm.b
	
	df.res<-data.frame(a=c(orbit.a), b=c(orbit.b), theta=a[["theta"]],
			sign=c(orbit.sign), degree=c(orbit.degree))
#	df.res<-unique(df.res)
	return(df.res)
}

GammaConst <- function(theta){
  res<-max(Abs(Unit(theta)), 1/Abs(Unit(theta)))
	return(res)
}

Dist <- function(x.a, x.b, x.theta, lower.a, upper.a, lower.b, upper.b){
	df.value<-expand.grid(a=seq(from=lower.a, to=upper.a, by=1), 
			b=seq(from=lower.b, to=upper.b, by=1))
	
	x<-list()
	x[["a"]]<-x.a
	x[["b"]]<-x.b
	x[["theta"]]<-x.theta
	
	y<-list()
	y[["a"]]<-df.value$a
	y[["b"]]<-df.value$b
	y[["theta"]]<-x.theta
	
	df.value$dist<-NormIdeal(sum(x, Inverse(y)))
	df.value<-df.value[df.value$dist==min(df.value$dist),]
	return(list(nearest.a=df.value$a[1], nearest.b=df.value$b[1], dist=df.value$dist[1]))
}

FindNearest <- function(a, b){
	eps<-10^{-6}
	orbit.a<-GetOrbit(a, b)
	max.dist<-sqrt(GammaConst(theta=a[["theta"]]))
	if (a[["theta"]] %% 4 != 1) {
  	orbit.a$x<-orbit.a$a+sqrt(a[["theta"]])*orbit.a$b
  	orbit.a$y<-orbit.a$a-sqrt(a[["theta"]])*orbit.a$b
	} else {
	  orbit.a$x<-orbit.a$a+(1+sqrt(a[["theta"]]))*orbit.a$b/2
	  orbit.a$y<-orbit.a$a+(1-sqrt(a[["theta"]]))*orbit.a$b/2
	}
	
	orbit.a$lower.x<-orbit.a$x-max.dist
	orbit.a$upper.x<-orbit.a$x+max.dist
	orbit.a$lower.y<-orbit.a$y-max.dist
	orbit.a$upper.y<-orbit.a$y+max.dist

	if (a[["theta"]] %% 4 != 1) {
  	orbit.a$lower.a<-ceiling((orbit.a$lower.x+orbit.a$lower.y)/2)
  	orbit.a$upper.a<-floor((orbit.a$upper.x+orbit.a$upper.y)/2)
  	orbit.a$lower.b<-ceiling((orbit.a$lower.x-orbit.a$upper.y)/2/sqrt(a[["theta"]]))
  	orbit.a$upper.b<-floor((orbit.a$upper.x-orbit.a$lower.y)/2/sqrt(a[["theta"]]))
	} else {
	  orbit.a$lower.a<-ceiling(
	    (orbit.a$lower.x+orbit.a$lower.y - (orbit.a$lower.x-orbit.a$lower.y)/sqrt(a[["theta"]]))/2)
	  orbit.a$upper.a<-floor(
	    (orbit.a$upper.x+orbit.a$upper.y - (orbit.a$upper.x-orbit.a$upper.y)/sqrt(a[["theta"]]))/2)
	  orbit.a$lower.b<-ceiling((orbit.a$lower.x-orbit.a$upper.y)/sqrt(a[["theta"]]))
	  orbit.a$upper.b<-floor((orbit.a$upper.x-orbit.a$lower.y)/sqrt(a[["theta"]]))
	}	
	FindNearestOne <- function(x){
		Z<-Dist(x.a=orbit.a$a[x], x.b=orbit.a$b[x], x.theta=a[["theta"]],
				lower.a=orbit.a$lower.a[x], upper.a=orbit.a$upper.a[x],
				lower.b=orbit.a$lower.b[x], upper.b=orbit.a$upper.b[x])
		return(Z)
	}
	
	QQ<-lapply(c(1:length(orbit.a[,1])), FUN=FindNearestOne)
	QQ<-rbindlist(QQ)
	orbit.a<-cbind(orbit.a, QQ)

	orbit.a<-orbit.a[orbit.a$dist<min(orbit.a$dist)+eps,]
	orbit.a<-orbit.a[orbit.a$degree==min(orbit.a$degree),]
	
	z<-list()
	z[["a"]]<-orbit.a$a[1]
	z[["b"]]<-orbit.a$b[1]
	z[["theta"]]<-a[["theta"]]
	
	Z<-list()
	Z[["a"]]<-orbit.a$nearest.a[1]
	Z[["b"]]<-orbit.a$nearest.b[1]
	Z[["theta"]]<-a[["theta"]]
	
	n<-orbit.a$degree[1]
	sgn<-orbit.a$sign[1]
	nearest.dist<-orbit.a$dist[1]
	
	return(list(z=z, Z=Z, nearest.dist=nearest.dist, n=n, sgn=sgn))
}

intToBin <- function(x){
	if (x == 1)
		1
	else if (x == 0)
		NULL
	else {
		mod <- x %% 2
		c(intToBin((x-mod) %/% 2), mod)
	}
}

Round <- function(a){
	res<-list()
	res[["a"]]<-round(a[["a"]])
	res[["b"]]<-round(a[["b"]])
	res[["theta"]]<-round(a[["theta"]])
	return(res)
}

DivMod <- function(a, b){
	base.unit<-Unit(a[["theta"]])
	a.b<-Div(a, b)
	res<-FindNearest(a=a, b=b)
	diff.Z<-sum(res[["Z"]], Inverse(res[["z"]]))
	n.abs<-abs(res[["n"]])
	for (i in 1:n.abs){
		if (n.abs>0){
			if (res[["n"]]<0){
				diff.Z<-Prod(diff.Z, base.unit)
			}else{
				diff.Z<-Div(diff.Z, base.unit)
			}
		}
	}
	if (res[["sgn"]]==-1){
		diff.Z<-Inverse(diff.Z)
	}
	q<-Round(sum(a.b, diff.Z))
	r<-sum(a, Inverse(Prod(b, q)))
	return(list(q=q, r=r))
}

Euclidean <- function(a, b){
	res<-list()
	eps<-10^{-6}
	Z<-DivMod(a,b)
	steps<-1
	a.new<-b
	b.new<-Z[["r"]]
	res[["1"]]<-a
	res[["2"]]<-b
	res[["3"]]<-b.new
	while (NormIdeal(b.new)>eps){
		Z<-DivMod(a.new, b.new)
		a.new<-b.new
		b.new<-Z[["r"]]
		steps<-steps+1
		res[[as.character(steps+2)]]<-b.new
	}
	res<-rbindlist(res)
	return(list(steps=steps, res=res))
}

###################
global_theta<-19
found = FALSE

for ( i in 70:100 ) {
  if (found) {
    break
  }
  for ( j in 20:49 ) {
    if ( i %% j != 0 ) {
      print(c(i, j))
      a<-list()
      a[["a"]]<-i
      a[["b"]]<-0
      a[["theta"]]<-global_theta
    
      b<-list()
      b[["a"]]<-j
      b[["b"]]<-0
      b[["theta"]]<-global_theta
      
      Z<-Euclidean(a=a, b=b)
      
      c<-list()
      c[["a"]]<-b[["a"]]
      c[["b"]]<-0
      c[["theta"]]<-global_theta

      d<-list()
      d[["a"]]<-a[["a"]] %% b[["a"]]
      d[["b"]]<-0
      d[["theta"]]<-global_theta
      
      Z1<-Euclidean(a=c, b=d)

      e<-list()
      e[["a"]]<-b[["a"]]
      e[["b"]]<-0
      e[["theta"]]<-global_theta
      
      f<-list()
      f[["a"]]<-a[["a"]] %% b[["a"]] - b[["a"]]
      f[["b"]]<-0
      f[["theta"]]<-global_theta
      
      Z2<-Euclidean(a=e, b=f)
      
      print(c("z z1", Z[["steps"]], Z1[["steps"]]))
      if (Z[["steps"]] - Z1[["steps"]] > 1) {
        print("found")
        print(Z)
        print(Z1)
        found = TRUE
        break
      }
      
      print(c("z z2", Z[["steps"]], Z2[["steps"]]))
      if (Z[["steps"]] - Z2[["steps"]] > 1) {
        print("found")
        print(Z)
        print(Z2)
        found = TRUE
        break
      }
    }
  }
}