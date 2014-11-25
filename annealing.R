require(ggmap)
require(mapproj) # for us.cities data
require(Imap) # for gdist
require(Matrix) # for sparse
require(maps)
library(geosphere)
require(googleVis)

#install.packages('googleVis')
#install.packages("geosphere")


data(us.cities)
head(us.cities)

n = nrow(us.cities) #1005
memo = Matrix(0,n,n)
city_dist = function (i,j){
  if(i==j){return(0)}  
  
  # enforce i<j
  if(i>j){
    t=i
    i=j
    j=t
  }  
  
  d=memo[i,j]
  if(d>0){return(d)}    
  d = gdist(us.cities[i,'long'], us.cities[i,'lat'], 
        us.cities[j,'long'], us.cities[j,'lat'])  
  #memo[i,j] = d  
  #memo[i,j] <- d  
  memo[i,j] <<- d  # I don't like this global assignment. I just want to assign to the next scope
  return(d)
}

path_dist = function(path){
  trips = embed(path, 2)
  sum(apply(trips, 1, function(row)city_dist(row[1], row[2])))
}

# proposal function to generate a "neighbor" path for MCMC walk of solution space
propose_path = function(last_path, k_swaps=1){
  path=last_path
  swaps = k_swaps*2
  ix = sample(length(last_path), swaps, replace=FALSE)
  if(k_swaps==1){
    path[ix[1]] = last_path[ix[2]]
    path[ix[2]] = last_path[ix[1]]    
  } else {
    m=matrix(ix, ncol=2)
    for(i in 1:k_swaps){
      path[m[i,1]] = last_path[m[i,2]]
      path[m[i,2]] = last_path[m[i,1]]
    }    
  }
  return(path)
}

annealing = function(x0,          
                     func=function(x){-path_dist(x)}, #SA will find a maxima, we're looking for a minima 
                     prop=propose_path,
                     cooling=.003,
                     iters=1e4, # This is really relative to the cooling rate
                     temp = 10000,
                     mintemp=1e-6,
                     #maxreject=1e3, # unimplemented. Lag change stop condition. Requires some historical tracking (could just use a deque)
                     history=FALSE,
                     verbose=TRUE
          ){
  hx=NULL
  dx=NULL
  if(history){
    hx=matrix(NA, iters, length(x0))
    dx=rep(NA, iters)
  }
  last = x0
  d = func(last)
  U = runif(iters)
  for(i in 1:iters){
    if(verbose & i%%100==0){print(c(i,d, temp))}
    next_state = prop(last)
    delta = func(next_state) - d
    R = exp(delta/temp)
    #last = ifelse(U[i]<R, last, next_state)
    if(U[i]<R){
      last=next_state
      d = func(last)
    }
    temp = temp*(1-cooling)
    if(history){
      hx[i,]=last
      dx[i]=d
    }
    if(temp<mintemp){
      print(i)
      if(history){
        hx=hx[1:i,]
        dx=dx[1:i]
      }
      break
      }
    }
  list(solution=last, history=hx, values=dx)
}

plotGPS <- function(path) {
  n = length(path)
  lat_long = paste(us.cities[path,]$lat, us.cities[path,]$long, sep=':')
  print(lat_long)
  nn <- list(LatLong = lat_long 
             #Tip = path$timestamp[1:v]
             ,Tip = seq_along(path)
  )
  
  m <- gvisMap(nn, 'LatLong' , 'Tip',
               options=list(showTip=TRUE, 
                            showLine=TRUE,
                            lineColor='#800000',
                            enableScrollWheel=TRUE,
                            mapType='hybrid', useMapTypeControl=TRUE,
                            width=800,height=400),
               chartid="TravelingSalesman")
  
  plot(m)
}


plot_path = function(path){    
  trips = embed(path, 2)
  print(trips)
  inter <- gcIntermediate(us.cities[trips[,2],c('long','lat')], us.cities[trips[,1],c('long','lat')], n=50, addStartEnd=TRUE)
  map('usa', col="#f3f3f3", fill=TRUE, border=FALSE)
  map('state', col="#999999", add=TRUE)
  map.cities(us.cities[path,], pch=19, cex=1.1, label=FALSE)
  lapply(inter, lines)  
}


x0=1:30
test = annealing(x0)
path_dist(x0)
path_dist(test$solution)
plot(test$values, type='l')
plot_path(test$solution) # not very special, need a slower cooling schedule
#plotGPS(test$solution)

xv=1:1e5
yv=1000*(1-.001)^xv
#plot(xv, yv)
yv[3e4] # stabilizes at temp = 0.081949
yv[4700] # stabilizes at temp = 0.081949

test.002 = annealing(x0, cooling = .002) # slower cooling schedule. Better, not great.
test = annealing(x0, cooling = .001, iters=2e4) # getting better, still not great. Halted early.
plot_path(test.002$solution)

# Let's try modifying the proposal slightly
test.001 = annealing(x0, cooling = .001, iters=3e4, prop=function(p){propose_path(p, k_swaps=2)}) #  8805.397
# just barely, but actually still improving up to iteration 22200. Maybe I'm setting the mintemp too high?
# also... this solution sucks. I'm calling it at stuck in local minima. Cooled too slow?

