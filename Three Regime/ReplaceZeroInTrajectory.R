
# Replace zero in the trajectory with the average of neighbors

replace.zero <- function(trajectory){
  index.zero <- which(trajectory == 0)
  for (ind in index.zero){
    if (ind == 1){
      trajectory[ind] <- trajectory[ind+1]/2
    }else if(ind == length(trajectory)){
      trajectory[ind] <- trajectory[ind-1]/2
    }else{
      trajectory[ind] <- (trajectory[ind-1] + trajectory[ind+1])/2
    }
  }
  return(trajectory)
}

  