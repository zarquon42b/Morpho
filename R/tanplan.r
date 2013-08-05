tanplan <- function(x)
{		if (sum(x^2)==0)
		{stop(cat("zero vector has no orthogonal subspace"))}
		
		
		
		if (0 %in% x)
			{
			y <- c(0,0,0)	
			y[which(x==0)] <- 1
                        y <- y/sqrt(sum(y^2))
			}
		else 
			{
			y <- c(1,1,-(x[1]+x[2])/x[3])
			y <- y/sqrt(sum(y^2))
			}
		z <- crossp(y,x)
		z <- z/sqrt(sum(z^2))
	return(list(z=z,y=y))
}
			
		
