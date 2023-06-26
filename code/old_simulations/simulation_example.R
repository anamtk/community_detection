#https://stackoverflow.com/questions/56010141/generating-random-0-1-matrix-with-fixed-features-in-r

M = 10 # rows
N = 10 # columns
m = 6 #ones are only in m of the M rows
n = 4 #ones are only in n of the N columns

ni = sample(1:N, n)
mi = sample(1:M, m)

expand.grid(N = 1:N, M = 1:M) %>% 
  mutate(value = ifelse(N %in% ni & M %in% mi, 1, 0)) %>% 
  .$value %>% 
  matrix(., nrow = M, byrow = TRUE)



comm <- c(0,0,0,1,0,1,1,1)
nsteps <- 10
thresh <- 0.2

matrix <- matrix(data = NA, nrow = nsteps, ncol = length(comm))
matrix[1,] <- comm

for(i in 1:nsteps){
  randdraw <- runif(length(comm))
  trans <- ifelse(randdraw > thresh, 1, 0)
  comm <- ifelse(trans == 1, abs(comm-1), comm)
  matrix[i,] <- comm
   }
