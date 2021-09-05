#Problem 1 
#a) 

generate_sigma = function(sd, ex){
  sigma_sq = log((sd/ex)^2 + 1)
  sigma_sq
} 

sds = c(1.0, 3.0, 5.0)

sigmas_sq = rep(0, length(sds))

for (i in 1:length(sds)){
  sigmas_sq[i] = generate_sigma(sds[i], 2)
} 
sigmas_sq

sigmas = sqrt(sigmas_sq)
sigmas

#Now we want to calculate the xi's: 
xi = log(2) - sigmas_sq/2
xi  

#plotting of the densities
x<- seq(0,10,length = 100)
a <- dlnorm(x, meanlog = xi[1], sdlog = sigmas[1] , log = FALSE)
b <- dlnorm(x, meanlog = xi[2], sdlog = sigmas[2], log = FALSE)
g <- dlnorm(x, meanlog = xi[3], sdlog = sigmas[3], log = FALSE)

matplot(x, cbind(a,b,g), type = "l", ylab = "density", main = "log-normal",
        col = 1:3, lty = 1:3)

legend("topright",
       legend = c("xi_1= 0.470 , sigma_1 =0.473 ", "xi_2 = -0.486, sigma_2 = 1.086", "xi_3 = -1.289, sigma_3 = 1.407"),
       col = 1:3,
       lty = 1:3)





#b) 
"""
J = 1000; mu = 0.01; T=1
lambda = J*mu*T 
m = 10000 
N = rpois(m,lambda)

X = 1:m*0 # 10.000 zeros 
a1 = xi[1]; b1 = sigmas[1]
for (i in 1:m){
  Z = rlnorm(N[i], a1, b1); 
  X[i] = sum(Z)
}
X
sort(X)[0.95*m]
sort(X)[0.99*m] 
"""
reserve = function(xi, sigma, lambda, m){
  N = rpois(m, lambda) 
  X = 1:m*0 #generate m zeros
  for (i in 1:m){
    Z = rlnorm(N[i], xi, sigma)
    X[i] =sum(Z)
  }
  lower_reserve = sort(X)[0.95*m] #95% reserve 
  upper_reserve = sort(X)[0.99*m] #99% reserve
  cbind(lower_reserve, upper_reserve)
}

J = 1000; mu = 0.01; T=1
lambda = J*mu*T 

for (i in 1:length(sigmas)){
  print(reserve(xi[i], sigmas[i], lambda, 10000))
}


#c)

reserve_d = function(xi, sigma, lambda, a, b, m){
  N = rpois(m, lambda) 
  X = 1:m*0 #generate m zeros
  for (i in 1:m){
    Z = rlnorm(N[i], xi, sigma) #the claim sizes 
    H = pmin(pmax(Z-a, 0), b)   #the claim sizes reduced according to contract
    X[i] =sum(H)
  }
  lower_reserve = sort(X)[0.95*m] #95% reserve 
  upper_reserve = sort(X)[0.99*m] #99% reserve
  cbind(lower_reserve, upper_reserve)
}

J = 1000; mu = 0.01; T=1
lambda = J*mu*T 

for (i in 1:3){
  print(reserve_d(xi[i], sigmas[i], lambda, 0.5, 3.0, 10000))
}


####Problem 2: 


#a) Computing price of put option

price_put = function(sigma, r, rg, v0, T=1){
  a = (log(1+rg) - r*T + (sigma**2)*(T/2))/sigma*sqrt(T) 
  price = ((1+rg)*exp(-r*T)*pnorm(a) - pnorm(a-sigma*sqrt(T)))*v0
  price
} 

price_put(sigma = 0.25, r = 0.04, rg = 0.06, v0=1, T=1) 



#b) Monte Carlo simulations for put prices: 

monte_carlo_price = function(sigma, r, rg, v0, T=1, m){
  eps = rnorm(m); #drawing m normal variables.
  R = exp(r*T-sigma**2*(T/2) + sigma*sqrt(T)*eps) - 1
  price = exp(-r*T)*mean(pmax(rg-R,0))*v0
  price
}  
sigma_p = c(0.25, 0.30, 0.35)

for (i in 1:length(sigma_p)){
  print(monte_carlo_price(sigma_p[i],r= 0.04, rg=0.06, v0=1,T= 1, m=100000))
}
monte_carlo_price(0.25, 0.04, 0.06, 1, 1, 100000) 

#c) 

monte_carlo_corr = function(sigma1, sigma2, r, rg, v0, T=1, rho, m){
  eps1 = rnorm(m); 
  R1 = exp(r*T-sigma1**2*(T/2) + sigma1*sqrt(T)*eps1) - 1 #first asset simulated. 
  
  eta2 = rnorm(m) #drawing random number
  eps2 = rho*eps1 + sqrt(1-rho^(2))*eta2  
  R2 = exp(r*T-sigma2**2*(T/2) + sigma2*sqrt(T)*eps2) - 1 
  
  R = 0.5*(R1 + R2) # the weights are equal 
  price = exp(-r*T)*mean(pmax(rg-R,0))*v0 
  price
}
rho = c(-0.9, -0.5, 0.0, 0.5, 0.9) 

monte_carlo_corr(sigma1 = 0.25, sigma2 = 0.35, r=0.04, rg=0.06, v0=1, T=1, rho[1], 10000)

for(i in 1:length(rho)){
  print(monte_carlo_corr(sigma1 = 0.25, sigma2 = 0.35, r=0.04, rg=0.06, v0=1, T=1, rho[i], 10000))
}

