library(GA)
library(rgl)

#### exemplo - otimizacao de uma variavel
f <- function(x) (x^2 + x) * cos(x)

lbound <- -10
ubound <- 10

curve(f, lbound, ubound, n = 1000)

GA <- ga(type = "real-valued", fitness = f, lower = c(th = lbound), upper = ubound)
summary(GA)
plot(GA)

curve(f, lbound, ubound, n = 1000)
points(GA@solution, GA@fitnessValue, col = 2, pch = 19)

### exemplo 2 - otimizacao de duas variaveis

Rastrigin <- function(w) {
  20 + w[1]^2 + w[2]^2 - 10 * ( cos(2 * pi * w[1]) + cos(2 * pi * w[2]) )
}

Rastrigin_wrapper <- function(x1, x2) {
  Rastrigin(c(x1, x2))
}

x1 <- x2 <- seq(-5.12, 5.12, by = 0.1)

f <- outer(x1, x2, Vectorize(Rastrigin_wrapper))

persp3D(x1, x2, f, theta = 50, phi = 20, col.palette = bl2gr.colors)

filled.contour(x1, x2, f, color.palette = bl2gr.colors)

GA <- ga(type = "real-valued", 
         fitness = function(x) -Rastrigin(x),
         lower = rep(-5.12, 2), upper = rep(5.12, 2),
         popSize = 50, maxiter = 1000, run = 100, seed = 42)

summary(GA)

filled.contour(x1, x2, f, color.palette = bl2gr.colors, 
               plot.axes = { axis(1); axis(2); 
                 points(GA@solution[,1], GA@solution[,2], 
                        pch = 3, cex = 2, col = "white", lwd = 2) }
)

### exemplo 3 - otimizacao de duas variaveis com constraints

x1 <- x2 <- seq(0, 1, by = .05)
fitness <- function(w) {
  pen <- sqrt(.Machine$double.xmax)
  penaltyGt1 <- max(sum(w) - 1, 0) * pen
  penaltyLt1 <- max(1 - sum(w), 0) * pen
  -(w[1]^2 + w[2]^2) - penaltyGt1 - penaltyLt1
}

fitness_wrap <- Vectorize(function(x1, x2) {fitness(c(x1, x2))})
f <- outer(x1, x2, fitness_wrap)
filled.contour(x1, x2, f, color.palette = bl2gr.colors)

GA <- ga(type = "real-valued", 
         fitness = function(x) fitness(x),
         lower = rep(0, 2), upper = rep(1, 2),
         popSize = 50, maxiter = 1000, run = 100, seed = 42)

summary(GA)
plot(GA)
filled.contour(x1, x2, f, color.palette = bl2gr.colors, 
               plot.axes = { axis(1); axis(2); 
                 points(GA@solution[,1], GA@solution[,2], 
                        pch = 3, cex = 2, col = "red", lwd = 2) })
