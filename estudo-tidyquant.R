library(tidyverse)
library(tidyquant)
library(corrr)
library(GA)

Rcpp::sourceCpp("fitness.cpp")


acoes <- c("COGN3", "EQTL3", "SBSP3", "IRBR3", "RADL3",
           "ENBR3", "ITUB4", "CYRE3", "BRDT3", "WEGE3")

qtde_xm <- c(1300, 400, 200, 800, 400, 500, 300, 400, 400, 200)

carteira_xm <- tibble(acoes, qtde_xm) %>%
  mutate(prop = qtde_xm / sum(qtde_xm)) %>%
  mutate(symb = paste0(acoes, ".SA")) %>%
  arrange(desc(prop))

carteira_xm

av_key <- "" # put your alpha vantage api key here
av_api_key(av_key)
av_api_rate_limit <- 5 # at most 5 calls per minute
ohlc <- list()

# Download data from alpha vantage, respecting API rate limits
# for (i in 1:nrow(carteira_xm)) {
#   cat(sprintf("Obtendo dados da ação %s... %d/%d\n", carteira_xm$symb[i], i, nrow(carteira_xm)))
#   ohlc[[i]] <- tq_get(carteira_xm$symb[i], get = "alphavantage",
#                                       av_fun = "TIME_SERIES_DAILY_ADJUSTED",
#                                       outputsize = 'full')
#   Sys.sleep(60 / 5 + 0.1)
# }

# cache the data
#saveRDS(ohlc, file = 'ohlc.RDS')
ohlc <- readRDS('ohlc.RDS')

quotes <- bind_rows(ohlc)

# make sure we dont have missing data on the close price
stopifnot(!any(is.na(quotes$adjusted_close)))

# determine the set of days we will use for the data points, ignoring days we
# have missing observations from any one of the stocks in the pool
dates_range <- quotes %>%
  pivot_wider(id_cols = timestamp, names_from = symbol, values_from = adjusted_close) %>%
  drop_na() %>%
  select(timestamp) %>%
  arrange(timestamp) %>%
  pull

filtered_quotes <- quotes %>%
  filter(timestamp %in% dates_range)

filtered_quotes %>%
  count(symbol)

filtered_quotes %>%
  filter(symbol == 'ITUB4.SA') %>%
  ggplot(aes(timestamp, adjusted_close)) +
  geom_line(col = 'blue') +
  scale_x_date(breaks = "6 months") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), labels = scales::dollar) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('Data') +
  ylab('Fechamento ajustado')

Ra <- filtered_quotes %>%
  group_by(symbol) %>%
  tq_transmute(select = adjusted_close,
               mutate_fun = periodReturn,
               period = 'monthly',
               col_rename = 'Ra') %>%
  ungroup()

Ra %>%
  ggplot(aes(timestamp, Ra)) +
  geom_line(aes(color = symbol), lwd = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), labels = scales::percent) +
  scale_x_date(breaks = "1 months") +
  scale_color_discrete(name = 'Ação') +
  ylab('Retorno mensalizado') +
  xlab('Data')

# mean and std returns: mu_i, sigma_i
returns <- Ra %>%
  group_by(symbol) %>%
  summarise(mu = mean(Ra), sigma = sd(Ra)) %>%
  ungroup()

returns

corrs <- Ra %>%
  pivot_wider(id_cols = timestamp, names_from = symbol, values_from = Ra) %>%
  arrange(timestamp) %>%
  select(-timestamp) %>%
  correlate(diagonal = 1)

idxs <- sapply(corrs$term, function(x) which(returns$symbol == x))

returns <- returns[idxs,]

term_seq <- corrs$term
corrs <- corrs[,2:11]
covs <- corrs
covs[,] <- 0
N <- nrow(covs) # number of assets

for(i in 1:(N-1)) {
  for(j in (i+1):N) {
    #cat(sprintf("i = %d j = %d\n", i, j))
    covs[i,j] <- corrs[i, j] * returns$sigma[i] * returns$sigma[j]
  }
}

covs <- as.matrix(covs)

generate_lambdas <- function(n = 10) {
  seq(0, n - 1, by = 1) / (n - 1)
}

fitness <- function(w) {
  pen <- sqrt(.Machine$double.xmax)
  # penalties for making the sum of allocations equal to 1
  penaltyGt1 <- max(sum(w) - 1, 0) * pen
  penaltyLt1 <- max(1 - sum(w), 0) * pen

  returnFactor <- - (1 - lambda) * sum(w * returns$mu)
  riskFactor <- 0

  for(i in 1:(N-1)) {
    for(j in (i+1):N) {
      #cat(sprintf("i = %d j = %d\n", i, j))
      riskFactor <- riskFactor + w[i] * w[j] * covs[i , j]
    }
  }

  riskFactor <- lambda * as.numeric(riskFactor)

  finalTerm <- riskFactor + returnFactor #+ penaltyGt1 + penaltyLt1

  - finalTerm
}

w <- rep(1/N, N)

lambdas_seq <- generate_lambdas(5)
solutions <- list()
#fitness(w)

for (i in 1:length(lambdas_seq)) {
  lambda <- lambdas_seq[i]

  cat(sprintf("***** Iteração %d / %d ******", i, length(lambdas_seq)))

  GA <- ga(type = "real-valued",
           fitness = function(x) fitnessCpp(x, lambda, returns$mu, covs),
           lower = rep(0, N), upper = rep(1, N),
           popSize = 50, maxiter = 1000, run = 100, seed = 42)

  normalized_solution <- as.numeric(GA@solution)
  normalized_solution <- normalized_solution / sum(normalized_solution)

  solutions[[i]] = list()
  solutions[[i]]$w <- normalized_solution
  solutions[[i]]$lambda <- lambda
  solutions[[i]]$return <- sum(normalized_solution * returns$mu)
}

#saveRDS(solutions, file = "solutions_backtest_lambda_5.RDS")
#saveRDS(solutions, file = "solutions_backtest_lambda_5_with_sum_penalization.RDS")

solutions <- readRDS(file = "solutions_backtest_lambda_5.RDS")

simulation <- tibble(lambda = sapply(solutions, function(x) x$lambda),
                      return = sapply(solutions, function(x) x$return))

simulation %>%
  ggplot(aes(1 - lambda, return)) +
  geom_line(lwd = 2) +
  xlab(expression(Fator~de~Risco~lambda)) +
  ylab('Retorno médio do portfólio') +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(text = element_text(size=20)) +
  ggtitle('Fronteira eficiente de investimento')

# carteira_xm %>%
#   mutate(ga = normalized_solution) %>%
#   arrange(desc(ga))
#
# returns %>% arrange(desc(mu))
