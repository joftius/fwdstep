source('sim_gen_data.R')

d = run_simulation(nsim=10, n=50, groups=1:10, k=3, upper=3, lower=2, max.steps=8)

groups = sort(rep(1:10, 5))
dg = run_simulation(nsim=10, n=50, groups=groups, k=3, upper=3, lower=2, max.steps=8)

gl = run_simulation(nsim=10, type="glint", n=100, groups=1:10, k=3, upper=3, lower=2, max.steps=8)

gs = run_simulation(nsim=10, type="gamsel", design="uniform", n=100, groups =1:10, k=3, upper=3, lower=2, max.steps=8)
