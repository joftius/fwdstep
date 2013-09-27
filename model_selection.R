
source('stop_rules.R')

# Input incl.list, p.list, true.active.set
# Use stopping rules stop_first_large, stop_forward, stop_hybrid
# (those input p.list and alpha, output index of last inclusion)

# from fwd_step_test.R results
# names = null.p, signal.p, active.set, true.active

# write a function taking  these inputs and calculating:
# Rob/Alex/Max/Stefan show FDR, V, R, and average power (fraction of non-null rejected)
