
library(DPQ)

#
p = 3
(mu = plogis(p))

# log(prob)
log(mu)
p - logspace.add(0, p)

# log(1 - prob)
log(1-mu)
0 - logspace.add(0, p)
