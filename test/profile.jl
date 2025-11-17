using Profile, PProf, Oscar, BuchbergerLearning

n = 4
R, _ = polynomial_ring(GF(32003), n; internal_ordering=:degrevlex)
K = katsura(R)

# profile cpu time
Profile.clear()
Profile.@profile buchberger(
  K,
  complete_reduction=true,
  elimination_strategy=ProductCriterion(),
  selection_strategy=NormalStrategy(),
)
PProf.pprof()

# # profile allocations
# Profile.Allocs.clear()
# @Profile.Allocs.profile dimension_qa(A, R)
# PProf.Allocs.pprof()
