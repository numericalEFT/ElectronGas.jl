"""
Template of convention.

    Stores conventions that will not change for all purposes.
"""
module Convention

# using StaticArrays

# const Weight = SVector{2,Float64}
# const Base.abs(w::Weight) = abs(w[1]) + abs(w[2]) # define abs(Weight)
const INL, OUTL, INR, OUTR = 1, 2, 3, 4
const EPS = 1e-16

# export all conventions
for n in names(@__MODULE__; all=true)
	  if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
		    @eval export $n
	  end
end

end
