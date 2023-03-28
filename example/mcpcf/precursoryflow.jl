module MCPCF

include("propagators.jl")
using .Propagators

end


if abspath(PROGRAM_FILE) == @__FILE__
    using .MCPCF

end

