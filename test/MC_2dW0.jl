using MCIntegration
using ElectronGas

const k0, p0 = 1.000000000017037, 1.0000000000422837
param = LegendreInteraction.Parameter.defaultUnit(1e-4, 1.0, dim = 2)
const e0 = param.e0
# Define the integrand 
function integrand(config)
    #config.var is a tuple of variable types specified in the second argument of `MCIntegration.Configuration(...)`
    X = config.var[1]
    x = cos(X[1])
    q2 = (k0 - p0)^2 + 2k0 * p0 * (1 - x)
    if x == 1 || q2 <= 0
        q2 = (k0 - p0)^2 + k0 * p0 * (X[1])^2
    end

    if config.curr == 1 #config.curr is the index of the currently sampled integral by MC
        # return 4π / sqrt(q2)
        return 4π * e0^2 / sqrt(q2)
    end
end

# Define how to measure the observable
function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    weight = integrand(config)
    config.observable[config.curr] += weight / abs(weight) * factor #note that config.observable is an array with two elements as discussed below
end

# MC step of each block
const blockStep = 1e7

# Define the types variables, the first argument sets the range, the second argument gives the largest change to the variable in one MC update. see the section [variable](#variable) for more details.
T = MCIntegration.Continuous([0.0, π], 0.5)

# Define how many (degrees of freedom) variables of each type. 
# For example, [[n1, n2], [m1, m2], ...] means the first integral involves n1 varibales of type 1, and n2 variables of type2, while the second integral involves m1 variables of type 1 and m2 variables of type 2. 
dof = [[1],]

# Define the container for the observable. It must be a number or an array-like object. In this case, the observable has two elements, corresponds to the results for the two integrals. 
obs = [0.0,]

# Define the configuration struct which is container of all kinds of internal data for MC,
# the second argument is a tuple listing all types of variables, one then specify the degrees of freedom of each variable type in the third argument.  
config = MCIntegration.Configuration(blockStep, (T,), dof, obs)

# perform MC integration. Nblock is the number of independent blocks to estimate the error bar. In MPI mode, the blocks will be sent to different workers. Set "print=n" to control the level of information to print.
avg, err = MCIntegration.sample(config, integrand, measure; Nblock = 64, print = 1)

#avg, err are the same shape as obs. In MPI mode, only the root node return meaningful estimates. All other workers simply return nothing
if isnothing(avg) == false
    println("bare kernel: $(avg[1]) +- $(err[1]) (exact: )")
end