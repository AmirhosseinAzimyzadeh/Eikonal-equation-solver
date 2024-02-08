using DelimitedFiles
using Plots


# Solves the Eikonal equation using the Fast Sweeping Method
function solveEikonal(ni, nj, h, source, tol = 1E-6)
    F = zeros(ni+2, nj+2) # Speed function
    U = fill(Inf, ni+2, nj+2) # Solution matrix
    iterations = 0

    initialize!(F, U, ni, nj, source) # Read input from file

    maxI = ni+1
    maxJ = nj+1

    maxerr = tol+1
    while maxerr > tol
        maxErrors = fill(0.0, 4)

        maxErrors[1] = sweep!(U, maxI, 2, 2, maxJ, F, h, maxerr)
        maxErrors[2] = sweep!(U, maxI, 2, maxJ, 2, F, h, maxerr)
        maxErrors[3] = sweep!(U, 2, maxI, maxJ, 2, F, h, maxerr)
        maxErrors[4] = sweep!(U, 2, maxI, 2, maxJ, F, h, maxerr)
        
        iterations += 1

        maxerr = maximum(maxErrors)
    end

    result = @view U[2:end-1, 2:end-1]

    println("Number of iterations: ", iterations, ", Max error: ", maxerr)
    # Uncomment the following line to show the result
    # @show result

    return result
end

# Perform one set of sweeps on matrix U using directions specified in 
# ia,ib,ja,jb, and speed function F, and current value of maximum
# absolute error (err = (U[i, j]-Unew)/U[i, j]), where Unew
# is new value of U[i,j] calculated using solveQuadratic()
# h is the spacing between points
# Returns updated value of maxerr (Max Error)
function sweep!(U, ia, ib, ja, jb, F, h, maxerr)
    stepi = ib < ia ? -1 : 1
    stepj = jb < ja ? -1 : 1

    localMaxerr::Float64 = 0.0

    for j in ja:stepj:jb
        for i in ia:stepi:ib
            Unew = solveQuadratic(U, i, j, F, h)
            if Unew < U[i, j]
                err::Float64 = (U[i, j]-Unew)/U[i, j]
                if err > localMaxerr
                    localMaxerr = err
                end
                U[i, j] = Unew
            end
        end
    end
    return localMaxerr
end

# Solve the discretized Eikonal equation at (i,j), given 
# speed function F, and spacing between points h
# Returns the new value of U[i,j]
function solveQuadratic(U, i, j, F, h)
    if U[i, j] <= 0
        return U[i, j]
    end
    a = min(U[i-1, j], U[i+1, j])
    b = min(U[i, j-1], U[i, j+1])
    if abs(a-b) >= h/F[i, j]
        return min(a, b) + h/F[i, j]
    else
        return (a+b+sqrt(2*((h/F[i, j])^2)-((a-b)^2)))/2
    end
end

# Reads from source a speed function F defined a ni x nj grid 
# as well as coordinates of boundary of front
# (input is terminated with a negative number on last line)
# Also initializes solution matrix U
# U and F are (ni+2) x (nj+2) 2D matrices
# Requires the DelimitedFiles package
function initialize!(F, U, ni, nj, source)
    temp = readdlm(source)
    for j in 1:nj, i in 1:ni
        F[i+1, j+1] = Float64(temp[i, j])
    end
    for i in ni+1:size(temp, 1)-1 # skip last line that has -1 terminator
        # +2: skip border and convert input from 0-based indexing
        U[temp[i,1]+2, temp[i,2]+2] = 0
    end

    nothing
end


width = 100;
filename = "test100_extra.txt"
h = 1
U = solveEikonal(width, width, h, filename)

# Use heatmap(U) "heatmap(U)"
