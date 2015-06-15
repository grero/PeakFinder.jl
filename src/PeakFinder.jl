module PeakFinder
using Docile
@docstrings

Docile.@doc meta("Find the set of contiguous intervals for which `X` exceeds `limits` for at least `minnbins` bins", return_type=Dict{Int64,Int64})->
function get_intervals(X::Array{Float64,1},limit::Real=3,minnbins::Integer=5)
    nbins = length(X)
    nsig = 0
    intervals = Dict{Int64,Int64}()
    for i=1:nbins
        x = X[i]
        if !isfinite(x)
            continue
        end
        if x > limit
            nsig += 1
        else
            if nsig >= minnbins
                intervals[i-nsig] = nsig
            end
            nsig = 0
        end
    end
	#make sure we pick up the last interval as well
	if nsig >= minnbins && X[end] > limit
		intervals[nbins-nsig+1] = nsig
	end
    return intervals
end

end #module
