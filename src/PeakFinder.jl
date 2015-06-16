module PeakFinder
using Docile
@docstrings

include("types.jl")

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

function get_peaks{T<:Real}(X::Array{T,1}, timepts::Array{Float64,1}, limit::T=0.0, minnbins::Int64=5)
	intervals = get_intervals(X,limit,minnbins)
	peaks = Array(Peak,length(intervals))
	i = 1
	for (k,v) in intervals
		_x = X[k:k-1+v]
		(mx,idx) = findmax(_x)
		_area = sum(_x)
		peaks[i] = Peak(k,v,mx,_area,idx)
		i +=1
	end
	return peaks
end

get_peaks(X, limit, minnbins) = get_peaks(X, [1.0:length(X)], limit,minnbins)

function group_peaks(peaks::Array{Peak,1})
	speaks = reverse(sort(peaks))
	newpeaks = Array(Peak,0)
	push!(newpeaks, speaks[1])
	i = 2
	while i <= length(speaks)
		addpeak = false
		j = 1
		while j <= length(newpeaks)
			if !overlaps(speaks[i], newpeaks[j])
				addpeak = true
			else
				addpeak = false
			end
			j += 1
		end
		if addpeak
			push!(newpeaks, speaks[i])
		end
		i += 1
	end
	newpeaks
end

end #module
