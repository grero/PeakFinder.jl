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

function get_contiguous!(counts::Dict{Int64,Int64},sidx::AbstractArray{Int64,1})
	nsig = length(sidx)
	nn = 0
	for j in 2:nsig
		if sidx[j] - sidx[j-1] == 1 #find neighburing bins
			nn += 1
		elseif nn >= 1
			nn += 1
			#println("Found group of length $(nn) starting at $(j-nn)")
			counts[nn] = get(counts, nn, 0) + 1
			nn = 0
		else
			nn = 0
		end
	end
	#handle edge
	if nn >= 1
		nn += 1
		#println("Found group of length $(nn) starting at $(nsig-nn)")
		counts[nn] = get(counts, nn, 0) + 1
	end
end

function get_contiguous(sidx::AbstractArray{Int64,1})
	counts = Dict{Int64,Int64}()
	get_contiguous!(counts, sidx)
	counts
end


function get_peaks{T<:Real}(X::Array{T,1}, timepts::Array{Float64,1}, limit::T=0.0, minnbins::Int64=5)
	intervals = get_intervals(X,limit,minnbins)
	peaks = Array(Peak,length(intervals))
	i = 1
	dt = diff(timepts)
	for (k,v) in intervals
		_x = X[k:k-1+v]
		(mx,idx) = findmax(_x)
		_area = sum(_x)
		peaks[i] = Peak(timepts[k],timepts[k+v-1]-timepts[k]+dt[k],mx,_area,timepts[k+idx-1])
		i +=1
	end
	return peaks
end

get_peaks(X, limit, minnbins) = get_peaks(X, [1.0:length(X)], limit,minnbins)

function get_peaks{T<:Real}(X::Array{T,2}, timepts::Array{Float64,1}, limit::T=0.0, minnbins::Union(Int64,Symbol)=5,pvalue::Float64=0.01)
    nbins, ncells = size(X)
    peaks = Array(Peak,0)
    cellidx = Array(Int64,0)
    nsigbins = Array(Int64,0)
    for i in 1:ncells
        if minnbins == :optimum
            nsig = sum(X[:,i].>limit)
            if nsig == 0
                continue #skip this cell
            end
            ngroups,PP = PeakFinder.check_random_groups(nbins, nsig)
            nn = ngroups[findfirst(1-PP.<pvalue)]
        elseif isa(minnbins, Symbol)
            nn = 5
        else
            nn = minnbins
        end
        _peaks = get_peaks(X[:,i], timepts, limit, nn)
        if !isempty(_peaks)
            append!(peaks, _peaks)
            append!(cellidx, fill(i,length(_peaks)))
            append!(nsigbins, fill(nn, length(_peaks)))
        end
    end
    if minnbins == :optimum
        return peaks, cellidx, nsigbins
    else
        return peaks, cellidx
    end
end

function group_peaks(peaks::Array{Peak,1})
	speaks = reverse(sort(peaks))
	newpeaks = Array(Peak,0)
	push!(newpeaks, speaks[1])
	i = 2
	while i <= length(speaks)
		addpeak = true
		j = 1
		while j <= length(newpeaks)
			if overlaps(speaks[i], newpeaks[j])
				addpeak = false
				break
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

function check_random_groups(nbins::Int64, nsig::Int64,nruns::Int64=10000)
    @assert nsig <= nbins
    counts = Dict{Int64, Int64}()
    idx = [1:nbins]
    sidx = zeros(Int64,nsig)
    for i in 1:nruns
        shuffle!(idx)
        for j in 1:nsig
            @inbounds _idx = idx[j]
            @inbounds sidx[j] = _idx
        end
        sort!(sidx)
        get_contiguous!(counts, sidx)
    end
    sidx = sortperm(collect(keys(counts)))
    collect(keys(counts))[sidx], cumsum(collect(values(counts))[sidx])./sum(values(counts))
end

end #module
