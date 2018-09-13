module PeakFinder

include("types.jl")

"""
Find the set of contiguous intervals for which `X` exceeds `limits` for at least `minnbins` bins.
"""
function get_intervals(X::Array{Float64,1},limits::Array{Float64,1},minnbins::Integer=5)
    nbins = length(X)
    nsig = 0
    intervals = Dict{Int64,Int64}()
    for i=1:nbins
        limit = limits[i]
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
    if nsig >= minnbins && X[end] > limits[end]
		intervals[nbins-nsig+1] = nsig
	end
    return intervals
end

function get_intervals(X::Array{Float64,1},limit::Real=3,minnbins::Integer=5)
    limits = fill(limit,length(X))
    get_intervals(X,limits,minnbins)
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

"""
Find peaks, i.e. periods in which the signal `X` exceeds `limit`. The optimum number of consecutive time points for which `X > limit` is determined by a cluster analysis, using the p-value `pvalue`.

    function get_peaks{T<:Real}(X::Array{T,1}, timepts::Array{Float64,1}, limit::Union{Array{T,1}, T}=0.0, minnbins::Symbol=:optimum,pvalue::Real=0.05)
"""
function get_peaks(X::Array{T,1}, timepts::AbstractVector{Float64}, limit::Union{Array{T,1}, T}=0.0, minnbins::Symbol=:optimum,pvalue::Real=0.05) where T<:Real
    nbins = length(X)
    nsig = sum(X.>limit)
    if nsig == 0
        return Peak[],0
    end
    ngroups,PP = PeakFinder.check_random_groups(nbins, nsig)
    idx = findfirst(1-PP.<pvalue)
    if idx == 0
        return Peak[],0
    end
    nn = ngroups[idx]
    _peaks = get_peaks(X, timepts, limit, nn)
    return _peaks,nn
end

function get_peaks(X::Array{T,1}, timepts::AbstractVector{Float64}, limit::Union{Array{T,1}, T}=0.0, minnbins::Int64=5) where T<:Real
	intervals = get_intervals(X,limit,minnbins)
    peaks = Array{Peak}(undef, length(intervals))
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

get_peaks(X, limit, minnbins) = get_peaks(X, [1.0:length(X);], limit,minnbins)

function get_peaks(X::Array{T,2}, timepts::Array{Float64,1}, limit::Union{Array{T,1}, T}=0.0, minnbins::Union{Int64,Symbol}=5,pvalue::Float64=0.01) where T<:Real
    nbins, ncells = size(X)
    peaks = Peak[]
    cellidx = Int64[]
    nsigbins = Int64[]
    for i in 1:ncells
        if minnbins == :optimum
            nsig = sum(X[:,i].>limit)
            if nsig == 0
                continue #skip this cell
            end
            ngroups,PP = PeakFinder.check_random_groups(nbins, nsig)
            idx = findfirst(1-PP.<pvalue)
            if idx == 0
                continue #ksip this cell
            end
            nn = ngroups[idx]
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

function group_peaks(peaks::Array{T,1}) where T<:AbstractPeak
	speaks = reverse(sort(peaks))
    newpeaks = T[]
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
    nsig <= nbins || ArgumentError("Number of significant bins cannot exceed total number of bins")
    counts = Dict{Int64, Int64}()
    idx = collect(1:nbins)
    sidx = fill(0, nsig)
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
