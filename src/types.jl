import JLD
import Base.==, Base.isless, Base.intersect, Base.in
abstract AbstractPeak

type Peak <: AbstractPeak
	start::Float64
	length::Float64
	height::Float64
	area::Float64
	pidx::Float64
end

function ==(p1::AbstractPeak, p2::AbstractPeak)
	return p1.start == p2.start && p1.length == p2.length && p1.height == p2.height && p1.pidx == p2.pidx && p1.area == p2.area
end

function isless(p1::AbstractPeak, p2::AbstractPeak)
	p1.area < p2.area
end

function overlaps(p1::AbstractPeak, p2::AbstractPeak)
	return (p1.start <= p2.start <= p1.start + p1.length-1) || (p2.start <= p1.start <= p2.start+p2.length-1)
end

type RandomContiguousDistribution
    n::Int64
    k::Int64
    P::Dict{Int64,Float64}
end

function RandomContiguousDistribution(n::Int64,k::Int64)
    #check for cached results
    _path = joinpath((homedir()),".cached")
    if !ispath(_path)
        mkpath(_path)
    end
    fname = joinpath(_path,"randomcontiguous.jd")
    kk = "$n-$k"
    if isfile(fname)
        jff = JLD.jldopen(fname,"r")
        if kk in names(jff)
            pp = read(jff, kk)
            close(jff)
            return RandomContiguousDistribution(n,k,pp)
        end
        close(jff)
        fmode = "r+"
    else
        fmode = "w"
    end
    ngroups,_pp = check_random_groups(n,k,100000)
    PP = Dict(ngroups, _pp)
    #append the results to the cache
    JLD.jldopen(fname,fmode) do ff
        write(ff, kk, PP)
    end
    RandomContiguousDistribution(n,k,PP)
end

function lowerlimit(rcg::RandomContiguousDistribution,pvalue::Float64=0.05)
    _qq = collect(keys(rcg.P))
    sidx = sortperm(_qq)
    #_qq = _qq[sidx
    _pp = collect(values(rcg.P))[sidx]
    _xi = findfirst(1-_pp .< 0.05)
    if _xi != 0
        _minnbins = _qq[_xi]
    else
        _minnbins = 2
    end
end




