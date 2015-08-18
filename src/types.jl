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





