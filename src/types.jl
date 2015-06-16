import Base.==, Base.isless

type Peak
	start::Float64
	length::Float64
	height::Float64
	area::Float64
	pidx::Float64
end

function ==(p1::Peak, p2::Peak)
	return p1.start == p2.start && p1.length == p2.length && p1.height == p2.height && p1.pidx == p1.pidx && p1.area == p2.area
end

function isless(p1::Peak, p2::Peak)
	p1.area < p2.area
end

function overlaps(p1::Peak, p2::Peak)
	return (p1.start <= p2.start <= p1.start + p1.length-1) || (p2.start <= p1.start <= p2.start+p2.length-1)
end



