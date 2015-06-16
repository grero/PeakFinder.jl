type Peak{T<:Real}
	start::Int64
	length::Int64
	height::T
	area::T
	pidx::Int64
end

function ==(p1::Peak, p2::Peak)
	return p1.start == p2.start && p1.length == p2.length && p1.height == p2.height && p1.pidx == p1.pidx && p1.area == p2.area
end

function isless(p1::Peak, p2::Peak)
	p1.area < p2.area
end



