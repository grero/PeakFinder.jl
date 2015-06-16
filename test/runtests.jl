import PeakFinder
using Base.Test

function test_sigint()
	X = zeros(50)
	X[21:30] = 1.0
	intervals = PeakFinder.get_intervals(X,0.0, 10)
	@test intervals == {21 => 10}
	fill!(X, 0.0)
	X[41:50] = 1.0
	intervals = PeakFinder.get_intervals(X,0.0, 10)
	@test intervals == {41 => 10}
	println("Sigint test passed")
end

function test_peaks()
	X = zeros(50)
	x = [-5:4]
	X[21:30] = exp(-x.^2/2) #gaussian peak
	peaks = PeakFinder.get_peaks(X, 0.0, 5)
	@test peaks[1] == PeakFinder.Peak(21, 10, 1.0, 2.506624530883954, 6)
	println("Peak finding test passed")
end

function test_peak_overlaps()
	p1 = PeakFinder.Peak(5,5,1.0,2.0,7)
	p2 = PeakFinder.Peak(6,5,1.0,2.0,10)
	@test PeakFinder.overlaps(p1,p2)
	println("Peak overlap test passed")
end

test_sigint()
test_peaks()
test_peak_overlaps()
