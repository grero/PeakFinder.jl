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
	println("sigint test passed")
end

test_sigint()
