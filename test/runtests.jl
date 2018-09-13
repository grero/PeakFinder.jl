import PeakFinder
using Random
using Test

function test_random_grouping()
    RNG = MersenneTwister(1234)
    ngroups, PP = PeakFinder.check_random_groups(250,20;RNG=RNG)
    @test ngroups == [2,3,4,5,6]
    @test PP[1] ≈ 0.9291383625312836
    @test PP[2] ≈ 0.9955666785841973
    @test PP[3] ≈ 0.9997139792634966
    @test PP[4] ≈ 0.9999284948158742
    @test PP[5] ≈ 1.0
    println("Random grouping test passed")
end

function test_contiguous()
	sidx = [1,2,3, 6,7, 10, 12,13,14]
	counts = PeakFinder.get_contiguous(sidx)
	@test counts[2] == 1 && counts[3] == 2
	println("Contigous bins test passsed")
end

function test_sigint()
	X = zeros(50)
	X[21:30] .= 1.0
	intervals = PeakFinder.get_intervals(X,0.0, 10)
	@test intervals == Dict((21 => 10))
	fill!(X, 0.0)
	X[41:50] .= 1.0
	intervals = PeakFinder.get_intervals(X,0.0, 10)
	@test intervals == Dict((41 => 10))
	println("Sigint test passed")
end

function test_peaks()
	X = zeros(50)
	x = collect(-5:4)
	X[21:30] = exp.(-x.^2/2) #gaussian peak
	peaks = PeakFinder.get_peaks(X, 0.0, 5)
	@test peaks[1] == PeakFinder.Peak(21, 10, 1.0, 2.506624530883954, 26.0)
	#test edge case using time point vector
	fill!(X, 0)
	X[41:50] = exp.(-x.^2/2) #gaussian peak
	t = cumsum(fill(5.0, 50))
	peaks = PeakFinder.get_peaks(X, t, 0.0, 5)
	@test peaks[1] == PeakFinder.Peak(205.0,50.0,1.0,2.506624530883954,230.0)
	println("Peak finding test passed")
    X = zeros(50,2)
	X[21:30,1] = exp.(-x.^2/2) #gaussian peak
	X[41:50,2] = exp.(-x.^2/2) #gaussian peak
	peaks,cellidx = PeakFinder.get_peaks(X, t, 0.0, 5)
	@test peaks[1] == PeakFinder.Peak(105, 50.0, 1.0, 2.506624530883954, 130.0)
    @test cellidx[1] == 1
	@test peaks[2] == PeakFinder.Peak(205.0,50.0,1.0,2.506624530883954,230.0)
    @test cellidx[2] == 2
    println("Finding peaks for multiple cells test passed")
    srand(1234)
	peaks,cellidx,minnbins = PeakFinder.get_peaks(X, t, 0.0, :optimum, 0.01)
	@test peaks[1] == PeakFinder.Peak(105, 50.0, 1.0, 2.506624530883954, 130.0)
    @test cellidx[1] == 1
    @test minnbins[1] == 4
	@test peaks[2] == PeakFinder.Peak(205.0,50.0,1.0,2.506624530883954,230.0)
    @test cellidx[2] == 2
    @test minnbins[2] == 4
    println("Finding peaks for multiple cells with optimum bins test passed")
end

function test_peak_overlaps()
	p1 = PeakFinder.Peak(5,5,1.0,2.0,7)
	p2 = PeakFinder.Peak(6,5,1.0,2.0,10)
	@test PeakFinder.overlaps(p1,p2)
	println("Peak overlap test passed")
end

function test_group_peaks()
    peaks = PeakFinder.Peak[]
	push!(peaks, PeakFinder.Peak(4.0, 5.0, 4.0, 10.0, 6.0))
	push!(peaks, PeakFinder.Peak(4.0, 6.0, 4.0, 8.0, 6.0))
	push!(peaks, PeakFinder.Peak(10, 4.0, 4.0, 12.0, 12.0))
	push!(peaks, PeakFinder.Peak(11, 3.0, 7.0, 13.0, 13.0))
	push!(peaks, PeakFinder.Peak(4.0, 6.0, 4.0, 3.0, 6.0))
	newpeaks = PeakFinder.group_peaks(peaks)
	@test newpeaks[1] == peaks[end-1] && newpeaks[2] == peaks[1]
	println("Peak grouping test passed")
end

#TODO: Add test for finding peaks using cluster analysis

test_contiguous()
test_sigint()
test_random_grouping()
test_peaks()
test_peak_overlaps()
test_group_peaks()
