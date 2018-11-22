function foxWolframMoments(particles)
    entries = length(particles)
    h0 = 0.0
    hd = 0.0
    # init
    h10 = 0.0
    h20 = 0.0
    h30 = 0.0
    h40 = 0.0
    
    if entries < 2
        println(stderr, "foxWolframMoments: Not enough entries")
        return [-1., -1., -1., -1.]
    end 
    particleArray = Vector{SVector{4, Float64}}(undef, entries)
	# copy all particles into an array and calculate h0 and hd
    @inbounds for i in 1:entries
        p = particles[i]
	    momentum = getMomentum(p)
	    magsquared = sum(momentum.^2)
        particleArray[i] = SVector{4, Float64}(momentum[1], momentum[2], momentum[3], sqrt(magsquared))
	    h0 += particleArray[i][4]
	    hd += magsquared
	end
	h0 = h0*h0 # square it
	# now do the moments
	@inbounds for i = 1:entries-1
	    @inbounds for j = i+1:entries
    		ptemp = particleArray[i][4]*particleArray[j][4]
	    	cosθ = (particleArray[i][1]*particleArray[j][1] +
	    	        particleArray[i][2]*particleArray[j][2] +
                    particleArray[i][3]*particleArray[j][3]
                   ) / ptemp
            h10 += ptemp*cosθ
		    h20 += ptemp*(1.5 * cosθ^2 - 0.5)
		    h30 += ptemp*(2.5 * cosθ^3 - 1.5*cosθ)
		    h40 += ptemp*(4.375*cosθ^4 - 3.75*cosθ^2 + 0.375)
        end
    end
	# normalize the moments
	h10 = (hd + 2*h10) / h0
	h20 = (hd + 2*h20) / h0
	h30 = (hd + 2*h30) / h0
    h40 = (hd + 2*h40) / h0
    return [h10, h20, h30, h40] 
end
