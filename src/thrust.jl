# a simple implementation avoiding vectors should make it simpler to port to Python and Julia
function CrossProduct(v, i, j, result)
    @inbounds result[1] = v[i][2]*v[j][3]-v[i][3]*v[j][2]
	@inbounds result[2] = v[i][3]*v[j][1]-v[i][1]*v[j][3]
    @inbounds result[3] = v[i][1]*v[j][2]-v[i][2]*v[j][1]
	@inbounds result[4] = 0
    return 0
end

function Thrust(particles)
    entries = length(particles)
    particleArray = Vector{SVector{4, Float64}}(undef, entries)
    PartVector = MVector{4, Float64}(0., 0., 0., 0.)
    RefVector = MVector{4, Float64}(0., 0., 0., 0.)
    SumVector = MVector{4, Float64}(0., 0., 0., 0.)
    maxVal = 0.0
       
    if (entries<2)
        println(stderr, "Thrust : Not enough entries")
        return -1.0
    end
    # copy all particles into an array and start calcualting the thrust
    @inbounds for i in 1:entries
        p = particles[i]
        momentum = getMomentum(p)
        magsquared = sum(momentum.^2)
        particleArray[i] = SVector{4, Float64}(momentum[1], momentum[2], momentum[3], sqrt(magsquared)) # store the abs(p)
        SumVector += particleArray[i]
    end
    
    # now do the thrust    
    @inbounds for i = 1:entries-1
        @inbounds for j = i+1:entries
            CrossProduct(particleArray, i, j, RefVector)
			RefVector ./= sqrt(sum(RefVector.^2))
			PartVector =0;
            # Add all momenta with sign; two choices for each reference particle.
            @inbounds for k = 1:entries 
                if k != i && k != j 
                    dotProduct = particleArray[k] ⋅ RefVector
                    PartVector += sign(dotProduct) * particleArray[k]
                end
            end
            c1 = 0
            c2 = 0
            for l = 1:4 
                if l == 0
                    c1 = 1
                    c2 = 1
                elseif l == 1
                    c1 = 1
                    c2 = -1
                elseif l == 2
                    c1 = -1
                    c2 = 1
                else
                    c1 = -1
                    c2 = -1                       
                end
                
				tempVal = 0.0
                @inbounds for n=1:3 # only copy three of the four as we overwrite the fourth element again ...
                    tempVal += (PartVector[n] + c1*particleArray[i][n] + c2*particleArray[j][n])^2
                end            
                tempVal = sqrt(tempVal)
                if (tempVal > maxVal) 
                    maxVal = tempVal
                end    
            end
        end
    end
    return maxVal / SumVector[4]
end
