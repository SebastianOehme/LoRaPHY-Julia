# Working with the Doppler Est Centroid Data to do linear regression
using GLM
using DataFrames

# x_1 = [1,2,3,4,5,6]
# x_2 = [2,3,4,5,6,15]
# dopplerEstCentroid1 = [2.2,2.4,2.5,2.7,2.8,3.0]
# dopplerEstCentroid2 = [2.4,2.5,2.7,2.8,3.0, 5.0]

# predictionIntervall = 1:1:15

# scatter!(x_2, dopplerEstCentroid2)

# plot!(linearDopplerEstimate)

# linearDopplerEstimate, m, n = linearRegression(x_1,dopplerEstCentroid1, predictionIntervall)
# linearDopplerEstimate, m, n = linearRegression(x_2,dopplerEstCentroid2, predictionIntervall)


# x = 1:15
# y = m .*x .+ n
# y_new = m_new .* x .+ n_new
# plot(linearDopplerEstimate)
# plot!(y_new)



function linearRegression(data_x::StepRange{Int64, Int64}, data_y::Vector{Float64}, predictionIntervall::StepRange{Int64, Int64})

    linearDopplerEstimate = Vector{Float64}[]

    x = data_x
    y = data_y
        
    # Put the data into a DataFrame
    df = DataFrame(x=x, y=y)
    # Fit a linear regression model: y ~ x (y is dependent, x is independent)
    model = lm(@formula(y ~ x), df)
    
    m = coef(model)[2]
    n = coef(model)[1]

    x_new = predictionIntervall
    linearDopplerEstimate = vcat(linearDopplerEstimate,  m .* x_new .+ n)
    
    return linearDopplerEstimate, m, n

end


#plot(dopplerEstSimple[1023:1033])
#plot!(dopplerEstCentroid[1023:1033])
#plot!(linearDopplerEstimate[1023:1033])

#plot(dopplerEstSimple)
#plot!(dopplerEstCentroid)
#plot!(linearDopplerEstimate)