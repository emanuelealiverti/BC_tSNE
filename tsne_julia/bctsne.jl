"""
BC-TSNE
Batch-corrected t-SNE via adjustment and projected gradient descent
- The underlying implementation of the standard t-SNE  has been forked by https://github.com/lejon/TSne.jl
"""

module TSne
using LinearAlgebra, Statistics, Distances, ProgressMeter, StatsBase, DataFrames, GLM
using Printf: @sprintf

export tsne
function Hbeta!(P::AbstractVector, D::AbstractVector, beta::Number)
	@inbounds P .= exp.(-beta .* D)
	sumP = sum(P)
	@assert (isfinite(sumP) && sumP > 0.0) "Degenerated P: sum=$sumP, beta=$beta"
	H = log(sumP) + beta * dot(D, P) / sumP
	@inbounds P .*= 1/sumP
	return H
end

function perplexities(D::AbstractMatrix{T}, tol::Number = 1e-5, perplexity::Number = 30.0;
		      max_iter::Integer = 50,
		      verbose::Bool=false, progress::Bool=true) where T<:Number
	(issymmetric(D) && all(x -> x >= 0, D)) ||
	throw(ArgumentError("Distance matrix D must be symmetric and positive"))

	# initialize
	n = size(D, 1)
	P = fill(zero(T), n, n) # perplexities matrix
	beta = fill(one(T), n)  # vector of Normal distribution precisions for each point
	Htarget = log(perplexity) # the expected entropy
	Di = fill(zero(T), n)
	Pcol = fill(zero(T), n)

	# Loop over all datapoints
	progress && (pb = Progress(n, "Computing point perplexities"))
	for i in 1:n
		progress && update!(pb, i)

		# Compute the Gaussian kernel and entropy for the current precision
		betai = 1.0
		betamin = 0.0
		betamax = Inf

		copyto!(Di, view(D, :, i))
		Di[i] = prevfloat(Inf) # exclude D[i,i] from minimum(), yet make it finite and exp(-D[i,i])==0.0
		minD = minimum(Di) # distance of i-th point to its closest neighbour
		@inbounds Di .-= minD # entropy is invariant to offsetting Di, which helps to avoid overflow

		H = Hbeta!(Pcol, Di, betai)
		Hdiff = H - Htarget

		# Evaluate whether the perplexity is within tolerance
		tries = 0
		while abs(Hdiff) > tol && tries < max_iter
			# If not, increase or decrease precision
			if Hdiff > 0.0
				betamin = betai
				betai = isfinite(betamax) ? (betai + betamax)/2 : betai*2
			else
				betamax = betai
				betai = (betai + betamin)/2
			end

			# Recompute the values
			H = Hbeta!(Pcol, Di, betai)
			Hdiff = H - Htarget
			tries += 1
		end
		verbose && abs(Hdiff) > tol && warn("P[$i]: perplexity error is above tolerance: $(Hdiff)")
		# Set the final column of P
		@assert Pcol[i] == 0.0 "Diagonal probability P[$i,$i]=$(Pcol[i]) not zero"
		@inbounds P[:, i] .= Pcol
		beta[i] = betai
	end
	progress && finish!(pb)
	# Return final P-matrix
	verbose && @info(@sprintf("Mean σ=%.4f", mean(sqrt.(1 ./ beta))))
	return P, beta
end

"""
pca(X::Matrix, ncols::Integer = 50)
Run PCA on `X` to reduce the number of its dimensions to `ndims`.
"""

function pca(X::AbstractMatrix, ndims::Integer = 50)
	(n, d) = size(X)
	(d <= ndims) && return X
	Y = X .- mean(X, dims=1)
	C = Symmetric((Y' * Y) ./ (n-1))
	Ceig = eigen(C, (d-ndims+1):d) # take eigvects for top ndims largest eigvals
	return Y * reverse(Ceig.vectors, dims=2)
end

"""
OG, remove linear dependenencies in the data
"""
function og(X::AbstractMatrix, 
	    Z::Union{AbstractMatrix,AbstractVector};
	    pca_init::Bool = false, 
	    red_dims::Integer = 50)
	str = "PCA"*ifelse(pca_init," ", " NOT ")*"provided"
	println(str)
	if pca_init
		Y = X
	else
		Y = pca(X, red_dims)
	end
	println("doing OG correction")
	Y_bar = copy(Y)
	for j in 1:size(Y_bar,2)
		tmp = residuals(lm(Z,Y_bar[:,j]))
		Y_bar[:,j] = copy(tmp) 
	end
	return Y_bar
end
# K-L divergence element
kldivel(p, q) = ifelse(p > zero(p) && q > zero(q), p*log(p/q), zero(p))


"""
tsne(X::Union{AbstractMatrix, AbstractVector}, ndims::Integer=2, reduce_dims::Integer=0,
max_iter::Integer=1000, perplexity::Number=30.0; [keyword arguments])
"""
function tsne(X::AbstractMatrix, 
	      Z::Union{AbstractMatrix,AbstractVector},
	      ndims::Integer = 2,
	      reduce_dims::Integer = 30,
	      max_iter::Integer = 2000,
	      perplexity::Number = 30.0;
	      og_corr::Bool = true,
	      grad_corr::Bool = true,
	      distance::Union{Bool, Function, SemiMetric} = false,
	      min_gain::Number = 0.01,
	      eta::Number = 200.0,
	      pca_init::Bool = false,
	      approx_dist::Bool = false, 
	      tol_dist::Number = 0.01,
	      initial_momentum::Number = 0.5,
	      final_momentum::Number = 0.8,
	      momentum_switch_iter::Integer = 250,
	      stop_cheat_iter::Integer = 250,
	      cheat_scale::Number = 12.0,
	      verbose::Bool = false,
	      progress::Bool=true,
	      eps::Number = 0.01,
	      extended_output = false)

	n = size(X, 1)
	Y = randn(n, ndims)


	dY = fill!(similar(Y), 0)     # gradient vector
	iY = fill!(similar(Y), 0)     # momentum vector
	gains = fill!(similar(Y), 1)  # how much momentum is affected by gradient


	if og_corr
		println("Using OG to reduce dimensions, and computing euclidean distance. Might take a while")
		str = "Using "*string(ifelse(pca_init,size(X,2),reduce_dims))*" components"
		println(str)
		X_bar = og(X,Z, pca_init = pca_init,red_dims = reduce_dims)
		X_bar .-= mean(X_bar, dims = 1)
		D = pairwise(Euclidean(),X_bar,X_bar,dims=1)
	else
		#We are not doing correction. We either have PCA as input, or compute them, or trow error
		X_bar = copy(X)
		if !pca_init
			X_bar = pca(X, reduce_dims) 
		end
		if size(X_bar,2) > 100
			println("X has too many columns. Are you sure you don't want to do pca?")
			return NaN
		end
		println("OG not estimated, and computing euclidean distance")
		X_bar .-= mean(X_bar, dims = 1)
		D = pairwise(Euclidean(),X_bar,X_bar,dims=1)
	end

	if approx_dist
		D[D .< tol_dist] .= 0
	end


	P, beta = perplexities(D, 1e-5, perplexity,
			       verbose=verbose, progress=progress)
	P .+= P' # make P symmetric
	P .*= cheat_scale/sum(P) # normalize + early exaggeration
	sum_P = cheat_scale

	# Run iterations
	progress && (pb = Progress(max_iter, "Computing t-SNE"))
	Q = fill!(similar(P), 0)     # temp upper-tri matrix with 1/(1 + (Y[i]-Y[j])²)
	Ymean = similar(Y, 1, ndims) # average for each embedded dimension
	sum_YY = similar(Y, n, 1)    # square norms of embedded points
	L = fill!(similar(P), 0)     # temp upper-tri matrix for KLdiv gradient calculation
	Lcolsums = similar(L, n, 1)  # sum(Symmetric(L), 2)
	last_kldiv = NaN

	# Structure for batch correcting

	old_kldiv = -9999
	for iter in 1:max_iter
		# Compute pairwise affinities
		BLAS.syrk!('U', 'N', 1.0, Y, 0.0, Q) # Q=YY', updates only the upper tri of Q
		@inbounds for i in 1:size(Q, 2)
			sum_YY[i] = Q[i, i]
		end
		sum_Q = 0.0
		@inbounds for j in 1:size(Q, 2)
			sum_YYj_p1 = 1.0 + sum_YY[j]
			Qj = view(Q, :, j)
			Qj[j] = 0.0
			for i in 1:(j-1)
				sqdist_p1 = sum_YYj_p1 - 2.0 * Qj[i] + sum_YY[i]
				@fastmath Qj[i] = ifelse(sqdist_p1 > 1.0, 1.0 / sqdist_p1, 1.0)
				sum_Q += Qj[i]
			end
		end
		sum_Q *= 2 # the diagonal and lower-tri part of Q is zero

		# Compute the gradient
		inv_sum_Q = 1.0 / sum_Q
		fill!(Lcolsums, 0.0) # column sums
		# fill the upper triangle of L (gradient)
		@inbounds for j in 1:size(L, 2)
			Lj = view(L, :, j)
			Pj = view(P, :, j)
			Qj = view(Q, :, j)
			Lsumj = 0.0
			for i in 1:(j-1)
				@fastmath Lj[i] = l = (Pj[i] - Qj[i]*inv_sum_Q) * Qj[i]
				Lcolsums[i] += l
				Lsumj += l
			end
			Lcolsums[j] += Lsumj
		end
		@inbounds for (i, ldiag) in enumerate(Lcolsums)
			L[i, i] = -ldiag
		end
		# dY = -4LY
		BLAS.symm!('L', 'U', -4.0, L, Y, 0.0, dY)

		# Perform the update
		momentum = iter <= momentum_switch_iter ? initial_momentum : final_momentum
		@inbounds for i in eachindex(gains)
			gains[i] = max(ifelse((dY[i] > 0) == (iY[i] > 0),
					      gains[i] * 0.8,
					      gains[i] + 0.2),
				       min_gain)
			iY[i] = momentum * iY[i] - eta * (gains[i] * dY[i])
			Y[i] += iY[i]
		end

		# Projected gradient. Do update and send back
		if grad_corr
			for j in 1:ndims
				tmp = residuals(lm(Z,Y[:,j]))
				Y[:,j] = copy(tmp) 
			end
		end
		Ymean = mean(Y,dims=1) 
		@inbounds Y .-= Ymean 

		# stop cheating with P-values
		if sum_P != 1.0 && iter >= min(max_iter, stop_cheat_iter)
			P .*= 1/sum_P
			sum_P = 1.0
		end
		# Compute the current value of cost function
		if !isfinite(last_kldiv) || iter == max_iter ||
			(progress && mod(iter, max(max_iter÷20, 10)) == 0)
			local kldiv = 0.0
			@inbounds for j in 1:size(P, 2)
				Pj = view(P, :, j)
				Qj = view(Q, :, j)
				kldiv_j = 0.0
				for i in 1:(j-1)
					# P and Q are symmetric (only the upper triangle used)
					@fastmath kldiv_j += kldivel(Pj[i], Qj[i])
				end
				kldiv += 2*kldiv_j + kldivel(Pj[j], Q[j])
			end
			old_kldiv = last_kldiv
			last_kldiv = kldiv/sum_P + log(sum_Q/sum_P) # adjust wrt P and Q scales
		end

		progress && update!(pb, iter)
	end
	progress && (finish!(pb))

	# Return solution
	if !extended_output
		return Y
	else
		return Y, beta, last_kldiv
	end
end

end
