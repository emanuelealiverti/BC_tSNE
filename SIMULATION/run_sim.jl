println("Loading module\n")
include("../tsne_julia/bctsne.jl")
using .TSne
using DelimitedFiles 

perplexity = 30.0
println("Reading the data\n")
X = readdlm(ARGS[1],Float64,header=false);
Z = readdlm(ARGS[2],Float64,header=false);

iterations = 5000 
using Random
Random.seed!(1)
println("Starting tsne\n")
println("Might take a while in time and memory \n")

Y = tsne(X,Z, 2, 30, iterations,pca_init=false,grad_corr=true)
println("saving")
open(pwd()*"/Y30C.txt", "w") do io
	writedlm(io, Y)
end

Y = tsne(X,Z, 2, 50, iterations,pca_init=false,grad_corr=true)
println("saving")
open(pwd()*"/Y50C.txt", "w") do io
	writedlm(io, Y)
end

Y = tsne(X,Z, 2, 10, iterations,pca_init=false,grad_corr=true)
println("saving")
open(pwd()*"/Y10C.txt", "w") do io
	writedlm(io, Y)
end


Y = tsne(X,Z, 2, 30, iterations,pca_init=false,og_corr=false, grad_corr=false)
println("saving")
open(pwd()*"/Y30U.txt", "w") do io
	writedlm(io, Y)
end

Y = tsne(X,Z, 2, 50, iterations,pca_init=false,og_corr=false, grad_corr=false)
println("saving")
open(pwd()*"/Y50U.txt", "w") do io
	writedlm(io, Y)
end

Y = tsne(X,Z, 2, 10, iterations,pca_init=false,grad_corr=false, og_corr=false)
println("saving")
open(pwd()*"/Y10U.txt", "w") do io
	writedlm(io, Y)
end

