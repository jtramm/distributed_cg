using Plots
matrix = readdlm("dense.out")
heatmap(matrix, aspect_ratio=1.0)
savefig("picture.png")
