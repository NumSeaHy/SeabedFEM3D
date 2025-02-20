using PyPlot
PyPlot.rc("text", usetex=true)
PyPlot.rc("font", family="serif")

# Define the mesh sizes (are arbitrary, but the point is that each mess size is half the previous one)
h = [0.5, 0.25, 0.125]

# This errors has been obtained from the CSV simulation file of the problem (TODO: Automating this process)
errors = [12.574723200826659, 4.0777691513143415, 1.247497041253719]

# Ellaborate the figure
figure()
loglog(h, errors, marker="o", linestyle="-", color="k")
# Create the O(h^2) triangle
# Plot the triangle for slope 2
plt.plot([0.25, 0.5, 0.5, 0.25],
         [2, 2, 7.5, 2],
         color="k", linestyle="-")

# Annotate the triangle with "O(h^2)"
plt.text(0.4, 3, L"$O(h^2)$", fontsize=12)
xlabel(L"$h$", fontsize=12)
ylabel(L"$e_h[\%]$", fontsize=12)
grid(true, which="both", linestyle="--", linewidth=0.7)  # Add grid lines
tight_layout()
savefig("./images/convergenceplot_marron_triangle.pdf", transparent=true)
display(gcf())
