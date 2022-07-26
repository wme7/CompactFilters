[![View Compact Filters with Holes on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://fr.mathworks.com/matlabcentral/fileexchange/115515-compact-filters-with-holes)

# Explicit Selective Filters
An implementation of the explicit Taylor-based filters presented by Bogey and Bailly (2004).
Using Matlab, these filters are here implemented as matrix operators suitable for using along with explicit/compact finite-difference schemes. 

|  Test  | True / Noised / Filtered f(x) | Error w.r.t. True f(x) |
|:------:|:-------------------------------:|:---------------------------:|
|   1D   | ![Test1D_compact](.figure/../figures/Test_compactFilters1d.png) | ![Test1D_Error](.figure/../figures/Test_compactFilters1d_error.png) |


|  Test  |  Direct Explicit Filter |  Explicit Filter Operator |
|:------:|:-----------------------:|:-------------------------:|
|   2D   | ![Test2D_explicit](.figure/../figures/Test_explicitFilters2d.png) | ![Test2D_compact](.figure/../figures/Test_compactFilters2d.png) |

## References
- Bogey, C., & Bailly, C. (2004). A family of low dispersive and low dissipative explicit schemes for flow and noise computations. Journal of Computational physics, 194(1), 194-214.
