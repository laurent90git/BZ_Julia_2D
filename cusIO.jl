function export_sol_2dcart(filename, sol, grid)
	# Export de la solution 2D au cours du temps
    io = open(filename, "w");
    nx, ny, dx, dy, xmin, ymin, xmax, ymax = grid 

    nvar = size(sol,1)/nx/ny
    @assert nvar == Int(nvar)
    nvar = Int(nvar)

    delim = " "
    for it in 1:size(sol.t,1)
        Y = reshape(sol[:,it], nvar, nx, ny)
        for iny in 1:ny
            yi = ymin + (iny-1) * dy
            for inx in 1:nx
                xi = xmin + (inx-1) *dx
                out = string(sol.t[it])*delim*string(xi)*delim*string(yi)
                for ivar in 1:nvar
                    out=out*delim*string(Y[ivar,inx,iny])
                end
                println(io, out)
            end
        end
    end
    close(io)
end
