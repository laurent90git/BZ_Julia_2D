using DifferentialEquations, BenchmarkTools
using Plots
using Logging
using TerminalLoggers, ProgressBars
include("cusIO.jl")

function laplacian!(u, lap, ncol, nlin,dx,dy)
	# Laplacian for constant coefficient diffusion
    #Flux nul
    lap[2:ncol-1,2:nlin-1] = (u[1:ncol-2,2:nlin-1]+u[3:ncol,2:nlin-1]-2*u[2:ncol-1,2:nlin-1])/dx^2 +
                            +(u[2:ncol-1,1:nlin-2]+u[2:ncol-1,3:nlin]-2*u[2:ncol-1,2:nlin-1])/dy^2
    
    lap[2:ncol-1,1] =        (u[1:ncol-2,1]+u[3:ncol,1]-2*u[2:ncol-1,1])/dx^2 +
                            +(u[2:ncol-1,2]-u[2:ncol-1,1])/dy^2

    lap[2:ncol-1,nlin] =    (u[1:ncol-2,nlin]+u[3:ncol,nlin]-2*u[2:ncol-1,nlin])/dx^2 +
                            +(u[2:ncol-1,nlin-1]-u[2:ncol-1,nlin])/dy^2
    
    lap[1,2:nlin-1]     = (u[2,2:nlin-1]-u[1,2:nlin-1])/dx^2 +
                            +(u[1,1:nlin-2]+u[1,3:nlin]-2*u[1,2:nlin-1])/dy^2
    
    lap[ncol,2:nlin-1]     = (u[ncol-1,2:nlin-1]-u[ncol,2:nlin-1])/dx^2 +
                            +(u[ncol,1:nlin-2]+u[ncol,3:nlin]-2*u[ncol,2:nlin-1])/dy^2
    
    lap[1,1] = (u[2,1]-u[1,1])/dx^2 + (u[1,2]-u[1,1])/dy^2
    lap[1,nlin] = (u[2,nlin]-u[1,nlin])/dx^2 + (u[1,nlin-1]-u[1,nlin])/dy^2
    lap[ncol,1] = (u[ncol-1,1]-u[ncol,1])/dx^2 + (u[ncol,2]-u[ncol,1])/dy^2
    lap[ncol,nlin] = (u[ncol-1,nlin]-u[ncol,nlin])/dx^2 + (u[ncol,nlin-1]-u[ncol,nlin])/dy^2
end

function ODEheat!(du,u,p,t)
		# Heat equation in 2D
        nvar , nlin , ncol, dx, dy, lap = p
        
        dv = reshape(du,nvar,ncol,nlin)
        v = reshape(u, nvar, ncol, nlin)
        T = view(v,1,:,:)
        dT = view(dv,1,:,:)
        
        laplacian!(T, lap, ncol, nlin,dx,dy)
        @. dT = lap
end

function ODEbz!(du,u,p,t)
		# complete BZ system (3 variables) in 2D
        nvar , nlin , ncol, dx, dy, Da, Db, Dc, q, f, eps, mu, lapa, lapb, lapc = p
        
        dv = reshape(du,nvar,ncol,nlin)
        v = reshape(u, nvar, ncol, nlin)
        
        a = view(v,1,:,:)
        da = view(dv,1,:,:)
        b = view(v,2,:,:)
        db = view(dv,2,:,:)
        c = view(v,3,:,:)
        dc = view(dv,3,:,:)
        
        laplacian!(a, lapa, ncol, nlin,dx,dy)
        laplacian!(b, lapb, ncol, nlin,dx,dy)
        laplacian!(c, lapc, ncol, nlin,dx,dy)
    
        @. da = Da*lapa + 1/mu * (-q*a - a*b + f*c)
        @. db = Db*lapb + 1/eps * (q*a - a*b + b*(1-b))
        @. dc = Dc*lapc + b - c 
end

function bz_3eq_init_sol!(u,p,ymin,xmin)
	# Initialise BZ solution
    nvar , nlin , ncol, dx, dy, Da, Db, Dc, q, f, eps, mu, lapa, lapb, lapc = p
    
    for iny in 1:nlin
      yi = ymin + (iny-1)*dy - 0.5
      for inx in  1:ncol
        xi = xmin + (inx-1)*dx
  
        inxy = inx + (iny-1)*ncol
        irow = 1 + (inxy-1)*nvar
        
        if (0.3*xi>=yi && yi>=0 && xi>=0)
          u[irow+1] = 0.8        
        else
          u[irow+1] = q*(f+1)/(f-1)
        end
  
        if (xi>0)
          if (yi>=0)
            u[irow+2] = q*(f+1)/(f-1) + (atan(yi/xi))/(8*pi*f)
          else
            u[irow+2] = q*(f+1)/(f-1) + (atan(yi/xi)+2*pi)/(8*pi*f)
          end
        else 
          if (xi<0)
            u[irow+2] = q*(f+1)/(f-1) + (atan(yi/xi)+pi)/(8*pi*f)
          else 
            if (yi>=0)
              u[irow+2] = q*(f+1)/(f-1)+1/16/f
            else
              u[irow+2]= q*(f+1)/(f-1)+3/16/f
            end
          end
        end
        u[irow]=(f*u[irow+2])/(q+u[irow+1])
      end
    end
end

function compute_timestep()
    return 0.0001
end

function solve_strang(u,tspan, dt ,p, int_diff, int_reac)
    tspan_fake = (tspan[1], tspan[2]*2) # pour être sûr que l'itnégration se poursuit --> peut-être problématique si la RHS dépend de t...
    prob_reac = ODEProblem(ODEbz_reac!,u,tspan_fake,p)
    prob_diff = ODEProblem(ODEbz_diff!,u,tspan_fake,p)
    
    solver_reac = init(prob_reac,int_reac, save_on = true, save_everystep = false, save_end = false, calck = true,
                       abstol = 1e-6, reltol = 1e-6)
    solver_diff = init(prob_diff,int_diff, save_on = false, save_everystep = false, save_end = false, calck = true,
                       abstol = 1e-6, reltol = 1e-6)    
    
    
    
    soltime = [copy(tspan[1])]
    solhist=[copy(u)]
    t=tspan[1]
    dt0 = dt
    while(t<tspan[2]*0.999999)
        #dt = compute_timestep()

        dt = min(dt0, tspan[2]-t)
        println(t, " ",dt)
        step!(solver_reac,dt*0.5,true)
        set_u!(solver_diff ,solver_reac.u )
        step!(solver_diff,dt,true)
        set_u!( solver_reac, solver_diff.u )
        step!(solver_reac,dt*0.5,true)
        savevalues!(solver_reac, true)
        t=t+dt
        push!(soltime, t)
        push!(solhist, copy(solver_reac.u))
    end
    return (soltime, solhist)
end

# Define mesh
nvar = 3
nlin = 60
ncol = 61

ymin = -0.5
ymax = 0.5
xmin = 0.
xmax = 1.

dx = (xmax-xmin)/(ncol-1)
dy = (ymax-ymin)/(nlin-1)

# Define model parameters
Da = 2.5e-3
Db = 2.5e-3
Dc = 1.5e-3
q = 2.e-4
f = 3
eps = 1e-2
mu = 1e-5

# Affect memory
u = zeros(nvar,ncol,nlin)
u = reshape(u,:)
lapa = zeros(ncol,nlin)
lapb = zeros(ncol,nlin)
lapc = zeros(ncol,nlin)
p = [ nvar , nlin , ncol, dx, dy, Da, Db, Dc, q, f, eps, mu, lapa, lapb, lapc ]

# Initialise solution
bz_3eq_init_sol!(u,p,xmin,ymin)

# Set up time integration
tspan = (0.0,1.)
prob = ODEProblem(ODEbz!,u,tspan,p)

# Solve the problem in a monolithic manner
sol=solve(prob,ROCK4(), progress = true, progress_steps = 10)
# export the solution (only works with monolithic right now :\)
export_sol_2dcart("save_test_singlestep.txt", sol, (ncol,nlin,dx,dy,xmin,ymin,xmax,ymax))

# Or solve the problem in a splitted manner
# soltime, solhist = solve_strang(u,tspan, 0.01 ,p, ROCK4(), ROCK4())


