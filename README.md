# unptrack

unptrack is an UNstructured mesh, lagrangian Particle-TRACKing model    
that simulates the transport pathways and dispersion/dispersal of    
pelagic biota or chemical contaminants using flow fields generated      
by unstructured mesh hydrodynamic models. The model runs offline using     
velocity fields defined at cell centres (e.g. FVCOM) or at element nodes.  
Advection can be treated using either a fourth-order Runge-Kutta algorithm    
or a simple Euler approach. A random walk model is used to simulate        
horizontal and vertical eddy diffusion. Various aspects of biological      
development and behaviour can be simulated. For chemical contaminants, a   
decay half-life can be included.                                           
                                                                           
Reference: Gillibrand, P.A. and K.J. Willis, 2007. Dispersal of sea lice   
larvae from salmon farms: A model study of the influence of environmental  
conditions and larval behaviour. Aquatic Biology, 1, 63 - 75.              
                                                                           
This version uses velocity fields from unstructured (flexible) mesh        
hydrodynamic models e.g. FVCOM, RiCOM.                                     
                                                                           
Created by: Philip Gillibrand (philip.gillibrand@mowi.com)                 