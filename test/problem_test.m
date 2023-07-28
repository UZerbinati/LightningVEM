function [f, g, u, deru_x, deru_y] = problem_test(k, matProps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: problem_test
%
% Created by : M. Trezzi
%
%---------------------------------------------------------------------------------------------------
% Purpose
% =======
% This function returns the function related to the equation that we want
% to solve
%
% Input
% =====
% k        : The index of the problem that we want to solve
% matProps : This struct contains the constants of the equation
%
% Output
% ======
% f      : Right-hand side of the problem
% g      : Boundary condition
% u      : Solution of the equation
% deru_x : x-derivative of u
% deru_y : y-derivative of u
%      
%---------------------------------------------------------------------------------------------------
% Function's updates history
% ==========================
% Mai  7, 2022: first realease (by M. Trezzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

switch(k)
        
    case(1)
                
        u = @(x,y) x + y;
                
        deru_x = @(x,y) 1 + 0.*x + 0.*y;
        deru_y = @(x,y) 1 + 0.*x + 0.*y;
                
        u_xx = @(x,y) 0 + 0.*x + 0.*y;
        u_yy = @(x,y) 0 + 0.*x + 0.*y;
                
    case(2)
                
        u = @(x,y) sin(pi.*x).*sin(pi.*y);
                 
        deru_x = @(x,y) pi.*cos(pi.*x).*sin(pi.*y);
        deru_y = @(x,y) pi.*cos(pi.*y).*sin(pi.*x);
                
        u_xx = @(x,y) -pi.*pi.*sin(pi.*x).*sin(pi.*y);
        u_yy = @(x,y) -pi.*pi.*sin(pi.*x).*sin(pi.*y);

    case(3)
                
        u = @(x,y) exp(x+3*y);
                 
        deru_x = @(x,y) exp(x+3*y);
        deru_y = @(x,y) 3*exp(x+3*y);
                
        u_xx = @(x,y) exp(x+3*y);
        u_yy = @(x,y) 9*exp(x+3*y);
        
    case(4)
                
        u = @(x,y) x.^3 + 0.*y;
                 
        deru_x = @(x,y) 3*x.^2;
        deru_y = @(x,y) 0 + 0.*x + 0.*y;
                
        u_xx = @(x,y) 6.*x;
        u_yy = @(x,y) 0 + 0.*x + 0.*y;
    
    case(5)
                
        u = @(x,y) x.^2 + 0.*y;
                 
        deru_x = @(x,y) 2*x;
        deru_y = @(x,y) 0 + 0.*x + 0.*y;
                
        u_xx = @(x,y) 2;
        u_yy = @(x,y) 0 + 0.*x + 0.*y;

    case(6)
                
        u = @(x,y) 20 + (x.^2 - 10*cos(pi*x*20)) + (y.^2 - 10*cos(pi*y*20));
                 
        deru_x = @(x,y) 2*x + 200*pi*sin(pi*20*x) + 0.*y;
        deru_y = @(x,y) 2*y + 200*pi*sin(pi*20*y) + 0.*x;
                
        u_xx = @(x,y) 2 + 4000*pi*pi*cos(pi*20*x) + 0.*y;
        u_yy = @(x,y) 2 + 4000*pi*pi*cos(pi*20*y) + 0.*x;

    case(7)
                
        u = @(x,y) 0 + 1.*x + 0.*y;
                
        deru_x = @(x,y) 1 + 0.*x + 0.*y;
        deru_y = @(x,y) 0 + 0.*x + 0.*y;
                
        u_xx = @(x,y) 0 + 0.*x + 0.*y;
        u_yy = @(x,y) 0 + 0.*x + 0.*y;

end
    
f = @(x,y) - matProps.epsilon      .* (u_xx(x,y) + u_yy(x,y)) ...
           + matProps.sigma        .* u(x,y);
   
g = @(x,y) u(x,y);

end