function [f1, f2, u1, u1_x, u1_y, u2, u2_x, u2_y] = problem_test_elasticity(k, matProps)

% problem_test_lighting: this function defines the problem that we want to
% solve
%
% Input parameters: 
% k       : The index of the problem that we want to solve;
% matProps: This struct contains the constants of the equation.
%
% Output parameters:
% f      : right-hand side of the problem;
% g      : boundary condition;
% u      : solution of the equation;
% deru_x : x-derivative of u;
% deru_y : y-derivative of u.

switch(k)
        
    case(1)
                
        u1 = @(x,y)  0.25 * (sin(2*pi*x)).^2 .* sin(4*pi*y) + log(1 + x.*y);
        u2 = @(x,y) -0.25 * (sin(2*pi*y)).^2 .* sin(4*pi*x);
                
        u1_x = @(x,y)  pi *  sin(4*pi*x)     .* sin(4*pi*y) / 2 + y ./ (1 + x.*y);
        u1_y = @(x,y)  pi * (sin(2*pi*x)).^2 .* cos(4*pi*y) + x ./ (1 + x.*y);
        
        u1_xx = @(x,y)   2*pi^2 *  cos(4*pi*x)     .* sin(4*pi*y) - y.^2 ./ ((1 + x.*y).^2) ;
        u1_xy = @(x,y)   2*pi^2 *  sin(4*pi*x)     .* cos(4*pi*y) -   ((1 + x.*y) - y.*x) ./ ((1 + x.*y).^2);
        u1_yy = @(x,y) - 4*pi^2 * (sin(2*pi*x)).^2 .* sin(4*pi*y) - x.^2 ./ ((1 + x.*y).^2);

        u2_x = @(x,y)  - pi * (sin(2*pi*y)).^2 .* cos(4*pi*x);
        u2_y = @(x,y)  - pi *  sin(4*pi*x)     .* sin(4*pi*y) / 2;
        
        u2_xx = @(x,y) + 4*pi^2 * (sin(2*pi*y)).^2 .* sin(4*pi*x);
        u2_xy = @(x,y) - 2*pi^2 *  sin(4*pi*y) .* cos(4*pi*x);
        u2_yy = @(x,y) - 2*pi^2 *  cos(4*pi*y) .* sin(4*pi*x);

    case(2)
                
        u1 = @(x,y) x.^2;
        u2 = @(x,y) y.^2;
                
        u1_x = @(x,y)  2.*x + 0.*y;
        u1_y = @(x,y)  0.*x + 0.*y;
        
        u1_xx = @(x,y)   2 + 0.*x + 0.*y;
        u1_xy = @(x,y)   0 + 0.*x + 0.*y;
        u1_yy = @(x,y)   0 + 0.*x + 0.*y;

        u2_x = @(x,y)  0 + 0.*x + 0.*y;
        u2_y = @(x,y)      0.*x + 2.*y;
        
        u2_xx = @(x,y) 0 + 0.*x + 0.*y;
        u2_xy = @(x,y) 0 + 0.*x + 0.*y;
        u2_yy = @(x,y) 2 + 0.*x + 0.*y;

case(3)
                
        u1 = @(x,y) x;
        u2 = @(x,y) y;
                
        u1_x = @(x,y)  1 + 0.*y;
        u1_y = @(x,y)  0 + 0.*y;
        
        u1_xx = @(x,y)   0 + 0.*x + 0.*y;
        u1_xy = @(x,y)   0 + 0.*x + 0.*y;
        u1_yy = @(x,y)   0 + 0.*x + 0.*y;

        u2_x = @(x,y)  0 + 0.*x + 0.*y;
        u2_y = @(x,y)  1 + 0.*x + 0.*y;
        
        u2_xx = @(x,y) 0 + 0.*x + 0.*y;
        u2_xy = @(x,y) 0 + 0.*x + 0.*y;
        u2_yy = @(x,y) 0 + 0.*x + 0.*y;

    case(4)
                
        u1 = @(x,y)  0.25 * (sin(2*pi*x)).^2 .* sin(4*pi*y);
        u2 = @(x,y) -0.25 * (sin(2*pi*y)).^2 .* sin(4*pi*x);
                
        u1_x = @(x,y)  pi *  sin(4*pi*x)     .* sin(4*pi*y) / 2;
        u1_y = @(x,y)  pi * (sin(2*pi*x)).^2 .* cos(4*pi*y) ;
        
        u1_xx = @(x,y)   2*pi^2 *  cos(4*pi*x)     .* sin(4*pi*y);
        u1_xy = @(x,y)   2*pi^2 *  sin(4*pi*x)     .* cos(4*pi*y);
        u1_yy = @(x,y) - 4*pi^2 * (sin(2*pi*x)).^2 .* sin(4*pi*y);

        u2_x = @(x,y)  - pi * (sin(2*pi*y)).^2 .* cos(4*pi*x);
        u2_y = @(x,y)  - pi *  sin(4*pi*x)     .* sin(4*pi*y) / 2;
        
        u2_xx = @(x,y) + 4*pi^2 * (sin(2*pi*y)).^2 .* sin(4*pi*x);
        u2_xy = @(x,y) - 2*pi^2 *  sin(4*pi*y) .* cos(4*pi*x);
        u2_yy = @(x,y) - 2*pi^2 *  cos(4*pi*y) .* sin(4*pi*x);
                

end
    
f1 = @(x,y) -2 * matProps.mu * (u1_xx(x,y) + 0.5 * (u1_yy(x,y) + u2_xy(x,y))) - matProps.lambda * (u1_xx(x,y) + u2_xy(x,y));
f2 = @(x,y) -2 * matProps.mu * (u2_yy(x,y) + 0.5 * (u2_xx(x,y) + u1_xy(x,y))) - matProps.lambda * (u2_yy(x,y) + u1_xy(x,y));
   

end