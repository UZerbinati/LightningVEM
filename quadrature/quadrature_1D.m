function [int] = quadrature_1D(e_quad, edge, f)

    int = 0;
    
    for i = 1:size(e_quad.x,1)
        for j = 1:size(e_quad.x,2)
            int = int + e_quad.w(j)*f(e_quad.x(i,j),e_quad.y(i,j))*edge(i)/2;
        end
    end
    
end