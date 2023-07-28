function [x] = sort_boundary_edges_right(coords, x)

    i = 1;
        
    while (i < numel(x))
            
        if (coords(x(i),2) > coords(x(i+1),2))

            temp = x(i);
            x(i) = x(i+1);
            x(i+1) = temp;
            
            if (i == 1)
                i = 0;
            else
                i = i - 2;
            end
            
        end

        i = i + 1;
        
    end
    
return
end