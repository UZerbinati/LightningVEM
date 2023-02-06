function vem_star_domain(n, perc)

%% PARAMETERS OF THE MESH
h    = 1/n;                                                                                          %Mesh size                                               
h2   = 1/(2*n);                                                                                      %Half of the mesh size
dist = h2 * sqrt(2) * perc / 100;                                                                    %Parameter of the star

%% MEMORY ALLOCATION
vertices = zeros(7*n^2 + 4*n + 1, 2);
boundary = zeros(n*8,1);
elements = cell(5*n*n,1);

%% GENERATE THE FIRST SQUARE LOCATED AT THE BOTTOM LEFT CORNER

vertices(1,:) = [0 0];                                                                               %Coordinates of the vertices
vertices(2,:) = [h 0];
vertices(3,:) = [h h];
vertices(4,:) = [0 h];

vertices(5,:) = [h2 0];                                                                              %Coordinates of the midpoints 
vertices(6,:) = [h h2];
vertices(7,:) = [h2 h];
vertices(8,:) = [0 h2];

xmid = h2;                                                                                           %Center of the square
ymid = h2;

vertices(9,:)  = [xmid - dist, ymid - dist];                                                         %Coordinates of the inner points of the star
vertices(10,:) = [xmid + dist, ymid - dist];
vertices(11,:) = [xmid + dist, ymid + dist];
vertices(12,:) = [xmid - dist, ymid + dist];

boundary(1:5) = [1; 5; 2; 4; 8];                                                                     %Boundary nodes

if (n == 1)                                                                                          %If we choose n=1, we have also this boundary nodes

    boundary(6:8) = [6; 3; 7];

end

elements{1} = [1; 5; 9; 8];                                                                          %Creation of the first five elements
elements{2} = [2; 6; 10; 5];
elements{3} = [3; 7; 11; 6];
elements{4} = [4; 8; 12; 7];
elements{5} = [5; 10; 6; 11; 7; 12; 8; 9];

n_elem  = 5;                                                                                         %Update the number of objects that we have inseterted
n_nodes = 12;
n_bound = 5;

%% CREATION OF THE FIRST ROW
for j = 2:n                                                                                          %From the second square to the last one

    jh  = j*h;                                                                                       %Utility variable
    mid = jh-h2;

    vertices(n_nodes + 1, :) = [jh 0];                                                               %Coordinate of the bottom right node
    vertices(n_nodes + 2, :) = [jh h];                                                               %Coordinate of the top    right node

    vertices(n_nodes + 3, :) = [mid 0];                                                              %Coordinate of the bottom mid node
    vertices(n_nodes + 4, :) = [jh h2];                                                              %Coordinate of the right  mid node
    vertices(n_nodes + 5, :) = [mid h];                                                              %Coordinate of the top    mid node

    vertices(n_nodes + 6, :) = [mid - dist, h2 - dist];                                              %Coordinates of the the star
    vertices(n_nodes + 7, :) = [mid + dist, h2 - dist];
    vertices(n_nodes + 8, :) = [mid + dist, h2 + dist];
    vertices(n_nodes + 9, :) = [mid - dist, h2 + dist];

    if j == 2                                                                                        %Check if we are adding the second square

        elements{n_elem + 1} = [ 2; 15; 18;  6];
        elements{n_elem + 2} = [13; 16; 19; 15];
        elements{n_elem + 3} = [14; 17; 20; 16];
        elements{n_elem + 4} = [ 3;  6; 21; 17];
        elements{n_elem + 5} = [15; 19; 16; 20; 17; 21; 6; 18];

        boundary(6:7) = [15; 13];                                                                    %Add the boundary nodes 

        n_nodes = n_nodes + 9;                                                                       %Update the number of objects that we have inserted
        n_elem  = n_elem  + 5;
        n_bound = n_bound + 2;
        
        if (j == n)
            
            boundary(n_bound + 1 : n_bound + 2) = [16; 14];                                          %If it is the last element of the row, we add 
            n_bound = n_bound + 2;                                                                   %the boundary nodes

        end

    else

        elements{n_elem + 1} = [n_nodes - 8; n_nodes + 3; n_nodes + 6; n_nodes - 5];
        elements{n_elem + 2} = [n_nodes + 1; n_nodes + 4; n_nodes + 7; n_nodes + 3];
        elements{n_elem + 3} = [n_nodes + 2; n_nodes + 5; n_nodes + 8; n_nodes + 4];
        elements{n_elem + 4} = [n_nodes - 7; n_nodes - 5; n_nodes + 9; n_nodes + 5];
        elements{n_elem + 5} = [n_nodes + 3; n_nodes + 7; n_nodes + 4; n_nodes + 8; ...
                                n_nodes + 5; n_nodes + 9; n_nodes - 5; n_nodes + 6];

        boundary(n_bound + 1 : n_bound + 2) = [n_nodes + 1; n_nodes + 3];                            %Add the boundary nodes
        n_bound = n_bound + 2;

        if (j == n)
            
            boundary (n_bound + 1 : n_bound + 2) = [n_nodes + 4; n_nodes + 2];                        %If it is the last element of the row, we add         
            n_bound = n_bound + 2;                                                                   %the boundary nodes

        end

        n_nodes = n_nodes + 9;                                                                       %Update the number of objects that we have inserted
        n_elem  =  n_elem + 5;
  
    end

end

%% GENERATE THE OTHER ROWS

for i = 2:n

    %Construct the first element of the row
    if i == 2                                                                                        %Check how many elements we have to subtract
        minus = 9*n + 3;                                                                             %to recover the indexes of the nodes
    else
        minus = 7*n + 2;
    end

    midx = 1*h-h2;
    midy = i*h-h2;

    vertices(n_nodes + 1, :) = [h, i*h];                                                             %Coordinate of the top right node
    vertices(n_nodes + 2, :) = [0, i*h];                                                             %Coordinate of the top left  node

    vertices(n_nodes + 3, :) = [h,  i*h-h2];                                                         %Coordinate of the right mid node
    vertices(n_nodes + 4, :) = [h2, i*h];                                                            %Coordinate of the top   mid node
    vertices(n_nodes + 5, :) = [0,  i*h-h2];                                                         %Coordinate of the left  mid node

    vertices(n_nodes + 6, :) = [midx - dist, midy - dist];                                           %Coordinates of the the star
    vertices(n_nodes + 7, :) = [midx + dist, midy - dist];
    vertices(n_nodes + 8, :) = [midx + dist, midy + dist];
    vertices(n_nodes + 9, :) = [midx - dist, midy + dist];

    if (i == 2)

        elements{n_elem + 1} = [          4;           7; n_nodes + 6; n_nodes + 5];
        elements{n_elem + 2} = [          3; n_nodes + 3; n_nodes + 7;           7];
        elements{n_elem + 3} = [n_nodes + 1; n_nodes + 4; n_nodes + 8; n_nodes + 3];
        elements{n_elem + 4} = [n_nodes + 2; n_nodes + 5; n_nodes + 9; n_nodes + 4];
        elements{n_elem + 5} = [          7; n_nodes + 7; n_nodes + 3; n_nodes + 8; ...
                                n_nodes + 4; n_nodes + 9; n_nodes + 5; n_nodes + 6];

    else

        elements{n_elem + 1} = [n_nodes-minus+2; n_nodes-minus+4; n_nodes + 6;     n_nodes + 5];
        elements{n_elem + 2} = [n_nodes-minus+1;       n_nodes+3; n_nodes + 7; n_nodes-minus+4];
        elements{n_elem + 3} = [    n_nodes + 1;     n_nodes + 4; n_nodes + 8;     n_nodes + 3];
        elements{n_elem + 4} = [    n_nodes + 2;     n_nodes + 5; n_nodes + 9;     n_nodes + 4];
        elements{n_elem + 5} = [n_nodes-minus+4;     n_nodes + 7; n_nodes + 3;      n_nodes + 8; ...
                                    n_nodes + 4;     n_nodes + 9; n_nodes + 5;     n_nodes + 6];
    
    end

    boundary (n_bound + 1 : n_bound + 2) = [n_nodes + 5; n_nodes + 2];                               %Add the left boundary nodes
    n_bound = n_bound + 2;

    if (i == n)                                                                                      %If it is the last row
                                                                                                     %We had the top boundary nodes
        boundary (n_bound + 1 : n_bound + 2) = [n_nodes + 1; n_nodes + 4];
        n_bound = n_bound + 2;

    end

    n_nodes = n_nodes + 9;                                                                           %Update
    n_elem  = n_elem  + 5;

    %Construct the rest of the row
    for j = 2:n

        midx = j*h-h2;
        midy = i*h-h2;
    
        vertices(n_nodes + 1, :) = [j*h, i*h];                                                       %Coordinate of the top right node                            
    
        vertices(n_nodes + 2, :) = [j*h,    i*h-h2];                                                 %Coordinate of the right mid node
        vertices(n_nodes + 3, :) = [j*h-h2,    i*h];                                                 %Coordinate of the top   mid node
    
        vertices(n_nodes + 4, :) = [midx - dist, midy - dist];                                       %Coordinates of the top star
        vertices(n_nodes + 5, :) = [midx + dist, midy - dist];
        vertices(n_nodes + 6, :) = [midx + dist, midy + dist];
        vertices(n_nodes + 7, :) = [midx - dist, midy + dist];

        if (i == 2)                                                                                  %Recover the indexes of the other nodes

            if (j == 2) 

                bottom_left  = 3;
                bottom_right = 14;
                bottom_mid   = 17;
                left_mid     = 9*n + 6;
                top_left     = left_mid - 2;
               
            else

                bottom_left  = 12 + (j-3)*9 + 2;
                bottom_mid   = bottom_left + 12;
                bottom_right = bottom_mid  -  3;
                left_mid     = n_nodes - 5;
                top_left     = n_nodes - 6;

            end
            
        else

            bottom_right = 9*n + 3 + (i-3)*(7*n + 2) + 9 + (j-2) * 7 + 1;
            bottom_mid   = bottom_right + 2;

            if (j == 2)

                bottom_left = bottom_right - 9;
                left_mid    = n_nodes - 6;
                top_left    = n_nodes - 8;

            else

                bottom_left = bottom_right - 7;
                left_mid    = n_nodes - 5;
                top_left    = n_nodes - 6;

            end

        end

        elements{n_elem + 1} = [ bottom_left;  bottom_mid; n_nodes + 4;    left_mid];                %Add the elements
        elements{n_elem + 2} = [bottom_right; n_nodes + 2; n_nodes + 5;  bottom_mid];
        elements{n_elem + 3} = [ n_nodes + 1; n_nodes + 3; n_nodes + 6; n_nodes + 2];
        elements{n_elem + 4} = [    top_left;    left_mid; n_nodes + 7; n_nodes + 3];
        elements{n_elem + 5} = [  bottom_mid; n_nodes + 5; n_nodes + 2; n_nodes + 6; ...
                                 n_nodes + 3; n_nodes + 7;    left_mid; n_nodes + 4];


        if (j == n) && (i ~= n)                                                                      %Update the boundary

            boundary(n_bound + 1 : n_bound + 2) = [n_nodes + 1; n_nodes + 2];
            n_bound = n_bound + 2;

        elseif (j == n) && (i == n)

            boundary(n_bound + 1 : n_bound + 3) = [n_nodes + 1; n_nodes + 2; n_nodes + 3];
            n_bound = n_bound + 3;

        elseif (j ~= n) && (i == n)

            boundary(n_bound + 1 : n_bound + 2) = [n_nodes + 1; n_nodes + 3];
            n_bound = n_bound + 2;

        end

        n_nodes  = n_nodes + 7;
        n_elem   = n_elem + 5;

    end

end

filename = 'star_' + string(n) + '_' + string(perc);
save(filename, 'elements', 'vertices', 'boundary');