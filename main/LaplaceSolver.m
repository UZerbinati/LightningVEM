function [f, g] = LaplaceSolver(polygon)

outward = zeros(polygon.nedges,1);

w = complex(polygon.vertex);                                                %Convert the vertex of the polygon to complex numbers

%%COMPUTE THE BISECTORS
for k = 1:polygon.nedges

    forward = w(k+1) - w(k);
    j = mod(k-2, polygon.nedges) + 1;
    backward = w(j) - w(k);
    tmp = 1i*backwar*sqrt(-forward/backward);
    outward(k) = tmp/abs(tmp);

end

%% MAIN LOOP, INCREASE THE NUMBER OF POLES
errc = ones(polygon.nedges,1);                                              %Error at each corner
np   = zeros(polygon.nedges,1);                                             %N of poles at each corner
max_step = 0; err0 = Inf;

sigma = ?

for step = 1:max_step

   Z = [];           % col vector of sample points on boundary
   G = [];           % col vector of boundary values at these points
   T = [];           % col vector of unit tangent vectors at these points
   pol = [];         % row vector of poles of the rational approximation
   J = [];           % row vector of indices of which corner each pole belongs to
   d = [];           % row vector of distances from poles to their corners
   tt = cell(nw,1);  % cell array of distances of sample points along each side

   for k = 1:polygon.nedges

       npk = np(k);                                                         %Number of points in this corner
       sk  = sqrt(1:npk) - sqrt(npk);
       dk = exp(sigma(k)*sk);
       dk = dk * polygon.diameter;
       dk = dk(dk > 1e-15*polygon.diameter);
       polk = w(k) + outward(k)*dk;
   end

end