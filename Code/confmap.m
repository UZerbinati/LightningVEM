function [f, finv, pol, polinv] = confmap(P, varargin)
%CONFMAP  Compute conformal map
%
%   F = CONFMAP(P) computes the conformal map to the unit disk 
%   from the simply-connected region Omega bounded by P, which
%   may be a polygon, circular polygon, or other region with corners.
%   P must enclose the point z=0, which is mapped to w=0.
%
%   [F,FINV] = CONFMAP(P) also computes the inverse map.
%
%   F and FINV are function handles that evaluate rational functions
%   obtained by AAA approximation after the lightning Laplace solver
%   LAPLACE is called to solved a Green's function problem.
%   [F,FINV,POL,POLINV] = CONFMAP(P) returns the poles of these
%   rational functions.
%  
%   This software requires Chebfun in the Matlab path: see www.chebfun.org.
%
%   v1, (c) Lloyd N. Trefethen, U. of Oxford, March 2020. 
%   For info see https://people.maths.ox.ac.uk/trefethen/lightning.html.
%   Please send corrections or suggestions to trefethen@maths.ox.ac.uk.
%
%  Inputs:
%     P = vector of corners as complex numbers z = x+iy in counterclockwise
%             order to specify a polygon
%         or cell array of corners v and pairs [v r] to specify a circular
%             polygon: r = radius of curvature of arc from this v to the next
%         or 'pent'[agon], 'snow'[flake], 'iso'[spectral], 'L', or 'circleL'
%         or integer >= 3, the number of corners of a random polygon
%         or integer <= -3, -1 x no. of corners of a random circular polygon
%
%  Further optional inputs:
%    'tol', tol = approximate tolerance (default 1e-6)
%    'numbers' to print certain numbers
%    'noplots' to suppress plots
%
%  For the standard software for conformal mapping of polygons, see Driscoll's
%  Schwarz-Christoffel Toolbox at math.udel.edu/~driscoll/software/.
%
% Reference: L. N. Trefethen, Numerical conformal mapping with rational 
%     functions, Computational Methods and Function Theory, to appear.
%
% Examples:
%
%   confmap([-1-1i 2-1i 2+2i -1+2i]);      % square
%   confmap('iso');                        % isospectral octagon
%   confmap('L','numbers');                % L shape
%   confmap(8);                            % random octagon
%   confmap(-8);                           % random circular octagon
%   confmap({[2-1i 1] 2+1i -2+1i -2-1i});  % bullet
%   confmap('circleL');                    % circular L-shape

tic
warnstate = warning; warning off
[tol, plots, numbers, P] = parseinputs(P, varargin{:});
[~, ~, w, Z, Zplot] = laplace(P,@(z) -log(abs(z)), 'tol', tol, 'noplots');
W = Z.*exp(w(Z));
[f, pol,~,~,~,~,~,evec] = ...
  aaa(W,Z,'tol', tol,'cleanup',1,'cleanuptol',tol/1,'mmax',140,'lawson',0);
[finv, polinv,~,~,~,~,~,evecinv] = ...
  aaa(Z,W,'tol', tol,'cleanup',1,'cleanuptol',tol/1,'mmax',140,'lawson',0);
tcomp = toc;

if plots 
    clf
    LW = 'linewidth'; PO = 'position'; MS = 'markersize';
    FW = 'fontweight'; NO = 'normal'; XT = 'xtick'; YT = 'ytick';
    circ = exp(2i*pi*(0:900)/900);

    h1 = axes(PO,[.04 .42 .45 .53]);   % plot P and poles of forward map
    xmin = min(real(Zplot)); xmax = max(real(Zplot));
    dx = xmax-xmin; xmid = mean([xmin xmax]);
    ymin = min(imag(Zplot)); ymax = max(imag(Zplot));
    dy = ymax-ymin; ymid = mean([ymin ymax]);
    dd = .65*max(dx,dy);
    plot(Zplot, 'b', LW, 1)
    axis([xmid-dd xmid+dd ymid-dd ymid+dd])
    axis square, hold on
    plot(pol, '.r', MS, 8)
    title([num2str(length(pol)) ' poles'], FW, NO)

    h2 = axes(PO,[.52 .42 .45 .53]);   % plot disk and poles of inverse map
    plot(circ, 'b', LW, 1), axis(1.7*[-1 1 -1 1])
    axis square, hold on, set(gca, XT,-1:1, YT,-1:1)
    plot(polinv, '.r', MS, 8)
    title([num2str(length(polinv)) ' poles'], FW, NO)

    ncirc = 8;                % plot concentric circles and their images
    for r = (1:ncirc-1)/ncirc
        axes(h1), plot(finv(r*circ), '-k', LW,.5)
        axes(h2), plot(r*circ, '-k', LW,.5)
    end
    
    nrad = 12;                % plot radii of unit circle and their images
    ray = chebpts(301);
    ray = ray(ray>=0);
    for k = 1:nrad
        axes(h1), plot(finv(ray*exp(2i*pi*k/nrad)), '-k', LW,.5)
        axes(h2), plot(ray*exp(2i*pi*k/nrad), '-k', LW,.5)
    end

    axes(h1), hold off
    axes(h2), hold off

end

% Print various quantities if requested
if numbers
    disp(' ')
    fprintf('                              solution time, in secs:%6.2f\n', tcomp)
    ns = 1e4; zz = .1*randn(ns,1) + .1i*randn(ns,1);
    tic, f(zz); tcompz = 1e6*toc/ns;
    tic, finv(zz); tcompw = 1e6*toc/ns;
    fprintf('eval time per point from Omega to disk, in microsecs:%6.2f\n', tcompz)
    fprintf('   eval time per point for inverse map, in microsecs:%6.2f\n', tcompw)
    fprintf('back-and-forth boundary error norm(Z-finv(f(Z)),inf):  %6.1e\n',...
        norm(Z-finv(f(Z)), inf))
    fprintf(' inverse back-and-forth error norm(W-f(finv(W)),inf):  %6.1e\n',...
        norm(W-f(finv(W)), inf))
    fprintf('                      error at points with |w| = 0.9:  %6.1e\n',...
        norm(.9*W-f(finv(.9*W)), inf))
    disp(' ')
end

warning(warnstate)

end   % end of conformal

function [tol, plots, numbers, P] = parseinputs(P, varargin)
if strcmp(P,'L'), P = [2 2+1i 1+1i 1+2i 2i 0] - (.75+.75i); end
if strcmp(P,'circleL'), P = {1-1i [1 -1] 1i -1+1i -1-1i}; end
tol = 1e-6; plots = 1; numbers = 0;     % defaults
j = 1;
while j < nargin
    j = j+1;
    v = varargin{j-1};
    if strcmp(v, 'tol')
        j = j+1;
        tol = varargin{j-1};
    elseif strcmp(v, 'noplots')
        plots = 0;
    elseif strcmp(v, 'numbers')
        numbers = 1;
    end
end
end   % end of parseinputs