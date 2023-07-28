function [y, x] = max(f, flag, dim)
%MAX   Maximum value of a CHEBFUN.
%   MAX(F) and MAX(F, 'global') return the maximum value of the CHEBFUN F.
%
%   [Y, X] = MAX(F) returns also a point X such that F(X) = Y.
%
%   [Y, X] = MAX(F, 'local') returns not just the global maximum value but all
%   of the local maxima.
%
%   If F is complex-valued, absolute values are taken to determine the maxima,
%   but the resulting values correspond to those of the original function.
%
%   If F is array-valued, then the columns of X and Y correspond to the columns
%   of F. NaNs are used to pad Y and X when the 'local' flag is used and the
%   columns are not of the same length.
%
%   H = MAX(F, G), where F and G are CHEBFUNs defined on the same domain,
%   returns a CHEBFUN H such that H(x) = max(F(x), G(x)) for all x in the
%   domain of F and G. Alternatively, either F or G may be a scalar.
%
%   MAX(F, [], DIM) computes the maximum of the CHEBFUN F in the dimension DIM.
%   If DIM = 1 and F is a column CHEBFUN or DIM = 2 and F is a row CHEBFUN, this
%   is equivalent to MAX(F). Otherwise, MAX(F, [], DIM) returns a CHEBFUN which
%   is the maximum across the discrete dimension of F. For example, if F is a
%   quasimatrix with two columns, MAX(F, [], 2) = MAX(F(:,1), F(:,2)).
%
% See also MIN, MINANDMAX, ROOTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Trivial empty case.
if ( isempty(f) )
    x = [];
    y = [];
    return
end

if ( nargin == 3 )
    if ( ~any(dim == [1 2]) )
        error('CHEBFUN:CHEBFUN:max:badDim', ...
            'DIM input to CHEBFUN MAX must be 1 or 2.');
    end

    if ( dim ~= 1 + f(1).isTransposed )
        % Take max across discrete dimension of a quasimatrix:
        f = cheb2cell(f);
        y = -realmax;
        for k = 1:numel(f)
            y = maxOfTwoChebfuns(y, f{k});
        end
        y = merge(y);
        return
    end
end

if ( (nargin > 2) && isempty(flag) )
    % MAX(F, [], 1).
    flag = 'global';
end

if ( (nargin == 1) || strcmp(flag, 'global') ) 
    % MAX(F) or MAX(F, 'global')
    [y, x] = globalMax(f);    
    
elseif ( isa(flag, 'chebfun') || isnumeric(flag) )
    % MAX(F, G)
    y = maxOfTwoChebfuns(f, flag);
    
elseif ( strcmp(flag, 'local') )
    % MAX(F, 'local')
    [y, x] = localMax(f);
    
else
    error('CHEBFUN:CHEBFUN:max:flag', 'Unrecognized flag.');
    
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL MAXIMUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, x] = globalMax(f)

% Call MINANDMAX():
[y, x] = minandmax(f);

% Extract the maximum:
y = y(2,:);
x = x(2,:);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL MAXIMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, x] = localMax(f)

% Call MINANDMAX():
[y, x] = minandmax(f, 'local');

% Determine which are maxima.

ends = f(1).domain([1, end]).'; % Endpoints of the domain are special.
f = mat2cell(f); % Convert f into a cell of scalar-valued CHEBFUNs.

% Loop over the FUNs:
for k = 1:numel(f)
    % Compute 1st and 2nd derivatives.
    dfk = diff(f{k});
    dfk2 = diff(dfk);

    % For interior extrema, look at 2nd derivative:
    maximaLoc = feval(diff(f{k}, 2), x(:,k)) < 0;

    % For end-points, look at 1st derivative:
    dfk_ends = feval(dfk, ends);
    endptMaxLoc = dfk_ends.*[1, -1]' < 0;

    % If 1st derivative is small at an endpoint, assume it's zero and try to
    % use 2nd derivative to determine if it's a minimum:
    %
    % [TODO]:  What if the 2nd derivative is zero, so that rounding error
    % precludes us from accurately determining the sign?
    smallEndDer = abs(dfk_ends) < 1e3*vscale(dfk)*eps;
    endptMaxLoc(smallEndDer) = feval(dfk2, ends(smallEndDer)) < 0;

    maximaLoc(1) = endptMaxLoc(1);
    maximaLoc(x(:,k) == ends(2)) = endptMaxLoc(2);

    % Set points corresponding to local maxima to NaN:
    y(~maximaLoc,k) = NaN;
    x(~maximaLoc,k) = NaN;

    % Sort the result
    [x(:,k), maximaLoc] = sort(x(:,k));
    y(:,k) = y(maximaLoc,k);
end

% Remove any rows which contain only NaNs.
x(all(isnan(x), 2),:) = []; 
y(all(isnan(y), 2),:) = [];

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAX(F, G) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = maxOfTwoChebfuns(f, g)
% Return the function h(x) = max(f(x), g(x)) for all x. 

% If one is complex, use abs(f) and abs(g) to determine which function values to
% keep. (experimental feature)
if ( isreal(f) && isreal(g) && (nargin < 3) )
	S = sign(f - g);
else
	S = sign(abs(f) - abs(g));
end

% Heaviside function (0 where f > g, 1 where f < g);
H = 0.5*(S + 1);
notH = 0.5*(1 - S); % ~H.

% Combine for output:
h = H.*f + notH.*g;

% [TODO]: Enforce continuity?

% TODO: Why simplify?
% Simplify:
h = simplify(h);

end
