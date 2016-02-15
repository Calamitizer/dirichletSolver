% dirichletSolver.m
% J. Alex Ruble
% MA 428 Project 1
% Viewable online at goo.gl/0nqhjA

function [U,Boundary,xP,yP,A] = dirichletSolver(region, domainDefFunc, nonhomFunc, res)
    % region: A (rectangular) region enveloping the domain Omega
    %   of the form {[x_min,x_max],[y_min,y_max]}
    % domainDefFunc: A function which returns whether an input
    %   point is inside Omega
    % nonhomFunc: The nonhomogeneous function in the differential equation
    % res: Resolution parameter (grid spacing = 2^(-res))
    % U: The solution matrix
    % Boundary: A 1/0/NaN matrix marking the edge of Omega (for visual plotting)
    % xP,yP: Coordinate matrices for plotting
    % A: The central difference operator matrix
    h = 2^(-res);
    xB = region{1};
    yB = region{2};
    [xP, yP] = meshgrid(xB(1):h:xB(2), yB(1):h:yB(2));
    [yN, xN] = size(xP); % Switched because Matlab counts row, column
    xScale = @(j) ((xN*xB(1)-xB(2))/(xN-1)) + j*((xB(2)-xB(1))/(xN-1));
    yScale = @(i) ((yN*yB(1)-yB(2))/(yN-1)) + i*((yB(2)-yB(1))/(yN-1));
    % These return the x-,y-coordinates of a matrix index
    Omega = arrayfun(domainDefFunc, xP, yP);
    Omega = double(Omega);
    Interior = getInterior(Omega);
    [to1D,to2D] = getLinearLookup(Interior);
    A = getAdjMatrix(to1D,to2D);
%     A = speye(size(to2D,1));
    b = getFuncVector(to2D,nonhomFunc,xScale,yScale,h);
    x = A\b;
    % Populate U:
    U = zeros(size(Omega));
    for i = 1:size(x,1)
        U(to2D(i,1),to2D(i,2)) = x(i);
    end
    % Do some stuff for nice plots
    Boundary = sparse(Omega - Interior);
    Boundary = trim(Boundary, Interior);
    Boundary = .25*max(U(:))*Boundary;
    U(Omega==0) = NaN; % instead of 0, so mesh() skips plotting them
end

function [Interior] = getInterior(D)
    Interior = zeros(size(D));
    % If a point is along the rectangular edge of definition for D, it is
    %   not an interior point, so we iterate over the non-ege elements.
    for i = 2:size(Interior,1)-1
        for j = 2:size(Interior,2)-1
            % If a point or any of its neighbors is empty
            %   in D, it is not an interior point. Otherwise, it is.
            Interior(i,j) = double(all(all(D(i-1:i+1,j-1:j+1))));
        end
    end
end

function [to1D,to2D] = getLinearLookup(D)
    % D: A binary matrix defining the inside of a region
    to1D = zeros(size(D));
    % A D-size matrix whose nonzero elements are the row indices of the corresponding
    %   points in to2D
    to2D = zeros(sum(D(:)),2);
    % A (# nonzeros in D)-by-2 matrix whose k-th row contains the i,j
    %   coordinates of point k in to1D
    hit = 0;
    for i = 1:size(D,1)
        for j = 1:size(D,2)
            if D(i,j)
                hit = hit + 1;
                to1D(i,j) = hit;
                to2D(hit,:) = [i,j];
            end
        end
    end
end

function [A] = getAdjMatrix(to1D,to2D)
    % This computes the (negative) central difference adjacency matrix for the
    %   interior points. We'll just do the upper half and exploit main diagonal symmetry.
    A = 2*speye(size(to2D,1));
    count = size(to1D,1);
    for i = 2:size(to1D,1)-1
        for j = 2:size(to1D,2)-1
            if to1D(i,j)
                if to1D(i,j+1)
                    A(to1D(i,j),to1D(i,j+1)) = -1;
                end
                if to1D(i+1,j)
                    A(to1D(i,j),to1D(i+1,j)) = -1;
                end
            end
        end
        disp([int2str(i),'/',int2str(count)]);
    end
    A = A + A.';
end

function [b] = getFuncVector(to2D,f,xS,yS,h)
    b = zeros(size(to2D,1),1);
    for i = 1:size(b)
        b(i) = h^2 * f(xS(to2D(i,2)),yS(to2D(i,1))); % Switched because row, col
    end
end

function [B] = trim(oldB, I)
    % Leaves a bit of zero space around nonzero elements so that mesh()
    %   will plot the outline, but renders most of them NaN
    B = oldB;
    B(oldB==0) = NaN;
    for i = 2:size(oldB,1)-1
        for j = 2:size(oldB,2)-1
            if I(i,j) && any(any(oldB(i-1:i+1,j-1:j+1)))
                B(i,j) = 0;
            end
        end
    end
end