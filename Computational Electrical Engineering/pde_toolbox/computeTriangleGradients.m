function [gradients, centroids] = computeTriangleGradients(p, t, u)
% computeTriangleGradients - Computes per-triangle gradients of a scalar field
%                             and the centroids of the triangles
%
% Syntax:
%   [gradients, centroids] = computeTriangleGradients(p, t, u)
%
% Inputs:
%   p - N x 2 array of node coordinates (x, y)
%   t - M x 3 array of triangle vertex indices
%   u - N x 1 vector of nodal scalar values
%
% Outputs:
%   gradients - M x 2 array of gradient vectors [∂u/∂x, ∂u/∂y] for each triangle
%   centroids - M x 2 array of triangle centroids (x, y)

    num_tri = size(t, 1);
    gradients = zeros(num_tri, 2);
    centroids = zeros(num_tri, 2);

    for i = 1:num_tri
        tri = t(i, :);              % Node indices of triangle
        coords = p(tri, :);         % Node coordinates (3x2)
        x = coords(:,1);
        y = coords(:,2);
        u_tri = u(tri);             % Nodal values in triangle

        % Solve for coefficients in u(x, y) = a + b*x + c*y
        A = [ones(3,1), x, y];
        coeff = A \ u_tri;

        % Gradient is [b; c]
        gradients(i, :) = coeff(2:3)';

        % Centroid is mean of the triangle's vertex coordinates
        centroids(i, :) = mean(coords, 1);
    end
end
