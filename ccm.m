%
% Implements Convergent Cross Mapping (CCM) (Sugihara et al. Science 2012)
% 20120408 - Leonardo Barbosa
% 
% [x_r, y_r, x_o, y_o] = ccm(x,y,tau)
%
% Parameters:
%
% x, y : original time series
% tau  : step between samples to construct the shadow manifold (default 1)
% E    : size of the window to construct the manifold (default 10% original size)
%
% Returns:
%
% x_r, y_r : reconstructed time series
% x_o, y_o : truncated original time series (so that first sample of x_o is
% contemporaneous of x_r, i.e. reconstructed time series are missing N
% samples in the beggining, where N = 1+(E-1)*tau)
%

function [x_r, y_r, x_o, y_o] = ccm(x, y, tau, E)

    if nargin < 3
        tau = 1;
    end

    if nargin < 4
        E = 2;
    end
    
    L = length(x);

    if L ~= length(y)
        error('series should have same size')
    end
    
    % construct the shadow manifolds
    M_x = shadow(x, tau, E);
    M_y = shadow(y, tau, E);

    t0 = 1+(E-1)*tau;
    
    % intialize reconstructed time series
    x_r = zeros(1, L-t0+1);
    y_r = zeros(1, L-t0+1);
    
    % shrink time series to remove first points and make t = t_x = t_y 
    % for (M(t), x(t_x), y(t_y))
    % i.e. we can only reconstruct from 1+(E-1)*tau and future samples
    x_o = x(t0:L);
    y_o = y(t0:L);
    for t=1:L-t0+1
        
        % reconstuct x
        x_r(t) = reconstruct(t, x_o, M_y);

        % reconstuct y
        y_r(t) = reconstruct(t, y_o, M_x);
        
    end
end

function M = shadow(v, tau, E)
    
    t0 = 1+(E-1)*tau;
    L = length(v);
    M = zeros(L-t0+1, floor(E/tau)); 
    
    for t = t0:L
        M(t-t0+1, :) = v(t-(E-1)*tau:tau:t);
    end
    
end

function r = reconstruct(t, v, M)

    % find nearest neighbors
    [d, n] = sort(sqrt(sum((M-repmat(M(t,:), [size(M,1) 1])).^2, 2)));
    
    E = size(M,2);
    
    n = n(2:E+2);
    d = d(2:E+2);

    % calculate weights
    u = exp(-d./d(1))';
    w = u ./ sum(u);

    % reconstuct
    r = sum(w .* v(n));

end