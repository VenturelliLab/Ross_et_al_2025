function params = getLHS(n, p, scale, seed)
    % getLHS generates n, p-dimensional parameter sets based on Latin
    % Hypercube Sampling.
    % Inputs:
    %   n - number of samples
    %   p - number of dimensions/parameters
    %   scale - if size(scale) = [1, p], each parameter is scaled between
    %      0.5*p and 1.5*p, otherwise each parameter is scaled between
    %      scale(1, :) and scale(2, :)
    %   seed - if defined, seed is a positive integer used to set the state
    %      for random number generation
    % Output:
    %   params - an n x p matrix of scaled parameter sets

    if exist('seed', 'var')
        rng(seed)
    else
        rng('default')
    end

    if size(scale, 1) == 1
        lower_B = scale*0.5;
        upper_B = scale*1.5;
    
    elseif size(scale, 1) == 2
        lower_B = scale(1, :);
        upper_B = scale(2, :);
    end    

    params = lhsdesign(n, p)';
    params = params.*repmat(upper_B' - lower_B', 1, n);
    params = params + repmat(lower_B', 1, n);

end