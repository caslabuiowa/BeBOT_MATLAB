function [tnodes, w, Dm] = PiecewiseLGL_PS_mod(N, tknots)
    % This function outputs the LGL nodes, LGL integration weigths, and LGL differentiation matrix for piecewise polynomials


    % Written by Venanzio Cichella

    % Preallocate arrays based on expected sizes
    M = length(tknots) - 1; % M is the number of piece-wise polynomials
    tnodes = zeros(1, M * (N + 1));
    w = zeros(M * N + 1 , 1);
    Dm = zeros(M * (N+1), M * (N+1)); % Adjust based on the size returned by LGL_PS

    for i = 1 : M
        [nodeSegment, wSegment, Diff] = LGL_PS(N, tknots(i + 1) - tknots(i));
        nodeSegment = nodeSegment + tknots(i); % Adjust nodes

        % Populate the preallocated arrays
        tnodes((i - 1) * (N + 1) + 1 : i * (N + 1)) = nodeSegment;
        w((i - 1) * (N + 1) + 1 : i * (N + 1)) = wSegment;
        
        % Efficient block diagonal construction for Dm
        % Efficient block diagonal construction for Dm
        startIndex = (i - 1) * (N + 1) + 1;
        endIndex = i * (N + 1);
        Dm(startIndex:endIndex, startIndex:endIndex) = Diff;
    end
    
    
    
end