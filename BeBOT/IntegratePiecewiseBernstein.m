function integratedPoly = IntegratePiecewiseBernstein(coeffs, tknots, initialValue)
    % coeffs: A cell array where each cell contains the coefficients of a Bernstein polynomial piece.
    % tknots: A vector specifying the end-points of the sub-intervals of the domain.
    % initialValue: The initial value of the integrated polynomial at the first point of the first interval.

%%    Example
%     for i = 1:K
%     coeffs{i} = rand(N+1,1)';
%     end
%     tknots = linspace(0,1,K+1);
%     integratedPoly = IntegratePiecewiseBernstein(coeffs, tknots, 2);
    
    nPieces = length(coeffs); % Number of polynomial pieces
    integratedPoly = cell(1, nPieces); % Initialize the cell array for the integrated pieces
    
    
    for i = 1:nPieces
        % Degree of the current polynomial piece
        N = numel(coeffs{i})-1;
        I = BernsteinIntegrationMatrix(N,tknots(i+1)-tknots(i));

        
        % Adjust the constant of integration for the initial condition or continuity
        if i == 1
            % For the first piece, adjust based on the provided initial condition
            integratedCoeffs = coeffs{i}*I + initialValue;
        else
            % Adjust based on continuity with the end of the previous piece
            integratedCoeffs = coeffs{i}*I + integratedPoly{i-1}(end);
        end
        
        integratedPoly{i} = integratedCoeffs;
    end
end

