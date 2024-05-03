function integratedPoly = IntegratePiecewiseBernstein(cp, tknots, initialValue)
   
    % Number of pieces and degree
    nPieces = length(tknots) - 1; % Number of polynomial pieces
    totalCoeffs = length(cp); % Total number of coefficients
    coeffsPerPiece = totalCoeffs / nPieces; % Coefficients per piece
    N = coeffsPerPiece - 1; % Degree of the polynomials is N-1 but each piece is represented by N coefficients

    % Initialize output vector, which is longer due to integration increasing degree by 1 for each piece
    integratedPoly = zeros(1, totalCoeffs + nPieces); % One more coefficient per piece

    % Initial index for input and output vectors
    inIndex = 1; 
    outIndex = 1;

    for i = 1:nPieces
        % Extract the coefficients for the current piece from the input vector
        currentCoeffs = cp((i-1)*N+i:i*N+i);
        
        % Integration matrix computation. This needs to match the dimensions of currentCoeffs.
        % Assuming a simple integration matrix here, but you'll replace this with BernsteinIntegrationMatrix as needed.
        I = BernsteinIntegrationMatrix(N,tknots(i+1)-tknots(i));

        % Compute the integrated coefficients
        % Note: The actual integration logic and matrix must be applied here.
        integratedCoeffs = currentCoeffs * I; % This operation may need revision based on actual I
        
        % Adjust the constant term for initial value or continuity
        if i == 1
            % Add the initial value for the first polynomial piece
            integratedCoeffs = integratedCoeffs + initialValue;
        else
            % Adjust for continuity; you might need to adjust based on actual integration results
            integratedCoeffs = integratedCoeffs + integratedPoly((i-1)*(N+1)+i-1);
        end
        
        % Assign the integrated coefficients to the output vector
        integratedPoly((i-1)*(N+1)+i:i*(N+1)+i) = integratedCoeffs;

    end

end

