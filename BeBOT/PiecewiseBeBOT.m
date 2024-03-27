function [tnodes, w, Dm] = PiecewiseBeBOT(N, tknots)
% This function calculates the Bernstein basis polynomial differentiation matrix
% for a piecewise polynomial representation over a given set of knots.
% N is the order of one Bernstein polynomial
% tknots include the first and final time

M = length(tknots) - 1; % M is the number of piece-wise polynomials
tnodes = zeros(1, M * (N + 1));

% Preallocate Dm with an estimated size
DmSize = N * M + M; % Adjust based on expected size
Dm = zeros(DmSize, DmSize);

% Compute time nodes and fill Dm
for i = 1:M
    % Compute the start and end indices for the current segment
    startIndex = (i - 1) * N + i;
    endIndex = i * N + i;
    
    % Adjust the linspace operation based on the segment
    if i == 1
        tnodes(startIndex:endIndex) = linspace(tknots(i), tknots(i + 1) - eps, N + 1);
    elseif i == M
        tnodes(startIndex:endIndex) = linspace(tknots(i) + eps, tknots(i + 1), N + 1);
    else
        tnodes(startIndex:endIndex) = linspace(tknots(i) + eps, tknots(i + 1) - eps, N + 1);
    end
    
    % Populate Dm using block diagonal
    tempDm = BernsteinDifferentiationMatrix(N, tknots(i + 1) - tknots(i)) * DegElevMatrix(N - 1, N);
    DmBlockStart = (i - 1) * N + (i - 1) + 1;
    DmBlockEnd = DmBlockStart + N;
    Dm(DmBlockStart:DmBlockEnd, DmBlockStart:DmBlockEnd) = tempDm;
end

% Preallocate and compute weights
w = repmat((tknots(2:end) - tknots(1:end-1))' / (N + 1), N + 1, 1);
w = reshape(w, M * (N + 1), 1);

end

