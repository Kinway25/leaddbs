function varargout = ea_explorer_stats_2sampleweightedlinreg(varargin)
% Function to estimate a certain mass univariate test (see below) as used
% in explorer apps.
% - Expects valsin as V x N matrix where 
%                                    V is the number of
%                                        voxels/streamlines/etc and 
%                                    N is the number of 
%                                        E-Fields/VTAs/Lesions/etc.
% - Expects outcomein as improvement values with dimension N x 1 or 2.
% - Outputs valsout (test results) and psout (p-values of test results) as V x 1 or 2.

if nargin > 0
    % Map inputs
    valsin = varargin{1};
    outcomein = varargin{2};
else
    varargout{1}.name = "2-Sample Weighted Regression";
    varargout{1}.file = mfilename;
    varargout{1}.type = "2-Sample Tests";
    varargout{1}.outcometype = {'gradual'};
    varargout{1}.compatibility = {'Electric Field', 'Sigmoid Field'}; 
    return;
end

% Check if weights are outside the required range [0, 1]
if any(valsin(:) < 0) || any(valsin(:) > 1)
    valsin = normalize(valsin, 'range', [0, 1]);
end


% Actual test:
outcomein = repmat(outcomein', size(valsin, 1), 1);
group1 = outcomein;

% Initialize output arrays
valsout = nan(size(valsin, 1), 1);
psout = nan(size(valsin, 1), 1);

% Check for Parallel Computing Toolbox
if license('test', 'Distrib_Computing_Toolbox')
    % Start parallel pool if not already started
    if isempty(gcp('nocreate'))
        parpool;
    end

    % Parallel processing
    parfor i = 1:size(valsin, 1)
        
        % *** MODIFICATION START ***
        % 1. Find non-NaN indices for the current voxel (row i)
        valid_idx = ~isnan(valsin(i, :)) & ~isnan(group1(i, :));
        
        % 2. Subset the data vectors using the valid indices
        current_valsin = valsin(i, valid_idx);
        current_group1 = group1(i, valid_idx);
        
        % Check if enough valid data points remain for calculation
        N_valid = length(current_valsin);
        % The regression requires at least two data points for the two-column X matrix
        % size(X, 1) - size(X, 2) is the degrees of freedom (df). df must be > 0.
        % size(X, 1) = 2 * N_valid. size(X, 2) = 2.
        % So, 2 * N_valid - 2 > 0 -> 2 * N_valid > 2 -> N_valid > 1.
        if N_valid <= 1
            continue; % Skip to the next iteration, results remain NaN
        end
        
        % Prepare data for regression using the subsetted vectors
        Y = [current_group1'; current_group1'];
        X = [ones(size(Y)), [zeros(size(current_valsin')); ones(size(current_valsin'))]];
        W = [1 - current_valsin'; current_valsin'];
        % *** MODIFICATION END ***

        % Perform weighted least squares regression
        % The rest of the logic is unchanged but now operates on the filtered data
        Wsqrt = diag(sqrt(W));
        Xw = Wsqrt * X;
        Yw = Wsqrt * Y;
        
        % Note: \ performs least squares, which is more robust than inv(Xw' * Xw)
        b = Xw \ Yw; 
        
        % The standard error calculation requires the residual degrees of freedom:
        df = size(X, 1) - size(X, 2); 
        
        residuals = Yw - Xw * b;
        sigma2 = (residuals' * residuals) / df;
        
        % Use pinv (pseudo-inverse) for robustness in case of rank-deficiency, 
        % though standard inv might be sufficient if data is good.
        C = sigma2 * pinv(Xw' * Xw); 
        se = sqrt(diag(C));

        % Store results
        tStat = b(2) / se(2); % t-statistic
        valsout(i) = tStat;
        psout(i) = 2 * (1 - tcdf(abs(tStat), df)); % p-value
    end
else
    % Standard processing
    for i = 1:size(valsin, 1)
        
        % *** MODIFICATION START ***
        % 1. Find non-NaN indices for the current voxel (row i)
        valid_idx = ~isnan(valsin(i, :));
        
        % 2. Subset the data vectors using the valid indices
        current_valsin = valsin(i, valid_idx);
        current_group1 = group1(i, valid_idx);
        
        % Check if enough valid data points remain for calculation
        N_valid = length(current_valsin);
        if N_valid <= 1
            continue; % Skip to the next iteration, results remain NaN
        end
        
        % Prepare data for regression using the subsetted vectors
        Y = [current_group1'; current_group1'];
        X = [ones(size(Y)), [zeros(size(current_valsin')); ones(size(current_valsin'))]];
        W = [1 - current_valsin'; current_valsin'];
        % *** MODIFICATION END ***

        % Perform weighted least squares regression
        Wsqrt = diag(sqrt(W));
        Xw = Wsqrt * X;
        Yw = Wsqrt * Y;
        b = Xw \ Yw;
        
        df = size(X, 1) - size(X, 2);
        
        residuals = Yw - Xw * b;
        sigma2 = (residuals' * residuals) / df;
        
        C = sigma2 * pinv(Xw' * Xw);
        se = sqrt(diag(C));

        % Store results
        tStat = b(2) / se(2); % t-statistic
        valsout(i) = tStat;
        psout(i) = 2 * (1 - tcdf(abs(tStat), df)); % p-value
    end
end

% Map outputs
varargout{1} = valsout;
varargout{2} = psout;
end