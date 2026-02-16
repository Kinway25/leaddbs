function varargout = ea_explorer_stats_1sampleweightedlinreg(varargin)
% Function to estimate a certain mass univariate test as used in explorer apps.

if nargin > 0
    % Map inputs
    valsin = varargin{1};
    outcomein = varargin{2};
    H0_input = varargin{3};
else
    varargout{1}.name = "1-Sample Weighted Regression";
    varargout{1}.file = mfilename;
    varargout{1}.type = "1-Sample Tests";
    varargout{1}.outcometype = {'gradual'};
    varargout{1}.compatibility = {'Electric Field', 'Sigmoid Field'};    
    return;
end

nan_mask = isnan(outcomein);
if any(nan_mask)
    outcomein(nan_mask) = [];
    valsin(:, nan_mask) = []; % Ensure dimensions stay aligned
end

% Check if weights are outside the required range [0, 1]
if any(valsin(:) < 0, 'all') || any(valsin(:) > 1, 'all')
    valsin = normalize(valsin, 'range', [0, 1]);
end

% Resolve H0
if ischar(H0_input)
    switch H0_input
        case 'Average'
            H0_val = mean(outcomein, 'all', 'omitnan');
        case 'Zero'
            H0_val = 0;
    end
else
    H0_val = H0_input;
end

% Initialize output arrays
valsout = nan(size(valsin, 1), 1);
psout = nan(size(valsin, 1), 1);

% Check for Parallel Computing Toolbox
if license('test', 'Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
    parpool;
end

% We'll use a local variable for the outcome to avoid overhead in parfor
local_outcome = outcomein(:)'; 

%ICC_table = readtable('/home/interscan/Documents/data/JS/ReFitCohort_Avg_ICC.csv');


if license('test', 'Distrib_Computing_Toolbox')
    parfor i = 1:size(valsin, 1)
        % --- NaN HANDLING ---
        valid_idx = ~isnan(valsin(i, :)) & ~isnan(local_outcome);
        curr_vals = valsin(i, valid_idx)'; % These are our Weights (W)
        curr_outcome = local_outcome(valid_idx)'; % This is our Data (Y)
        
        N_valid = sum(valid_idx);
        if N_valid <= 1
            continue; 
        end
        
        % --- STATISTICAL CORRECTION ---
        % 1. Center the outcome around the Null Hypothesis
        Y = curr_outcome - H0_val;
        
        % 2. Predictor is just a constant (intercept) for 1-sample test
        X = ones(N_valid, 1);
        
        % 3. Apply weights
        W = curr_vals;
        Wsqrt = sqrt(W); % Using element-wise sqrt since W is a vector here
        Xw = X .* Wsqrt; 
        Yw = Y .* Wsqrt;
        
        % 4. Weighted Least Squares Solve
        b = Xw \ Yw; % This is the weighted mean difference from H0
        
        % 5. Degrees of Freedom (Corrected: N - 1)
        df = N_valid - 1;
        
        % 6. Standard Error Calculation
        residuals = Yw - Xw * b;
        sigma2 = (residuals' * residuals) / df;
        
        % Standard error of the weighted mean
        % se = sqrt(sigma2 * inv(Xw' * Xw))
        se = sqrt(sigma2 / (Xw' * Xw)); 
        
        % 7. Statistics
        tStat = full(b / se);
        valsout(i) = tStat;
        psout(i) = 2 * (1 - tcdf(abs(tStat), df));
    end
else
    % Standard processing (Logic identical to parfor block)
    for i = 1:size(valsin, 1)
        % --- NaN HANDLING ---
        valid_idx = ~isnan(valsin(i, :)) & ~isnan(local_outcome);

        % if full(sum(valid_idx)) ~= size(valsin,2)
        %     disp("NaNs detected")
        % end

        curr_vals = valsin(i, valid_idx)'; % These are our Weights (W)
        curr_outcome = local_outcome(valid_idx)'; % This is our Data (Y)
        
        N_valid = sum(valid_idx);
        if N_valid <= 1
            continue; 
        end
        
        % --- STATISTICAL CORRECTION ---
        % 1. Center the outcome around the Null Hypothesis
        Y = curr_outcome - H0_val;
        
        % 2. Predictor is just a constant (intercept) for 1-sample test
        X = ones(N_valid, 1);
        
        % 3. Apply weights
        W = curr_vals;
        Wsqrt = sqrt(W); % Using element-wise sqrt since W is a vector here
        Xw = X .* Wsqrt; 
        Yw = Y .* Wsqrt;
        
        % 4. Weighted Least Squares Solve
        b = Xw \ Yw; % This is the weighted mean difference from H0
        
        % 5. Degrees of Freedom (Corrected: N - 1)
        df = N_valid - 1;
        
        % 6. Standard Error Calculation
        residuals = Yw - Xw * b;
        sigma2 = (residuals' * residuals) / df;
        
        % Standard error of the weighted mean
        % se = sqrt(sigma2 * inv(Xw' * Xw))
        se = sqrt(sigma2 / (Xw' * Xw)); 
        
        % 7. Statistics
        tStat = full(b / se);
        valsout(i) = tStat;
        psout(i) = 2 * (1 - tcdf(abs(tStat), df));
    end
end

% Map outputs
varargout{1} = valsout;
varargout{2} = psout;
end