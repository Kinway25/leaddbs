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
        % --- NaN HANDLING START ---
        % Identify non-NaN weights for this specific voxel
        valid_idx = ~isnan(valsin(i, :)) & ~isnan(local_outcome);
        
        % Subset data
        curr_vals = valsin(i, valid_idx)';
        curr_outcome = local_outcome(valid_idx)';
        
        %ICCs = ICC_table.ICC_a_hemisphere(valid_idx);
        %curr_vals = curr_vals.*ICCs;

        % Degrees of freedom check: 
        % We have 2*N observations and 2 parameters. Need 2*N - 2 > 0.
        if sum(valid_idx) <= 1
            continue; 
        end
        % --- NaN HANDLING END ---

        % Prepare data (Group 2 is H0, Group 1 is the actual outcome)
        Y = [repmat(H0_val, size(curr_outcome)); curr_outcome];
        X = [ones(size(Y)), [zeros(size(curr_vals)); ones(size(curr_vals))]];
        W = [curr_vals; curr_vals];

        % Perform weighted least squares regression
        Wsqrt = diag(sqrt(W));
        Xw = Wsqrt * X;
        Yw = Wsqrt * Y;
        
        % Solve and calculate stats
        b = Xw \ Yw;
        df = size(X, 1) - size(X, 2);
        residuals = Yw - Xw * b;
        sigma2 = (residuals' * residuals) / df;
        
        % Using pinv for better stability with small weights
        C = sigma2 * pinv(Xw' * Xw); 
        se = sqrt(diag(C));

        tStat = b(2) / se(2);
        valsout(i) = tStat;
        psout(i) = 2 * (1 - tcdf(abs(tStat), df));
    end
else
    % Standard processing (Logic identical to parfor block)
    for i = 1:size(valsin, 1)
        valid_idx = ~isnan(valsin(i, :)) & ~isnan(local_outcome);
        curr_vals = valsin(i, valid_idx)';
        curr_outcome = local_outcome(valid_idx)';

        % ICCs = ICC_table.ICC_a_hemisphere(valid_idx);
        % curr_vals = curr_vals.*ICCs;
        % 
        if sum(valid_idx) <= 1, continue; end

        Y = [repmat(H0_val, size(curr_outcome)); curr_outcome];
        X = [ones(size(Y)), [zeros(size(curr_vals)); ones(size(curr_vals))]];
        W = [curr_vals; curr_vals];

        Wsqrt = diag(sqrt(W));
        Xw = Wsqrt * X;
        Yw = Wsqrt * Y;
        b = Xw \ Yw;
        df = size(X, 1) - size(X, 2);
        residuals = Yw - Xw * b;
        sigma2 = (residuals' * residuals) / df;
        C = sigma2 * pinv(Xw' * Xw);
        se = sqrt(diag(C));

        tStat = b(2) / se(2);
        valsout(i) = tStat;
        psout(i) = 2 * (1 - tcdf(abs(tStat), df));
    end
end

% Map outputs
varargout{1} = valsout;
varargout{2} = psout;
end