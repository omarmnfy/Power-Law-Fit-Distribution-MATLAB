% plfit: Fit a Power-Law Distribution to Data
% PLFIT(x) estimates x_min and alpha according to the goodness-of-fit based method described in Clauset, Shalizi, Newman (2007). x is a  vector of observations of some quantity to which we wish to fit the  power-law distribution p(x) ~ x^-alpha for x >= xmin. PLFIT automatically detects whether x is composed of real or integer values, and applies the appropriate method. For discrete data, if min(x) > 1000, PLFIT uses the continuous approximation, which is  a reliable in this regime.
% Original Code Source: http://www.santafe.edu/~aaronc/powerlaws/
% This function implements both the discrete and continuous maximum likelihood estimators for fitting the power-law distribution to data, along with the goodness-of-fit based approach to estimating the lower cutoff for the scaling region. Usage information is included in the file; type 'help plfit' at the Matlab prompt for more information.

function [alpha, xmin, L] = plfit(x, varargin)

% Varagin is an input variable in a function definition statement that enables the function to accept any number of input arguments.

vec     = [];
sample  = [];
xminx   = [];
limit   = [];
finite  = false;
nosmall = false;
nowarn  = false;

% Parse Command-line Parameters; Trap for Bad Input

i = 1;

% This initializes a counter i to 1, which will be used to iterate through the varargin cell array containing optional parameters passed to the function.

while i <= length(varargin)

% This starts a while loop that will run as long as i is less than or equal to the number of elements in varargin.

  argok = 1;

% argok is a flag variable set to 1 (true) initially. It will be used to check if the current argument is valid.
% A flag variable is a variable used to signal whether a particular condition has been met or not.
% As the loop processes each argument in varargin, it uses a switch statement to check if the argument matches any of the expected parameter names ('range', 'sample', 'limit', 'xmin', 'finite', 'nowarn', 'nosmall'). If the argument matches one of these cases, the corresponding action is taken (e.g., assigning a value or setting a boolean flag).
% If the argument does not match any of the expected cases (handled by the otherwise part of the switch statement), argok is set to 0. This indicates that the current argument is invalid.

  if ischar(varargin{i}) 
    switch varargin{i}
        case 'range',        vec     = varargin{i+1}; i = i + 1;
        case 'sample',       sample  = varargin{i+1}; i = i + 1;
        case 'limit',        limit   = varargin{i+1}; i = i + 1;
        case 'xmin',         xminx   = varargin{i+1}; i = i + 1;
        case 'finite',       finite  = true;
        case 'nowarn',       nowarn  = false;
        case 'nosmall',      nosmall = true;
        otherwise, argok=0; 
    end
  end
 
  % This switch statement checks the value of varargin{i} and assigns the appropriate value to the corresponding parameter:
    % 'range': Sets vec to the next element varargin{i+1} and increments i by 1.
	% 'sample': Sets sample to the next element varargin{i+1} and increments i by 1.
	% 'limit': Sets limit to the next element varargin{i+1} and increments i by 1.
	% 'xmin': Sets xminx to the next element varargin{i+1} and increments i by 1.
	% 'finite': Sets the boolean finite to true.
	% 'nowarn':
        % The nowarn parameter is a flag that, when set to true, suppresses warning messages during the execution of the plfit function. Specifically, it prevents the function from displaying a warning about potential finite-size bias when the number of data points is small. When nowarn is not set or set to false, the function will display a warning message if the number of data points used to estimate alpha is less than 50 and the finite correction is not applied.
	% 'nosmall':
        % The nosmall parameter is a flag that, when set to true, prevents the function from using values of xmin that would result in a large finite-size bias. This is particularly useful to avoid overestimating alpha when the number of observations above xmin is too small. When nosmall is not set or set to false, the function will consider all possible values of xmin during the fitting process. When nosmall is set to true, the function will skip values of xmin for which the finite-size bias exceeds a certain threshold (e.g., 0.1).
	% otherwise: If the argument does not match any of the specified cases, sets argok to 0 (false), indicating an invalid argument.

  if ~argok 
    disp(['(PLFIT) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end
if ~isempty(vec) && (~isvector(vec) || min(vec)<=1)
	fprintf('(PLFIT) Error: ''range'' argument must contain a vector; using default.\n');
    vec = [];
end
if ~isempty(sample) && (~isscalar(sample) || sample<2)
	fprintf('(PLFIT) Error: ''sample'' argument must be a positive integer > 1; using default.\n');
    sample = [];
end
if ~isempty(limit) && (~isscalar(limit) || limit<min(x))
	fprintf('(PLFIT) Error: ''limit'' argument must be a positive value >= 1; using default.\n');
    limit = [];
end
if ~isempty(xminx) && (~isscalar(xminx) || xminx>=max(x))
	fprintf('(PLFIT) Error: ''xmin'' argument must be a positive value < max(x); using default behavior.\n');
    xminx = [];
end

% Reshape Input Vector

x = reshape(x,numel(x),1);

% Reshapes x to a column vector. For example, if x [1.2, 2.3, 3.7, 4.1, 5.9, 6.4, 7.2, 8.8, 9.1]; then afrer reshaping x = [1.2; 2.3; 3.7; 4.1; 5.9; 6.4; 7.2; 8.8; 9.1]; 

% Select Method (Discrete or Continuous) for Fitting

if isempty(setdiff(x, floor(x)))
    f_dattype = 'INTS';  % Integer Data
elseif isreal(x)
    f_dattype = 'REAL';  % Real (Continuous) Data
else
    f_dattype = 'UNKN';  % Unknown Data Type
end
if strcmp(f_dattype, 'INTS') && min(x) > 1000 && length(x) > 100
    f_dattype = 'REAL';
end

% setdiff(x, floor(x)) checks if x contains non-integer values.
% isreal(x) checks if x contains real numbers.

% Estimate of xmin and alpha

switch f_dattype

% Continuous Data Case (REAL)

    case 'REAL'
        xmins = unique(x);
        xmins = xmins(1:end-1);

% unique(x) returns unique values of x sorted in ascending order.
% xmins = xmins(1:end-1) removes the last element from xmins to avoid including the maximum value as a candidate for xmin.

% Apply Optional Parameters

        if ~isempty(xminx)
            xmins = xmins(find(xmins>=xminx,1,'first'));
        end
        if ~isempty(limit)
            xmins(xmins>limit) = [];
        end
        if ~isempty(sample)
            xmins = xmins(unique(round(linspace(1,length(xmins),sample))));
        end

% These conditions modify xmins based on optional parameters (xminx, limit, sample). 

% Initialize Data Array

        dat   = zeros(size(xmins));

% dat is initialized to store the KS statistics for each candidate xmin.
% dat is an array of zeros with the same size as xmins.

% Sort Data

        z     = sort(x);

% Sorts x in ascending order.

        for xm = 1:length(xmins)
            xmin = xmins(xm);
            z    = z(z>=xmin); 
            n    = length(z);

% The loop iterates over each value in xmins.
% xmin is set to the current candidate value from xmins.
% z is filtered to include only values greater than or equal to the current xmin.
% n is the number of values in z.

% Estimate alpha Using Direct MLE

            a    = n ./ sum( log(z./xmin) );

% Optional Finite-Size Bias Correction

            if nosmall
                if (a-1)/sqrt(n) > 0.1
                    dat(xm:end) = [];
                    xm = length(xmins)+1;
                    break;
                end
            end

% If nosmall is true, the code checks if the finite-size bias is greater than 0.1.
% If the bias is too large, it removes the remaining xmin candidates and breaks out of the loop.

% Compute KS Statistic

            cx   = (0:n-1)'./n;
            cf   = 1-(xmin./z).^a;
            dat(xm) = max( abs(cf-cx) );
        end

% cx is the empirical cumulative distribution function (CDF) of the data.
% cf is the CDF of the fitted power-law model.
% The KS statistic is the maximum absolute difference between cf and cx and is stored in dat.

% Select the Best xmin

        D     = min(dat);
        xmin  = xmins(find(dat<=D,1,'first'));

% D is the minimum KS statistic value from dat.
% xmin is set to the corresponding value from xmins.

% Recalculate alpha and Calculate L (Log-Likelihood) with the Best xmin

        z     = x(x>=xmin);
        n     = length(z); 
        alpha = 1 + n ./ sum( log(z./xmin) );
        if finite, alpha = alpha*(n-1)/n+1/n; end % Finite-size correction
        if n < 50 && ~finite && ~nowarn
            fprintf('(PLFIT) Warning: finite-size bias may be present.\n');
        end
        L = n*log((alpha-1)/xmin) - alpha.*sum(log(z./xmin));

% z is filtered to include only values greater than or equal to the best xmin.
% n is the number of values in the new z.
% alpha is recalculated using the MLE formula.
% If finite is true, a finite-size correction is applied to alpha.
% If there are fewer than 50 data points and finite-size correction is not applied, a warning is displayed (unless nowarn is true).
% L is calculated as the log-likelihood of the data above xmin under the fitted power-law model.

% Discrete Data Case (INTS)

    case 'INTS'

% Scaling Parameters and Zeta Function Calculation
        
        if isempty(vec)
            vec  = (1.50:0.01:3.50);

% If the vector of scaling parameters (vec) is empty, set it to the range 1.50 to 3.50 in increments of 0.01.
            
        end
        zvec = zeta(vec);

% Calculates the Riemann zeta function for each value in vec.

        xmins = unique(x);
        xmins = xmins(1:end-1);

% unique(x) returns unique values of x sorted in ascending order.
% xmins = xmins(1:end-1) removes the last element from xmins to avoid including the maximum value as a candidate for xmin.

% Apply Optional Parameters

        if ~isempty(xminx)
            xmins = xmins(find(xmins>=xminx,1,'first'));
        end
        if ~isempty(limit)
            limit = round(limit);
            xmins(xmins>limit) = [];
        end
        if ~isempty(sample)
            xmins = xmins(unique(round(linspace(1,length(xmins),sample))));
        end
        if isempty(xmins)
            fprintf('(PLFIT) Error: x must contain at least two unique values.\n');
            alpha = NaN; xmin = x(1); D = NaN;
            return;
        end

% These conditions modify xmins based on optional parameters (xminx, limit, sample). 

        xmax   = max(x);
 
% Determine the maximum value of x.

        dat    = zeros(length(xmins),2);

% Initialize the data array to store the KS statistics and corresponding alpha values for each candidate xmin.

        z      = x;
        fcatch = 0;

% Sort the data.

        for xm = 1:length(xmins)
            xmin = xmins(xm);
            z    = z(z>=xmin);
            n    = length(z);

 % Loop through each candidate xmin.
    % For each candidate xmin, data points greater than or equal to xmin are selected.
    % The number of data points n is determined.
    % Alpha (alpha) is estimated by maximizing the likelihood function using either a vectorized or iterative calculation method.
    % The maximum likelihood value and corresponding alpha are found.
    % The KS statistic is computed by comparing the empirical cumulative distribution function (CDF) and the fitted CDF.
    % The KS statistic and corresponding alpha are stored in dat.

            if fcatch == 0
                try

% Estimate alpha via direct maximization of likelihood function.

                    zdiff = sum( repmat((1:xmin-1)',1,length(vec)).^-repmat(vec,xmin-1,1) ,1);
                    L = -vec.*sum(log(z)) - n.*log(zvec - zdiff);

% Vectorized version of numerical calculation.

                catch

% Catch: force loop to default to iterative version for remainder of the search.

                    fcatch = 1;
                end
            end

            if fcatch == 1

% Force iterative calculation (more memory efficient, but can be slower)

                L       = -Inf*ones(size(vec));
                slogz   = sum(log(z));
                xminvec = (1:xmin-1);
                for k = 1:length(vec)
                    L(k) = -vec(k)*slogz - n*log(zvec(k) - sum(xminvec.^-vec(k)));
                end
            end

% If vectorized calculation fails, switch to iterative calculation.

            [Y,I] = max(L);

% Find the maximum likelihood value and corresponding alpha.

            fit = cumsum((((xmin:xmax).^-vec(I)))./ (zvec(I) - sum((1:xmin-1).^-vec(I))));
            cdi = cumsum(hist(z,xmin:xmax)./n);
            dat(xm,:) = [max(abs( fit - cdi )) vec(I)];
        end

% Computing KS statistic.

        [D,I] = min(dat(:,1));
        xmin  = xmins(I);
        z     = x(x>=xmin);
        n     = length(z);

% Select the index for the minimum value of the KS statistic (D).

        alpha = dat(I,2);

% Set alpha to the corresponding value for the selected xmin.

        if finite, alpha = alpha*(n-1)/n+1/n; end

% Apply finite-size correction if specified.

        if n < 50 && ~finite && ~nowarn
            fprintf('(PLFIT) Warning: finite-size bias may be present.\n');
        end

% Print a warning if finite-size bias may be present and finite correction is not applied.

        L     = -alpha*sum(log(z)) - n*log(zvec(find(vec<=alpha,1,'last')) - sum((1:xmin-1).^-alpha));

% Calculate the log-likelihood of the power-law fit.

    otherwise
        fprintf('(PLFIT) Error: x must contain only reals or only integers.\n');
        alpha = [];
        xmin  = [];
        L     = [];
        return;


end
end

% otherwise Clause: This clause handles different data types (INTS for integers and REAL for continuous data). The otherwise clause catches any cases where x contains data that is neither entirely real numbers nor integers.
% Error Message (fprintf('(PLFIT) Error: x must contain only reals or only integers.\n');): This line prints an error message to the console, indicating that the input data x must consist solely of either real numbers or integers. If the data does not meet this criterion, the function cannot proceed with fitting a power-law distribution.
% Set Output Variables to Empty:
    % alpha = [];
    % xmin  = [];
    % L     = [];
        % These lines set the output variables alpha, xmin, and L to empty arrays. This indicates that the function could not produce valid results due to the incorrect data type.


