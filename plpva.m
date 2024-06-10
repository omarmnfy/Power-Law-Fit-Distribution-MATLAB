% plpva: P-value and Goodness of Fit Calculator the Given Power-Law Fit
% PLPVA(x, xmin) takes data x and given lower cutoff for the power-law behavior xmin and computes the corresponding p-value for the Kolmogorov-Smirnov test, according to the method described in  Clauset, Shalizi, Newman (2007). PLPVA automatically detects whether x is composed of real or integer values, and applies the appropriate method. For discrete data, if min(x) > 1000, PLPVA uses the continuous approximation, which is  a reliable in this regime.
% Original Code Source: http://www.santafe.edu/~aaronc/powerlaws/


    % Varagin is an input variable in a function definition statement that enables the function to accept any number of input arguments.

function [p, gof]=plpva(x, xmin, varargin)

% Initialization of Variables

vec    = [];
sample = [];
limit  = [];
xminx  = [];
Bt     = [];
quiet  = false;

persistent rand_state;

    % Initializes variables that will be used later in the function. rand_state is declared as persistent to retain its value across function calls.
    % The persistent keyword in MATLAB is used to define variables that retain their values between function calls. Unlike regular variables that are reinitialized each time a function is called, persistent variables maintain their state.
    % rand_state is a persistent variable used to store the state of the random number generator across multiple calls to the plpva function. This ensures that the random number generator can continue from where it left off, providing continuity and consistency in the random numbers generated during bootstrapping. 

% Parsing Optional Arguments

   % Initializes a counter i to 1, which will be used to iterate through the varargin cell array containing optional parameters passed to the function.

i=1; 

   % Starts a while loop that will run as long as i is less than or equal to the number of elements in varargin.

while i<=length(varargin), 

   % argok is a flag variable set to 1 (true) initially. It will be used to check if the current argument is valid.
   % A flag variable is a variable used to signal whether a particular condition has been met or not.
   % As the loop processes each argument in varargin, it uses a switch statement to check if the argument matches any of the expected parameter names ('range', 'sample', 'limit', 'xmin', 'reps', 'silent'). If the argument matches one of these cases, the corresponding action is taken (e.g., assigning a value or setting a boolean flag).
   % If the argument does not match any of the expected cases (handled by the otherwise part of the switch statement), argok is set to 0. This indicates that the current argument is invalid.

  argok = 1;

   % This switch statement checks the value of varargin{i} and assigns the appropriate value to the corresponding parameter.
        % 'range': Sets vec to the next element varargin{i+1} and increments i by 1.
        % 'sample': Sets sample to the next element varargin{i+1} and increments i by 1.
        % 'limit': Sets limit to the next element varargin{i+1} and increments i by 1.
        % 'xmin': Sets xminx to the next element varargin{i+1} and increments i by 1.
        % 'reps': Sets Bt to the next element varargin{i+1} and increments i by 1.
        % 'silent': Sets the boolean quiet to true, which will suppress output messages.
        % otherwise: If the argument does not match any of the specified cases, sets argok to 0 (false), indicating an invalid argument.

  if ischar(varargin{i}), 
    switch varargin{i},
        case 'range',        vec    = varargin{i+1}; i = i + 1;
        case 'sample',       sample = varargin{i+1}; i = i + 1;
        case 'limit',        limit  = varargin{i+1}; i = i + 1;
        case 'xmin',         xminx  = varargin{i+1}; i = i + 1;
        case 'reps',         Bt     = varargin{i+1}; i = i + 1;
        case 'silent',       quiet  = true;
        otherwise, argok=0; 
    end
  end

    % If the argument is invalid (argok is 0), it displays a message and skips the invalid argument.

  if ~argok, 
    disp(['(PLPVA) Ignoring invalid argument #' num2str(i+1)]); 
  end

    % Increments i by 1 to process the next argument in varargin.

  i = i+1; 
end

% Argument Validations

if ~isempty(vec) && (~isvector(vec) || min(vec)<=1),
	fprintf('(PLPVA) Error: ''range'' argument must contain a vector; using default.\n');
    vec = [];
end;
if ~isempty(sample) && (~isscalar(sample) || sample<2),
	fprintf('(PLPVA) Error: ''sample'' argument must be a positive integer > 1; using default.\n');
    sample = [];
end;
if ~isempty(limit) && (~isscalar(limit) || limit<1),
	fprintf('(PLPVA) Error: ''limit'' argument must be a positive value >= 1; using default.\n');
    limit = [];
end;
if ~isempty(Bt) && (~isscalar(Bt) || Bt<2),
	fprintf('(PLPVA) Error: ''reps'' argument must be a positive value > 1; using default.\n');
    Bt = [];
end;
if ~isempty(xminx) && (~isscalar(xminx) || xminx>=max(x)),
	fprintf('(PLPVA) Error: ''xmin'' argument must be a positive value < max(x); using default behavior.\n');
    xminx = [];
end;

% Reshape Input Vector

    % Ensures x is a column vector.

x = reshape(x,numel(x),1);

% Data Type Detection (Continuious (REAL) or Discrete(INTS))

    % If all elements in x are integers (i.e., there are no differences between x and its floor values), set f_dattype to 'INTS'.

if     isempty(setdiff(x,floor(x))), f_dattype = 'INTS';

    % If x contains real numbers (not just integers), set f_dattype to 'REAL'.

elseif isreal(x),    f_dattype = 'REAL';

    % If x contains any other type of data, set f_dattype to 'UNKN'.
    
else                 f_dattype = 'UNKN';
end;

    % Further check to see if the data in x, initially identified as integers, should instead be treated as real numbers:
        % If f_dattype is 'INTS', and the minimum value in x is greater than 1000, and the length of x is greater than 100, 
        % set f_dattype to 'REAL'. This is based on the assumption that for large integer values, the continuous approximation is reliable.

if strcmp(f_dattype,'INTS') && min(x) > 1000 && length(x)>100,
    f_dattype = 'REAL';
end;

    % N is the number of elements in x.
    % x is reshaped to be a column vector (if it isn't already).

N = length(x);
x = reshape(x,N,1); % guarantees x is a column vector

% Random State Initialization
    
    % Initializes the random number generator if it hasn't been initialized already. This ensures reproducibility in random number generation.

if isempty(rand_state)
    rand_state = cputime;
    rand('twister',sum(100*clock));
end;

% Set Default Number of Bootstraps (Bt)

    % If Bt is not provided, set it to 1000.
    % nof is an array to store the goodness-of-fit statistics from the bootstrap samples.

if isempty(Bt), Bt = 1000; end;
nof = zeros(Bt,1);

% Print Initial Information

    % Prints information about the calculation if not in quiet mode.

if ~quiet,
    fprintf('Power-law Distribution, p-value calculation\n');
    fprintf('   Copyright 2007-2010 Aaron Clauset\n');
    fprintf('   Warning: This can be a slow calculation; please be patient.\n');
    fprintf('   n    = %i\n   xmin = %6.4f\n   reps = %i\n',length(x),xmin,length(nof));
end;
tic;

% Switch Based on Data Type to Estimate xmin and alpha Accordingly

    % The next steps are executed differently based on whether f_dattype is 'REAL' or 'INTS'.

switch f_dattype,

% Empirical Distribution

    case 'REAL',
        
% Compute D for the Empirical Distribution

    % z = x(x >= xmin); selects data points greater than or equal to xmin.
    % nz = length(z); counts these data points.
    % y = x(x < xmin); selects data points less than xmin.
    % ny = length(y); counts these data points.
    % alpha = 1 + nz / sum(log(z / xmin)); calculates the scaling parameter.
    % cz = (0:nz-1)' / nz; calculates the empirical cumulative distribution function (CDF).
    % cf = 1 - (xmin ./ sort(z)).^(alpha - 1); calculates the fitted CDF.
    % gof = max(abs(cz - cf)); calculates the goodness-of-fit statistic.
    % pz = nz / N; calculates the proportion of data points greater than or equal to xmin.

        z     = x(x>=xmin);	nz   = length(z);
        y     = x(x<xmin); 	ny   = length(y);
        alpha = 1 + nz ./ sum( log(z./xmin) );
        cz    = (0:nz-1)'./nz;
        cf    = 1-(xmin./sort(z)).^(alpha-1);
        gof   = max( abs(cz - cf) );
        pz    = nz/N;

% Compute Distribution of Goodness of Fits from Semi-Parametric Bootstrap

    % Repeats the following steps Bt times (1000 times by default).

        % n1 = sum(rand(N, 1) > pz); calculates the number of points to sample from y.
        % q1 = y(ceil(ny * rand(n1, 1))); samples n1 points from y.
        % n2 = N - n1; calculates the number of points to generate from the power-law distribution.
        % q2 = xmin * (1 - rand(n2, 1)).^(-1 / (alpha - 1)); generates n2 points from the power-law distribution.
        % q = sort([q1; q2]); combines and sorts these points.
        % qmins = unique(q); identifies unique potential values for xmin.
        % Filters qmins based on xminx, limit, and sample parameters.
        
        for B=1:length(nof) % length of entire data set with fit

% Semi-Parametric Bootstrap of Data

    % n1 = sum(rand(N, 1) > pz); calculates the number of points to sample from y.
    % q1 = y(ceil(ny * rand(n1, 1))); samples n1 points from y.
    % n2 = N - n1; calculates the number of points to generate from the power-law distribution.
    % Generates n2 points from the power-law distribution using a discrete zeta generator.
    % q = [q1; q2]; combines these points.

            n1 = sum(rand(N,1)>pz);
            q1 = y(ceil(ny.*rand(n1,1)));
            n2 = N-n1;
            q2 = xmin*(1-rand(n2,1)).^(-1/(alpha-1));
            q  = sort([q1; q2]);

 % Estimate xmin and alpha via GoF-Method

    % qmins = unique(q); identifies unique potential values for xmin.
    % Filters qmins based on xminx, limit, and sample parameters.

            qmins = unique(q);
            qmins = qmins(1:end-1);

            if ~isempty(xminx),
                qmins = qmins(find(qmins>=xminx,1,'first'));
            end;

            if ~isempty(limit),
                qmins(qmins>limit) = [];
                if isempty(qmins), qmins = min(q); end;
            end;

            if ~isempty(sample),
                qmins = qmins(unique(round(linspace(1,length(qmins),sample))));
            end;

    % For each potential xmin in qmins:

        % zq = q(q >= qmin); selects points greater than or equal to qmin.
        % nq = length(zq); counts these points.
        % Estimates alpha and calculates the fitted CDF.
        % dat(qm) = max(abs(fit - cdi)); calculates the GoF statistic. 

            dat   = zeros(size(qmins));
            for qm=1:length(qmins)
                  qmin = qmins(qm);
                  zq   = q(q>=qmin);
                  nq   = length(zq);
                  a    = nq ./ sum( log(zq./qmin) );
                  cq   = (0:nq-1)'./nq;
                  cf   = 1-(qmin./zq).^a;
                  dat(qm) = max( abs(cq - cf) );
            end;

            if ~quiet,
                fprintf('[%i]\tp = %6.4f\t[%4.2fm]\n',B,sum(nof(1:B)>=gof)./B,toc/60);
            end;

% Store Distribution of Estimated Goodness of Fit Values

    % nof(B) = min(dat); stores the minimum GoF statistic.
    % p = sum(nof >= gof) / length(nof); calculates the p-value as the proportion of bootstrap GoF statistics that are greater than or equal to the empirical GoF statistic.

            nof(B) = min(dat);
        end;
        p = sum(nof>=gof)./length(nof);

    case 'INTS',

        if isempty(vec),
            vec  = (1.50:0.01:3.50);    % covers range of most practical 
        end;                            % scaling parameters

% Parameter Initialization

    % vec is initialized if empty, covering a range of scaling parameters.
    % zvec is the zeta function evaluated at each value in vec.

% Empirical Distribution

        % z = x(x >= xmin); selects data points greater than or equal to xmin.
        % nz = length(z); counts these data points.
        % y = x(x < xmin); selects data points less than xmin.
        % ny = length(y); counts these data points.
        % L is initialized to negative infinity for each value in vec.
        % For each value in vec, the log-likelihood is computed, and the alpha corresponding to the maximum log-likelihood is selected.
        % fit is the cumulative distribution function for the fitted power-law.
        % cdi is the cumulative distribution of the empirical data.
        % gof is the goodness-of-fit statistic.
        % pz is the proportion of data points greater than or equal to xmin.

        zvec = zeta(vec);

 % Compute D for the Empirical Distribution

        z     = x(x>=xmin);	nz   = length(z);	xmax = max(z);
        y     = x(x<xmin); 	ny   = length(y);
        L  = -Inf*ones(size(vec));

        for k=1:length(vec)
            L(k) = -vec(k)*sum(log(z)) - nz*log(zvec(k) - sum((1:xmin-1).^-vec(k)));
        end

        [Y,I] = max(L);
        alpha = vec(I);

        fit = cumsum((((xmin:xmax).^-alpha))./ (zvec(I) - sum((1:xmin-1).^-alpha)));
        cdi = cumsum(hist(z,(xmin:xmax))./nz);
        gof = max(abs( fit - cdi ));
        pz  = nz/N;

        mmax = 20*xmax;
        pdf = [zeros(xmin-1,1); (((xmin:mmax).^-alpha))'./ (zvec(I) - sum((1:xmin-1).^-alpha))];
        cdf = [(1:mmax+1)' [cumsum(pdf); 1]];

% Compute Distribution of Goodness of Fits from Semi-Parametric Bootstrap

    % Repeats the following steps Bt times (1000 times by default).

        % n1 = sum(rand(N, 1) > pz); calculates the number of points to sample from y.
        % q1 = y(ceil(ny * rand(n1, 1))); samples n1 points from y.
        % n2 = N - n1; calculates the number of points to generate from the power-law distribution.
        % Generates n2 points from the power-law distribution using a discrete zeta generator.
        % q = [q1; q2]; combines these points.
        % qmins = unique(q); identifies unique potential values for xmin.
        % Filters qmins based on xminx, limit, and sample parameters.

        for B=1:length(nof) % length of entire data set with fit
 
% Semi-Parametric Bootstrap of Data

            n1 = sum(rand(N,1)>pz);
            q1 = y(ceil(ny.*rand(n1,1)));
            n2 = N-n1;

 % Simple Discrete Zeta Generator

            r2 = sort(rand(n2,1));  c = 1;
            q2 = zeros(n2,1);	    k = 1;

            for i=xmin:mmax+1
                while c<=length(r2) && r2(c)<=cdf(i,2), c=c+1; end;
                q2(k:c-1) = i;
                k = c;
                if k>n2, break; end;
            end;

            q = [q1; q2];

% Estimate xmin and alpha via GoF-Method

            qmins = unique(q);
            qmins = qmins(1:end-1);

            if ~isempty(xminx),
                qmins = qmins(find(qmins>=xminx,1,'first'));
            end;

            if ~isempty(limit),
                qmins(qmins>limit) = [];
                if isempty(qmins), qmins = min(q); end;
            end;

            if ~isempty(sample),
                qmins = qmins(unique(round(linspace(1,length(qmins),sample))));
            end;

% For each potential xmin in qmins:

        % zq = q(q >= qmin); selects points greater than or equal to qmin.
        % nq = length(zq); counts these points.
        % Estimates alpha and calculates the fitted CDF.
        % dat(qm) = max(abs(fit - cdi)); calculates the GoF statistic.

            dat   = zeros(size(qmins));
            qmax  = max(q); zq = q;

            for qm=1:length(qmins)
                qmin = qmins(qm);
                zq   = zq(zq>=qmin);
                nq   = length(zq);

                if nq>1
                    try
                        % vectorized version of numerical calculation
                        zdiff = sum( repmat((1:qmin-1)',1,length(vec)).^-repmat(vec,qmin-1,1) ,1);
                        L = -vec.*sum(log(zq)) - nq.*log(zvec - zdiff);
                    catch
                       % iterative version (more memory efficient, but slower)
                       L       = -Inf*ones(size(vec));
                       slogzq  = sum(log(zq));
                       qminvec = (1:qmin-1);
                       for k=1:length(vec)
                           L(k) = -vec(k)*slogzq - nq*log(zvec(k) - sum(qminvec.^-vec(k)));
                       end;
                    end;
                    [Y,I] = max(L);

                    fit = cumsum((((qmin:qmax).^-vec(I)))./ (zvec(I) - sum((1:qmin-1).^-vec(I))));
                    cdi = cumsum(hist(zq,(qmin:qmax))./nq);
                    dat(qm) = max(abs( fit - cdi ));
                else
                    dat(qm) = -Inf;
                end;

            end
            if ~quiet,
                fprintf('[%i]\tp = %6.4f\t[%4.2fm]\n',B,sum(nof(1:B)>=gof)./B,toc/60);
            end;

% Store Distribution of Estimated Goodness of Fit Values

        % nof(B) = min(dat); stores the minimum GoF statistic.

            nof(B) = min(dat);
        end;

        % p = sum(nof >= gof) / length(nof); calculates the p-value as the proportion of bootstrap GoF statistics that are greater than or equal to the empirical GoF statistic.

        p = sum(nof>=gof)./length(nof);

    % otherwise Clause: This clause handles different data types (INTS for integers and REAL for continuous data). The otherwise clause catches any cases where x contains data that is neither entirely real numbers nor integers.
    % Error Message (fprintf('(PLFIT) Error: x must contain only reals or only integers.\n');): This line prints an error message to the console, indicating that the input data x must consist solely of either real numbers or integers. If the data does not meet this criterion, the function cannot proceed with fitting a power-law distribution.
    % Set Output Variables to Empty:
        % p   = [];
        % gof = [];
            % These lines set the output variables p, and gof to empty arrays. This indicates that the function could not produce valid results due to the incorrect data type.

    otherwise,
        fprintf('(PLPVA) Error: x must contain only reals or only integers.\n');
        p   = [];
        gof = [];
        return;
end;