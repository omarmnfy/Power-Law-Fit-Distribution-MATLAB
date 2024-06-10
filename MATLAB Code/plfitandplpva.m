%%

filePath = "";
data = load(filePath); % load the data from the .mat file
x = data.cropareas; % access the data points

%%

[alpha, xmin, L] = plfit(x); % call the plfit function to estimate the alpha, xmin, and L
[p, gof] = plpva(x, xmin); % call the plpva function to estimate the p-value and goo dness-of-fit
 
%%

disp(['Alpha: ', num2str(alpha)]); % disp to display answers/output.
disp(['Xmin: ', num2str(xmin)]); % disp to display answers/output.
disp(['Log-Likelihood: ', num2str(L)]); % disp to display answers/output.

disp(['p-value: ', num2str(p)]); % disp to display answers/output.
disp(['Goodness-of-fit: ', num2str(gof)]); % disp to display answers/output.


