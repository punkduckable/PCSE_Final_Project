% Read in results as comma separated values. We start on the 3rd line (the
% 2 argument skips the first two lines of the file).
U_k_with_BC = csvread('Results.txt',2);

% We only want to use interior elements (No BC's)
U_k = U_k_with_BC(2:(end-1), 2:(end-1));

% Now, we need the dimensions of the matrix.
[Rows, Cols] = size(U_k);

% Set up upper and lower bounds for x, y (this uses paramaters set in the
% main file).
x_max = 1;
x_min = -1;
x = linspace(x_min, x_max, Cols);

y_max = 1;
y_min = -1;
y = linspace(y_min, y_max, Rows);

% Create a meshgrid. This will be used to plot the results
[X,Y] = meshgrid(x,y);

% Plot as a surface
figure(1);
clf;
grid off;
p = surf(X,Y,U_k);
    set(p,'linestyle','none');
    xlabel('X');
    ylabel('Y');
    title(['Laplace solver with a ',num2str(Rows),'x',num2str(Cols),' grid.'])