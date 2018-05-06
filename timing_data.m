% n = 500 timing data
n_threads = [1,2,4,8,16,24,48,72,96];
t_alloc = [3.86e-04, 1.47e-04, 4.34e-04, 2.63e-04, 3.14e-04, 5.54e-04, 3.19e-05, 3.41e-05, 4.68e-04];
t_IC = [3.79e-04, 2.67e-04, 1.92e-04, 1.10e-04, 2.72e-04, 9.73e-04, 3.92e-03, 2.41e-02, 3.63e-02];
t_BC = [4.90e-04, 4.42e-04, 3.17e-04, 3.00e-04, 1.59e-04, 3.26e-04, 5.44e-04, 1.38e-03, 6.67e-04];
t_Iter = [7.04e+01, 3.70e+01,  1.68e+01, 9.14e+00, 6.28e+00, 5.80e+00, 8.02e+00, 1.06e+01, 1.58e+01];
Runtime = [7.04e+01, 3.71e+01, 1.68e+01, 9.20e+00, 6.34e+00, 5.86e+00, 8.10e+00, 1.08e+01, 1.60e+01];

% Number of data points per vector
N = numel(n_threads);

% Strong scaling speedup. Speedup = ser time / par time
speedup_alloc = t_alloc(1)./t_alloc;
speedup_IC = t_IC(1)./t_IC;
speedup_BC = t_BC(1)./t_BC;
speedup_Iter = t_Iter(1)./t_Iter;
speedup_Runtime = Runtime(1)./Runtime;

% Plot results
figure(1);
clf;
hold on;
grid on;
plot(n_threads, speedup_BC, 'bo-');  % speedup plot
plot(n_threads,n_threads,'r');          % ideal speedup plot
    xlabel('threads');
    ylabel('speedup');
    title('Strong scaling speed up of BC kernel');
    legend('Measured speedup','Ideal speedup');

% Calculate relative difference 
rel_dif_IC = (n_threads-speedup_IC)./(n_threads)
rel_dif_BC = (n_threads-speedup_BC)./(n_threads)
rel_dif_Iter = (n_threads-speedup_Iter)./(n_threads)

