
rng(42);

% Number of samples
num_samples = 101; %  101 samples for example

% Generate random samples for t in the range [0, 2]
t = linspace(0, 2, num_samples);

% Generate random samples for A from a normal distribution N(-5, 5)
X = -5 + 10 * randn(100, num_samples);

save('C:/Users/HELAL/Downloads/Documents/random_process_data_n.mat', 't','X');
