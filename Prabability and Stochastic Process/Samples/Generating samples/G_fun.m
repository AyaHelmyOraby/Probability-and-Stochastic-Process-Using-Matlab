
rng(42);

% Number of samples
num_samples = 101; 

% Generate random samples for t in the range [0, 2]
t = linspace(0, 2, num_samples);

% Generate random samples for theta in the range [0, pi]
X = rand(100, num_samples) * pi;


save('C:/Users/HELAL/Downloads/Documents/random_process1_data.mat', 't', 'X');
