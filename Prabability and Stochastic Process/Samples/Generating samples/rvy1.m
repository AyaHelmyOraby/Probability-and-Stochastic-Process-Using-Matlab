

% Number of samples
num_samples = 10000;


folder_path = 'C:\Users\HELAL\Documents\rv\';

% Generate random variable Y following N(-8, 4)
mu_Y = -8;
sigma_Y = 2;  % Standard deviation, not variance
Y = mu_Y + sigma_Y * randn(1, num_samples);




data_file_path_Y = fullfile(folder_path, 'random_variable_data_Y.mat');
save(data_file_path_Y, 'Y');

