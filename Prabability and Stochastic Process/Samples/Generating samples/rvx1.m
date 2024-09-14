

% Number of samples
num_samples = 10000;

% Specify custom folder location
folder_path = 'C:\Users\HELAL\Documents\rv\';  %my PC :)

% Generate random variable X following U(-3, 5)
X = -3 + 8 * rand(1, num_samples);





data_file_path = fullfile(folder_path, 'random_variable_datax1.mat');
save(data_file_path, 'X');

