classdef appprob_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        RandomVariableTab              matlab.ui.container.Tab
        EditField_12                   matlab.ui.control.EditField
        EditField_11                   matlab.ui.control.EditField
        SecondMomentatt0Button         matlab.ui.control.Button
        FirstMomentatt0Button          matlab.ui.control.Button
        PlotthefirstandsecondderivitiveButton  matlab.ui.control.Button
        PlotMGFButton                  matlab.ui.control.Button
        EditField_4                    matlab.ui.control.EditField
        EditField_3                    matlab.ui.control.EditField
        TheThirdmomentButton           matlab.ui.control.Button
        VarianceButton                 matlab.ui.control.Button
        EditField_2                    matlab.ui.control.EditField
        MeanButton                     matlab.ui.control.Button
        ImportButton                   matlab.ui.control.Button
        FileNameEditField              matlab.ui.control.EditField
        FileNameEditFieldLabel         matlab.ui.control.Label
        RandomProcessTab_2             matlab.ui.container.Tab
        ComparingButton                matlab.ui.control.Button
        TimeACFEditField_2             matlab.ui.control.EditField
        TimeACFEditField_2Label        matlab.ui.control.Label
        StatisticalACFEditField        matlab.ui.control.EditField
        StatisticalACFEditFieldLabel   matlab.ui.control.Label
        TimeMeanEditField              matlab.ui.control.EditField
        TimeMeanEditFieldLabel         matlab.ui.control.Label
        StatisticalMeanEditField       matlab.ui.control.EditField
        StatisticalMeanEditFieldLabel  matlab.ui.control.Label
        DACFButton                     matlab.ui.control.Button
        EditField_13                   matlab.ui.control.EditField
        ensamplemeanofMButton          matlab.ui.control.Button
        EmsampleplotButton             matlab.ui.control.Button
        EditField_10                   matlab.ui.control.EditField
        AveragePowerButton             matlab.ui.control.Button
        EditField_9                    matlab.ui.control.EditField
        BSDValueButton                 matlab.ui.control.Button
        EditField_8                    matlab.ui.control.EditField
        ACFValueButton                 matlab.ui.control.Button
        EditField_7                    matlab.ui.control.EditField
        emsampemeanvalueButton         matlab.ui.control.Button
        BSDPlotButton                  matlab.ui.control.Button
        EditField_6                    matlab.ui.control.EditField
        TofACFButton                   matlab.ui.control.Button
        EditField_5                    matlab.ui.control.EditField
        TimeofNButton                  matlab.ui.control.Button
        ACFButton                      matlab.ui.control.Button
        TheensamplemeanoftheprocessButton  matlab.ui.control.Button
        PlotMButton                    matlab.ui.control.Button
        ImportButton_2                 matlab.ui.control.Button
        nthoftheProcessEditField       matlab.ui.control.EditField
        nthoftheProcessEditFieldLabel  matlab.ui.control.Label
        FileEditField                  matlab.ui.control.EditField
        FileEditFieldLabel             matlab.ui.control.Label
        EnterNumberofsamplesEditField  matlab.ui.control.EditField
        EnterNumberofsamplesEditFieldLabel  matlab.ui.control.Label
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Callback function
        function UploadthefileButtonPushed(app, event)
       
        
        end

        % Button pushed function: ImportButton
        function ImportButtonPushed(app, event)

[filename, path] = uigetfile();
if isequal(filename, 0)
    % User canceled the operation
    return;
end
fullFilePath = fullfile(path, filename);
app.FileNameEditField.Value = fullFilePath;
figure(app.UIFigure);
try
    % Load all variables from the MAT file
    loadedData = load(fullFilePath);
    
    % Check if the loadedData structure is not empty
    if ~isempty(loadedData)
        % Display the loaded variables
        disp('Variables loaded successfully:');
        variableNames = fieldnames(loadedData);
        disp(variableNames);
        
        % Assign the variables to the base workspace
        for i = 1:numel(variableNames)
            assignin('base', variableNames{i}, loadedData.(variableNames{i}));
        end
        
        disp('Variables assigned to the workspace.');
    else
        disp('Error: No variables found in the loaded data.');
    end

catch
    disp('Error loading data from the MAT file.');
end



        end

        % Value changed function: FileNameEditField
        function FileNameEditFieldValueChanged(app, event)
            value = app.FileNameEditField.Value;
            
        end

        % Button pushed function: MeanButton
        function MeanButtonPushed(app, event)
% Button pushed function: MeanButton

    %{
        % Check if there's data loaded in the workspace
    if ~exist('base', 'var')
        disp('Error: No data loaded. Load data before calculating mean.');
        return;
    end
    
    % Calculate the mean of the loaded variables
    meanValue = mean(base(:));
    
    % Update the EditField_4 with the calculated mean
    app.EditField_4.Value = num2str(meanValue);
            %}

    % Check if there are any variables loaded in the workspace
  
    % Check if there are any variables loaded in the workspace
    workspaceVariables = evalin('base', 'who');
    
    if isempty(workspaceVariables)
        disp('Error: No data loaded. Load data before calculating mean.');
        return;
    end
    
    % Calculate the mean of all numeric variables in the workspace
    allData = [];
    
    for i = 1:numel(workspaceVariables)
        currentVar = evalin('base', workspaceVariables{i});
        
        % Check if the variable is numeric
        if isnumeric(currentVar)
            allData = [allData; currentVar(:)];
        end
    end
    
    if ~isempty(allData)
        meanValue = mean(allData);
        app.EditField_4.Value = num2str(meanValue);
    else
        disp('Error: No numeric data loaded. Load numeric data before calculating mean.');
    end


    
        end

        % Value changed function: EditField_4
        function EditField_4ValueChanged(app, event)
            value = app.EditField_4.Value;
            
        end

        % Button pushed function: VarianceButton
        function VarianceButtonPushed(app, event)

    % Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating variance.');
    return;
end

% Calculate the mean of all numeric variables in the workspace
allData = [];

for i = 1:numel(workspaceVariables)
    currentVar = evalin('base', workspaceVariables{i});

    % Check if the variable is numeric
    if isnumeric(currentVar)
        allData = [allData; currentVar(:)];
    end
end

if ~isempty(allData)
    % Calculate mean 
    meanValue = sum(allData) / numel(allData);

    % Update the EditField_4 with the calculated mean
    app.EditField_4.Value = num2str(meanValue);

    % Calculate the variance using the mean
    varianceValue = sum((allData - meanValue).^2) / numel(allData);
    app.EditField_2.Value = num2str(varianceValue);
else
    disp('Error: No numeric data loaded. Load numeric data before calculating variance.');
end


        end

        % Value changed function: EditField_2
        function EditField_2ValueChanged(app, event)
            value = app.EditField_2.Value;
            
        end

        % Button pushed function: TheThirdmomentButton
        function TheThirdmomentButtonPushed(app, event)
             % Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating the third moment.');
    return;
end

% Calculate the third moment of all numeric variables in the workspace
allData = [];

for i = 1:numel(workspaceVariables)
    currentVar = evalin('base', workspaceVariables{i});
    
    % Check if the variable is numeric
    if isnumeric(currentVar)
        allData = [allData; currentVar(:)];
    end
end

if ~isempty(allData)
    % Calculate the third moment
    meanValue = mean(allData);
    thirdMomentValue = (sum((allData - meanValue).^3) / numel(allData));  % third moment is Ex^3

    % Display the result
     app.EditField_3.Value = num2str(thirdMomentValue);

  
else
    disp('Error: No numeric data loaded. Load numeric data before calculating the third moment.');
end


        end

        % Value changed function: EditField_3
        function EditField_3ValueChanged(app, event)
            value = app.EditField_3.Value;
            
        end

        % Button pushed function: PlotMGFButton
        function PlotMGFButtonPushed(app, event)
% Get the mean and variance from the GUI 
meanValue = str2double(get(app.EditField_4, 'Value'));  
varianceValue = str2double(get(app.EditField_2, 'Value'));  

% Check if the inputs are valid
if isnan(meanValue) || isnan(varianceValue)
    disp('Error: Please enter valid mean and variance values.');
    return;
end

% Calculate the MGF over the specified range of t
t_values = 0:0.01:2;  % Define the range of t values
MGF_values = exp(t_values .* meanValue + 0.5 * (t_values.^2) .* varianceValue);

% Plot the MGF
figure;
plot(t_values, MGF_values, 'LineWidth', 2);
title('Moment Generating Function (MGF) vs t');
xlabel('t');
ylabel('M(t)');
grid on;


        end

        % Button pushed function: PlotthefirstandsecondderivitiveButton
        function PlotthefirstandsecondderivitiveButtonPushed(app, event)

meanValue = str2double(app.EditField_4.Value);
varianceValue = str2double(app.EditField_2.Value);

% Check if the inputs are valid
if isnan(meanValue) || isnan(varianceValue)
    disp('Error: Please enter valid mean and variance values.');
    return;
end

% Symbolic variables
syms t;

% Define the MGF
MGF = exp(t * meanValue + 0.5 * t^2 * varianceValue);

% Calculate the first and second derivatives
firstDerivative = diff(MGF, t);
secondDerivative = diff(MGF, t, 2);  % Calculate the second derivative directly

% Calculate the derivatives over the specified range of t
t_values = 0:0.01:2;  % Define the range of t values
firstDerivativeValues = double(subs(firstDerivative, t, t_values));
secondDerivativeValues = double(subs(secondDerivative, t, t_values));

% Plot the first and second derivatives
figure;
subplot(2, 1, 1);
plot(t_values, firstDerivativeValues, 'LineWidth', 2);
title('First Derivative of MGF vs t');
xlabel('t');
ylabel('M''(t)');
grid on;

subplot(2, 1, 2);
plot(t_values, secondDerivativeValues, 'LineWidth', 2);
title('Second Derivative of MGF vs t');
xlabel('t');
ylabel('M''''(t)');
grid on;

        end

        % Callback function
        function ImportButton_2Pushed(app, event)
            
[filenameP, path] = uigetfile();
if isequal(filenameP, 0)
    % User canceled the operation
    return;
end

fullFilePathP = fullfile(path, filenameP);
app.FileEditField.Value = fullFilePathP;
figure(app.UIFigure);

try
    % Load all variables from the MAT file
    loadedData = load(fullFilePathP);
    
    % Check if the loadedData structure is not empty
    if ~isempty(loadedData)
        % Display the loaded variables
        disp('Variables loaded successfully:');
        variableNames = fieldnames(loadedData);
        disp(variableNames);
        
        % Assign the variables to the base workspace
        for i = 1:numel(variableNames)
            assignin('base', variableNames{i}, loadedData.(variableNames{i}));
        end
        
        disp('Variables assigned to the workspace.');
    else
        disp('Error: No variables found in the loaded data.');
    end

catch
    disp('Error loading data from the MAT file.');
end

        end

        % Callback function
        function FileEditFieldValueChanged(app, event)
            value = app.FileEditField.Value;
            
        end

        % Callback function
        function PlotMButtonPushed(app, event)

% Get the number of samples from the edit field
M = str2double(app.EnterNumberofsamplesEditField.Value);

% Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before plotting sample functions.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');
    
    % Ensure the number of samples 'M' is not greater than the number of columns in 'X'
    if M > size(X, 2)
        disp('Error: Not enough samples in the data.');
        return;
    end
    
    % Randomly select 'M' unique indices
    selectedIndices = randperm(size(X, 2), M);
    
    % Plot selected sample functions
    figure;
    hold on;
    for i = 1:M
        % Get the current sample
        currentSample = X(:, selectedIndices(i));
        
        % Ensure that the sample length matches the length of 't'
        if length(currentSample) ~= length(t)
           
            currentSample = interp1(linspace(0, 1, length(currentSample)), currentSample, linspace(0, 1, length(t)));
        end
        
        % Plot the adjusted sample function against the common time vector 't'
        plot(t, currentSample, 'DisplayName', ['Sample Function ', num2str(selectedIndices(i))]);
    end
    
    title('Randomly Selected Sample Functions');
    xlabel('Time (t)');
    ylabel('Amplitude (X)');
    legend('show');
    hold off;
else
    disp('Error: The loaded data does not contain ''t'' and ''X''.');
end


        end

        % Callback function
        function TheensamplemeanoftheprocessButtonPushed(app, event)

    % Check if there are any variables loaded in the workspace
    workspaceVariables = evalin('base', 'who');

    if isempty(workspaceVariables)
        disp('Error: No data loaded. Load data before calculating the ensemble mean.');
        return;
    end

    % Check if the loaded data contains 't' and 'X'
    if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
        % Retrieve 't' and 'X' from the base workspace
        t = evalin('base', 't');
        X = evalin('base', 'X');

        % Calculate the ensemble mean
        ensembleMean = mean(X, 1);  % Calculate mean along columns

        % Plot the ensemble mean
        figure;
        plot(t, ensembleMean, 'LineWidth', 2);
        title('Ensemble Mean of the Process');
        xlabel('Time (t)');
        ylabel('Ensemble Mean Amplitude');
        grid on;
    else
        disp('Error: The loaded data does not contain ''t'' and ''X''.');
    end



        end

        % Callback function
        function ACFButtonPushed(app, event)

    

        end

        % Callback function
        function TimeofNButtonPushed(app, event)


        end

        % Callback function
        function TIMEOFACFEditFieldValueChanged(app, event)
            value = app.TIMEOFACFEditField.Value;
            
        end

        % Callback function
        function TofACFButtonPushed(app, event)

n = str2double(app.nthoftheProcessEditField.Value);

% Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating the time auto-correlation function.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');
    
    % Ensure n is a valid index
    if isnan(n) || n < 1 || n > size(X, 2)
        disp('Error: Invalid value for n.');
        return;
    end
    
    % Get the n-th sample function
    nthSample = X(:, n);
    
    % Calculate the time auto-correlation function
    acf = xcorr(nthSample, nthSample, 'coeff');
    
    % Extract the relevant information for display
    acfValue = acf(length(t));  % Assuming length(t) corresponds to the relevant value
    
    % Display the time auto-correlation function 
    app.EditField_6.Value = num2str(acfValue);
else
    disp('Error: The loaded data does not contain ''t'' and ''X''.');
end


        end

        % Callback function
        function BSDButtonPushed(app, event)
  % Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating the power spectral density.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');

    % Compute the power spectral density using FFT
    fs = 1 / (t(2) - t(1));  % Sampling frequency (assuming uniform sampling)
    n = length(X);
    f = (0:n-1) * fs / n;
    psd = abs(fft(X)).^2 / (n * fs);

    % Plot the power spectral density
    figure;
    plot(f, 10*log10(psd), 'LineWidth', 2);
    title('Power Spectral Density of the Process');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    grid on;
else
    disp('Error: The loaded data does not contain ''t'' and ''X''.');
end


        end

        % Callback function
        function emsampemeanvalueButtonPushed(app, event)

    % Check if there are any variables loaded in the workspace
    workspaceVariables = evalin('base', 'who');

    if isempty(workspaceVariables)
        disp('Error: No data loaded. Load data before calculating the ensemble mean.');
        return;
    end

    % Check if the loaded data contains 't' and 'X'
    if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
        % Retrieve 't' and 'X' from the base workspace
        t = evalin('base', 't');
        X = evalin('base', 'X');

        % Calculate the ensemble mean
        ensembleMean = mean(X, 1);  % Calculate mean along columns

        % Calculate the mean value to display
        meanValue = mean(ensembleMean);

        % Display the mean value 
        app.EditField_7.Value = num2str(meanValue);
    else
        disp('Error: The loaded data does not contain ''t'' and ''X''.');
    end
        
        end

        % Button pushed function: ACFValueButton
        function ACFValueButtonPushed(app, event)

    % Check if there are any variables loaded in the workspace
    workspaceVariables = evalin('base', 'who');

    if isempty(workspaceVariables)
        disp('Error: No data loaded. Load data before calculating the ACF.');
        return;
    end

    % Check if the loaded data contains 't' and 'X'
    if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
        % Retrieve 't' and 'X' from the base workspace
        t = evalin('base', 't');
        X = evalin('base', 'X');

        % Calculate the statistical auto-correlation function (ACF)
        acf = xcorr(X, 'coeff');

        % Find the midpoint of the ACF
        midPoint = ceil(length(acf) / 2);

        % Display the ACF value at the midpoint 
        app.EditField_8.Value = num2str(acf(midPoint));
    else
        disp('Error: The loaded data does not contain ''t'' and ''X''.');
    end


        end

        % Value changed function: EditField_8
        function EditField_8ValueChanged(app, event)
            value = app.EditField_8.Value;
            
        end

        % Button pushed function: ImportButton_2
        function ImportButton_2Pushed2(app, event)
            [filename, path] = uigetfile();
if isequal(filename, 0)
    % User canceled the operation
    return;
end

fullFilePath = fullfile(path, filename);
app.FileEditField.Value = fullFilePath;
figure(app.UIFigure);

try
    % Load all variables from the MAT file
    loadedData = load(fullFilePath);
    
    % Check if the loadedData structure is not empty
    if ~isempty(loadedData)
        % Display the loaded variables
        disp('Variables loaded successfully:');
        variableNames = fieldnames(loadedData);
        disp(variableNames);
        
        % Assign the variables to the base workspace
        for i = 1:numel(variableNames)
            assignin('base', variableNames{i}, loadedData.(variableNames{i}));
        end
        
        disp('Variables assigned to the workspace.');
    else
        disp('Error: No variables found in the loaded data.');
    end

catch
    disp('Error loading data from the MAT file.');
end

        end

        % Button pushed function: BSDPlotButton
        function BSDPlotButtonPushed2(app, event)
% Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating and plotting the power spectral density.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');

    % Calculate the power spectral density (PSD)
    Fs = 1 / (t(2) - t(1));  % Sampling frequency
    N = length(X);
    f = linspace(0, Fs/2, N/2 + 1);  % Frequency vector up to Nyquist frequency
    Y = fft(X);
    P = abs(Y(1:N/2 + 1)).^2 / (Fs * N);

    % Plot the power spectral density
    figure;
    plot(f, sqrt(P), 'LineWidth', 2);  % Square root to get amplitude
    title('Power Spectral Density of the Process');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    grid on;
else
    disp('Error: The loaded data does not contain ''t'' and ''X''.');
end



        end

        % Button pushed function: BSDValueButton
        function BSDValueButtonPushed(app, event)
% Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating the power spectral density.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');

    % Calculate the power spectral density (PSD)
    Fs = 1 / (t(2) - t(1));  % Sampling frequency
    N = length(X);
    Y = fft(X);   
    P = abs(Y(1:N/2 + 1)).^2 / (Fs * N);

    % Calculate the mean value of the PSD
    meanPSD = mean(P);

    % Display the mean PSD value 
    app.EditField_9.Value = num2str(meanPSD);
else
    disp('Error: The loaded data does not contain ''t'' and ''X''.');
end




        end

        % Value changed function: EditField_9
        function EditField_9ValueChanged(app, event)
            value = app.EditField_9.Value;
            
        end

        % Button pushed function: emsampemeanvalueButton
        function emsampemeanvalueButtonPushed2(app, event)
            % Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating the ensemble mean.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');

    % Calculate the ensemble mean
    ensembleMean = sum(X, 1) / size(X, 1);  % Calculate mean along columns

    % Calculate the mean value to display
    meanValue = sum(ensembleMean) / numel(ensembleMean);

    % Display the mean value 
    app.EditField_7.Value = num2str(meanValue);
else
    disp('Error: The loaded data does not contain ''t'' and ''X''.');
end

        end

        % Value changed function: EditField_10
        function EditField_10ValueChanged(app, event)
            value = app.EditField_10.Value;
            
        end

        % Button pushed function: AveragePowerButton
        function AveragePowerButtonPushed(app, event)
      % Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating the total average power.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');

    % Calculate the power spectral density (PSD)
    Fs = 1 / (t(2) - t(1));  % Sampling frequency
    N = length(X);
    f = linspace(0, Fs/2, N/2 + 1);  % Frequency vector up to Nyquist frequency
    Y = fft(X);
    P = abs(Y(1:N/2 + 1)).^2 / (Fs * N);

    % Calculate the total average power using numerical integration
    totalAveragePower = trapz(f, P);

    % Display the total average power value 
    app.EditField_10.Value = num2str(totalAveragePower);
else
    disp('Error: The loaded data does not contain ''t'' and ''X''.');
end

        end

        % Value changed function: EditField_11
        function EditField_11ValueChanged(app, event)
            value = app.EditField_11.Value;
            
        end

        % Button pushed function: FirstMomentatt0Button
        function FirstMomentatt0ButtonPushed(app, event)
   
    % Get the mean and variance from the GUI (assuming 'EditField_4' and 'EditField_2' are handles to UI components)
    meanValue = str2double(app.EditField_4.Value);
    varianceValue = str2double(app.EditField_2.Value);

    % Check if the inputs are valid
    if isnan(meanValue) || isnan(varianceValue)
        disp('Error: Please enter valid mean and variance values.');
        return;
    end

    % Symbolic variables
    syms t;

    % Define the MGF
    MGF = exp(t * meanValue + 0.5 * t^2 * varianceValue);  %  t^2 * varianceValue using series i find it good equation to be calculate with

    % Calculate the first derivative (first moment) at t = 0
    firstMomentAtZero = double(subs(diff(MGF, t), t, 0));

    % Display the value in the GUI EditField_11
    app.EditField_11.Value = num2str(firstMomentAtZero);
        
        end

        % Button pushed function: SecondMomentatt0Button
        function SecondMomentatt0ButtonPushed(app, event)

  
    meanValue = str2double(app.EditField_4.Value);
    varianceValue = str2double(app.EditField_2.Value);

    % Check if the inputs are valid
    if isnan(meanValue) || isnan(varianceValue)
        disp('Error: Please enter valid mean and variance values.');
        return;
    end

    % Symbolic variables
    syms t;

    % Define the MGF
    MGF = exp(t * meanValue + 0.5 * t^2 * varianceValue);

    % Calculate the second derivative (second moment) at t = 0
    secondMomentAtZero = double(subs(diff(MGF, t, 2), t, 0));

    % Display the value in the GUI EditField_11
    app.EditField_12.Value = num2str(secondMomentAtZero);

        end

        % Value changed function: EditField_12
        function EditField_12ValueChanged(app, event)
            value = app.EditField_12.Value;
            
        end

        % Value changed function: FileEditField
        function FileEditFieldValueChanged2(app, event)
            value = app.FileEditField.Value;
            
        end

        % Button pushed function: PlotMButton
        function PlotMButtonPushed2(app, event)
 M = str2double(app.EnterNumberofsamplesEditField.Value);

% Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before plotting sample functions.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');
    
    % Ensure the number of samples 'M' is not greater than the number of columns in 'X'
    if M > size(X, 2)
        disp('Error: Not enough samples in the data.');
        return;
    end
    
    % Create a subplot grid
    rows = ceil(sqrt(M));
    cols = ceil(M / rows);

    % Randomly select 'M' unique indices
    selectedIndices = randperm(size(X, 2), M);
    
    % Plot selected sample functions in subplots
    figure;
    for i = 1:M
        % Get the current sample
        currentSample = X(:, selectedIndices(i));

        % Interpolate the current sample
        currentSample = interp1(1:size(currentSample, 1), currentSample, linspace(1, size(currentSample, 1), length(t)), 'linear', 'extrap');

        % Create subplot
        subplot(rows, cols, i);
        plot(t, currentSample, 'DisplayName', ['Sample Function ', num2str(selectedIndices(i))]);
        title(['Sample Function ', num2str(selectedIndices(i))]);
    end
    sgtitle('Sample Functions');
end



        end

        % Value changed function: EnterNumberofsamplesEditField
        function EnterNumberofsamplesEditFieldValueChanged(app, event)
            value = app.EnterNumberofsamplesEditField.Value;
            
        end

        % Button pushed function: TheensamplemeanoftheprocessButton
        function TheensamplemeanoftheprocessButtonPushed2(app, event)
            % Check if there are any variables loaded in the workspace
  workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating the ensemble mean.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');

    % Calculate the ensemble mean 
    ensembleMean = sum(X, 1) / size(X, 1);  % Calculate mean along columns

    % Plot the ensemble mean
    figure;
    plot(t, ensembleMean, 'LineWidth', 2);
    title('Ensemble Mean of the Process');
    xlabel('Time (t)');
    ylabel('Ensemble Mean Amplitude');
    grid on;
else
    disp('Error: The loaded data does not contain ''t'' and ''X''.');
end


        end

        % Value changed function: EditField_7
        function EditField_7ValueChanged(app, event)
            value = app.EditField_7.Value;
            
        end

        % Value changed function: nthoftheProcessEditField
        function nthoftheProcessEditFieldValueChanged(app, event)
            value = app.nthoftheProcessEditField.Value;
            
        end

        % Button pushed function: TimeofNButton
        function TimeofNButtonPushed2(app, event)
% Get the number of samples from the edit field
N = str2double(app.nthoftheProcessEditField.Value);

% Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating the ensemble mean.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');
    
    % Ensure the number of samples 'N' is not greater than the number of columns in 'X'
    if N > size(X, 2)
        disp('Error: Not enough samples in the data.');
        return;
    end
    
    % Randomly select 'N' unique indices
    selectedIndices = randperm(size(X, 2), N);
    
    % Initialize an array to store the time values of the selected samples
    selectedSampleTimes = zeros(N, 1);
    
    % Record the time values of the selected samples
    for i = 1:N
        % Get the current sample
        currentSample = X(:, selectedIndices(i));
        
        % Ensure that the sample length matches the length of 't'
        if length(currentSample) ~= length(t)
            %  truncate the sample to match the length of 't'
            currentSample = interp1(linspace(0, 1, length(currentSample)), currentSample, linspace(0, 1, length(t)));
        end
        
        % Randomly select a range of time for each amplitude
        startTime = t(randi(length(t)));
        endTime = t(randi(length(t)));
        
        % Calculate the time taken for the current sample
        timeTaken = abs(endTime - startTime);
        
        % Store the time value of the selected sample
        selectedSampleTimes(i) = timeTaken;
    end
    
    % Calculate the mean time of the selected samples
    meanTime = mean(selectedSampleTimes);
    
    % Display the mean time in 
    app.EditField_5.Value = num2str(meanTime);
else
    disp('Error: The loaded data does not contain ''t'' and ''X''.');
end

        end

        % Button pushed function: TofACFButton
        function TofACFButtonPushed2(app, event)
% Get the value of n from the edit field
n = str2double(app.nthoftheProcessEditField.Value);

% Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating the time auto-correlation function.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');
    
    % Ensure n is a valid index
    if isnan(n) || n < 1 || n > size(X, 2)
        disp('Error: Invalid value for n.');
        return;
    end
    
    % Get the n-th sample function
    nthSample = X(:, n);
    
    % Calculate the time auto-correlation function
    acf = xcorr(nthSample, nthSample, 'coeff');
    
    % Extract the relevant information for display
    acfValue = acf(length(t));  % Assuming length(t) corresponds to the relevant value
    
    % Display the time auto-correlation function in 
    app.EditField_6.Value = num2str(acfValue);
else
    disp('Error: The loaded data does not contain ''t'' and ''X''.');
end

        end

        % Value changed function: EditField_5
        function EditField_5ValueChanged(app, event)
            value = app.EditField_5.Value;
            
        end

        % Value changed function: EditField_6
        function EditField_6ValueChanged(app, event)
            value = app.EditField_6.Value;
            
        end

        % Button pushed function: ACFButton
        function ACFButtonPushed2(app, event)
            % Check if there are any variables loaded in the workspace
    workspaceVariables = evalin('base', 'who');

    if isempty(workspaceVariables)
        disp('Error: No data loaded. Load data before calculating the ACF.');
        return;
    end

    % Check if the loaded data contains 't' and 'X'      %  I dentify X and
    % t to be more flexible and for other datafiles
    if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
        % Retrieve 't' and 'X' from the base workspace
        t = evalin('base', 't');
        X = evalin('base', 'X');

        % Calculate and plot the auto-correlation function (ACF)
        figure;
        [acf, lags] = xcorr(mean(X, 1), 'coeff');  % Calculate mean along columns
        plot(lags, acf, 'LineWidth', 2);
        title('Auto-correlation Function (ACF)');
        xlabel('Lags');
        ylabel('Correlation Coefficient');
        grid on;
    else
        disp('Error: The loaded data does not contain ''t'' and ''X''.');
    end
        end

        % Button pushed function: EmsampleplotButton
        function EmsampleplotButtonPushed(app, event)
           
    M = str2double(app.EnterNumberofsamplesEditField.Value);

    % Check if there are any variables loaded in the workspace
    workspaceVariables = evalin('base', 'who');

    if isempty(workspaceVariables)
        disp('Error: No data loaded. Load data before plotting sample functions.');
        return;
    end

    % Check if the loaded data contains 't' and 'X'
    if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
        % Retrieve 't' and 'X' from the base workspace
        t = evalin('base', 't');
        X = evalin('base', 'X');

        % Ensure the number of samples 'M' is not greater than the number of columns in 'X'
        if M > size(X, 2)
            disp('Error: Not enough samples in the data.');
            return;
        end

        % Calculate and plot ensemble mean
        ensembleMean = mean(X(:, 1:M), 2);

        % Interpolate the ensemble mean
        ensembleMean = interp1(1:size(ensembleMean, 1), ensembleMean, linspace(1, size(ensembleMean, 1), length(t)), 'linear', 'extrap');

        % Plot ensemble mean
        figure;
        plot(t, ensembleMean, 'DisplayName', 'Ensemble Mean');
        title('Ensemble Mean');

        % Optionally add other plot settings (xlabel, ylabel, legend, etc.) if needed
    end


        end

        % Value changed function: EditField_13
        function EditField_13ValueChanged(app, event)
            value = app.EditField_13.Value;
            
        end

        % Button pushed function: ensamplemeanofMButton
        function ensamplemeanofMButtonPushed(app, event)
         % Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating the ensemble mean.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');

    % Specify the number of samples 'M'
    M = str2double(app.EnterNumberofsamplesEditField.Value);

    % Ensure the number of samples 'M' is not greater than the number of columns in 'X'
    if M > size(X, 2)
        disp('Error: Not enough samples in the data.');
        return;
    end

    % Calculate the ensemble mean of the first M samples
    ensembleMean = sum(X(:, 1:M), 2) / M;  

    % Display the mean value 
    app.EditField_13.Value = num2str(mean(ensembleMean));
else
    disp('Error: The loaded data does not contain ''t'' and ''X''.');
end

        end

        % Button pushed function: DACFButton
        function DACFButtonPushed(app, event)

    % Check if there are any variables loaded in the workspace
    workspaceVariables = evalin('base', 'who');

    if isempty(workspaceVariables)
        disp('Error: No data loaded. Load data before calculating the ACF.');
        return;
    end

    % Check if the loaded data contains 't' and 'X'
    if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
        % Retrieve 't' and 'X' from the base workspace
        t = evalin('base', 't');
        X = evalin('base', 'X');

        % Number of sample functions
        numSamples = size(X, 2);

        % Initialize a cell array to store ACF values
        acfCellArray = cell(numSamples, numSamples);

        % Calculate ACF for each pair of sample functions (i, j)
        for i = 1:numSamples
            for j = 1:numSamples
                % Calculate ACF for the pair (i, j)
                [acf, lags] = xcorr(X(:, i), X(:, j), 'coeff');

                % Store ACF values in the cell array
                acfCellArray{i, j} = acf;
            end
        end

        % Extract the ACF values for the 3D plot
        acfValues = cellfun(@(x) x(1), acfCellArray);

        % Create a 3D plot of ACF values
        figure;
        [I, J] = meshgrid(1:numSamples, 1:numSamples);
        surf(I, J, acfValues, 'EdgeColor', 'none');
        title('3D Plot of ACF between Sample Functions');
        xlabel('Sample Function i');
        ylabel('Sample Function j');
        zlabel('ACF Value');
        view(3); % 3D view
        colorbar;
    else
        disp('Error: The loaded data does not contain ''t'' and ''X''.');
    
end

        end

        % Callback function: not associated with a component
        function etimeaverageandthetimeACFofarandomsamplefunctionButtonPushed(app, event)

    % Check if there are any variables loaded in the workspace
    workspaceVariables = evalin('base', 'who');

    if isempty(workspaceVariables)
        disp('Error: No data loaded. Load data before calculating time properties.');
        return;
    end

    % Check if the loaded data contains 't' and 'X'
    if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
        % Retrieve 't' and 'X' from the base workspace
        t = evalin('base', 't');
        X = evalin('base', 'X');

        % Get the number of samples from the EnterNumberofsamplesEditField
        numSamples = str2double(app.EnterNumberofsamplesEditField.Value);

        % Check if the specified number of samples is valid
        if isnan(numSamples) || numSamples <= 0
            disp('Error: Please enter a valid number of samples.');
            return;
        end

        % Choose random sample indices
        sampleIndices = randi(size(X, 2), 1, numSamples);

        % Initialize an array to store autocorrelation values
        autocorrelationAtLag = zeros(1, numSamples);

        % Specify the lag for ACF calculation 
        lag = 1;

        % Calculate the autocorrelation at the specified lag for each random sample
        for i = 1:numSamples
            % Extract the random sample
            randomSample = X(:, sampleIndices(i));

            % Calculate the mean
            mean_x = mean(randomSample);

            % Calculate the autocorrelation at the specified lag manually
            shiftedSample = circshift(randomSample - mean_x, lag);
            autocorrelationAtLag(i) = sum((randomSample - mean_x) .* shiftedSample) / sum((randomSample - mean_x).^2);
        end

        % Display the mean autocorrelation 
        app.TimeACFEditField.Value = num2str(mean(autocorrelationAtLag));
    else
        disp('Error: The loaded data does not contain ''t'' and ''X''.');
    end


        
        end

        % Value changed function: StatisticalMeanEditField
        function StatisticalMeanEditFieldValueChanged(app, event)
            value = app.StatisticalMeanEditField.Value;
            
        end

        % Value changed function: TimeMeanEditField
        function TimeMeanEditFieldValueChanged(app, event)
            value = app.TimeMeanEditField.Value;
            
        end

        % Value changed function: StatisticalACFEditField
        function StatisticalACFEditFieldValueChanged(app, event)
            value = app.StatisticalACFEditField.Value;
            
        end

        % Value changed function: TimeACFEditField_2
        function TimeACFEditField_2ValueChanged(app, event)
            value = app.TimeACFEditField_2.Value;
            
        end

        % Button pushed function: ComparingButton
        function ComparingButtonPushed(app, event)
% Check if there are any variables loaded in the workspace
workspaceVariables = evalin('base', 'who');

if isempty(workspaceVariables)
    disp('Error: No data loaded. Load data before calculating the ensemble statistics.');
    return;
end

% Check if the loaded data contains 't' and 'X'
if ismember('t', workspaceVariables) && ismember('X', workspaceVariables)
    % Retrieve 't' and 'X' from the base workspace
    t = evalin('base', 't');
    X = evalin('base', 'X');

    % Ensure there is at least one column in 'X'
    if isempty(X)
        disp('Error: No sample functions in the data.');
        return;
    end

    % Choose a random sample function index for each operation
    randomIndex = randi([1, size(X, 2)]);

    % Ensure that the randomIndex does not exceed the number of columns in X
    if randomIndex > size(X, 2)
        disp('Error: Random index exceeds the number of sample functions.');
        return;
    end

    % Perform trapezoidal numerical integration for the time mean
    dt = diff(t);
    timeMean = sum((X(1:end-1, randomIndex) + X(2:end, randomIndex)) .* dt) / 2;
    timeMean = mean(timeMean); % Take the mean of the resulting vector
    app.TimeMeanEditField.Value = num2str(timeMean);

    % Calculate Statistical Mean (Expected Value) for the random sample function
    values = X(:, randomIndex);
    pdf_values = calculateProbabilityDensityFunction(values);
    statisticalMean = sum(values .* pdf_values) / length(t);
    app.StatisticalMeanEditField.Value = num2str(statisticalMean);

    % Calculate Statistical ACF for the random sample function
    N = length(t);
    acf_values_stat = zeros(2*N-1, 1);

    for lag = 1:2*N-1
        acf_values_stat(lag) = sum(X(1:end-lag+1, randomIndex) .* X(lag:end, randomIndex));
    end

    statisticalACF = mean(acf_values_stat);
    app.StatisticalACFEditField.Value = num2str(statisticalACF);

    % Manually calculate the mean absolute autocorrelation value
    acf_values = zeros(2*N-1, 1);

    for lag = 1:2*N-1
        acf_values(lag) = sum(abs(X(1:end-lag+1, randomIndex) .* X(lag:end, randomIndex)));
    end

    timeACF = mean(acf_values);
    app.TimeACFEditField_2.Value = num2str(timeACF);
else
    disp('Error: The loaded data does not contain ''t'' and ''X''.');
end

% Function to calculate the probability density function
function pdf_values = calculateProbabilityDensityFunction(values)
   
    pdf_values = ones(size(values)) / (max(values) - min(values));
end





        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0 0.4471 0.7412];
            app.UIFigure.Position = [100 100 1907 730];
            app.UIFigure.Name = 'MATLAB App';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 1 1907 730];

            % Create RandomVariableTab
            app.RandomVariableTab = uitab(app.TabGroup);
            app.RandomVariableTab.Title = 'Random Variable';

            % Create FileNameEditFieldLabel
            app.FileNameEditFieldLabel = uilabel(app.RandomVariableTab);
            app.FileNameEditFieldLabel.HorizontalAlignment = 'right';
            app.FileNameEditFieldLabel.Position = [54 595 60 22];
            app.FileNameEditFieldLabel.Text = 'File Name';

            % Create FileNameEditField
            app.FileNameEditField = uieditfield(app.RandomVariableTab, 'text');
            app.FileNameEditField.ValueChangedFcn = createCallbackFcn(app, @FileNameEditFieldValueChanged, true);
            app.FileNameEditField.Position = [129 595 170 22];

            % Create ImportButton
            app.ImportButton = uibutton(app.RandomVariableTab, 'push');
            app.ImportButton.ButtonPushedFcn = createCallbackFcn(app, @ImportButtonPushed, true);
            app.ImportButton.Position = [153 551 100 23];
            app.ImportButton.Text = 'Import';

            % Create MeanButton
            app.MeanButton = uibutton(app.RandomVariableTab, 'push');
            app.MeanButton.ButtonPushedFcn = createCallbackFcn(app, @MeanButtonPushed, true);
            app.MeanButton.Position = [64 444 100 23];
            app.MeanButton.Text = 'Mean';

            % Create EditField_2
            app.EditField_2 = uieditfield(app.RandomVariableTab, 'text');
            app.EditField_2.ValueChangedFcn = createCallbackFcn(app, @EditField_2ValueChanged, true);
            app.EditField_2.Position = [199 398 100 22];

            % Create VarianceButton
            app.VarianceButton = uibutton(app.RandomVariableTab, 'push');
            app.VarianceButton.ButtonPushedFcn = createCallbackFcn(app, @VarianceButtonPushed, true);
            app.VarianceButton.Position = [67 398 100 23];
            app.VarianceButton.Text = 'Variance';

            % Create TheThirdmomentButton
            app.TheThirdmomentButton = uibutton(app.RandomVariableTab, 'push');
            app.TheThirdmomentButton.ButtonPushedFcn = createCallbackFcn(app, @TheThirdmomentButtonPushed, true);
            app.TheThirdmomentButton.Position = [67 350 113 23];
            app.TheThirdmomentButton.Text = 'The Third moment';

            % Create EditField_3
            app.EditField_3 = uieditfield(app.RandomVariableTab, 'text');
            app.EditField_3.ValueChangedFcn = createCallbackFcn(app, @EditField_3ValueChanged, true);
            app.EditField_3.Position = [199 350 100 22];

            % Create EditField_4
            app.EditField_4 = uieditfield(app.RandomVariableTab, 'text');
            app.EditField_4.ValueChangedFcn = createCallbackFcn(app, @EditField_4ValueChanged, true);
            app.EditField_4.Position = [199 444 100 22];

            % Create PlotMGFButton
            app.PlotMGFButton = uibutton(app.RandomVariableTab, 'push');
            app.PlotMGFButton.ButtonPushedFcn = createCallbackFcn(app, @PlotMGFButtonPushed, true);
            app.PlotMGFButton.Position = [67 300 100 23];
            app.PlotMGFButton.Text = 'Plot MGF';

            % Create PlotthefirstandsecondderivitiveButton
            app.PlotthefirstandsecondderivitiveButton = uibutton(app.RandomVariableTab, 'push');
            app.PlotthefirstandsecondderivitiveButton.ButtonPushedFcn = createCallbackFcn(app, @PlotthefirstandsecondderivitiveButtonPushed, true);
            app.PlotthefirstandsecondderivitiveButton.Position = [26 235 195 23];
            app.PlotthefirstandsecondderivitiveButton.Text = 'Plot the first and second derivitive';

            % Create FirstMomentatt0Button
            app.FirstMomentatt0Button = uibutton(app.RandomVariableTab, 'push');
            app.FirstMomentatt0Button.ButtonPushedFcn = createCallbackFcn(app, @FirstMomentatt0ButtonPushed, true);
            app.FirstMomentatt0Button.Position = [368 517 122 23];
            app.FirstMomentatt0Button.Text = 'First Moment at t =0';

            % Create SecondMomentatt0Button
            app.SecondMomentatt0Button = uibutton(app.RandomVariableTab, 'push');
            app.SecondMomentatt0Button.ButtonPushedFcn = createCallbackFcn(app, @SecondMomentatt0ButtonPushed, true);
            app.SecondMomentatt0Button.Position = [358 466 143 23];
            app.SecondMomentatt0Button.Text = 'Second Moment at t = 0';

            % Create EditField_11
            app.EditField_11 = uieditfield(app.RandomVariableTab, 'text');
            app.EditField_11.ValueChangedFcn = createCallbackFcn(app, @EditField_11ValueChanged, true);
            app.EditField_11.Position = [558 517 100 22];

            % Create EditField_12
            app.EditField_12 = uieditfield(app.RandomVariableTab, 'text');
            app.EditField_12.ValueChangedFcn = createCallbackFcn(app, @EditField_12ValueChanged, true);
            app.EditField_12.Position = [560 466 100 22];

            % Create RandomProcessTab_2
            app.RandomProcessTab_2 = uitab(app.TabGroup);
            app.RandomProcessTab_2.Title = 'Random Process';

            % Create EnterNumberofsamplesEditFieldLabel
            app.EnterNumberofsamplesEditFieldLabel = uilabel(app.RandomProcessTab_2);
            app.EnterNumberofsamplesEditFieldLabel.HorizontalAlignment = 'right';
            app.EnterNumberofsamplesEditFieldLabel.Position = [56 497 144 22];
            app.EnterNumberofsamplesEditFieldLabel.Text = 'Enter Number of  samples';

            % Create EnterNumberofsamplesEditField
            app.EnterNumberofsamplesEditField = uieditfield(app.RandomProcessTab_2, 'text');
            app.EnterNumberofsamplesEditField.ValueChangedFcn = createCallbackFcn(app, @EnterNumberofsamplesEditFieldValueChanged, true);
            app.EnterNumberofsamplesEditField.Position = [219 497 130 22];

            % Create FileEditFieldLabel
            app.FileEditFieldLabel = uilabel(app.RandomProcessTab_2);
            app.FileEditFieldLabel.HorizontalAlignment = 'right';
            app.FileEditFieldLabel.Position = [84 616 34 22];
            app.FileEditFieldLabel.Text = 'File   ';

            % Create FileEditField
            app.FileEditField = uieditfield(app.RandomProcessTab_2, 'text');
            app.FileEditField.ValueChangedFcn = createCallbackFcn(app, @FileEditFieldValueChanged2, true);
            app.FileEditField.Position = [133 616 183 22];

            % Create nthoftheProcessEditFieldLabel
            app.nthoftheProcessEditFieldLabel = uilabel(app.RandomProcessTab_2);
            app.nthoftheProcessEditFieldLabel.HorizontalAlignment = 'right';
            app.nthoftheProcessEditFieldLabel.Position = [459 636 102 22];
            app.nthoftheProcessEditFieldLabel.Text = 'nth of the Process';

            % Create nthoftheProcessEditField
            app.nthoftheProcessEditField = uieditfield(app.RandomProcessTab_2, 'text');
            app.nthoftheProcessEditField.ValueChangedFcn = createCallbackFcn(app, @nthoftheProcessEditFieldValueChanged, true);
            app.nthoftheProcessEditField.Position = [576 636 187 22];

            % Create ImportButton_2
            app.ImportButton_2 = uibutton(app.RandomProcessTab_2, 'push');
            app.ImportButton_2.ButtonPushedFcn = createCallbackFcn(app, @ImportButton_2Pushed2, true);
            app.ImportButton_2.Position = [164 580 100 23];
            app.ImportButton_2.Text = 'Import';

            % Create PlotMButton
            app.PlotMButton = uibutton(app.RandomProcessTab_2, 'push');
            app.PlotMButton.ButtonPushedFcn = createCallbackFcn(app, @PlotMButtonPushed2, true);
            app.PlotMButton.Position = [234 454 100 23];
            app.PlotMButton.Text = 'Plot M';

            % Create TheensamplemeanoftheprocessButton
            app.TheensamplemeanoftheprocessButton = uibutton(app.RandomProcessTab_2, 'push');
            app.TheensamplemeanoftheprocessButton.ButtonPushedFcn = createCallbackFcn(app, @TheensamplemeanoftheprocessButtonPushed2, true);
            app.TheensamplemeanoftheprocessButton.Position = [54 385 204 23];
            app.TheensamplemeanoftheprocessButton.Text = 'The ensample mean of the process';

            % Create ACFButton
            app.ACFButton = uibutton(app.RandomProcessTab_2, 'push');
            app.ACFButton.ButtonPushedFcn = createCallbackFcn(app, @ACFButtonPushed2, true);
            app.ACFButton.Position = [461 444 100 23];
            app.ACFButton.Text = 'ACF';

            % Create TimeofNButton
            app.TimeofNButton = uibutton(app.RandomProcessTab_2, 'push');
            app.TimeofNButton.ButtonPushedFcn = createCallbackFcn(app, @TimeofNButtonPushed2, true);
            app.TimeofNButton.Position = [819 635 100 23];
            app.TimeofNButton.Text = 'Time of N';

            % Create EditField_5
            app.EditField_5 = uieditfield(app.RandomProcessTab_2, 'text');
            app.EditField_5.ValueChangedFcn = createCallbackFcn(app, @EditField_5ValueChanged, true);
            app.EditField_5.Position = [968 636 100 22];

            % Create TofACFButton
            app.TofACFButton = uibutton(app.RandomProcessTab_2, 'push');
            app.TofACFButton.ButtonPushedFcn = createCallbackFcn(app, @TofACFButtonPushed2, true);
            app.TofACFButton.Position = [819 595 100 23];
            app.TofACFButton.Text = 'T of ACF';

            % Create EditField_6
            app.EditField_6 = uieditfield(app.RandomProcessTab_2, 'text');
            app.EditField_6.ValueChangedFcn = createCallbackFcn(app, @EditField_6ValueChanged, true);
            app.EditField_6.Position = [968 595 100 22];

            % Create BSDPlotButton
            app.BSDPlotButton = uibutton(app.RandomProcessTab_2, 'push');
            app.BSDPlotButton.ButtonPushedFcn = createCallbackFcn(app, @BSDPlotButtonPushed2, true);
            app.BSDPlotButton.Position = [461 497 100 23];
            app.BSDPlotButton.Text = 'BSD Plot';

            % Create emsampemeanvalueButton
            app.emsampemeanvalueButton = uibutton(app.RandomProcessTab_2, 'push');
            app.emsampemeanvalueButton.ButtonPushedFcn = createCallbackFcn(app, @emsampemeanvalueButtonPushed2, true);
            app.emsampemeanvalueButton.Position = [54 341 133 23];
            app.emsampemeanvalueButton.Text = 'emsampe mean value';

            % Create EditField_7
            app.EditField_7 = uieditfield(app.RandomProcessTab_2, 'text');
            app.EditField_7.ValueChangedFcn = createCallbackFcn(app, @EditField_7ValueChanged, true);
            app.EditField_7.Position = [214 341 100 22];

            % Create ACFValueButton
            app.ACFValueButton = uibutton(app.RandomProcessTab_2, 'push');
            app.ACFValueButton.ButtonPushedFcn = createCallbackFcn(app, @ACFValueButtonPushed, true);
            app.ACFValueButton.Position = [58 300 131 23];
            app.ACFValueButton.Text = 'ACF Value';

            % Create EditField_8
            app.EditField_8 = uieditfield(app.RandomProcessTab_2, 'text');
            app.EditField_8.ValueChangedFcn = createCallbackFcn(app, @EditField_8ValueChanged, true);
            app.EditField_8.Position = [214 300 100 22];

            % Create BSDValueButton
            app.BSDValueButton = uibutton(app.RandomProcessTab_2, 'push');
            app.BSDValueButton.ButtonPushedFcn = createCallbackFcn(app, @BSDValueButtonPushed, true);
            app.BSDValueButton.Position = [459 385 100 23];
            app.BSDValueButton.Text = 'BSD Value';

            % Create EditField_9
            app.EditField_9 = uieditfield(app.RandomProcessTab_2, 'text');
            app.EditField_9.ValueChangedFcn = createCallbackFcn(app, @EditField_9ValueChanged, true);
            app.EditField_9.Position = [576 385 180 33];

            % Create AveragePowerButton
            app.AveragePowerButton = uibutton(app.RandomProcessTab_2, 'push');
            app.AveragePowerButton.ButtonPushedFcn = createCallbackFcn(app, @AveragePowerButtonPushed, true);
            app.AveragePowerButton.Position = [459 307 100 23];
            app.AveragePowerButton.Text = 'Average Power';

            % Create EditField_10
            app.EditField_10 = uieditfield(app.RandomProcessTab_2, 'text');
            app.EditField_10.ValueChangedFcn = createCallbackFcn(app, @EditField_10ValueChanged, true);
            app.EditField_10.Position = [576 307 100 22];

            % Create EmsampleplotButton
            app.EmsampleplotButton = uibutton(app.RandomProcessTab_2, 'push');
            app.EmsampleplotButton.ButtonPushedFcn = createCallbackFcn(app, @EmsampleplotButtonPushed, true);
            app.EmsampleplotButton.Position = [234 420 100 23];
            app.EmsampleplotButton.Text = 'Emsample plot';

            % Create ensamplemeanofMButton
            app.ensamplemeanofMButton = uibutton(app.RandomProcessTab_2, 'push');
            app.ensamplemeanofMButton.ButtonPushedFcn = createCallbackFcn(app, @ensamplemeanofMButtonPushed, true);
            app.ensamplemeanofMButton.Position = [58 454 127 23];
            app.ensamplemeanofMButton.Text = 'ensample mean of M';

            % Create EditField_13
            app.EditField_13 = uieditfield(app.RandomProcessTab_2, 'text');
            app.EditField_13.ValueChangedFcn = createCallbackFcn(app, @EditField_13ValueChanged, true);
            app.EditField_13.Position = [74 420 100 22];

            % Create DACFButton
            app.DACFButton = uibutton(app.RandomProcessTab_2, 'push');
            app.DACFButton.ButtonPushedFcn = createCallbackFcn(app, @DACFButtonPushed, true);
            app.DACFButton.Position = [890 395 100 23];
            app.DACFButton.Text = '3D ACF';

            % Create StatisticalMeanEditFieldLabel
            app.StatisticalMeanEditFieldLabel = uilabel(app.RandomProcessTab_2);
            app.StatisticalMeanEditFieldLabel.HorizontalAlignment = 'right';
            app.StatisticalMeanEditFieldLabel.Position = [46 111 93 22];
            app.StatisticalMeanEditFieldLabel.Text = {'Statistical Mean:'; ''};

            % Create StatisticalMeanEditField
            app.StatisticalMeanEditField = uieditfield(app.RandomProcessTab_2, 'text');
            app.StatisticalMeanEditField.ValueChangedFcn = createCallbackFcn(app, @StatisticalMeanEditFieldValueChanged, true);
            app.StatisticalMeanEditField.Position = [154 111 100 22];

            % Create TimeMeanEditFieldLabel
            app.TimeMeanEditFieldLabel = uilabel(app.RandomProcessTab_2);
            app.TimeMeanEditFieldLabel.HorizontalAlignment = 'right';
            app.TimeMeanEditFieldLabel.Position = [75 48 68 22];
            app.TimeMeanEditFieldLabel.Text = {'Time Mean:'; ''};

            % Create TimeMeanEditField
            app.TimeMeanEditField = uieditfield(app.RandomProcessTab_2, 'text');
            app.TimeMeanEditField.ValueChangedFcn = createCallbackFcn(app, @TimeMeanEditFieldValueChanged, true);
            app.TimeMeanEditField.Position = [158 48 100 22];

            % Create StatisticalACFEditFieldLabel
            app.StatisticalACFEditFieldLabel = uilabel(app.RandomProcessTab_2);
            app.StatisticalACFEditFieldLabel.HorizontalAlignment = 'right';
            app.StatisticalACFEditFieldLabel.Position = [317 111 86 22];
            app.StatisticalACFEditFieldLabel.Text = 'Statistical ACF:';

            % Create StatisticalACFEditField
            app.StatisticalACFEditField = uieditfield(app.RandomProcessTab_2, 'text');
            app.StatisticalACFEditField.ValueChangedFcn = createCallbackFcn(app, @StatisticalACFEditFieldValueChanged, true);
            app.StatisticalACFEditField.Position = [418 111 100 22];

            % Create TimeACFEditField_2Label
            app.TimeACFEditField_2Label = uilabel(app.RandomProcessTab_2);
            app.TimeACFEditField_2Label.HorizontalAlignment = 'right';
            app.TimeACFEditField_2Label.Position = [344 48 58 22];
            app.TimeACFEditField_2Label.Text = 'Time ACF';

            % Create TimeACFEditField_2
            app.TimeACFEditField_2 = uieditfield(app.RandomProcessTab_2, 'text');
            app.TimeACFEditField_2.ValueChangedFcn = createCallbackFcn(app, @TimeACFEditField_2ValueChanged, true);
            app.TimeACFEditField_2.Position = [417 48 100 22];

            % Create ComparingButton
            app.ComparingButton = uibutton(app.RandomProcessTab_2, 'push');
            app.ComparingButton.ButtonPushedFcn = createCallbackFcn(app, @ComparingButtonPushed, true);
            app.ComparingButton.Position = [216 158 100 23];
            app.ComparingButton.Text = 'Comparing';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = appprob_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end