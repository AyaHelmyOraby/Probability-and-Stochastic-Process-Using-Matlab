# Probability-and-Stochastic-Process-Using-Matlab


a GUI-based tool that allows a user to:
1) Enter the values of random variable values and results in the statistics of such variable.
2) Enter any stochastic process and results in the ensemble and the time statistics of such process.
The GUI can be built using Matlab or any other software package.
GUI Description
The GUI should do the following:
1) Section 1: Random Variables
• Allow the user enter a random variable in the form of its sample space.
An example .m file of the sample space is attached.
• Display the mean, the variance and the third moment of the random variable
• Plot the MGF M(t) vs 0 < t < 2
• Plot the first and the second derivatives of M(t), and calculate their values at t = 0
2) Section 2: Random Processes
• Allow the user enter a random process in the form of the ensemble, i.e. all the sample functions,
each defined by two vectors; time and amplitude. Note that the time vector can be common to
all the sample functions.
An example .m file of the ensemble is attached.
• Allow the user to perform and display the following:
– Plot M sample functions of the ensemble of the process, where M is entered by the user
– Calculate and plot the ensemble mean of the process
– Calculate and plot the statistical auto-correlation function
– Calculate the time mean of the n

th sample function of the process, where n is entered by

the user
– Calculate the time auto-correlation function of the n

th sample function of the process, where

n is entered by the user
– Calculate and plot the power spectral density of the process
– Calculate the total average power of the process


Test  GUI for the following:
1) The RV in the sample file.
2) X is a RV, where X ∼ U(−3, 5).
3) Y is a RV, where Y ∼ N (−8, 4).
4) The RP in the sample file.
5) Z(t) is a RP, where Z(t) = cos(4πt + θ), where 0 ≤ t ≤ 2, θ ∼ U(0, π).
6) W(t) is a RP, where W(t) = A cos(4πt), where 0 ≤ t ≤ 2, A ∼ N (−5, 5).

