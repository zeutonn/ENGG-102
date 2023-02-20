%Question 1 
% Define constants
x0 = 10; % Initial position
v0 = 15; % Initial velocity
g = -9.81; % Acceleration due to gravity

% Define time interval
t = 0:0.01:3; % Time in seconds

% Calculate distance travelled
x = x0 + v0*t + (1/2)*g*t.^2;

% Plot distance versus time
plot(t,x)
xlabel('Time (s)')
ylabel('Distance (m)')
title('Distance Travelled by a Falling Ball')

% Add text to figure window
gtext('Ball')

%---------END
%Question 2 
% Define time interval
t = -2*pi:0.01:2*pi;

% Calculate y values
y1 = sin(t);
y2 = cos(t);

% Plot y1 versus t
plot(t,y1,'r','LineWidth',3)
hold on

% Plot y2 versus t
plot(t,y2,'b--','LineWidth',3)

% Add labels and title
xlabel('Time (s)')
ylabel('y(t)')
title('Plot of sin(t) and cos(t)')

% Add legend
legend('sin(t)','cos(t)')

% Release hold
hold off

%---------END
%Question 3 
% Define number of points
n1 = 10;
n2 = 50;
n3 = 100;

% Define x values
x1 = linspace(0,4*pi,n1);
x2 = linspace(0,4*pi,n2);
x3 = linspace(0,4*pi,n3);

% Calculate y values
y1 = exp(-0.4*x1).*sin(x1);
y2 = exp(-0.4*x2).*sin(x2);
y3 = exp(-0.4*x3).*sin(x3);

% Plot y1 versus x1
plot(x1,y1)
hold on

% Plot y2 versus x2
plot(x2,y2)

% Plot y3 versus x3
plot(x3,y3)

% Add labels and title
xlabel('x')
ylabel('y')
title('Plot of y = e^{-0.4x}*sin(x)')

% Release hold
hold off

%---------END
%Question 4 
% Define variables
N0 = 3.049; % Initial population in billions
K = 12.77593983; % Carrying capacity
r = 0.0270767630; % Intrinsic growth rate

% Define time points
t = [1960 1965 1970 1975 1980 1985 1990 1995 2000 2005 2010];

% Define population data
N = [3.049 3.358 3.721 4.103 4.475 4.882 5.249 5.679 6.127 6.514 6.900];

% Calculate population according to logistic model
N_logistic = N0*K./(N0 + (K - N0)*exp(-r*t));

% Plot discrete data points
plot(t,N,'o')
hold on

% Plot logistic model
plot(t,N_logistic)

% Add labels and title
xlabel('Year')
ylabel('Population (billions)')
title('Population according to logistic model')

% Release hold
hold off

%---------END
%Question 5 
% Define time points
t = 0:0.1:10;

% Calculate v(t)
v = 10*exp(-0.2+pi*1i)*t;

% Plot v(t) using plot(t,v)
figure
plot(t,v)
xlabel('t')
ylabel('v(t)')
title('Plot of v(t) using plot(t,v)')

% Plot v(t) using plot(v)
figure
plot(v)
xlabel('Real part')
ylabel('Imaginary part')
title('Plot of v(t) using plot(v)')

% Plot real and imaginary parts of v(t)
figure
plot(t,real(v),'b')
hold on
plot(t,imag(v),'r--')
xlabel('t')
ylabel('v(t)')
title('Plot of real and imaginary parts of v(t)')
legend('Real part','Imaginary part')
hold off

% Create polar plot of v(t)
figure
polarplot(angle(v),abs(v))
title('Polar plot of v(t)')

% Plot v(t) using plot3
figure
plot3(real(v),imag(v),t)
xlabel('Real part')
ylabel('Imaginary part')
zlabel('t')
title('Plot of v(t) using plot3')

%---------END
%Question 6 
x = -3:0.1:3;
y = x;
[X,Y] = meshgrid(x,y);
Z = X.*Y.*(X.^2 - Y.^2)./(X.^2 + Y.^2);
figure(1)
mesh(X,Y,Z);
h = meshc(X,Y,Z);
set(h,'LineWidth',2);
x = x + (x==0)*eps;
y = y + (y==0)*eps;

x = -3:0.1:3; % create an array of x values from -3 to 3 with a step size of 0.1
y = -3:0.1:3; % create an array of y values from -3 to 3 with a step size of 0.1
[X,Y] = meshgrid(x,y); % create meshgrid using x and y arrays
Z = sin(X) + sin(Y); % calculate the function z = sin x + sin y for all values of x and y
figure; % create a new figure
mesh(X,Y,Z); % plot the mesh of the function

x = -3:0.1:3; % create an array of x values from -3 to 3 with a step size of 0.1
y = -3:0.1:3; % create an array of y values from -3 to 3 with a step size of 0.1
[X,Y] = meshgrid(x,y); % create meshgrid using x and y arrays
Z = sin(X) + sin(Y); % calculate the function z = sin x + sin y for all values of x and y
figure; % create a new figure
meshc(X,Y,Z); % plot the meshc of the function

x = -3:0.1:3; % create an array of x values from -3 to 3 with a step size of 0.1
y = -3:0.1:3; % create an array of y values from -3 to 3 with a step size of 0.1
[X,Y] = meshgrid(x,y);


% Define the range for x and y
x = -5:0.1:5;
y = -5:0.1:5;

% Create a meshgrid using x and y
[X,Y] = meshgrid(x,y);

% Calculate the function z
Z = cos(X) .* cos(Y) .* exp(-sqrt(X.^2 + Y.^2)./4);

% Plot the surf plot
figure(1)
surf(X,Y,Z)

% Plot the surfc plot
figure(2)
surfc(X,Y,Z)


% Define the range for t
t = 0:0.1:pi/2;

% Calculate x and y
x = t;
y = exp(t.^2) .* sin(t);

% Plot using ezplot
figure(1)
ezplot(x,y)

% Plot using fplot
figure(2)
fplot(@(t) t, @(t) exp(t.^2) .* sin(t), [0,pi/2])

x = linspace(-pi/2, pi/2, 100);
y = linspace(-pi/2, pi/2, 100);
[X,Y] = meshgrid(x,y);
Z = X.^2 + sin(X.*Y) + Y.^2;
figure
ezcontour(X,Y,Z)

x = linspace(-pi, pi, 100);
y = linspace(-2, 2, 100);
[X,Y] = meshgrid(x,y);
Z = X.^2 + (1-cos(Y));
figure
surf(X,Y,Z)
hold on
contour(X,Y,Z)


%---------END
%Question 7 
x = linspace(-2, 2, 100);
y = linspace(-2, 2, 100);
[X, Y] = meshgrid(x, y);
Z = -X.*Y.*exp(-X.^2 - Y.^2);
contour(X, Y, Z, 20);
c = contour(X, Y, Z, 20);
clabel(c, 'manual');

contour3(X, Y, Z, 20);
grid off;
meshc(X, Y, Z);
hidden off;


%---------END
%Question 8 

%---------END
%Question 9 
% Define the time vector
t = 0 : pi/20 : 4*pi;

% Define the function y(t)
y = exp(-0.2t).(cos(t) + 1i*sin(t));

% Plot y(t) using the conventional plot command
figure;
plot(t, y, 'LineWidth', 3);
xlabel('Time');
ylabel('y(t)');
title('Plot of y(t)');

% Define the time vector
t = 0 : pi/20 : 4*pi;

% Define the function y(t)
y = exp(-0.2t).(cos(t) + 1i*sin(t));

% Plot the real and imaginary parts of y(t)
figure;
plot(t, real(y), '--r', 'LineWidth', 3);
hold on;
plot(t, imag(y), 'LineWidth', 3);
xlabel('Time');
ylabel('Real and Imaginary Parts of y(t)');
legend('Real Part', 'Imaginary Part');
title('Plot of Real and Imaginary Parts of y(t)');

% Define the time vector
t = 0 : pi/20 : 4*pi;

% Define the function y(t)
y = exp(-0.2t).(cos(t) + 1i*sin(t));

% Plot the real versus imaginary parts of y(t)
figure;
plot(real(y), imag(y), 'LineWidth', 3);
xlabel('Real Part of y(t)');
ylabel('Imaginary Part of y(t)');
title('Plot of Real vs Imaginary Parts of y(t)');

% Define the time vector
t = 0 : pi/20 : 4*pi;

% Define the function y(t)
y = exp(-0.2t).(cos(t) + 1i*sin(t));

% Plot the polar plot of the angle of y(t) versus the absolute value of y(t)
figure;
polarplot(angle(y), abs(y), 'LineWidth', 3);
title('Polar Plot of Angle of y(t) vs Absolute Value of y(t)');


%Question 10 
% Define the number of moles of the gas
n = 1;

% Define the universal gas constant
R = 8.314;

% Define the temperature of the gas (in K)
T = 273;

% Define the range of pressures to consider (in Kpa)
P = 1:1000;

% Calculate the volume of the gas for each pressure
V = nRT./P;

% Plot the volume versus pressure
figure;
plot(P, V, 'r', 'LineWidth', 3);
xlabel('Pressure (Kpa)');
ylabel('Volume (L)');
title('Volume vs Pressure for Ideal Gas');

% Define the number of moles of the gas
n = 1;

% Define the universal gas constant
R = 8.314;

% Define the range of pressures to consider (in Kpa)
P = 1:1000;

% Calculate the volume of the gas for each pressure at T = 273 K
V1 = nR273./P;

% Calculate the volume of the gas for each pressure at T = 373 K
V2 = nR373./P;

% Plot the volume versus pressure for both temperatures
figure;
plot(P, V1, 'r', 'LineWidth', 3);
hold on;
plot(P, V2, 'b--', 'LineWidth', 2);
xlabel('Pressure (Kpa)');
ylabel('Volume (L)');
legend('T = 273 K', 'T = 373 K');
title('Volume vs Pressure for Ideal Gas at Different Temperatures');

% Define the number of moles of the gas
n = 1;

% Define the universal gas constant
R = 8.314;

% Define the volume of the gas (in L)
V = 101;

% Define the range of temperatures to consider (in K)
T = 250:450;

% Calculate the pressure of the gas for each temperature
P = nRT./V;

% Plot the pressure versus temperature
figure;
plot(T, P);
xlabel('Temperature (K)');
ylabel('Pressure (Kpa)');
title('Pressure vs Temperature for Ideal Gas');
%---------END
%Question 11
% Define the time vector
t = 0:20;

% Define the incoming concentration function
C = 2 + sin(2*t);

% Plot the incoming concentration
figure;
plot(t, C, 'LineWidth', 3);
xlabel('Time');
ylabel('Incoming Concentration (gram/gallon)');
title('Variation in Incoming Concentration');

% Define the time vector
t = 0:20;

% Define the constants in the model
A0 = 0;
k = 0.5;

% Calculate the steady solution
A_steady = 20 - 40/17cos(2t) + 10/17sin(2t);

% Calculate the unsteady solution
A_unsteady = -300/17*exp(-t/2);

% Plot the solution
figure;
plot(t, A_steady + A_unsteady, 'LineWidth', 3);
xlabel('Time');
ylabel('Concentration (gallon)');
title('Variation in Concentration in Pond');


%---------END


