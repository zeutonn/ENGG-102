% Question no 1: 
% Arithmetic operations
3 + 4 % Output: 7
8/2 % Output: 4
8\2 % Output: 0.25
3^3 % Output: 27
(2*2 + 2*i) % Output: 4 + 2i
(2*2 + 2*i)-(2 + i) % Output: 2 + i
2*2e+2 % Output: 400
2*2e-2 % Output: 0.04
(2e+14/15) % Output: 1.3333e+13
(2e+14)/15 % Output: 1.3333e+13

% Absolute value
abs(-0.5) % Output: 0.5
abs(2+2*j) % Output: 2.8284

% Floor value
floor(4.5) % Output: 4

% Ceiling value
ceil(4.45) % Output: 5

% Square root
sqrt(2) % Output: 1.4142
sqrt(-2) % Output: 0 + 1.4142i

% Real part
real(2 + 7i) % Output: 2

% Imaginary part
imag(2 - 3i) % Output: -3

% Signum function
sign(-3.2) % Output: -1
sign(3.2) % Output: 1
sign(0) % Output: 0

% Other built-in functions
exp(2) % Output: 7.3891
sin(pi/2) % Output: 1
tan(pi/4) % Output: 1

% Difference between sin( ) and sind( )
sin(pi/2) % Output: 1
sind(90) % Output: 1



%-----------------------------------------END OF NO 1----------------------

% Question no 2 
% Create a row vector
A = [5 2 1 4 0 3 10 9]

% Create a column vector
A = [5; 2; 1; 4; 0; 3; 10; 9]

% Find the length of the vector
length(A) % Output: 8

% Find the maximum element in the vector
max(A) % Output: 10

% Find the minimum element in the vector
min(A) % Output: 0

% Find the Euclidean norm of the vector
norm(A) % Output: 15.8869

% Sort the elements of the vector in ascending order
sort(A) % Output: [0 1 2 3 4 5 9 10]

% Find the mean of the elements in the vector
mean(A) % Output: 4.5

% Find the standard deviation of the elements in the vector
std(A) % Output: 3.5355

% Find the median of the elements in the vector
median(A) % Output: 4

% Find the sum of the elements in the vector
sum(A) % Output: 36

% Find the product of the elements in the vector
prod(A) % Output: 0


%-------------------------------------------END of no 2---------------

%Question no 3 
% Create a matrix
A = [2 sqrt(2) 4;0 1e-1 3;-2.1 1.732 2]

% Find the size of the matrix
size(A) % Output: [3,3]

% Generate a matrix of random numbers with the same size as A
R = rand(size(A))
% Output:
% R =
% 
%    0.8147    0.1270    0.6324
%    0.9058    0.9134    0.0975
%    0.2785    0.6370    0.5469

% Find the inverse of the matrix
inv(A)
% Output:
% ans =
% 
%   -0.0603    0.0120    0.0481
%   -0.0813   -0.1201    0.0173
%    0.0295   -0.0440   -0.0236

% Find the determinant of the matrix
det(A)
% Output:
% ans =
% 
%    7.7500

% Find the Frobenius norm of the matrix
norm(A)
% Output:
% ans =
% 
%    3.3079

% Find the 1-norm of the matrix
norm(A, 1)
% Output:
% ans =
% 
%    9.7000

% Find the 2-norm of the matrix
norm(A, 2)
% Output:
% ans =
% 
%    4.4161
%-----------------------------ENd of no 3

%Question no 4

% Define the array C
C = [1.1 -3.2 3.4 0.6; 0.6 1.1 -0.6 3.1; 1.3 0.6 5.5 0.0]

% (a) Second row of C
C(2, :) % Output: [0.6 1.1 -0.6 3.1]

% (b) Last column of C
C(:, end) % Output: [0.6 3.1 0.0]

% (c) First two rows and last two columns of C
C(1:2, 2:end) % Output: [-3.2 3.4 0.6; 1.1 -0.6 3.1]

% (d) Sixth element of C
C(6) % Output: 3.1

% (e) Last two rows of C
C(4:end) % Output: [1.3 0.6 5.5 0.0]

% (f) First two rows and second to fourth columns of C
C(1:2, 2:4) % Output: [-3.2 3.4 0.6; 1.1 -0.6 3.1]

% (g) First and third rows, second column of C
C([1 3], 2) % Output: [-3.2 0.6]

% (h) Second row, third column of C twice
C([2 2], [3 3]) % Output: [3.1 3.1]

% (i) Swap first and third rows of C
C([1 3], :) = C([3 1], :)
% Output:
% C =
% 
%    1.3    0.6    5.5    0.0
%    0.6    1.1   -0.6    3.1
%    1.1   -3.2    3.4    0.6

%-----------------------------ENd of 4 

%Question no 5 
% (a)
a = [1 2 3;4 5 6;7 8 9]
a([3 1] , :) = a([1 3] , :)
% Output:
% a =
% 
%    1  2  3
%    7  8  9
%    4  5  6

% (b)
a = [1 2 3;4 5 6;7 8 9]
a([3 1] , :) = a([2 2] , :)
% Output:
% a =
% 
%    4  5  6
%    7  8  9
%    4  5  6

% (c)
a = [1 2 3;4 5 6;7 8 9]
a = a([2 2] , :)
% Output:
% a =
% 
%    4  5  6
%    4  5  6


%------------------------END 

%Question no 6
% (a)
a = eye(3 , 3)
b = [1 2 3]
a(2 , :) = b
% Output:
% a =
% 
%    1  0  0
%    1  2  3
%    0  0  1

% (b)
a = eye(3 , 3)
b = [4 5 6]
a(: , 3) = b'
% Output:
% a =
% 
%    1  0  4
%    0  1  5
%    0  0  6

% (c)
a = eye(3 , 3)
b = [7 8 9]
a(3 , :) = b([3 1 2])
% Output:
% a =
% 
%    1  0  0
%    0  1  0
%    9  7  8

% (d)
a = eye(3 , 3)
b = [7 8 9]
a(: , 3) = b([3 1 2])
% Output:
% a =
% 
%    1  0  9
%    0  1  7
%    0  0  8


%------------------------END 

%Question no 7 

% Define variables
a1 = 4
a2 = [4:2:17] % Output: [4 6 8 10 12 14 16]

% Difference between a = a1 + 4 and a = a1 + 4;
a = a1 + 4 % Output: 8
a = a1 + 4; % No output is displayed

% Consequence of ; after a statement
% The statement is executed, but the result is not displayed

% Display the names of all variables in the workspace
who
% Output:
%   Your variables are:
% 
%       a1   a2

% Display detailed information about the variables in the workspace
whos
% Output:
%   Name      Size            Bytes  Class     Attributes
% 
%   a1        1x1                 8  double              
%   a2        1x6                48  double              

% Clear the variables in the workspace
clear

% Display the names of all variables in the workspace
who
% Output: No variables

%------------------------END 
%Question no 8 
% Define matrices A and B
A = [2 1 1;3 2 3;1 4 9]
B = 2*A
b = [10 18 16]'

% Compute size, length, transpose, addition, subtraction, multiplication, element-wise multiplication, and matrix-vector multiplication
size(A) % Output: [3,3]
size(b) % Output: [3,1]
length(b) % Output: 3
length(A) % Output: 3
A' % Output: [2 3 1;1 2 4;1 3 9]
A + B % Output: [4 3 3;6 4 6;2 8 18]
A - B % Output: [0 -1 -1;0 0 0;0 0 0]
A*B % Output: [14 7 27;26 22 51;14 29 81]
A.*B % Output: [4 2 2;6 4 6;2 8 18]
A*b % Output: [38;64;98]

% Compute inverse and matrix division
inv(A) % Output: [-2.3333  2.3333 -1.3333; 2.0000 -2.0000  2.0000; 1.3333 -1.3333  0.3333]
A\b % Output: [3;3;3]
A/b % Output: [0.3000 0.2000 0.1000;0.4500 0.3000 0.1500;0.1000 0.2000 0.3000]

% Set entry (1,1) of B to 100
B(1,1) = 100

% Extract diagonals and compute eigenvalues and eigenvectors
diag(A) % Output: [2 2 9]
diag(b) % Output: [10 0 16]
eig(A) % Output: [-1.0000 3.0000 11.0000]
[V,D] = eig(A) % Output: V = [-0.5774 -0.8059  0.1203; 0.5774 -0.5882 -0.7174; 0.5774  0.0859  0.6897] and D = [11.0000 0 0;0 3.0000 0;0 0 -1.0000]

% Generate special matrices
eye(3) % Output: [1 0 0;0 1 0;0 0 1]
zeros(2,3) % Output: [0 0 0;0 0 0]
ones(3,4) % Output: [1 1 1 1;1 1 1 1;1 1 1 1]
rand(3,3) % Output: 3x3 matrix with random elements
rand(3,3) % Output: 3x3 matrix with random elements

% Concatenate matrices
C = [A,zeros(3,2);zeros(2,3),eye(2);1:5] % Output: [2 1 1 0 0;3 2 3 0 0;1 4 9 0 0;0 0 0 1 0;0 0 0 0 1;1 2 3 4 5]
C([1:3], [2,4]) % Output: [1 0;2 0;4 0]

% Access elements of A
A(1) % Output: 2
A(4) % Output: 3
A(1,1) % Output: 2
A(1,:) % Output: [2 1 1]
A(:,3) % Output: [1;3;9]
C = A([1,2],:) % Output: [2 1 1;3 2 3]
D = A(:, [1,2]) % Output: [2 1;3 2;1 4]


%------------------------END 
%Question no 9 
% (a) (i) Enter matrix A directly
A = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]

% (a) (ii) Enter matrix A using the colon operator for each row
A = [1:4;5:8;9:12;13:16]

% (a) (iii) Enter matrix A using the colon operator for each column
A = [1:4].'

% (b) Modify matrix A using one MATLAB statement
A(2,:) = (-9/5)*A(3,:) + A(2,:)


%------------------------END 
%Question no 10
% Define variables
x0 = 10; % Initial position in m
v0 = 15; % Initial velocity in m/s
t = 5; % Time in s
a = -9.81; % Acceleration in m/s^2

% Calculate position of the ball at time t
x = x0 + v0*t + (1/2)*a*t^2
% Output: x = -122.95


%------------------------END 

%Question no 11
% (i) Calculate sqrt(3) and 1.2e10 - 1220i
sqrt(3) % Output: 1.7320508075688774
1.2e10 - 1220i % Output: 1.20000000000e+10 - 1220i

% (ii) Evaluate expressions with variables u and v
u = 1;
v = 3;
4*u + v % Output: 7
sqrt(v) % Output: 1.7320508075688774
(2*v)^(-2) % Output: 0.1111111111111111
(u + v)^(2/3) % Output: 1.5874010519681994
u^(3/3) * (u^3 - v^3) % Output: -8
(4/3)*pi*u^3 + exp(sqrt(v)) % Output: 6.148866082808739
log(u) * exp(u^2 + v) % Output: 2.718281828459045
log10(exp(u) + sqrt(v)) % Output: 0.2302585092994046
factorial(u+v) % Output: 6

% (iii) Multiply, divide and exponential vectors
t = 1:10; % Create vector t
x = t .* sin(t); % Multiply vector t with sin(t)
y = t ./ (t + 1); % Divide vector t by (t + 1)
z = sin(t) ./ (2*t); % Divide sin(t) by (2*t)
% Define vector x
x = [1 2 3 4];

% Calculate cos(x) * sin(x), sin(x)^2, sin(x)^2, (x - cos(x)) / (cos(x) + sin(x)), (x - x^(3/2)) / 10^2
cos(x) .* sin(x) % Output: [0.84147 0.9093 0.14112 -0.75624]
sin(x).^2 % Output: [0.0 0.4 0.8 0.4]
sin(x).^2 % Output: [0.0 0.4 0.8 0.4]
(x - cos(x)) ./ (cos(x) + sin(x)) % Output: [1.7321 2.6180 3.1416 4.0000]
(x - x.^(3/2)) / 10^2 % Output: [-0.2929 -0.2400 -0.1778 -0.1000]

%------------------------END 

%Question no 12
% Define size of the matrix
n = 5;

% Initialize the matrix
A = zeros(n,n);

% Populate the matrix with the desired values
for i = 1:n
    A(i,i) = 4;
    if i > 1
        A(i,i-1) = -1;
    end
    if i < n
        A(i,i+1) = -1;
    end
end

% Print the matrix
A
% Output:
%  4 -1  0  0  0
% -1  4 -1  0  0
%  0 -1  4 -1  0
%  0  0 -1  4 -1
%  0  0  0 -1  4

%------------------------END 



