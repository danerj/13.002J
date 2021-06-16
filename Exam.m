%% Exam Problem 1

A = [10,7,8,7; 7,5,6,5; 8,6,10,9; 7,5,9,10]
Ainv = inv(A)
b = [4,3,3,1].'
x = Ainv*b
cond(A,inf)

% bhat = b + 0.009
% db = b-bhat
% norm(db,inf)
% 
% xhat = Ainv*bhat
% dx = x - xhat
% norm(dx,inf)

%% Exam Problem 2

clear all
close all
clc

% a)

f = @(x) x.*tan(x) - 2;
fprime = @(x) tan(x) + x.*(sec(x)).^2;
xspace = linspace(0,pi/2-.1,1000);

figure(1)
figure1 = plot(xspace,f(xspace));
title("$f(x) = x\tan x - 2$", 'interpreter', 'latex')
xlabel("$x$", 'interpreter', 'latex')
ylabel("$f(x)$", 'interpreter', 'latex')
movegui(figure1,'center')

tolerance = 10^-9;
p(1) = 1;
errors_NRI(1) = abs(f(p(1)));
k = 1;
j = 1;

while errors_NRI(end) > tolerance
    p(k+1) = p(1,k) - f(p(k)) / fprime(p(k));
    errors_NRI(k+1) = abs(f(p(k+1)));
    k = k+1;
end

q(1) = 1;
q(2) = 1.05;
errors_secant(1) = abs(f(q(2)));

j = 2;

while errors_secant(end) > tolerance
    q(j+1) = q(j) - f(q(j))*(q(j) - q(j-1)) / (f(q(j)) - f(q(j-1)));
    errors_secant(j) = abs(f(q(j+1)));
    j = j+1;
end

figure(2)
figure2 = plot(1:k-1, errors_NRI(2:end), '*', 1:j-2,...
    errors_secant(2:end), 'o', 'linewidth', 1);
title("Error at each iteration", 'interpreter', 'latex')
xlabel("$k$", 'interpreter','latex')
ylabel("Error")
legend('Newton-Raphson', 'Secant')
movegui('center')
