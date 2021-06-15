%% Problem 1 Part 5
clear all
close all
clc

N=50;
y = zeros(1,N+1);
y(1) = log(6/5);

for n = 2:N+1
    y(n) = 1/(n-1) - 5*y(n-1);
end

plot(0:N,y,'*','linewidth',1)
title('Part 5','interpreter','latex')
xlabel('$n$', 'interpreter','latex')
ylabel('$y_n$','interpreter','latex')

%% Problem 1 Part 6
clear all
close all
clc

N=50;
y = zeros(1,N+1);
y(end) = 0;

for n = N+1:-1:2
    y(n-1) = 1/(5*(n-1))-y(n)/5;
end

plot(0:N,y,'*','linewidth',1)
title('Part 6','interpreter','latex')
xlabel('$n$', 'interpreter','latex')
ylabel('$y_n$','interpreter','latex')