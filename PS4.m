%% Problem Set 4 Parts 2 and 3

clear all
close all
clc

% Since elimination without pivoting was stable for me up to alpha = 40, I
% tested out larger values of alpha and found that it was stable up to
% alpha = 709 and unstable for alpha >= 710;
alphas = [0,5,10,20,40,709,710];
k = length(alphas);
absolute_errors = zeros(1,k);
for iteration = 1:k
    alpha = alphas(iteration);
    A = [exp(-alpha) 1 0; -1 exp(-alpha) -1; 1 exp(-alpha) -2];
    n = length(A);
    A_copy = A; % To check if Ax = b for original A.
    L = eye(n,n);

    for jj = 1:n-1
        pivot = A(jj,jj);
        for ii = jj+1:n
            c = -A(ii,jj)/pivot;
            L(ii,jj) = -c;
            A(ii,:) = A(ii,:) + c*A(jj,:);
        end
    end

    U = A;
    %A_copy - L*U
    b = [1 exp(-alpha) exp(-alpha)].';

    y = zeros(n,1);
    y(1) = b(1);
    for ii = 2:n
        y(ii) = b(ii);
        for jj = 1:ii-1
            y(ii) = y(ii) - L(ii,jj)*y(jj);
        end
    end

    x = zeros(n,1);
    x(n) = y(n) / U(n,n);
    for ii = n-1:-1:1
        x(ii) = y(ii);
        for jj = ii+1:n
            x(ii) = x(ii) - U(ii,jj)*x(jj);
        end
        x(ii) = x(ii) / U(ii,ii);
    end

    absolute_errors(iteration) = norm(A_copy*x - b,1);
end

absolute_errors

%% Problem 1 Part 4

clear all
close all
clc

% Elimination without pivoting became unstable for alpha >= 710. Here we
% see that even for alpha as large as 10^6, elimination with partial
% pivoting remains stable (larger alpha works too of course!).
alphas = [0,5,10,20,40,709,710,10^3,10^6];
k = length(alphas);
absolute_errors = zeros(1,k);
for iteration = 1:k
    alpha = alphas(iteration);
    A = [exp(-alpha) 1 0; -1 exp(-alpha) -1; 1 exp(-alpha) -2];
    n = length(A);
    A_copy = A; % To check if Ax = b for original A.
    L = eye(n,n);
    
    b = [1 exp(-alpha) exp(-alpha)].';
    b_copy = b; % To check if Ax = b for original A and b.

    for jj = 1:n
        % partial pivoting by largest entry at or below current pivot.
        [max_pivot, max_index] = max(abs(A(jj:n,jj)));
        max_pivot_row = max_index + jj - 1;
        A_temp = A(jj,:);
        b_temp = b(jj);
        A(jj,:) = A(max_pivot_row,:);
        b(jj) = b(max_pivot_row);
        A(max_pivot_row,:) = A_temp;
        b(max_pivot_row) = b_temp;
        pivot = A(jj,jj);
        % partial pivoting complete - back to elimination as before.
        for ii = jj+1:n
            c = -A(ii,jj)/pivot;
            L(ii,jj) = -c;
            A(ii,:) = A(ii,:) + c*A(jj,:);
        end
    end

    U = A;
    %A_copy - L*U
    y = zeros(n,1);
    y(1) = b(1);
    for ii = 2:n
        y(ii) = b(ii);
        for jj = 1:ii-1
            y(ii) = y(ii) - L(ii,jj)*y(jj);
        end
    end

    x = zeros(n,1);
    x(n) = y(n) / U(n,n);
    for ii = n-1:-1:1
        x(ii) = y(ii);
        for jj = ii+1:n
            x(ii) = x(ii) - U(ii,jj)*x(jj);
        end
        x(ii) = x(ii) / U(ii,ii);
    end
    
    % partial pivoting does not change (correct) solution for x but does
    % change b, so check results using the original b (before row swaps).
    absolute_errors(iteration) = norm(A_copy*x - b_copy,1);
end

absolute_errors