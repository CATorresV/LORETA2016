function X_s = myLoreta(M, N, y_obs)
% Computing the inverse problem solution following the Tikhonov
% Regularization method with LORETA for estimating the solution.
% ||Ax-y||^2 + ||Gamma x||^2
% y_obs size N x M
% X_s   size P x M
% A     size N x P

%Gamma in Tikhonov regularization is defined as Gamma = alpha * I
% alpha = 0.5; %Regularization parameter
% Gamma = alpha * I;
% 
% X_s = inv(A' * A + Gamma' * Gamma) * A' * y_obs ; 


%% LEad Field MAtrix
% Points coordinates of the LEad Field
%M = 1000 ; %Number of points in the discrete space
Mn = M^(1/3);
[x, y, z] = ndgrid(linspace(-0.8,0.8,round(Mn)));

V = [x(:),y(:),z(:)];

%N = 132; %number of electrodes
[x,y,z]=sphere;
x=(x(:));
y=(y(:));
z=(z(:)); %Centered at (0,0,0)

x=x((z>0));
y=y((z>0));
z=z((z>0));

x=x(1:N);
y=y(1:N);
z=z(1:N);

S = [x,y,z];   %Coordinates of the measurement electrodes

K = zeros(N,M,3); %Initialization of the lead field matrix
sigma = 0.3; %Conductivity of the medium

for i=1:N
    for j=1:M
        term1 = (1/(4*pi*sigma))*((S(i,:)-V(j,:))/(norm(S(i,:)-V(j,:),3)));
        term2 = (1/(4*pi*sigma))*((S(i,:)-V(j,:))/(norm(S(1,:)-V(1,:),3)));
        K(i,j,:) = term1 - term2;
    end
end

%% Laplacian matrix operator construction

A_1= zeros(M,M);

d =1; % reference distance

for i=1:M
    for j=1:M
        if (norm(V(i,:)-V(j,:))==d)
            A_1(i,j)=1/6;
        end
    end
end

A_0 = (1/2)*( eye(M) + inv( eye(M) * diag(A_1 * ones(M))))*A_1;

A = kron(A_0,eye(3));

B = (6/d^2) * (A - eye(3*M));

%% LORETA 
Omega = zeros(M,M);
for i=1:M
    sum=0;
    for j=1:N
        sum = sum + K(j,i)*K(j,i);
    end
    Omega(i,i) = sqrt(sum);
end


W_loreta = ( kron(Omega,eye(3)) ) * B' * B ( kron(Omega,eye(3))  );


%% Inverse problem solution

T = inv(W_loreta) * K' * inv(K* W_loreta * K') ;

X_s = T * y_obs;

%% From Eduardo Giraldo Thesis
%W_loreta2 = 

% M = A; %Lead field matrix
% 
% W = eye(size(M)); %Identity Matrix derives into the Tikhonov regularization
% 
% 
% X_s = inv(M' * inv(C_e) * M + lambda_k^2 * W) * (M' * inv(C_e) * y_obs + lambda_k^2 * X_prior);