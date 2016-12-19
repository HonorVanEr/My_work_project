%% Part 1

% Choose domain size.
N = 128;

% Task 1, generate and show a phantom
f = phantom(N); % f \in R^{n x n}

% Visualize the phantom
figure()
imagesc(f);
axis equal; colormap gray; axis off;

% Select projection angles
% ...
num_angles = 158;
theta = linspace(0,360,num_angles+1);
theta = theta(1:end-1);

% Create system matrix A and projection data b;
% Use create_projections( ... )
[W, P] = create_projections(f,theta);

% Verify that b = A*x
x = W*reshape(f, 128*128,1); %?
norm(P-x,2); %?

% Dummy variables for dummy art_solver( ... )
A = 0; b = 1;
I = 10;
% Task 2, implement ART solver
X = art_solver(A, b, I); % X \in R^{n*n x 1}

% Show the reconstructed image


% Plot the error (Mean squared error)


%% Part 2


% Task 1, Construct system matrix A and projection data b, based on 
% real data, use the function built_system_matrix( ... ) 
% - Edit the geometric details in the function to match what is found 
% in the name.xtekct, other things to play around with is how b/sino is
% downsampled (see line 99-103).
N = 128;
[A, sino] = built_system_matrix(N, 10, 500, 10);

% Perform reconstruction using ART
I = 60000;
X = art_solver(A, sino, I);

% Show the reconstructed image
figure();
imagesc(reshape(X, N, N)); colormap gray; axis off; axis equal;

%% Filtering
B = X;
p = quantile(X, 0.50);
truth = X<=p;
B(truth) = p;

figure();
B = reshape(B, N, N);
%k = round(N/4);
k = 1;
B(1:k,1:k) = p;
B(1:k, (end-k):end) = p;
B((end-k):end, 1:k) = p;
B((end-k):end, (end-k):end) = p;


imagesc(B); colormap gray; axis off; axis equal;
