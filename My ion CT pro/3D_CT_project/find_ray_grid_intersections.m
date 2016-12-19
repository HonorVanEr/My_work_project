function [Segment_Lengths, Cells, Valid_Intersections, Intersections] = find_ray_grid_intersections(Rays, varargin)
% FIND_RAY_GRID_INTERSECTIONS Calculate the length of ray segments within an n-dimensional grid
%   Calculate the length of ray segments within an n-dimensional irregular grid 
%   and also intersection locations for multiple rays
%
% INPUTS
%   Rays contains all rays that may intersect the grid
%   Rays is an No_Rays by 2 by No_Dims 3-dimensional array
%   No_Rays is the number of rays
%   No_Dims is the number of dimensions
%   Rays(i, 1, :) contains the starting coordinates for ray i
%   Rays(i, 2, :) contains the ending   coordinates for ray i
%   varargin contains No_Dims vectors. Each vector contains the grid
%   locations for that coordinate.
%
% OUTPUTS
%   Segment_Lengths contains the length of segments within each cell for each ray
%   Segment_Lengths is a 2-dimensional array
%   Segment_Lengths(i, j) is the segment length for ray i in cell j
%   Cells(i, j, k) identifies the grid cell number for ray i, cell j in coordinate k
%   Valid_Intersections is a logical array that identifies the valid
%   ray intersections for ray i and cell j
%   Intersections(i, j, k) contains the value of the intersection for ray i and coordinate k
%
% Commentary
%   The output results are contained in matrices. There is 1 row for
%   each ray. The number of segments within each ray will vary
%   depending on the number of cells that are intersected from start
%   to finish. Valid_Intersections identifies which intersection
%   locations are valid. The valid entries in the segment matrices
%   are greater than zero. These identify the segment lengths for
%   each ray in each grid cell.
%   Each ray is modelled in parametric form as a line of the form (3D example)
%   x = x0 + u*dx, where dx = x1-x0
%   y = y0 + u*dy, where dy = y1-y0
%   z = z0 + u*dy, where zy = z1-z0
%   and 0 <= u <= 1

% Form a matrix of u values that describes all intersections
% There is a row for each ray
% Each column represents intersections

No_Rays = size(Rays, 1);
No_Dims = size(Rays, 3);
% Calcluate the maximum no of intersections for each ray
Ni_Max = 2;
for dim = 1:No_Dims
    Ni_Max = Ni_Max + length(varargin{dim});
end
u = zeros(No_Rays, Ni_Max);
% The 1st column is 0 for the start of each ray
u(:, 2) = 1;    % This is the u value for the end of each ray
Grid_Idx = 3;   % Starting index location for coordinate dim
for dim = 1:No_Dims
    Grid_Values = reshape(varargin{dim}, 1, []);   % row vector
    Rs = Rays(:, 1, dim);       % column vector containing start of ray for coordinate dim
    Re = Rays(:, 2, dim);       % column vector containing end of ray for coordinate dim
    dR = Re - Rs;               % ray displacements for coordnate dim
    % Find all intersections with the x-axis grid lines
    GV_Rep = repmat(Grid_Values, No_Rays, 1);
    No_GV = length(Grid_Values);
    Rs_Rep = repmat(Rs, 1, No_GV);      % Replicated Ray starts
    dr_Rep = repmat(dR, 1, No_GV);      % Replicated Ray displacements
    u_Grid = (GV_Rep - Rs_Rep)./dr_Rep; % Find all u values that intersect grids
    u(:, Grid_Idx : Grid_Idx + No_GV-1) = u_Grid;
    Grid_Idx = Grid_Idx + No_GV;        % calculate grid index for next dimension
end

u = sort(u, 2); % sorted u values for all intersections
Valid_Intersections = u >= 0 & u <= 1;   % valid intersections
Ii = not(Valid_Intersections);           % Invalid intersections
Intersections = zeros(No_Rays, Ni_Max, No_Dims);
for dim = 1:No_Dims
    Rs = Rays(:, 1, dim);   % start of rays for coordinate dim
    Re = Rays(:, 2, dim);   % end of rays for coordinate dim
    dR = Re - Rs;           % ray displacements for dimension dim
    Ri = repmat(Rs, 1, Ni_Max) + u.*repmat(dR, 1, Ni_Max);
    Ri(Ii) = 0;
    Intersections(:,:,dim) = Ri;
end

% Calculate segment lengths
Ns_Max = Ni_Max - 1;    % No of segments
ssi = 1:Ns_Max;         % Segment start index
sei = 1 + ssi;          % Segment end index;
dRi = diff(Intersections, 1, 2);
Segment_Lengths = sqrt(sum(dRi.^2, 3));
Vs =  u(:, ssi) >=0 & u(:, sei) <= 1; % Valid segments
Is = not(Vs);                       % Invalid Segments
Segment_Lengths(Is) = 0;

Imin = min(Intersections(:, ssi, :), Intersections(:, sei, :));
Cells = zeros(No_Rays, Ns_Max, No_Dims);
for dim = 1:No_Dims
    Grid_Values = reshape(varargin{dim}, 1, []);   % row vector
    Grid_Rep = reshape(repmat(Grid_Values, Ns_Max*No_Rays, 1), No_Rays, Ns_Max, []);
    No_GV = length(Grid_Values);
    Int_Rep = repmat(Imin(:, :, dim), 1, 1, No_GV);
    %Int_Rep = repmat(Imin(:, :, dim), [1, 1, No_GV]);
    % Calculate the cell index value for this dimension
    Cell_Idx = sum(Int_Rep >= Grid_Rep, 3);
    Cell_Idx(Is) = 0;
    Cells(:,: , dim) = Cell_Idx;
end
end