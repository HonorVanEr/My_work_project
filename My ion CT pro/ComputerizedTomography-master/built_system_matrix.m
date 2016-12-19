function [A, sino] = built_system_matrix(no_pixels, sub_sample, slice_num, down_sample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Constructs a system matrix for fan-beam recon of the center slice of 
% micro ct data.
% Unit: millimeters.
%
% Based on the test04_fanbeamtomolinear.m by 
% Jakob Sauer Joergensen, DTU Compute, 2014-04-06
% jakj@dtu.dk
% Adapted to 02526 - Matematisk modellering by
% Michael Andersen, DTU Compute, 2015-04-20
% miand@dtu.dk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify parameters
% default parameters
% Slice number in the Z direction, for projections
slicenum   = 523;
% Sub sampling of projection angles
a_sampling = 10;
% Number of pixel in our domain
N          = 128;
% Down sampling of the sensor panel
downsample = 10;
if nargin < 4
    downsample = 10;
else
    downsample = down_sample;
end
if nargin < 3
    slicenum = 523;
    downsample = 10;
else
    slicenum = slice_num;
end

if nargin < 2
    % Subsampling factor for projections,
    % 1 = no subsampling, 10 = use only every 10th proj.
    a_sampling = 10;
    slicenum   = 523;
else
    a_sampling = sub_sample;
end

if nargin < 1
    N          = 128;
    a_sampling = 10;
    slicenum   = 523;
else
    N = no_pixels;
end

%% Paths and data specific geometry parameters

% MUST SET DATAPATH.
homedir = getenv('HOME');                        % Get the users home directory
datapath  = [homedir '/Documents/MYgithub/My_idlwork_pro/My ion CT pro/ComputerizedTomography-master/'];   % Remember folder separator after (i.e. '/' in this case)
groupname = 'fox100kV20W.to-short500';                       % Groupname aka filename prefix
total_number  = 1572;                            % Number of projections
num_bins      = 1000;                            % Number of detector sensors
source_origin = 459.251480102539;                % Distance from source to origin/object
source_detector = 962.7192;                      % Distance from source to detector
origin_det      = source_detector-source_origin; % Distance from origin/object to detector
voxel_side_1024 = 0.190814309219281;             % Voxel side length.

%% Pick out projections and set geometry parameters.

% List of projections to use
alist = 1:a_sampling:total_number;  
num_angles = size(alist,2);

% Physical length of one detector pixel
detector_pixel_side = ...
    voxel_side_1024 * (source_origin + origin_det)/source_origin;

% And relative to the object voxel size
detector_pixel_side_relative = ...
    detector_pixel_side/voxel_side_1024/num_bins*N;

% The physical length of the detector
object_physical_length = 1024*voxel_side_1024;

% The physical side length of the object
detector_physical_length = num_bins*detector_pixel_side;


%% Load projections and compute sinogram data

% Initizalize sinogram
sino = zeros(num_bins,num_angles); 

% Loop over projections, load, extract central line, which is a 1D
% projection of the 2D central slice of the total 3D object.
for k = 1:length(alist)
    fprintf('Loading projection %d of %d (%d of %d)...\n',alist(k),...
        total_number,k,num_angles)
    filename = sprintf('%s_%04d.tif',groupname,alist(k));
    P =  imread(fullfile(datapath,filename));
    sino(:,k) = P(slicenum,:);
end

% The loaded projection data are given as transmission percentages of the
% initial (background) intensity. To obtain the projection data, divide by
% 100 and do minus logarithm.
sino = -log(sino/100);

sub_idx = 1:downsample:size(sino,1);
sino = sino(sub_idx, :);
sino = sino(:);
%% Set up fanbeamtomolinear projection matrix

% Angles: num_angles equally distributed on 360 degrees.
theta = linspace(0,360,num_angles+1);
theta = theta(1:end-1);

% Source to center of rotation distance in units of the object side length.
R = source_origin / object_physical_length;

% Detector's width in units of the object side length.
dw = detector_physical_length / object_physical_length;

% Source to detector distance in units of the object side length.
sd = (source_origin + origin_det) / object_physical_length;

isDisp = 0.0;

% Set up the system matrix.
disp('Begin building system matrix...');
tic
A = fanbeamtomolinear(N,theta,num_bins/down_sample,R,dw,sd,isDisp);
t = toc();
disp(['System matrix built in ' num2str(t) ' seconds.']);

end