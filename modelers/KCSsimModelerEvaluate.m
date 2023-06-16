%% Read design of KCSsimModeler and map parametric points on surface.

%design must be written to filename_sample as a column

%parametric points to be mapped onto surface must be written to
%   filename_planar with each point occupying a line.
%   coordinates are separated by spaces

%on-surface points are written with each point occupying a line
%   and coordinates being separated by a space in filename_surface

addpath('nurbs-toolbox');

filename_design = "../../design.txt";

filename_planar = "../../planar_points.txt";

filename_surface = "../../surface_points.txt";

surface = KCSsimModeler(readmatrix(filename_design));

planar_points = transpose(readmatrix(filename_planar));

surface_points = transpose(nrbeval(surface, planar_points));

writematrix(surface_points,filename_surface,'Delimiter'," ");

exit;

