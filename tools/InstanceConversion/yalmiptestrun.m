function [] = yalmiptestrun(varargin)
%yalmiptestrun loads a file listing .dat-s files and invokes sdpaintoyalmip
%for each of them and creates an output file with solutions and solving
%times
%
%    yalmiptestrun('filename') solves all problems listed in 'filename'

addpath(genpath('/sw/cplex/cplex1261'))
addpath(genpath('/local/gally/programs/MOSEK/8.1.0.54'))
addpath(genpath('/home/gally/MATLAB_libraries/SeDuMi_1_3/conversion'))
addpath(genpath('/home/gally/MATLAB_libraries/YALMIP-R20180926-master'))

filename = varargin{1};

% Does the file exist
if ~exist(filename)
    error(['No such file.']);
end
    
filestring = fileread(filename);
instances = strread(filestring, '%s', 'delimiter', sprintf('\n'));
resultfile = fopen(strcat(filename, '.results'), 'a');
time = 0;
obj = 0;
fprintf(resultfile, 'instance                                                                                   objective       time            nodes\n');
global yalmiptestnnodes;
global yalmiptesterror;

for i=1:numel(instances)
    char(instances(i))
    yalmiptesterror = 0;
    [time, obj, nodes] = sdpainttoyalmip(char(instances(i)))
    yalmiptesterror
    if yalmiptesterror == yalmiptestnnodes
        fprintf(resultfile, '%-90.89s error           %9.7f      %9.0d \n', char(instances(i)), 3600, yalmiptestnnodes);
    else
        if obj < 0
            fprintf(resultfile, '%-90.89s %9.7f      %9.7f      %9.0d \n', char(instances(i)), obj, time, nodes);
        else
            fprintf(resultfile, '%-90.89s  %9.7f      %9.7f      %9.0d \n', char(instances(i)), obj, time, nodes);
        end
    end
end

fclose(resultfile);
