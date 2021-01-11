function [A] = readMatrix(file)
% liest Matrix aus 'file' und speichert sie in A

fid = fopen(file,'r');

% first line: type of matrix
line = fgetl(fid);

% second line: dimension of matrix
line = fgetl(fid);
dim = sscanf(line,"m = %d n = %d k = %d");
m = dim(1);
n = dim(2);
k = dim(3);
A = zeros(m,n);

% remaining lines: actual matrix
for i = 1:m
    line = fgetl(fid);
    entries = sscanf(line,"%f");
    A(i,:) = entries';
end

end