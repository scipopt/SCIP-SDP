function [time,obj,nodes] = sdpainttoyalmip(varargin)
%sdpainttoyalmip loads a problem definition in the SDPA with additional
%integrality constraints and inserts it to yalmip for solving it with bnb
%
%    [F,h] = sdpainttoyalmip('filename')  Loads the problem min(h(x)), F(x)>0, x(i) from file 'filename'
%    [F,h] = sdpainttoyalmip         A "Open" - box will be opened

% the following code is copied from loadsdpafile.m in yalmip with only the 
% integrality conditions newly added

filename = strtrim(varargin{1});

% Does the file exist
if ~exist(filename)
    %filename = [filename '.dat-s'];
    if ~exist(filename)
        error(['No such file.']);
    end
end
    
% Load using SeDuMi
[At,b,c,K] = fromsdpa(filename);

nvars = length(b);
x = sdpvar(nvars,1);

F = ([]);
top = 1;
if isfield(K,'l')
    if K.l(1)>0
        X = c(top:top+K.l-1)-At(top:top+K.l-1,:)*x;
        F = F + (X(:)>=0);
        top = top + K.l;
    end
end

if isfield(K,'s')
    if K.s(1)>0
        for i = 1:length(K.s)
            [ix,iy,iv] = find([c(top:top+K.s(i)^2-1) At(top:top+K.s(i)^2-1,:)]);
            off = (ix-1)/(K.s(i)+1);
            if all(off == round(off))           
                X = c(top:top+K.s(i)^2-1)-At(top:top+K.s(i)^2-1,:)*x;
                F = F + (diag(reshape(X,K.s(i),K.s(i))) >= 0);
                top = top + K.s(i)^2;
            else
                X = c(top:top+K.s(i)^2-1)-At(top:top+K.s(i)^2-1,:)*x;
                F = F + (reshape(X,K.s(i),K.s(i)) >= 0);
                top = top + K.s(i)^2;
            end
        end
    end
end

h = -b'*x;

str = fileread(filename);
intpos = strfind(str, '*INTEGER');
str = str(intpos + 10 : length(str)); % the +10 drops the *INTEGERS*\n as well
while true
   [intpos,str] = strtok(str, '*');
   if isempty(intpos),  break;  end
   intpos = str2num(intpos);
   F = F + integer(x(intpos));
end
global yalmiptestnnodes;
%diagnostics = optimize(F,h,sdpsettings('solver', 'cutsdp', 'cutsdp.maxtime', 3600, 'cutsdp.feastol', 1e-6, 'cplex.threads',1))
diagnostics = optimize(F,h,sdpsettings('bnb.maxtime', 3600, 'bnb.inttol', 1e-6, 'bnb.feastol', 1e-6, 'bnb.gaptol', 1e-6, 'bnb.prunetol', 1e-9, 'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS', 1e-7, 'mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS', 1e-7, 'mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS', 1e-7, 'mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED', 1e-5, 'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1e-5, 'mosek.MSK_IPAR_NUM_THREADS', 1, 'mosek.MSK_IPAR_INTPNT_MULTI_THREAD', 0))
%value(x)
time = diagnostics.yalmiptime + diagnostics.solvertime;
obj = value(h);
nodes = yalmiptestnnodes
yalmip('clear')