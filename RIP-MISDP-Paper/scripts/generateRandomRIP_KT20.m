% Generate random RIP instances as used by Kobayashi and Takano (2020): A
% branch-and-cut algorithm for solving mixed-integersemidefinite
% optimization problems

valsm = [10];
valsn = [10,15,20,25,30,35];
valsk = [3,5];


seed = 1234;
instances = 5; % number of instances per (type,m,n,k)
types = ['norm';];

mkdir Matrices_KT20
mkdir MISDP_KT20

% initialize random number generator with given seed
rng('default')
rng(seed)

m = 10;
type = 'norm';
for n=valsn
    for k=valsk
        for instance=1:instances
            file = sprintf('%s%d%d%d%s',type,m,n,k,char(instance+64));
            fid = fopen(strcat('Matrices_KT20/',file),'w');

            % generate random matrix A depending on type
            fprintf(fid,'randomization = N(0,1), seed = %d\n',seed);
            A = randn(m,n);
    
            % write matrix A to file
            fprintf(fid,'m = %d  n = %d  k = %d\n',m,n,k);
            for i=1:m
                for j=1:n-1
                    fprintf(fid,"%.15g ",A(i,j));
                end
                fprintf(fid,"%.15g\n",A(i,n));
            end
            fclose(fid);
    
            % 2. CBF
            name = strcat('MISDP_KT20/',file,'_MISDPl_KT20.cbf');
            RIPCBFdual(A,k,'l',name,0,0,0,0,2);
            name = strcat('MISDP_KT20/',file,'_MISDPr_KT20.cbf');
            RIPCBFdual(A,k,'r',name,0,0,0,0,2);
        end
    end        
end



