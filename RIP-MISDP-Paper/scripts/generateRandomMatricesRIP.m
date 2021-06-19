valsm = [15,25,30];
valsn = [30,35,40];
valsk = [5,4,3];

seed = 1234;
instances = 3; % number of instances per (type,m,n,k)
types = ['0+-1';'band';'bern';'bina';'norm';'rnk1';'wish'];

% mkdir Asp07

% initialize random number generator with given seed
rng('default')
rng(seed)

matdir = 'Matrices_testset';
datadir = 'MISDP_testset';

for t=1:length(types(:,1))
    type = types(t,:);
    for i=1:length(valsm)
        if strcmp(type,'band')
            m = valsn(i);
            bandwidth = 1 + 2*i;
        elseif strcmp(type,'rnk1')
            m = valsn(i);
            bandwidth = 0;
        else
            m = valsm(i);
            bandwidth = 0;
        end
        n = valsn(i);
        k = valsk(i);
        generateRIPA(m,n,k,seed,instances,type,bandwidth,matdir,datadir);
    end
end
