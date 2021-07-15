valsm = [15,25,30];
valsn = [30,35,40];
valsk = [5,4,3];

% valsm = [40,50,60];
% valsn = [60,70,80];
% valsk = [5,4,3];

seed = 1234;
instances = 5; % number of instances per (type,m,n,k)
types = ['0+-1';'band';'bern';'bina';'norm';'rnk1';'wish'];

% mkdir Asp07

% initialize random number generator with given seed
rng('default')
rng(seed)

matdir = 'RIP_EURO21_Standard_Matrices';
datadir = 'RIP_EURO21_Standard_Instances';

% Should pre-generated matrices be read?
readmatrices = false;

cnt = 0;
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
        tmp = generateRIPA(m,n,k,seed,instances,type,bandwidth,matdir,datadir,readmatrices);
        cnt = cnt + tmp;
%         fid1 = fopen(strcat(datadir,'/',testsetname,upper(type(1)),type(2:end),'L.test'),'w');
%         fid2 = fopen(strcat(datadir,'/',testsetname,upper(type(1)),type(2:end),'R.test'),'w');
%         str1 = strcat(datadir,'/',type,'*MISDPl*.cbf*');
%         str2 = strcat(datadir,'/',type,'*MISDPr*.cbf*');
%         content1 = dir(str1);
%         content2 = dir(str2);
%         cnt1 = 0;
%         cnt2 = 0;
%         for idx = 1:length(content1)
%           fprintf(fid1, strcat('RIPtest/',content1(idx).name,'\n'));
%           cnt1 = cnt1 + 1;
%         end
%         for idx = 1:length(content2)
%           fprintf(fid2, strcat('RIPtest/',content2(idx).name,'\n'));
%           cnt2 = cnt2 + 1;
%         end
%         fclose(fid1);
%         fclose(fid2);
    end
end

fprintf("Number of created instances: %d\n", cnt);

% fid3 = fopen(strcat(datadir,'/',testsetname,'All.test'),'w');
% str3 = strcat(datadir,'/*.cbf*');
% content3 = dir(str3);
% cnt3 = 0;
% for idx = 1:length(content3)
%   fprintf(fid3, strcat('RIPtest/',content3(idx).name,'\n'));
%   cnt1 = cnt1 + 1;
% end
% fclose(fid3);
% 
% fid4 = fopen(strcat(datadir,'/',testsetnameSDPA,'All.test'),'w');
% str4 = strcat(datadir,'/*.dat-s*');
% content4 = dir(str4);
% cnt4 = 0;
% for idx = 1:length(content4)
%   fprintf(fid4, strcat('RIPtest/',content4(idx).name,'\n'));
%   cnt1 = cnt1 + 1;
% end
% fclose(fid4);

% TODO: gzip files!
