function [] = generateRIPA(m,n,k,seed,instances,type,bandwidth)
%generiert zufällige Matrizen zum Berechnen der RIP
%Optionen:
%m = Zeilen der zufälligen Matrix
%n = Spalten der zufälligen Matrix
%k = Ordnung der RIC
%seed = random seed zur Reproduzierbarkeit der Matrizen
%instances = Anzahl der zu erzeugenden Matrizen
%type = '0+-1' : Einträge in 0 (prob 2/3), +/-sqrt(3/m) (prob 1/6)
%type = 'band' : Band-Matrix mit bandwidth, Einträge Unif({0,1})
%type = 'bern' : Einträge Unif({+/-1/sqrt(m)})
%type = 'bina' : Einträge Unif({0,1})
%type = 'norm' : Einträge N(0,1)
%type = 'rnk1' : rank1-Matrix A = aa^T, Einträge von a in N(0,1)
%type = 'wish' : Einträge in N(0,1/m)
%bandwidth = Bandbreite, falls type = 'band'

% initialize random number generator with given seed
rng('default')
rng(seed)

for instance=1:instances
    file = sprintf('%s%d%d%d%s',type,m,n,k,char(instance+64));
    fid = fopen(strcat('Matrices/',file),'w');
    
    % generate random matrix A depending on type
    switch type
        case '0+-1'
            fprintf(fid,'randomization = P(+sqrt(3/m))=1/6, P(0)=2/3, P(-sqrt(3/m))=1/6\n');
            A = randi([0 5],m,n);
            for i=1:m
                for j=1:n
                    if A(i,j) < 4
                        A(i,j) = 0.0;
                    elseif A(i,j) == 4
                        A(i,j) = sqrt(3/m);
                    elseif A(i,j) == 5
                        A(i,j) = -sqrt(3/m);
                    else
                        fclose(fid);
                        error('Error: Unexpected Case!')
                    end
                end
            end
        case 'band'
            if m~=n
                fclose(fid);
                error('Error: "type = band", but "m ~= n"')
            end
            if ~exist('bandwidth','var')
                fclose(fid);
                error('Error: "type=band", but no bandwidth specified!')
            end
            if mod(bandwidth,2) ~= 1
                fclose(fid);
                error('Error: "type=bandwidth", but bandwidth is no odd integer!')
            end
            fprintf(fid,'mxm band matrix with entries uniformly in {0,1}, bandwith= %d\n',bandwidth);
            A = randi([0 1],m,n);
            for i=1:m
                for j=1:n
                    if j < i - (bandwidth - 1)/2 || j > i + (bandwidth - 1)/2
                        A(i,j) = 0;
                    end
                end
            end
        case 'bern'
            fprintf(fid,'randomization = uniformly in +/- 1/sqrt(m)\n');
            A = randi([0 1],m,n) * 2/sqrt(m) - 1/sqrt(m);
        case 'bina'
            fprintf(fid,'randomization = uniformly in 0/1\n');
            A = randi([0 1],m,n);
        case 'norm'
            fprintf(fid,'randomization = N(0,1)\n');
            A = randn(m,n);
        case 'rnk1'
            if m~=n
                fclose(fid);
                error('Error: "type = rnk1", but "m ~= n"')
            end
            fprintf(fid,'rank 1 Matrix aa^T, randomization for a = N(0,1)\n');
            a = randn(m,1);
            A = a*a';
        case 'wish'
            fprintf(fid,'randomization = N(0,1/m)\n');
            A = 1/sqrt(m).*randn(m,n);
        otherwise
            fclose(fid);
            error('Error: Input for "type" not recognized!')
    end
    
    % write matrix A to file
    fprintf(fid,'m = %d  n = %d  k = %d\n',m,n,k);
    for i=1:m
        for j=1:n-1
            fprintf(fid,"%.15g ",A(i,j));
        end
        fprintf(fid,"%.15g\n",A(i,n));
    end
    fclose(fid);
    
    % generate SDPA-file of MISDP-formulation
    RIPSDPA(A,k,'l',strcat('MISDP/',file,'_MISDPl.dat-s'),0);
    RIPSDPA(A,k,'r',strcat('MISDP/',file,'_MISDPr.dat-s'),0);
    
    % generate SDPA-file of SDP-Relaxation from d'Aspremont 2007
    Aspr07SDPA(A,k,'l',strcat('Asp07/',file,'_Asp07l.dat-s'),0);
    Aspr07SDPA(A,k,'r',strcat('Asp07/',file,'_Asp07r.dat-s'),0);
    
end
end
        
        