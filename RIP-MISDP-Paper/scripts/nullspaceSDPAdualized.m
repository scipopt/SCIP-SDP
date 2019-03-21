function [] = nullspaceSDPAdualized(A, k, file, Integer, Rank) 
%schreibt SDP-File für NSP-SDP-Relaxierung nach d'Aspremont 11 für gegebene Matrix A, Ordnung k, schreibt in 'file'
%Optionen:
%Integer = 1: alle Y Variablen müssen ganzzahlig sein
%Rank = 1: zusätzliche Rangbedingung
fid = fopen(file, 'w');
P=null(A);
m=length(A(:,1));
n=length(A(1,:));
p=length(P(1,:));
N=p+n; %Matrix(X,Z^T,Z,Y) hat Dimension (p+n)x(p+n)=NxN
    t= [m n p k];
    s= '"nullspace-SDP-Relaxation, Matrix-Dimensions: m=%-4.0f, n=%-4.0f, p=%-3.0f, order k=%-4.0f"';
    s= [s '\n'];
        fprintf(fid, s, t);
T=[(N+1)*N*0.5+n*(n+1)*0.5+n*(n+1)*0.5+n^2 2 N -1*(1+(n+1)*n*0.5+2+2*(n*(n+1)*0.5)+2*(n*(n+1)*0.5)+2*n^2)]; %Anzahl Variablen, Blöcke und Blockgrößen
    %%% Variablen: (N+1)*N*0.5 für Matrix an sich, dann n*(n+1)*0.5 für Beträge
    %%% von PXP^T, dann n*(n+1)*0.5 für Beträge von Y, dann n² für Beträge von PZ^T
S= '%-4.0f\n%-4.0f\n%-4.0f%-4.0f\n';
P=null(A);
%% Zielfunktionswerte
Z=zeros(1,(N+1)*N*0.5);
for l=1:1:((N+1)*N*0.5) 
    [i,j]=MatrixPos(N,l-1); %% -1 um MatrixPos benutzen zu können, das für SCIP-Nummerierung 0, 1, ... gedacht ist.
    if j > p && i <= p %%falls in Z^T
        Z(1,l) = -1*P((j-p),i); %%% j-p um auf Position in Z statt kompletter Matrix zu kommen
    else
        Z(1,l)=0;
    end
   T = [T Z(1,l)];
   S= [S '%-4.4g '];
end
for l=1:1:n*(n+1)*0.5+n*(n+1)*0.5+n^2
    T = [T 0];
     S= [S '%-4.2g '];
end
%% Matrix positiv semidefinit
l=1;
for i=1:1:N 
    for j=i:1:N
    T=[T l 1 i j 1];
    S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];
    l=l+1;
    end
end
%% |PXP^T|_1 <= 1
T=[T 0 2 1 1 -1];
S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];
i=1;
for l=1:1:(n*(n+1)*0.5)
    if l==RIPPos(n,i,i)
        T=[T 0.5*(N+1)*N+l 2 1 1 -1];
        i=i+1;
    else
        T=[T 0.5*(N+1)*N+l 2 1 1 -2];
    end
    S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];    
end

%% |Y|_\infty <= 1
z=2;
for l=1:1:(n*(n+1)*0.5)
        T=[T 0 2 z z -1];
        S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];
        T=[T (N+1)*N*0.5+n*(n+1)*0.5+l 2 z z -1];
        S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];
        z=z+1;
end


%% |Y|_1 <= k²
T=[T 0 2 (n*(n+1)*0.5)+2 (n*(n+1)*0.5)+2 -1*(k^2)];
S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];
i=1;
for l=1:1:(n*(n+1)*0.5)
    if l==RIPPos(n,i,i)
        T=[T (N+1)*N*0.5+n*(n+1)*0.5+l 2 (n*(n+1)*0.5)+2 (n*(n+1)*0.5)+2 -1];
        i=i+1;
    else
        T=[T (N+1)*N*0.5+n*(n+1)*0.5++l 2 (n*(n+1)*0.5)+2 (n*(n+1)*0.5)+2 -2];
    end
    S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];    
end

%% |PZ^T|_1 <= k
T=[T 0 2 (n*(n+1)*0.5)+3 (n*(n+1)*0.5)+3 -1*k];
S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];
for l=1:1:n^2
            T=[T (N+1)*N*0.5+n*(n+1)*0.5+n*(n+1)*0.5+l 2 (n*(n+1)*0.5)+3 (n*(n+1)*0.5)+3 -1];
            S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];
end

%% zweite Hälfte der Variablen im X-Teil Beträge von PXP^T_{ij}
z=1;
coef=0;
for i=1:1:n 
    for j=i:1:n
        T=[T (N+1)*N*0.5+RIPPos(n,i,j) 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z 1];
        S= [S '\n%8.8g %8.8g %8.8g %8.8g %8.8g'];           
        for s=1:1:p % Betrag-Eintrag >=0
            for t=s:1:p
                if s==t
                    coef = -1*P(j,t)*P(i,s);
                    if (coef <-0.0001 | coef >0.0001) % don't write if 0 would be written
                        T=[T nullspacePosdualized(p,n,s,t,'X') 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z coef];
                        S= [S '\n%8.8g %8.8g %8.8g %8.8g %8.8g'];
                    end
                else
                    coef = -1*(P(j,t)*P(i,s)+P(i,t)*P(j,s));
                    if (coef <-0.0001 | coef >0.0001) % don't write if 0 would be written
                        T=[T nullspacePosdualized(p,n,s,t,'X') 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z coef];
                        S= [S '\n%8.8g %8.8g %8.8g %8.8g %8.8g'];
                    end
                end
            end
        end
        z=z+1;
        T=[T (N+1)*N*0.5+RIPPos(n,i,j) 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z 1];
        S= [S '\n%8.8g %8.8g %8.8g %8.8g %8.8g'];           
        for s=1:1:p  % Betrag+Eintrag >=0, also Betrag >= -Eintrag
            for t=s:1:p
                if s==t
                    coef = P(j,t)*P(i,s);
                    if (coef <-0.0001 | coef >0.0001)  % don't write if 0 would be written
                        T=[T nullspacePosdualized(p,n,s,t,'X') 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z coef];
                        S= [S '\n%8.8g %8.8g %8.8g %8.8g %8.8g'];
                    end
                else
                    coef = (P(j,t)*P(i,s)+P(j,s)*P(i,t));
                    if (coef <-0.0001 | coef >0.0001)   % don't write if 0 would be written                
                        T=[T nullspacePosdualized(p,n,s,t,'X') 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z coef];
                        S= [S '\n%8.8g %8.8g %8.8g %8.8g %8.8g'];
                    end
                end
            end
        end
        z=z+1;
    end
end


%% zweite Hälfte der Variablen im Y-Teil Beträge von Y
for i=1:1:n
    for j=i:1:n
        T=[T nullspacePosdualized(p,n,i,j,'Y') 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z -1]; % Betrag-Eintrag >=0
        S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];
        T=[T (N+1)*N*0.5+n*(n+1)*0.5+RIPPos(n,i,j) 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z 1];
        S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];
        z=z+1;
        T=[T nullspacePosdualized(p,n,i,j,'Y') 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z 1]; %Betrag+Eintrag >=0, also Betrag >= -Eintrag
        S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];
        T=[T (N+1)*N*0.5+n*(n+1)*0.5+RIPPos(n,i,j) 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z 1];
        S= [S '\n%8.4g %8.4g %8.4g %8.4g %8.4g'];
        z=z+1;
    end
end


%% zweite Hälfte der Variablen im Z-Teil Beträge von PZ^T
for i=1:1:n 
    for j=1:1:n
        T=[T (N+1)*N*0.5+n*(n+1)*0.5+n*(n+1)*0.5+(i-1)*n+j 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z 1];
        S= [S '\n%8.8g %8.8g %8.8g %8.8g %8.8g'];  
        for l=1:1:p % Betrag-Eintrag >=0
            coef = -1*P(i,l);
            if (coef <-0.0001 | coef >0.0001)  % don't write if 0 would be written
                T=[T nullspacePosdualized(p,n,j,l,'Z') 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z coef];   
                S= [S '\n%8.8g %8.8g %8.8g %8.8g %8.8g'];
            end
        end
        z=z+1;
        T=[T (N+1)*N*0.5+n*(n+1)*0.5+n*(n+1)*0.5+(i-1)*n+j 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z 1];
        S= [S '\n%8.8g %8.8g %8.8g %8.8g %8.8g'];  
        for l=1:1:p % Betrag-Eintrag >=0
            coef = P(i,l);
            if (coef <-0.0001 | coef >0.0001)      % don't write if 0 would be written       
                T=[T nullspacePosdualized(p,n,j,l,'Z') 2 1+(n+1)*n*0.5+2+z 1+(n+1)*n*0.5+2+z coef];   
                S= [S '\n%8.8g %8.8g %8.8g %8.8g %8.8g'];
            end
        end
        z=z+1;
    end
end

        
%% Ganzzahligkeit
if Integer == true
    S=[S '\n*INTEGER'];
    for i=1:1:n
        for j=i:1:n
            T=[T nullspacePosdualized(p,n,i,j,'Y')];
            S=[S '\n*%-4.0f'];
        end
    end
end

%% Rank 1
if Rank == true
    S = [S '\n*Add RANKONE_REAL*\n*R1'];
end

%% Schreiben
fprintf(fid, S, T);
fclose(fid);
end
