function [] = RIPSDPA(A, k, side, file, Rank)
%schreibt SDP-File für ganzzahlige RIP-SDP-Relaxierung für Matrix A, Ordnung k, schreibt in 'file'
%side ='l' für linke Seite/alpha_k, side='r' für rechte Seite/beta_k
%Rank=1 für zusätzliche Rang-NB, sonst Rank=0
fid = fopen(file, 'w');
m=length(A(:,1));
n=length(A(1,:));
t= [m n k];
if side=='l'
    s= '"RIP-SDP-Relaxation for lower bound, Matrix-Dimensions: m= %-4.0f, n= %-4.0f, order k= %-4.0f, left-hand side/lower bound ';
else
    s= '"RIP-SDP-Relaxation for lower bound, Matrix-Dimensions: m= %-4.0f, n= %-4.0f, order k= %-4.0f, right-hand side/upper bound ';
end
s= [s '\n'];
fprintf(fid, s, t);
T=[(n+1)*n*0.5+n 2 n -1*(2+2*n^2+2*n)]; %Anzahl Variablen, Blöcke und Blockgrößen
S= '%-9.0f\n%-9.0f\n%-9.0f%-9.0f\n';
%% Zielfunktionswerte
Z=zeros(n);
B=transpose(A)*A
if side=='l'
    for i=1:1:n
        for j=i:1:n
            if i==j
                T=[T B(j,i)];
                S= [S '%-9.5g '];  
            else
                T=[T B(j,i)+B(i,j)];
                S= [S '%-9.5g '];  
            end
        end
    end
else
    for i=1:1:n
        for j=i:1:n
            if i==j
                T=[T -1*B(j,i)];
                S= [S '%-9.5g '];  
            else
                T=[T -1*(B(j,i)+B(i,j))];
                S= [S '%-9.5g '];  
            end
        end
    end
end
for j=1:1:n
    T = [T 0];
     S= [S '%-9.2g '];
end
%% X positiv semidefinit
for i=1:1:n 
    for j=i:1:n
    T=[T RIPPos(n,i,j) 1 i j 1];
    S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
    end
end
%% Spur von X >=1
if side=='l' %sonst durch Maximierung sichergestellt
T=[T 0 2 1 1 1];
S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
for i=1:1:n
    T=[T RIPPos(n,i,i) 2 1 1 1];
    S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
end
end
%% Spur von -X >=-1
if side=='r' %sonst durch Minimierung sichergestellt
T=[T 0 2 1 1 -1];
S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
for i=1:1:n
    T=[T RIPPos(n,i,i) 2 1 1 -1];
    S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
end
end
%% -Summe z_j >= -k
T=[T 0 2 2 2 -k];
S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
for j=1:1:n
    T=[T (n+1)*n*0.5+j 2 2 2 -1];
    S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
end
%% X(i,j) >= -z_j & -X(i,j) >= -z_j
z=1;
for i=1:1:n
    for j=i:1:n
        if i==j
            T=[T RIPPos(n,i,j) 2 2+z 2+z 1];
            S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
            T=[T (n+1)*n*0.5+j 2 2+z 2+z 1];
            S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
            z=z+1;
            T=[T RIPPos(n,i,j) 2 2+z 2+z -1];               
            S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
            T=[T (n+1)*n*0.5+j 2 2+z 2+z 1];
            S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];    
            z=z+1;
        else
            T=[T RIPPos(n,i,j) 2 2+z 2+z 1];
            S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
            T=[T (n+1)*n*0.5+j 2 2+z 2+z 1];
            S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
            z=z+1;
            T=[T RIPPos(n,i,j) 2 2+z 2+z -1];
            S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
            T=[T (n+1)*n*0.5+j 2 2+z 2+z 1];
            S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];    
            z=z+1;
            T=[T RIPPos(n,i,j) 2 2+z 2+z 1];                 %%% Bedingungen für Einträge unterhalb der Diagonalen
            S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
            T=[T (n+1)*n*0.5+i 2 2+z 2+z 1];
            S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
            z=z+1;
            T=[T RIPPos(n,i,j) 2 2+z 2+z -1];
            S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
            T=[T (n+1)*n*0.5+i 2 2+z 2+z 1];
            S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];    
            z=z+1;
        end
    end
end
%% 0 <=z_j<=1
for j=1:1:n
    T=[T (n+1)*n*0.5+j 2 2+z 2+z 1];
    S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
    z=z+1;
    T=[T 0 2 2+z 2+z -1];
    S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];
    T=[T (n+1)*n*0.5+j 2 2+z 2+z -1];
    S= [S '\n%9.8g %9.8g %9.8g %9.8g %9.8g'];    
    z=z+1;
end
%% Ganzzahligkeit
S=[S '\n*INTEGER'];
for j=1:1:n
    T=[T (n+1)*n*0.5+j];
    S=[S '\n*%-9.0f'];
end

%% Rank 1
if Rank == true
    S = [S '\n*Add RANKONE_REAL*\n*R1'];
end

%% Schreiben
fprintf(fid, S, T);
fclose(fid);
end

