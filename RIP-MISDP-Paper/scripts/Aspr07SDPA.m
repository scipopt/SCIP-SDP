function [] = Aspr07SDPA(A, k, side, file, Rank)
%schreibt SDP-File der RIP-SDP-Relaxierung nach d'Aspremont 2007 für obere/untere Schranke, gegebene Matrix A, Ordnung k, schreibt in 'file'
%Optionen:
%side='l': linke Seite/untere Schranke
%side='r': rechte Seite/obere Schranke
%Rank = 1 fügt Rangebedingung hinzu

m=length(A(:,1));
n=length(A(1,:));
fid = fopen(file, 'w');
    t= [m n k];
    if side=='l'
        s= '"RIP-SDP-Relaxation by dAspremont 2007 for lower bound, Matrix-Dimensions: m= %-4.0f, n= %-4.0f, order k= %-4.0f';
    end
    if side=='r'
        s='"RIP-SDP-Relaxation by dAspremont 2007 for upper bound, Matrix-Dimensions: m= %-4.0f, n= %-4.0f, order k= %-4.0f';
    end
    s= [s '\n'];
    fprintf(fid, s, t);
    
Sigma = transpose(A)*A;

T=[(n+1)*n 2 n -1*(2+n*(n+1))]; %Anzahl Variablen, Blöcke und Blockgrößen
S= '%-4.0f\n%-4.0f\n%-4.0f%-4.0f\n';

%% Zielfunktionswerte
Z=zeros(1,n^2);
if side =='l'
for i=1:1:n
    for j=i:1:n 
        if i==j
            T=[T Sigma(j,i)];
            S= [S '%-.15g '];   
        else
            T=[T Sigma(j,i) + Sigma(i,j)];
            S= [S '%-.15g '];  
        end
    end
end
end
if side =='r'
for i=1:1:n
    for j=i:1:n 
        if i==j
            T=[T -1*(Sigma(j,i))];
            S= [S '%-.15g '];   
        else
            T=[T -1*(Sigma(j,i) + Sigma(i,j))];
            S= [S '%-.15g '];  
        end
    end
end
end
for i=1:1:(n+1)*n*0.5
        T=[T 0];
        S= [S '%-.15g '];     
end

%% X positiv semidefinit
for i=1:1:n 
    for j=i:1:n
    T=[T RIPPos(n,i,j) 1 i j 1];
    S= [S '\n%.15g %.15g %.15g %.15g %.15g'];
    end
end

%% Spur von X >=1 nur falls minimiert wird
if side=='l'
T=[T 0 2 1 1 1];
S= [S '\n%.15g %.15g %.15g %.15g %.15g'];
for i=1:1:n
    T=[T RIPPos(n,i,i) 2 1 1 1];
    S= [S '\n%.15g %.15g %.15g %.15g %.15g'];
end
end
%% Spur von -X >=-1 nur falls maximiert wird
if side=='r'
T=[T 0 2 1 1 -1];
S= [S '\n%.15g %.15g %.15g %.15g %.15g'];
for i=1:1:n
    T=[T RIPPos(n,i,i) 2 1 1 -1];
    S= [S '\n%.15g %.15g %.15g %.15g %.15g'];
end
end
%% 1^T*|X|*1 <= k <=> Summe aller Beträge <= k
T=[T 0 2 2 2 -k];
S= [S '\n%.15g %.15g %.15g %.15g %.15g'];
for i=1:1:n*(n+1)*0.5
    T=[T n*(n+1)*0.5+i 2 2 2 -1];
    S= [S '\n%.15g %.15g %.15g %.15g %.15g'];
end

%% zweite Hälfte der Variablen Beträge der ersten Hälfte
z=1;
for i=1:1:n*(n+1)*0.5
    T=[T i 2 2+z 2+z -1];
    S= [S '\n%.15g %.15g %.15g %.15g %.15g'];
    T=[T n*(n+1)*0.5+i 2 2+z 2+z 1];                  % Betrag-Eintrag >=0
    S= [S '\n%.15g %.15g %.15g %.15g %.15g']; 
    z=z+1;
    T=[T i 2 2+z 2+z 1];
    S= [S '\n%.15g %.15g %.15g %.15g %.15g'];
    T=[T n*(n+1)*0.5+i 2 2+z 2+z 1];                  % Betrag+Eintrag >=0, also Betrag >= -Eintrag
    S= [S '\n%.15g %.15g %.15g %.15g %.15g']; 
    z=z+1;
end

%% Rank 1
if Rank == true
    S = [S '\n*Add RANKONE_REAL*\n*R1'];
end


fprintf(fid, S, T);
fclose(fid);
end

