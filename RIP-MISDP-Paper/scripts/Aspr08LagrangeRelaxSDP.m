function [ X,L] = Aspr08LagrangeRelaxSDP( A,k,side,rho,alpha,outputfile )
%berechnet einzelne SDP-Relaxierung für Lagrange Relaxierung der RIP wie in d'Aspremont 2008 
%A,k,rho,alpha wie dort, side <- {'l','r'}: linke oder rechte Seite der RIP RIP

m=length(A(:,1));
n=length(A(1,:));

if side=='l'
    C=sqrtm(alpha*eye(n)-transpose(A)*A);
else
    C=sqrtm(transpose(A)*A);
end
B=zeros(n,n,n);
if side=='l'
    for i=1:1:n
       B(i,:,:)=C(:,i)*transpose(C(:,i))-rho*eye(n,n);
    end
else
    for i=1:1:n
       B(i,:,:)=C(:,i)*transpose(C(:,i))-rho*eye(n,n);
    end
end

%% berechne Schranke für Lagrange-Relaxierung
% X=(x_11, x_12, ..., x_nn, p1_11, p1_12, ..., p1_nn, p2_11, ..., pn_nn, (x-p1)_11, ..., (x-pn)_nn)
A=zeros(1+n^3,(2*n+1)*n^2);% 1. Zeile Spur=1, restliche Zeilen sorgen dafür, dass hintere Variablen tatsächlich (x-p_l)_ij
b=zeros(1+n^3,1);
C=zeros((2*n+1)*n^2,1);
for i=1:1:n
    A(1,(i-1)*n+i)=1; % SpurX=1
end
b(1)=1;
z=1;
for l=1:1:n;
    for i=1:1:n
        for j=1:1:n % [Variable für (x-p_l)_ij]-x_ij +(p_l)_ij = 0 
            z=z+1;
            A(z,n^3+n^2+(l-1)*n^2+(i-1)*n+j)=1; % Variable für (x-p_l)_ij
            A(z,(i-1)*n+j)=-1; %Variable für x_ij
            A(z,l*n^2+(i-1)*n+j)=1; %Variable für (p_l)_ij            
        end
    end
end
for l=1:1:n
    for i=1:1:n
        for j=1:1:n
            C(l*n^2+(i-1)*n+j) = -1*B(l,j,i); % ZF: sum Tr(BiPi), also ZF-Vorfaktor von (Pl)_ij = (Bl)_ji, -1 da SeDuMi minimiert statt maximiert
        end
    end
end
T=[n];
for l=1:1:n % ein SDP-Block für X, dann für l=1:1:n zwei weitere für P_l und X-P_l
    T=[T n n];
end
K.s = T;
pars.fid= fopen(outputfile, 'a');
[X,y,info]=sedumi(A,b,C,K, pars);
fclose(pars.fid);
pars.fid=1;
if side=='l'
    L=transpose(C)*X+alpha-rho*k; % transpose(C)*X ist größter (sparse) Eigenwert von -(alpha ID - AA^T), deshalb+ alpha für größten EW von AA^T, - rho*k für eigentliche Schranke
end
if side=='r'
    L=-1*transpose(C)*X+rho*k;
end
