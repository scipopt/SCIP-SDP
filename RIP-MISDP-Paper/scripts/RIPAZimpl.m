function [] = RIPAZimpl(A, k, side, file)
%schreibt zpl-File exaktes um alpha_k/beta_k der RIP für gegebene Matrix A, Ordnung k als MINLP zu berchnen, schreibt in 'file'
%side = 'l' : linke Seite, also alpha, side = 'r': rechte Seite, also beta
fid = fopen(file, 'w');
m=length(A(:,1));
n=length(A(1,:));
t= [m n k];
if side=='l'
    s= '#RIP, compute (1-delta_k) for : m= %-4.0f, n= %-4.0f, order k= %-4.0f\n';
end
if side=='r'
    s= '#RIP, compute (1+delta_k) for : m= %-4.0f, n= %-4.0f, order k= %-4.0f\n';
end
fprintf(fid, s, t);

%% Definitionen
T=[m n];
S=['set M := {1 .. %-4.0f };\n'];
S=[S 'set N := {1 .. %-4.0f };\n'];
S=[S 'set MN := M * N;\n'];
S=[S 'param A[MN] := \n      |'];
for j=1:1:n
    T=[T j];
    S=[S '     %-4.0f'];
    if j<n
        S=[S ','];
    end
end
S=[S ' |\n'];
for i=1:1:m
    T=[T i];
    S=[S ' |%-4.0f| '];
    for j=1:1:n
        T=[T A(i,j)];
        S=[S '%8.4g'];
        if j < n
            S=[S ', '];
        else if i<m
            S=[S ' |\n'];
            else
                S=[S ' |;\n'];   
            end
        end
    end
end
S=[S 'var x[N] real >= -infinity; \n'];
S=[S 'var z[N] binary;\n'];
if side=='l'
    S=[S 'var alpha real; \n'];
else
    S=[S 'var beta real; \n'];
end

%% Zielfunktion
if side=='l'
    S=[S 'minimize ZF: alpha;\n']; %Zimple erlaubt Nichtlinearität nur in Nebenbedingungen, nicht in Zielfunktion
end
if side=='r'
    S=[S 'maximize ZF: beta;\n']; %Zimple erlaubt Nichtlinearität nur in Nebenbedingungen, nicht in Zielfunktion
end

%% Nebenbedinungen
if side=='l'
    S=[S 'subto helpvar: sum<i> in M: abs(sum<j> in N: A[i,j]*x[j])^2 == alpha;\n'];
else
    S=[S 'subto helpvar: sum<i> in M: abs(sum<j> in N: A[i,j]*x[j])^2 == beta;\n'];
end
S=[S 'subto lb: forall<j> in N: x[j] >= -1*z[j];\n'];
S=[S 'subto ub: forall<j> in N: x[j] <= z[j];\n'];
S=[S 'subto norm: sum<j> in N: abs(x[j])^2 == 1;\n'];
T=[T k];
S=[S 'subto cardinality: sum <j> in N: z[j] <= %-4.0f;'];  

%% Schreiben
fprintf(fid, S, T);
fclose(fid);
end
