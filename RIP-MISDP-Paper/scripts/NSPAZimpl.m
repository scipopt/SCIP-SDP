function [] = NSPAZimpl(A, k, file)
%schreibt zpl-File um exaktes alpha_k für NSP für gegebene Matrix A und Ordnung k als MINLP zu berechnen
fid = fopen(file, 'w');
m=length(A(:,1));
n=length(A(1,:));
t= [m n k];
s= '#nullspace-property, compute alpha_k for : m= %-4.0f, n= %-4.0f, order k= %-4.0f';
s= [s '\n'];
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
        S=[S '%8.8g'];
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
S=[S 'var y[N] binary;\n'];
S=[S 'var alpha real; \n'];

%% Zielfunktion
S=[S 'maximize ZF: alpha;\n']; %Zimpl erlaubt Nichtlinearität nur in Nebenbedingungen, nicht in Zielfunktion

%% Nebenbedinungen
S=[S 'subto helpvar: sum<j> in N: x[j]*y[j]==alpha;\n'];
S=[S 'subto nullspace: forall <i> in M: sum <j> in N: A[i,j]*x[j] == 0;\n'];
S=[S 'subto norm: sum <j> in N: abs(x[j]) == 1;\n'];
T=[T k];
S=[S 'subto cardinality: sum <j> in N: y[j] <= %-4.0f;'];  

%% Schreiben
fprintf(fid, S, T);
fclose(fid);
end
