function [ k ] = MatrixCard( A, tol )
% berechnet Anzahl der nicht-null (bzw. abs(A_ij) > tol) Einträge der Matrix A
k=0;
for i=1:1:length(A(:,1))
    for j=1:1:length(A(1,:))
        if abs(A(i,j))>tol
            k=k+1;
        end
    end
end
end

