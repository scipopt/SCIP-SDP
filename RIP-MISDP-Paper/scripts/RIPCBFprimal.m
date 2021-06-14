function [] = RIPCBFprimal(A, order, side, file, Rank, socp, strgbnds, trineq, boundver, sumineq)
% schreibt SDP-File für ganzzahlige RIP-SDP-Relaxierung in primaler Form
% (mit PSD-Variablen) für Matrix A, Ordnung k, schreibt in 'file' 
% side ='l' für linke Seite/alpha_k, side='r' für rechte Seite/beta_k
% Rank = 1 für zusätzliche Rang-NB, sonst Rank = 0
% socp = Variante der zusätzlichen gültigen SOCP-Ungleichungen von Li/Xie
%        0: keine zusätzlichen gültigen SOCP-Ungleichungen
%        1: gültige SOCP-Ungleichungen als SDP-Constraint
%        2: gültige SOCP-Ungleichungen als SOCP-Constraint
% strgbnds = 1 falls für nichtnegative Matrizen A die Schranke 0 <= X_{ij}
% statt -z_j <= X_{ij} für alle i,j benutzt werden soll, sonst = 0
% trineq = 1 falls trace(X) <= 1 (rechte Seite) bzw trace(X) >= 1 (linke
% Seite) statt trace(X) == 1 geschrieben werden soll, sonst 0.
% boundver = Version der Constraints -z_j <= X_{ij} <= z_j, i,j = 1,...,n
%            0: -z_i <= X_{ii} <= z_i, i = 1,...,n (schwächste Variante)
%            1: Standard, -z_j <= X_{ij} <= z_j, i,j = 1,...,n
%            2: -0.5*z_j <= X_{ij} <= 0.5*z_j, i ~= j (stärkste Variante)
% sumineq = 1 falls die gültige Ungleichung sum_{i\neq j} X_{ij} \leq k-1
% hinzugefügt werden soll, sonst 0.
% ACHTUNG: schreibt untere Dreiecksmatrizen!

    fid = fopen(file, 'w');
    m=length(A(:,1));
    n=length(A(1,:));
    
    %% check input parameters
    if socp ~= 0 && socp ~= 1 && socp ~= 2
        error("Error: Option <%s> for parameter <socp> not valid!\n", socp);
    end
    if strgbnds ~= 0 && strgbnds ~= 1
        error("Error: Option <%s> for parameter <strgbnds> not valid!\n", strgbnds);
    end
    if trineq ~= 0 && trineq ~= 1
        error("Error: Option <%s> for parameter <trineq> not valid!\n", trineq);
    end
    if boundver ~= 0 && boundver ~= 1 && boundver ~= 2
        error("Error: Option <%s> for parameter <boundver> not valid!\n", boundver);
    end
    if sumineq ~= 0 && sumineq ~= 1
        error("Error: Option <%s> for parameter <strgbnds> not valid!\n", sumineq);
    end


    %% check if A is entrywise nonnegative:
    if strgbnds == 1 && ~all(A(:) >= 0)
        strgbnds = 0;
        fprintf("Setting strgbnds = 0, since matrix A is not nonnegative!\n");
    elseif strgbnds == 1 && ~side == 'r'
        strgbnds = 0;
        fprintf("Setting strgbnds = 0, since this only works for right side of RIP!\n");
    end
    
    %% SOCP-inequality is only valid for right side of the RIP
    if socp >= 1 && ~side == 'r'
        socp = 0;
        fprintf("Setting socp = 0, since this only works for right side of RIP!\n");
    end
    
    %% prepare output
    % compute B = A^T A
    B = transpose(A)*A;

    % output file
    fprintf("Writing problem to file <%s> ...\n", file);

    fprintf(fid, "VER\n1\n\n");

    if side == 'l'
        fprintf(fid, "OBJSENSE\nMIN\n\n");
    elseif side == 'r'
        fprintf(fid, "OBJSENSE\nMAX\n\n");
    else
        error("Error: Option <%s> for parameter side not valid!\n", side);
    end

    %% add scalar variables z
    fprintf(fid, "VAR\n");
    fprintf(fid, "%d 1\n", n);
    fprintf(fid, "L+ %d\n", n);
    fprintf(fid, "\n");

    %% add Matrix variables X
    fprintf(fid, "PSDVAR\n");
    if socp == 0 || socp == 2
        fprintf(fid, "1\n");
        fprintf(fid, "%d\n", n);
    elseif socp == 1
        fprintf(fid, "%d\n", n+1);
        fprintf(fid, "%d\n", n);
        for j = 0:n-1
            fprintf(fid, "%d\n", n+2);
        end
    else
        error("Something went wrong!\n");
    end
    fprintf(fid, "\n");

    %% objective
    fprintf(fid, "OBJFCOORD\n");
    fprintf(fid, "%d\n", 0.5*n*(n+1));
    for i = 0:n-1
        for j = 0:i
            fprintf(fid, "0 %d %d %.15g\n", i, j, B(i+1,j+1));
            % check that matrix is symmetric
            if ( abs(B(i+1,j+1) - B(j+1,i+1)) > 1e-6 )
                error("Error: B matrix not symmetric!\n");
            end
        end
    end

    fprintf(fid, "\n");

    %% all scalar variables z are binary
    fprintf(fid, "INT\n");
    fprintf(fid, "%d\n", n);
    for j = 0:n-1
        fprintf(fid, "%d\n", j);
    end
    fprintf(fid, "\n");

    %% add rank1-constraint
    if ( Rank == 1 )
        fprintf(fid, "PSDVARRANK1\n");
        fprintf(fid, "1\n");
        fprintf(fid, "0\n");
        fprintf(fid, "\n");
    end

    %% constraints
    fprintf(fid, "CON\n");
    nsocpcons = 0.5*n*(n+3)*(n+2);
    if socp == 2
        ncons = n*(n+2);
        ncones = n;
    elseif socp == 1
        ncons = nsocpcons;
        ncones = 1;
    else
        ncons = 0;
        ncones = 0;
    end
    if boundver == 0
        ncons = ncons + 3*n+2;
    else
        ncons = ncons + 2*n^2+n+2 - strgbnds*0.5*n*(n-1);
    end
    ncons = ncons + sumineq;
    ncones = ncones + 5 + sumineq;
    fprintf(fid, "%d %d\n", ncons, ncones);
    if socp == 1
        % SOCP constraint (as PSD cons) is the first constraint
        fprintf(fid, "L= %d\n", nsocpcons);    % \sum_j X_{ij}^2 <= X_{ii}z_i
    end
    fprintf(fid, "L+ %d\n", n);                % -z_j + 1 >= 0

    if boundver == 0
        fprintf(fid, "L- %d\n", n);            % -z_j - X_{ij} <= 0
        fprintf(fid, "L+ %d\n", n);            % z_j - X_{ij} >= 0
    else
        fprintf(fid, "L- %d\n", n^2 - strgbnds*0.5*n*(n-1));    % -z_j - X_{ij} <= 0
        fprintf(fid, "L+ %d\n", n^2);          % z_j - X_{ij} >= 0
    end
    if trineq == 0
        fprintf(fid, "L= 1\n");                % \sum_j X_{jj} == 1
    elseif trineq == 1
        if side == 'l'
            fprintf(fid, "L+ 1\n");            % \sum_j X_{jj} >= 1
        else
            fprintf(fid, "L- 1\n");            % \sum_j X_{jj} <= 1
        end
    else
        error("Error: Option <%s> for parameter trineq not valid!\n", trineq);
    end              
    fprintf(fid, "L- 1\n");                    % \sum_j z_j <= k
    if sumineq == 1
        fprintf(fid, "L- 1\n");                % \sum_{i\neq j} X_{ij} <= k-1
    end
    % SOCP constraints (as SOCP cons) are the last constraints
    if socp == 2
        for i = 0:n-1
            fprintf(fid, "Q %d\n", n+2);       % \sum_j x_{ij}^2 <= X_{ii}z_i
        end
    end
    fprintf(fid, "\n");

    %% write FCOORD
    fprintf(fid, "FCOORD\n");
    if socp == 2
        nFCOORD = n*(n+2);
    elseif socp == 1
        nFCOORD = n*(3*n+5+0.5*n*(n+1));
    else
        nFCOORD = 0;
    end
    if boundver == 0
        nFCOORD = nFCOORD + 3*n;
    else
        nFCOORD = nFCOORD + 2*n^2+n - strgbnds*0.5*n*(n-1);
    end
    nFCOORD = nFCOORD + sumineq * (n*(n-1)/2);
    fprintf(fid, "%d\n", nFCOORD);
    cnt = 0;
    conscnt = 0;
    
    % SOCP constraint
    if socp == 1
        for i = 0:n-1
            for j = 0:n
                fprintf(fid, "%d %d %d %d 1.0\n",conscnt,i+1,j+1,j+1);
                cnt = cnt + 1;
                conscnt = conscnt + 1;
            end
            fprintf(fid, "%d %d 0 0 1.0\n",conscnt,i+1);
            fprintf(fid, "%d 0 %d %d -0.5\n",conscnt,i,i);
            conscnt = conscnt + 1;
            cnt = cnt + 2;
            for j = 0:n-1
                if j == i
                    fprintf(fid, "%d %d %d 0 0.5\n",conscnt,i+1,j+1);
                else
                    fprintf(fid, "%d %d %d 0 1.0\n",conscnt,i+1,j+1); 
                end
                if j <= i
                    fprintf(fid, "%d 0 %d %d -1.0\n",conscnt,i,j);
                else
                    fprintf(fid, "%d 0 %d %d -1.0\n",conscnt,j,i);
                end
                conscnt = conscnt + 1;
                cnt = cnt + 2;
            end
            fprintf(fid, "%d %d %d 0 0.5\n",conscnt,i+1,n+1);
            fprintf(fid, "%d 0 %d %d -0.5\n",conscnt,i,i);
            conscnt = conscnt + 1;
            cnt = cnt + 2;
            for k = 0:n
                for j = 0:k-1
                    fprintf(fid, "%d %d %d %d 1.0\n",conscnt,i+1,k+1,j+1);
                    conscnt = conscnt + 1;
                    cnt = cnt + 1;
                end
            end
        end
        if ( conscnt ~= nsocpcons || cnt ~= n*(3*n+5+0.5*n*(n+1)) )
            fprintf("conscnt is = %d, should be = %d\n", conscnt, nsocpcons);
            fprintf("cnt is = %d, should be = %d\n", cnt, n*(3*n+5+0.5*n*(n+1)) );
            error("Error while writing FCOORD for SOCP constraint!\n");
        end
    end

    % no FCOORD for -z_j + 1 >= 0
    conscnt = conscnt + n;

    % add coupling constraints -z_j - X_{ij} <= 0 (if strgbnds = 1, only
    % half of them are needed for -X_{ij} <= 0)
    for i = 0:n-1
        if boundver == 0
            fprintf(fid, "%d 0 %d %d -1.0\n", conscnt, i, i);
            cnt = cnt + 1;
            conscnt = conscnt + 1;
        else
            % case i > j
            for j = 0:i-1
                fprintf(fid, "%d 0 %d %d -0.5\n", conscnt, i, j);
                cnt = cnt + 1;
                conscnt = conscnt + 1;
            end
            % case i == j
            fprintf(fid, "%d 0 %d %d -1.0\n", conscnt, i, i);
            cnt = cnt + 1;
            conscnt = conscnt + 1;
            % case i < j (only if strgbnds = 0)
            if strgbnds == 0
                for j = i+1:n-1
                    fprintf(fid, "%d 0 %d %d -0.5\n", conscnt, j, i);
                    cnt = cnt + 1;
                    conscnt = conscnt + 1;
                end
            end
        end
    end

    % add coupling constraints z_j - X_{ij} >= 0
    for i = 0:n-1
        if boundver == 0
            fprintf(fid, "%d 0 %d %d -1.0\n", conscnt, i, i);
            cnt = cnt + 1;
            conscnt = conscnt + 1;
        else
            % case i > j
            for j = 0:i-1
                fprintf(fid, "%d 0 %d %d -0.5\n", conscnt, i, j);
                cnt = cnt + 1;
                conscnt = conscnt + 1;
            end
            % case i == j
            fprintf(fid, "%d 0 %d %d -1.0\n", conscnt, i, i);
            cnt = cnt + 1;
            conscnt = conscnt + 1;
            % case i < j
            for j = i+1:n-1
                fprintf(fid, "%d 0 %d %d -0.5\n", conscnt, j, i);
                cnt = cnt + 1;
                conscnt = conscnt + 1;
            end
        end
    end

    % trace constraint
    for i = 0:n-1
        fprintf(fid, "%d 0 %d %d 1.0\n", conscnt, i, i);
        cnt = cnt + 1;
    end
    conscnt = conscnt + 1;

    % no FCOORD for sparsity constraint
    conscnt = conscnt + 1;
    
    % off-diagonal constraint
    if sumineq == 1
        for i = 0:n-1
            for j = 0:i-1
                fprintf(fid, "%d 0 %d %d 0.5\n", conscnt, i, j);
                cnt = cnt + 1;
            end
        end
        conscnt = conscnt + 1;
    end

    
    % SOCP constraints
    if socp == 2
        for i = 0:n-1
            fprintf(fid, "%d 0 %d %d 0.5\n", conscnt, i, i);
            conscnt = conscnt + 1;
            cnt = cnt + 1;
            for j = 0:n-1
                if i >= j
                    fprintf(fid, "%d 0 %d %d 1.0\n", conscnt, i, j);
                else
                    fprintf(fid, "%d 0 %d %d 1.0\n", conscnt, j, i);
                end
                conscnt = conscnt + 1;
                cnt = cnt + 1;
            end
            fprintf(fid, "%d 0 %d %d 0.5\n", conscnt, i, i);
            conscnt = conscnt + 1;
            cnt = cnt + 1;
        end
    end

    fprintf(fid, "\n");

    if ( conscnt ~= ncons || cnt ~= nFCOORD )
        fprintf("conscnt is = %d, should be = %d\n",conscnt, ncons);
        fprintf("cnt is = %d, should be = %d\n",cnt,nFCOORD);
        error("Error: Something went wrong when writing FCOORD!\n");
    end

    %% ACOORD
    fprintf(fid, "ACOORD\n");
    if socp == 1 || socp == 2
        nACOORD = 2*n;
    else
        nACOORD = 0;
    end
    if boundver == 0
        nACOORD = nACOORD + 4*n - strgbnds*n;
    else
        nACOORD = nACOORD + 2*n^2+2*n - strgbnds*(n^2);
    end
    fprintf(fid, "%d\n", nACOORD);
    cnt = 0;

    % SOCP constraint
    if socp == 1
        for i = 0:n-1
            pos = i*(0.5*(n+3)*(n+2));
            fprintf(fid, "%d %d -0.5\n",pos+n+1,i);
            cnt = cnt + 1;
            fprintf(fid, "%d %d 0.5\n",pos+2*n+2,i);
            cnt = cnt + 1;
        end
        if cnt ~= 2*n
            error("Error while writing ACOORD for SOCP constraint!\n");
        end
    end

    if socp == 1
        conscnt = nsocpcons;
    else
        conscnt = 0;
    end
    
    % add upper strgbnds on z
    for j = 0:n-1
        fprintf(fid, "%d %d -1.0\n", conscnt, j);
        conscnt = conscnt + 1;
        cnt = cnt + 1;
    end

    % add coupling constraints -z_j - X_{ij} <= 0 (if strgbnds = 1, this
    % changes to -X_{ij} <= 0, so that nothing needs to be specified here)
    if strgbnds == 0
        for i = 0:n-1
            if boundver == 0
                fprintf(fid, "%d %d -1.0\n", conscnt, i);
                cnt = cnt + 1;
                conscnt = conscnt + 1;
            else
                for j = 0:n-1
                    if ( i ~= j && boundver == 2 )
                        fprintf(fid, "%d %d -0.5\n", conscnt, j);
                    else
                        fprintf(fid, "%d %d -1.0\n", conscnt, j);
                    end
                    conscnt = conscnt + 1;
                    cnt = cnt + 1;
                end
            end
        end
    else
        if boundver == 0
            conscnt = conscnt + n;
        else
            conscnt = conscnt + 0.5*n*(n+1);
        end
    end
    
    % add coupling constraints z_j - X_{ij} >= 0
    for i = 0:n-1
        if boundver == 0
            fprintf(fid, "%d %d 1.0\n", conscnt, i);
            cnt = cnt + 1;
            conscnt = conscnt + 1;
        else
            for j = 0:n-1
                if ( i ~= j && boundver == 2 )
                    fprintf(fid, "%d %d 0.5\n", conscnt, j);
                else
                    fprintf(fid, "%d %d 1.0\n", conscnt, j);
                end
                conscnt = conscnt + 1;
                cnt = cnt + 1;
            end
        end
    end

    % no ACOORD for trace constraint
    conscnt = conscnt + 1;

    % sparsity constraint
    for j = 0:n-1   
        fprintf(fid, "%d %d 1.0\n", conscnt, j);
        cnt = cnt + 1;
    end
    conscnt = conscnt + 1;
    
    % no ACOORD for off-diagonal constraint
    if sumineq == 1
        conscnt = conscnt + 1;
    end
    
    % SOCP constraints
    if socp == 2
        for i = 0:n-1
            fprintf(fid, "%d %d 0.5\n", conscnt, i);
            conscnt = conscnt +1;
            for j = 0:n-1
                conscnt = conscnt + 1;
            end
            fprintf(fid, "%d %d -0.5\n", conscnt, i);
            cnt = cnt + 2;
            conscnt = conscnt + 1;
        end
    end

    fprintf(fid, "\n");
    if ( conscnt ~= ncons || cnt ~= nACOORD )
        fprintf("conscnt is = %d, should be = %d\n", conscnt, ncons);
        fprintf("cnt is = %d, should be = %d\n", cnt, nACOORD);
        error("Error: Something went wrong when writing ACOORD!\n");
    end

    %% write BCOORD
    fprintf(fid, "BCOORD\n");
    if socp == 1
        nBCOORD = n+2 + n*(n+1);
    else
        nBCOORD = n+2;
    end
    nBCOORD = nBCOORD + sumineq;
    fprintf(fid, "%d\n", nBCOORD);
    cnt = 0;

    % SOCP constraint
    if socp == 1
        for i = 0:n-1
            for j = 0:n
                pos = i*(0.5*(n+3)*(n+2));
                fprintf(fid, "%d -1.0\n",pos+j);
                cnt = cnt + 1;
            end
        end
        if cnt ~= n*(n+1)
            error("Error while writing BCOORD for SOCP constraint!\n");        
        end
    end
    if socp == 1
        conscnt = nsocpcons;
    else
        conscnt = 0;
    end
    
    % z < = 1
    for j = 0:n-1
        fprintf(fid, "%d 1.0\n", conscnt);
        cnt = cnt + 1;
        conscnt = conscnt + 1;
    end

    % nothing to add for coupling constraints -z_j - X_{ij} <= 0
    for i = 0:n-1
        if boundver == 0
            conscnt = conscnt + 1;
        else
            if strgbnds == 0
                for j = 0:n-1
                    conscnt = conscnt + 1;
                end
            else   
                for j = 0:i
                    conscnt = conscnt + 1;
                end
            end
        end
    end

        
    % nothing to add for coupling constraints z_j - X_{ij} >= 0
    for i = 0:n-1
        if boundver == 0
            conscnt = conscnt + 1;
        else
            for j = 0:n-1
                conscnt = conscnt + 1;
            end
        end
    end


    % trace constraint
    fprintf(fid, "%d -1.0\n", conscnt);
    conscnt = conscnt + 1;
    cnt = cnt + 1;

    % sparsity constraint
    fprintf(fid, "%d %d\n", conscnt, -order);
    conscnt = conscnt + 1;
    cnt = cnt + 1;
    
    % off-diagonal constraint
    if sumineq == 1
       fprintf(fid, "%d %d\n", conscnt, -order + 1);
       cnt = cnt + 1;
       conscnt = conscnt + 1;
    end

    % nothing to add for SOCP constraint
    if socp == 2
        conscnt = conscnt + n*(n+2);
    end
    
    fprintf(fid, "\n");
    if ( conscnt ~= ncons || cnt ~= nBCOORD )
        fprintf("conscnt is = %d, should be = %d\n", conscnt, ncons);
        fprintf("cnt is = %d, should be = %d\n", cnt, nBCOORD);
        error("Error: Something went wrong when writing BCOORD!\n");
    end

    %% close file
    fclose(fid);
end
