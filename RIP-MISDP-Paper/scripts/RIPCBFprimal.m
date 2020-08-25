function [] = RIPCBFprimal(A, k, side, file, Rank, socp, bounds)
% schreibt SDP-File für ganzzahlige RIP-SDP-Relaxierung in primaler Form
% (mit PSD-Variablen) für Matrix A, Ordnung k, schreibt in 'file' 
% side ='l' für linke Seite/alpha_k, side='r' für rechte Seite/beta_k
% Rank = 1 für zusätzliche Rang-NB, sonst Rank = 0
% socp = 1 für zusätzliche gültige SOCP-Ungleichung von Li/Xie, sonst = 0
% bounds = 1 falls für nichtnegative Matrizen A die Schranke 0 <= X_{ij}
% statt -z_j <= X_{ij} für alle i,j benutzt werden soll, sonst = 0
% ACHTUNG: schreibt untere Dreiecksmatrizen!

    fid = fopen(file, 'w');
    m=length(A(:,1));
    n=length(A(1,:));

    % check if A is entrywise nonnegative:
    if bounds == 1 && ~all(A >= 0, 'all')
        bounds = 0;
        fprintf("Setting bounds = 0, since matrix A is not nonnegative!\n");
    elseif bounds == 1 && ~side == 'r'
        bounds = 0;
        fprintf("Setting bounds = 0, since  only works for right side of RIP!\n");
    end

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

    % add scalar variables z
    fprintf(fid, "VAR\n");
    fprintf(fid, "%d 1\n", n);
    fprintf(fid, "L+ %d\n", n);
    fprintf(fid, "\n");

    % add Matrix variables X
    fprintf(fid, "PSDVAR\n");
    if socp == 0
        fprintf(fid, "1\n");
        fprintf(fid, "%d\n", n);
    elseif socp == 1
        fprintf(fid, "2\n");
        fprintf(fid, "%d\n", n);
        fprintf(fid, "%d\n", n+2);
    else
        error("Error: Option <%s> for parameter socp not valid!\n", socp);
    end
    fprintf(fid, "\n");

    % objective
    fprintf(fid, "OBJFCOORD\n");
    fprintf(fid, "%d\n", n * (n + 1)/2);
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

    % all scalar variables z are binary
    fprintf(fid, "INT\n");
    fprintf(fid, "%d\n", n);
    for j = 0:n-1
        fprintf(fid, "%d\n", j);
    end
    fprintf(fid, "\n");

    % add rank1-constraint
    if ( Rank == 1 )
        fprintf(fid, "PSDVARRANK1\n");
        fprintf(fid, "1\n");
        fprintf(fid, "0\n");
        fprintf(fid, "\n");
    end

    % ------------------------------------------
    % constraints
    fprintf(fid, "CON\n");
    ncons = 2*n^2+n+2;
    if socp == 0
        fprintf(fid, "%d 5\n", ncons);
    else
        % SOCP constraint is the first constraint
        nsocpcons = n*(2*n+3+0.5*n*(n+1));
        fprintf(fid, "%d 6\n", ncons + nsocpcons);
        fprintf(fid, "L= %d\n", nsocpcons);
    end
    fprintf(fid, "L+ %d\n", n);         % -z_j + 1 >= 0
    fprintf(fid, "L- %d\n", n * n);     % -z_j - X_{ij} <= 0
    fprintf(fid, "L+ %d\n", n * n);     % z_j - X_{ij} >= 0
    fprintf(fid, "L= 1\n");             % \sum_j X_{jj} == 1
    fprintf(fid, "L- 1\n");             % \sum_j z_j <= k
    if socp == 1                        % \sum_j X_{ij}^2 <= X_{ii}z_i
    end
    fprintf(fid, "\n");


    % ------------------------------------------
    % precompute stuff for valid socp constraint
    FCOORDcnt = 0;
    ACOORDcnt = 0;
    BCOORDcnt = 0;
    socpcnt = 0;
    if socp == 1
        % FCOORD
        socpFCOORD = "";
        for i = 0:n-1
            for j = 0:n
                socpFCOORD = socpFCOORD + sprintf("%d 1 %d %d 1.0\n",socpcnt,j+1,j+1);
                socpcnt = socpcnt + 1;
                FCOORDcnt = FCOORDcnt + 1;
            end
            socpFCOORD = socpFCOORD + sprintf("%d 1 0 0 1.0\n",socpcnt);
            socpFCOORD = socpFCOORD + sprintf("%d 0 %d %d -0.5\n",socpcnt,i,i);
            socpcnt = socpcnt + 1;
            FCOORDcnt = FCOORDcnt + 2;
            for j = 0:n-1
                socpFCOORD = socpFCOORD + sprintf("%d 1 %d 0 1.0\n",socpcnt,j+1);
                if j <= i
                    socpFCOORD = socpFCOORD + sprintf("%d 0 %d %d -1.0\n",socpcnt,i,j);
                else
                    socpFCOORD = socpFCOORD + sprintf("%d 0 %d %d -1.0\n",socpcnt,j,i);
                end
                socpcnt = socpcnt + 1;
                FCOORDcnt = FCOORDcnt + 2;
            end
            socpFCOORD = socpFCOORD + sprintf("%d 1 %d 0 0.5\n",socpcnt,n+1);
            socpFCOORD = socpFCOORD + sprintf("%d 0 %d %d -0.5\n",socpcnt,i,i);
            socpcnt = socpcnt + 1;
            FCOORDcnt = FCOORDcnt + 2;
            for k = 0:n
                for j = 0:k-1
                    socpFCOORD = socpFCOORD + sprintf("%d 1 %d %d 1.0\n",socpcnt,k+1,j+1);
                    socpcnt = socpcnt + 1;
                    FCOORDcnt = FCOORDcnt + 1;
                end
            end
        end
        if ( socpcnt ~= nsocpcons || FCOORDcnt ~= n*(3*n+5+0.5*n*(n+1)) )
            fprintf("socpcnt is = %d, should be = %d\n",socpcnt, nsocpcons);
            fprintf("FCOORDcnt is = %d, should be = %d\n",FCOORDcnt,n*(3*n+5+0.5*n*(n+1)));
            error("Error!\n");
        else
            socpFCOORD = socpFCOORD + "\n";
        end
        
        %ACOORD
        socpACOORD = "";
        for i = 0:n-1
            pos = i*(2*n+3+0.5*n*(n-1));
            socpACOORD = socpACOORD + sprintf("%d %d -0.5\n",pos+n+1,i);
            ACOORDcnt = ACOORDcnt + 1;
            socpACOORD = socpACOORD + sprintf("%d %d -0.5\n",pos+2*n+2,i);
            ACOORDcnt = ACOORDcnt + 1;
        end
        if ACOORDcnt ~= 2*n
            error("Error!\n");
        else
            socpACOORD = socpACOORD + "\n";
        end
        
        % BCOORD
        socpBCOORD = "";
        for i = 0:n-1
            for j = 0:n-1
                pos = i*(2*n+3+0.5*n*(n-1));
                socpBCOORD = socpBCOORD + sprintf("%d -1.0\n",pos+j);
                BCOORDcnt = BCOORDcnt + 1;
            end
        end
        if BCOORDcnt ~= n^2
            error("Error!\n");
        else
            socpBCOORD = socpBCOORD + "\n";
        end        
    end
    
    
    % write ACOORD

    fprintf(fid, "ACOORD\n");
    if bounds == 0
        fprintf(fid, "%d\n", 2*n^2+2*n + ACOORDcnt);
    elseif bounds == 1
        fprintf(fid, "%d\n", n^2+2*n + ACOORDcnt);
    else
        error("Error: Option <%s> for parameter bounds not valid!\n", bounds);
    end
    cnt = 0;
    conscnt = socpcnt;       % SOCP constraint comes first

    % SOCP constraint
    if socp == 1
        fprintf(fid, socpACOORD);
    end

    % add upper bounds on z
    for j = 0:n-1
        fprintf(fid, "%d %d -1.0\n", conscnt, j);
        conscnt = conscnt + 1;
        cnt = cnt + 1;
    end

    % add coupling constraints -z_j - X_{ij} <= 0 (if bounds = 1, this
    % changes to -X_{ij} <= 0, so that nothing needs to be specified here
    if bounds == 0
        for j = 0:n-1
            for i = 0:n-1
                if ( i ~= j )
                    fprintf(fid, "%d %d -0.5\n", conscnt, j);
                else
                    fprintf(fid, "%d %d -1.0\n", conscnt, j);
                end
                conscnt = conscnt + 1;
                cnt = cnt + 1;
            end
        end
    else
        conscnt = conscnt + n^2;
    end
    
    % add coupling constraints z_j - X_{ij} >= 0
    for j = 0:n-1
        for i = 0:n-1
            if ( i ~= j )
                fprintf(fid, "%d %d 0.5\n", conscnt, j);
            else
                fprintf(fid, "%d %d 1.0\n", conscnt, j);
            end
            conscnt = conscnt + 1;
            cnt = cnt + 1;
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

    fprintf(fid, "\n");
    if ( conscnt ~= ncons + socpcnt || cnt ~= 2*n^2+2*n - bounds*n^2)
        error("Error: Something went wrong when writing ACOORD!\n");
    end


    % write FCOORD
    fprintf(fid, "FCOORD\n");
    fprintf(fid, "%d\n", 2*n^2+n + FCOORDcnt);
    cnt = 0;
    conscnt = socpcnt;      % SOCP constraint comes first

    % SOCP constraint
    if socp == 1
        fprintf(fid, socpFCOORD);
    end

    % no FCOORD for -z_j + 1 >= 0
    conscnt = conscnt + n;

    % add coupling constraints -z_j - X_{ij} <= 0
    for j = 0:n-1
        for i = 0:n-1
            fprintf(fid, "%d 0 %d %d -1.0\n", conscnt, i, j);
            cnt = cnt + 1;
            conscnt = conscnt + 1;
        end
    end

    % add coupling constraints z_j - X_{ij} >= 0
    for j = 0:n-1
        for i = 0:n-1
            fprintf(fid, "%d 0 %d %d -1.0\n", conscnt, i, j);
            cnt = cnt + 1;
            conscnt = conscnt + 1;
        end
    end

    % trace constraint
    for i = 0:n-1
        fprintf(fid, "%d 0 %d %d 1.0\n", cnt, i, i);
        cnt = cnt + 1;
    end
    conscnt = conscnt + 1;

    % no FCOORD for sparsity constraint
    conscnt = conscnt + 1;

    fprintf(fid, "\n");

    if ( conscnt ~= ncons + socpcnt || cnt ~= 2*n^2+n )
        error("Error: Something went wrong when writing FCOORD!\n");
    end

    % write BCOORD
    fprintf(fid, "BCOORD\n");
    fprintf(fid, "%d\n", n+2 + BCOORDcnt);
    conscnt = socpcnt;      % SOCP constraint comes first
    cnt = 0;

    % SOCP constraint
    if socp == 1
        fprintf(fid, socpBCOORD);
    end

    % z < = 1
    for j = 0:n-1
        fprintf(fid, "%d 1.0\n", conscnt);
        cnt = cnt + 1;
        conscnt = conscnt + 1;
    end

    % nothing to add for coupling constraints -z_j - X_{ij} <= 0
    for j = 0:n-1
        for i = 0:n-1
            conscnt = conscnt + 1;
        end
    end

        
    % nothing to add for coupling constraints z_j - X_{ij} >= 0
    for j = 0:n-1
        for i = 0:n-1
            conscnt = conscnt + 1;
        end
    end

    % trace constraint
    fprintf(fid, "%d -1.0\n", conscnt);
    conscnt = conscnt + 1;
    cnt = cnt + 1;

    % sparsity constraint
    fprintf(fid, "%d %d\n", conscnt, -k);
    conscnt = conscnt + 1;
    cnt = cnt + 1;

    fprintf(fid, "\n");
    if ( conscnt ~= ncons + socpcnt || cnt ~= n+2 )
        error("Error: Something went wrong when writing BCOORD!\n");
    end

    % close file
    fclose(fid);
end
