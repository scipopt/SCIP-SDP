if [ ! -f MISDP_backup/RIPSDPA.solu ]; then
    echo "ERROR: File /MISDP_backup/RIPSDPA.solu is missing!";
    return
fi


# For large instances:
if [ ! -d MISDP ]; then
    echo "ERROR: No instances available in folder MISDP/ !";
    return
fi

cd MISDP/

fileSDPA=RIPSDPA.test
for x in `ls *.dat-s`;
do
    printf "RIPtest/%s\n" $x >> $fileSDPA;
done;

cnt_type=0;
cnt_instance=0;
testfileAll=RIPCBFall.test;
solufileAll=RIPCBFall.solu;

# options:
arr_rank1=(0 1);
arr_trineq=(0 1);
arr_boundver=(0 1 2);
arr_pd=("p" "d");
arr_socp=(0 1 2);
arr_strgbnds=(0 1);

# rank1-constraint
for r in "${arr_rank1[@]}";
do
    # trace-inequality
    for t in "${arr_trineq[@]}";
    do
	# boundver
	for b in "${arr_boundver[@]}";
	do
	    # primal/dual
	    for pd in "${arr_pd[@]}";
	    do
		# socp-inequality
		for so in "${arr_socp[@]}";
		do
		    # strong bounds
		    for st in "${arr_strgbnds[@]}";
		    do
			cmd="ls *\$pd\$r\$so\$st\$t\$b.cbf"
			testfile=RIPCBF$pd$r$so$st$t$b.test
			solufile=RIPCBF$pd$r$so$st$t$b.solu
			type=$pd$r$so$st$t$b;
			for x in `eval $cmd`
			do
			    oldname="${x/.cbf/}"
			    newname="${x/$type.cbf/}"
			    tmp=$(grep $newname ../MISDP_backup/RIPSDPA.solu)
			    if [[ $tmp =~ "MISDPr" ]]; then
				arr=($tmp)
				result=$(bc -l <<<"${arr[2]}*(-1)")
				printf "%s %s %s\n" ${arr[0]} ${arr[1]/$newname/$oldname} $result >> $solufile;
				printf "%s %s %s\n" ${arr[0]} ${arr[1]/$newname/$oldname} $result >> $solufileAll;
			    else
				echo "${tmp/$newname/$oldname}" >> $solufile;
				echo "${tmp/$newname/$oldname}" >> $solufileAll;
			    fi
			    printf "RIPtest/%s\n" $x >> $testfile;
			    printf "RIPtest/%s\n" $x >> $testfileAll;
			    ((cnt_instance+=1));
			done
			((cnt_type+=1))
		    done
		done
	    done
	done
    done
done

printf "LARGE INSTANCES: Created overall testset with %d instances and %d testsets for the different RIP formulations.\n" $cnt_instance $cnt_type;

cd ..


# For small testing instances:
if [ ! -d MISDP_small ]; then
    echo "No small instances available!";
    return
fi

if [ ! -f MISDP_backup/RIPSDPASMALL.solu ]; then
    echo "ERROR: File /MISDP_backup/RIPSDPASMALL.solu is missing!";
    return
fi


cd MISDP_small/

fileSDPA=RIPSDPASMALL.test
for x in `ls *.dat-s`;
do
    printf "RIPtest/%s\n" $x >> $fileSDPA;
done;

cnt_type=0;
cnt_instance=0;
testfileAll=RIPCBFSMALLall.test;
solufileAll=RIPCBFSMALLall.solu;

# options:
arr_rank1=(0 1);
arr_trineq=(0 1);
arr_boundver=(0 1 2);
arr_pd=("p" "d");
arr_socp=(0 1 2);
arr_strgbnds=(0 1);

# rank1-constraint
for r in "${arr_rank1[@]}";
do
    # trace-inequality
    for t in "${arr_trineq[@]}";
    do
	# boundver
	for b in "${arr_boundver[@]}";
	do
	    # primal/dual
	    for pd in "${arr_pd[@]}";
	    do
		# socp-inequality
		for so in "${arr_socp[@]}";
		do
		    # strong bounds
		    for st in "${arr_strgbnds[@]}";
		    do
			cmd="ls *\$pd\$r\$so\$st\$t\$b.cbf"
			testfile=RIPCBFSMALL$pd$r$so$st$t$b.test
			solufile=RIPCBFSMALL$pd$r$so$st$t$b.solu
			type=$pd$r$so$st$t$b;
			for x in `eval $cmd`
			do
			    oldname="${x/.cbf/}"
			    newname="${x/$type.cbf/}"
			    tmp=$(grep $newname ../MISDP_backup/RIPSDPASMALL.solu)
			    if [[ $tmp =~ "MISDPr" ]]; then
				arr=($tmp)
				result=$(bc -l <<<"${arr[2]}*(-1)")
				printf "%s %s %s\n" ${arr[0]} ${arr[1]/$newname/$oldname} $result >> $solufile;
				printf "%s %s %s\n" ${arr[0]} ${arr[1]/$newname/$oldname} $result >> $solufileAll;
			    else
				echo "${tmp/$newname/$oldname}" >> $solufile;
				echo "${tmp/$newname/$oldname}" >> $solufileAll;
			    fi
			    printf "RIPtest/%s\n" $x >> $testfile;
			    printf "RIPtest/%s\n" $x >> $testfileAll;
			    ((cnt_instance+=1));
			done
			((cnt_type+=1))
		    done
		done
	    done
	done
    done
done

printf "SMALL INSTANCES: Created overall testset with %d instances and %d testsets for the different RIP formulations.\n" $cnt_instance $cnt_type;

cd ..

