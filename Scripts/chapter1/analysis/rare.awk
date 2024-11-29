{
	if($1!="site" && $4+$5+$6+$7>0) {
		if($6/($4+$5+$6+$7)>0.001) nix++;
		total+=$4+$5+$6+$7;
		CtoT+=$6;
	}
	prefix=$1 "\t" $2;
}
END {
	if(nix==0) print prefix "\t" total "\t" CtoT/total;
}
