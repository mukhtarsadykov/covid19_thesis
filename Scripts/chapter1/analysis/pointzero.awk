{
	if($4!=$3 && substr($4,1,1)==substr($3,1,1) && substr($4,3)==substr($3,3) && substr($4,2,1)~/[ACGT]/) {
		if($6<min[substr($4,2,1)] || min[substr($4,2,1)]=="") min[substr($4,2,1)]=$6;
	}
	pos=$2;
	ref=$3;
}
END {
	printf("%s\t%s", pos, ref);
	if(min["T"]!="") printf("\t%s", min["T"]);
	else printf("\tNA");
	if(min["A"]!="") printf("\t%s", min["A"]);
	else printf("\tNA");
	if(min["G"]!="") printf("\t%s", min["G"]);
	else printf("\tNA");
	printf("\n");
}
