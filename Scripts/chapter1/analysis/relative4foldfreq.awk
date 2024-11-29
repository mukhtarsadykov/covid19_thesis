BEGIN {
	FS="\t";
}
{
	if($1!="site" && $4+$5+$6+$7>0) {
		if($3==0) {
			a++;
			site[a]=$1;
		}
		f[$1, $3]=$6/($4+$5+$6+$7);
		if(f[$1, $3]>max[$1]) max[$1]=f[$1, $3];
	}
}
END {
	for(i=1;i<=a;i++) {
		if(max[site[i]]>=minfreq) printf("\t%s", site[i]);
	}
	printf("\n");
	for(j=0;j<=30;j++) {
		printf("%s", j);
		for(i=1;i<=a;i++) {
			if(max[site[i]]>=minfreq) printf("\t%s", f[site[i], j]/max[site[i]]);
		}
		printf("\n");
	}
}