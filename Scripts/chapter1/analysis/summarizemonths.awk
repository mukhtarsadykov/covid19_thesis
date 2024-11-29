BEGIN {
	codon[1]="TCT"
	codon[2]="TCC"
	codon[3]="TCA"
	codon[4]="TCG"
	codon[5]="CCT"
	codon[6]="CCC"
	codon[7]="CCA"
	codon[8]="CCG"
	codon[9]="ACT"
	codon[10]="ACC"
	codon[11]="ACA"
	codon[12]="ACG"
	codon[13]="GCT"
	codon[14]="GCC"
	codon[15]="GCA"
	codon[16]="GCG"
	FS="\t";
}
{
	if($4+$5+$6+$7>0 && $1!="site") {
		count[$2, $3]++;
		accu[$2, $3]+=$6/($4+$5+$6+$7);
	}
}
END {
	for(i=1;i<=16;i++) printf("\t%s", codon[i]);
	printf("\n");
	for(j=0;j<=30;j++) {
		printf("%s", j);
		for(i=1;i<=16;i++) {
			if(count[codon[i], j]>0) printf("\t%s", accu[codon[i], j]/count[codon[i], j]);
			else printf("\t");
		}
		printf("\n");
	}
}