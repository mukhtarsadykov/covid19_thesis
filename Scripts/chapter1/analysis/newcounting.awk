BEGIN {
	triplet[1]="TCT";
	triplet[2]="TCC";
	triplet[3]="TCA";
	triplet[4]="TCG";
	triplet[5]="CCT";
	triplet[6]="CCC";
	triplet[7]="CCA";
	triplet[8]="CCG";
	triplet[9]="ACT";
	triplet[10]="ACC";
	triplet[11]="ACA";
	triplet[12]="ACG";
	triplet[13]="GCT";
	triplet[14]="GCC";
	triplet[15]="GCA";
	triplet[16]="GCG";
}
{
	if(file==0) {
		get[$2]=1;
		site[$2]=$3;
	}
	else {
		if($1~/>/) {
			split($1,t,"|");
			printf("%s", substr(t[1],2));
		}
		else {
			for(k in all) delete all[k];
			for(k in mut) delete mut[k];
			for(i=1;i<=length($0);i++) {
				if(get[i]==1) {
					if(substr($0,i-1,3)!~/[BDEFHIJKLMNOPQRTUVXYZ-]/) {
						all[site[i]]++;
						if(substr($0,i-1,3)==substr(site[i],1,1) "T" substr(site[i],3)) mut[site[i]]++;
					}
				}
			}
			for(j=1;j<=16;j++) printf("\t%s;%s", mut[triplet[j]]*1, all[triplet[j]]*1);
			printf("\n");
		}
	}
}