{
	if(file==0) {
		months="";
		if(length($3)>=7 && $4=="Human") {
			split($3,t,"-");
			if(t[1]=="2019") months=0;
			if(t[1]=="2020") months=(t[2]*1);
			if(t[1]=="2021") months=(t[2]*1)+12;
			if(t[1]=="2022") months=(t[2]*1)+24;
		}
		age[$1]=months;
	}
	else {
		split($1,s,"|");
		if(age[s[1]]!="" && $3==substr(ref,1,1) "T" substr(ref,3)) mut[age[s[1]]]++;
		if(age[s[1]]!="" && $3!~/BDEFHIJKLMNOPQRSUVXYZ-/) all[age[s[1]]]++;
		pos=$2;
	}
}
END {
	printf("%s\t%s", pos, ref);
	for(i=0;i<=30;i++) printf("\t%s;%s", mut[i]*1, all[i]*1);
	printf("\n");
}