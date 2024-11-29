{
	if(file==0) {
		if(length($3)>=7 && $4=="Human") {
			split($3,t,"-");
			if(t[1]=="2019") months[$1]=0;
			if(t[1]=="2020") months[$1]=(t[2]*1);
			if(t[1]=="2021") months[$1]=(t[2]*1)+12;
			if(t[1]=="2022") months[$1]=(t[2]*1)+24;
		}
	}
	else {
		split($1,s,"|");
		if(months[s[1]]!="") {
			if($3=="TCG") ref[months[s[1]]]++;
			else if($3=="TTG") mut[months[s[1]]]++;
			else if($3!~/N/ && $3!~/\-/) oth[months[s[1]]]++;
			else non[months[s[1]]]++;
			site=$2;
		}
	}
}
END {
	print "site\tmonths\treference\tUCG>UUG\tother\tmissing";
	for(i=0;i<=30;i++) print site "\t" i "\t" ref[i] "\t" mut[i] "\t" oth[i] "\t" non[i];
}
