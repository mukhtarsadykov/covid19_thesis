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
			if(print_all==1) printf("%s\t%s\t%s", s[1], months[s[1]], $3);
			if($3==ref) { # original/reference
				ori[months[s[1]]]++; 
				if(print_all==1) printf("\tORIGINAL\n");
			}
			else if($3~/\-/ || $3~/[BDEFHIJKLMNOPQRSUVWXYZ]/) { # other, gaps or Ns
				oth[months[s[1]]]++;
				if(print_all==1) printf("\tOTHER\n");
			}
			else if(substr(ref,1,1)!=substr($3,1,1) || substr(ref,3)!=substr($3,3)) { # neighbor change(s)
				nei[months[s[1]]]++;
				if(print_all==1) printf("\tNEIGHBOR\n");
			}
			else {
				mut["C>" substr($3,2,1), months[s[1]]]++;
				if(print_all==1) printf("\tC>%s\n", substr($3,2,1));
			}
			site=$2;
		}
	}
}
END {
	print "site\tref\tmonths\treference\tC>A\tC>U\tC>G\tneighborchange\tother"
	for(i=0;i<=30;i++) print site "\t" ref "\t" i "\t" ori[i] "\t" mut["C>A", i] "\t" mut["C>T", i] "\t" mut["C>G", i] "\t" nei[i] "\t" oth[i];;
}