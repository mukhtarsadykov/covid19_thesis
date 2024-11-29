{
	months="";
	if(length($3)>=7) {
		split($3,t,"-");
		if(t[1]=="2019") months=0;
		if(t[1]=="2020") months=(t[2]*1);
		if(t[1]=="2021") months=(t[2]*1)+12;
		if(t[1]=="2022") months=(t[2]*1)+24;
	}
	print $3 "\t" months;
}
