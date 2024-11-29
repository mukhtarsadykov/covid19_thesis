BEGIN {
	add[1]=31;
	add[2]=28;
	add[3]=31;
	add[4]=30;
	add[5]=31;
	add[6]=30;
	add[7]=31;
	add[8]=31;
	add[9]=30;
	add[10]=31;
	add[11]=30;
}
{
	if(file==0) {
		if(length($3)==10 && $4=="Human") {
			date[$1]=$3;
			split($3,t,"-");
			year=t[1];
			month=t[2]*1;
			day=t[3]*1;
			if(year=="2019") time[$1]=t[3]-24;
			if(year=="2020") {
				time[$1]=7;
				for(i=1;i<=month-1;i++) time[$1]+=add[i];
				time[$1]+=day;
			}
			if(year=="2021") {
				time[$1]=7+365;
				for(i=1;i<=month-1;i++) time[$1]+=add[i];
				time[$1]+=day;
			}
			if(year=="2022") {
				time[$1]=7+365+365;
				for(i=1;i<=month-1;i++) time[$1]+=add[i];
				time[$1]+=day;
			}
			#print $1 "\t" $3 "\t" year "\t" month "\t" day "\t" time[$1];
		}
	}
	else {
		split($1,s,"|");
		if(time[s[1]]!="") {
			print s[1] "\t" $2 "\t" ref "\t" $3 "\t" date[s[1]] "\t" time[s[1]];
		}
	}
}