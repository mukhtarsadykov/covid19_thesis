{
	if(file==0) done[$2]=1;
	else {
		if(done[$NF]==1) print "#" $0;
		else print $0;
	}
}
