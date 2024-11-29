{
	if($1~/>/) name=substr($1,2);
	else print name "\t" pos "\t" substr($0,pos-1,3);
}
