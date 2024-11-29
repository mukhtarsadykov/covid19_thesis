{
	site=$2;
	if($3=="TCG") ref++;
	if($3=="TTG") mut++;
}
END {
	if(ref+mut>0) print site "\t" ref "\t" mut "\t" mut/(mut+ref);
}
