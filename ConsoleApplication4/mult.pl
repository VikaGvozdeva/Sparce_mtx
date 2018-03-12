use File::Basename;
$path_to_mat = "../matrix";
$prog_dir = "/common/home/gvozdeva_v/test/bin";

@mat = <$path_to_mat/*.mtx>;
@format = ("COO");

for($j = 0; $j < $#format + 1; $j++)
{
	mkdir $format[$j];
}

for($j = 0; $j < $#format + 1; $j++)
{
	for($i = 0; $i < $#mat + 1; $i++)
	{
	    $fr = $format[$j];
		$file_name = basename($mat[$i]);
		($file_name) = split(/\./, $file_name);
		#print $fr." ".$file_name."\n";
		$dir = $fr."/".$file_name;
		$mt = $prog_dir."/".$mat[$i];
		mkdir $dir;
		#print "(cd $dir; sh $prog_dir/run_convert.sh $prog_dir $mt $fr)\n";
		system "(cd $dir; sbatch $prog_dir/run_mult.sh $prog_dir $mt $fr 5)";
	}
}

