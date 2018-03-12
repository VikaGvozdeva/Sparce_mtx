use File::Basename;
$path_to_mat = "../matrix";
$prog_dir = "/common/home/gvozdeva_v/sequential/bin";
$cnt_repit = 10;

@mat = <$path_to_mat/*.mtx>;
#@format = ("COO","CRS","CCS","JD");
@format = ("MKL"); 

$task = "r";

if($task eq "wr")
{
	for($j = 0; $j < $#format + 1; $j++)
	{
		for($i = 0; $i < $#mat + 1; $i++)
		{
			$fr = $format[$j];
			$file_name = basename($mat[$i]);
			($file_name) = split(/\./, $file_name);
			#print $fr." ".$file_name."\n";
			$dir = $file_name;
			$mt = $prog_dir."/".$mat[$i];
			mkdir $dir;
#			print "(cd $dir; sh $prog_dir/run_convert.sh $prog_dir $mt $fr)\n";
			system "(cd $dir; sbatch $prog_dir/run_convert.sh $prog_dir $mt $fr $task)";
		}
	}
}

if($task eq "r")
{
	for($j = 0; $j < $#format + 1; $j++)
	{
		for($i = 0; $i < $#mat + 1; $i++)
		{
			$fr = $format[$j];
			$file_name = basename($mat[$i]);
			($file_name) = split(/\./, $file_name);
			$dir = $file_name;
			$mt = $prog_dir."/".$mat[$i];
			#print $dir."\n";
			system "(cd $dir; sbatch $prog_dir/run_mult.sh $prog_dir $mt $fr $cnt_repit)";
#			system "(cd $dir; sh $prog_dir/run_mult.sh $prog_dir $mt $fr $cnt_repit)";
#			print "(cd $dir; sbatch $prog_dir/run_mult.sh $prog_dir $mt $fr $cnt_repit)\n";
		}
	}
}

if($task eq "mr")
{
	open(W, ">all_avg_time.txt");
	for($j = 0; $j < $#format + 1; $j++)
	{
		for($i = 0; $i < $#mat + 1; $i++)
		{
			$file_name = basename($mat[$i]);
			($file_name) = split(/\./, $file_name);
			$dir = $fr."/".$file_name;
			open(R, "./$dir/avg_time.txt");
			$time = <R>;
			print W "$file_name $time \n";
			close(R);
		}
	}
	close(W);
}