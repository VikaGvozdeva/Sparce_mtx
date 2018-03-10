$prog = shift(@ARGV);
$mat  = shift(@ARGV);
$fr   = shift(@ARGV);
$cnt_repit= shift(@ARGV);

$run = $prog." ".$mat." r ".$fr." 1";
print $run."\n";

$time = 10000000.0;
for($i = 0; $i < $cnt_repit; $i++)
{
	$res = `$run`;
	print $res;
	if($res =~ /time: 		([-+]?[0-9]*\.?[0-9]+)/)
	{
		if($time > $1)
		{
			$time = $1;
		}
	}
}
#$time /= $cnt_repit;
open(W, ">avg_time.txt");

print W $time;

close(W);