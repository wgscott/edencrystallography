#!/usr/bin/perl -w

@inf = glob("*inf");
@hlp = glob("*hlp");

$ii = 0;
$hh = 0;

print "nhlp = $#hlp, ninf = $#inf\n";

sub pr_head {
    my $FP = shift;
    my $ff = shift;
    
    printf ($FP "<title>EDEN manual for $ff</title>\n");
}

sub pr_foot {
    my $FP = shift;
    my $ff = shift;
    
    printf ($FP "<a href=\"index.html\">Index for EDEN manual</a>\n<hr>\n");
}

sub pr_inf {
    my $FP = shift;
    my $fn = shift;
    
    printf ($FP "<h2>Synopsis for $fn</h2>\n");
    open (FN,"$fn.inf") || die("$fn");
    while ($line = <FN>) {
	chomp($line);
	if (index($line,"******") < 0) {
	    print ($FP "$line<br>\n");
	}
    }
    close FN;
    printf ($FP "<hr>\n");
}

sub pr_hlp {
    my $FP = shift;
    my $fn = shift;
    
    printf ($FP "<h2>Background for $fn</h2>\n");
    open (FN,"$fn.hlp") || die("$fn");
    while ($line = <FN>) {
	chomp($line);
	if (index($line,"------") < 0) {
	    print ($FP "$line<br>");
	}
    }
    close FN;
    printf ($FP "<hr>\n");
}

sub name {
    my $str = shift;
    my @ttt = split('\.',$str);
    print "Processing $str, nttt = $#ttt, @ttt\n";    
    return $ttt[0];
}

open(NDX,">index.html") || die("index.html");
printf (NDX "<title>EDEN online manual</title>\n");
printf (NDX "<h2>EDEN online manual</h2>\n");
printf (NDX "<ol>\n");
while (($ii <= $#inf) || ($hh <= $#hlp)) {
    $i = name($inf[$ii]);
    $h = name($hlp[$hh]);
    print "i = $i, h = $h\n";
    if ($i eq $h) {
	print (NDX "<li><a href=\"$i.html\">$i</a></li>\n");
	open(FP,">$i.html");
	pr_head(FP,$i);
	pr_inf(FP,$i);
	pr_hlp(FP,$h);
	pr_foot(FP,$i);
	close FP;
	$ii++;
	$hh++;
    }
    elsif ($i lt $h) {
	open(FP,">$i.html");
	pr_head(FP,$i);
	pr_inf(FP,$i);
	pr_foot(FP,$i);
	close FP;
	$ii++;
    }
    else {
	open(FP,">$h.html");
	pr_head(FP,$h);
	pr_hlp(FP,$h);
	pr_foot(FP,$h);
	close FP;
	$hh++;
    }
}
printf (NDX "</ol>\n");
printf (NDX "<hr><a href=\"http://www.edencrystallography.org\">EDEN homepage</a>\n<hr>\n"); 
close NDX;
