#!/bin/env perl
use strict;
use warnings;

my $usage = "Usage:\n \t$0 <dir>/<card name>\n";
die $usage unless @ARGV == 1;
my %cards = ();
my $inputCard = $ARGV[0];
my ($dir,$name) = ($inputCard =~ /^(.+?)\/([^\/]+)$/);

print $dir, "\n";
print $name, "\n";



die $usage if (!defined $dir || !defined $name); 
$name =~ s/.*?([^\/]+)$/$1/;
$name =~ s/\.[^\.]+$//;
my $nCardsPerMassPoint = 0;
foreach my $line (`cat $inputCard`){
    $line =~ s/\n//;
    next if ($line =~ /^\s*\#/);
    next if ($line !~ /\S/);
    my @elements = split(/\s+/,$line);
    die "Wrong format of input file. Each line show be: <mass> <card1> [<card2>...]\n" 
	unless @elements>1;
    my $mass = shift @elements;
    $nCardsPerMassPoint = @elements if ($nCardsPerMassPoint == 0);
    die "Inconsistent number of cards per mass point. It should be the same number to parse it right\n"
	unless $nCardsPerMassPoint == @elements;
    $cards{$mass} = \@elements;
}
die "0 cards is found\n" unless $nCardsPerMassPoint>0;

sub FilterColumns{
    my @input = @_;
    if ($input[1] =~ /lnN/){
#	@input = @input[0,5..$#input];
	@input = @input[0,4..$#input];
    } else {
#	@input = @input[0,4..$#input];
	@input = @input[0,3..$#input];
    }
    return @input;
}


my $tableBegin = << 'EOF';
\begin{table}
{%\footnotesize
 \small
 \begin{center}
 \begin{tabular}{SEPARATORS}
 \hline
 TITLE
 \hline
EOF

my $tableEnd = << 'EOF';
\hline
\end{tabular}
\end{center}
}
\caption{TEXT}
\end{table}
EOF

open(OUT,">$dir/$name-report.tex");
print OUT << 'EOF';
\documentclass{report}
\usepackage{fullpage}
\begin{document}
EOF

# print "Number of mass points to process: ", scalar keys %cards,"\n";
foreach my $card(0..$nCardsPerMassPoint-1){
    my $title;
    my @names = ("") x $nCardsPerMassPoint;
    foreach my $mass( sort {$a<=>$b} keys %cards ){
	my $cardname = $cards{$mass}->[$card];
	my $file = "$dir/$cardname";
	$cardname =~ s/(\d+\/)//;
	$names[$card] = $cardname;
	# print "$file\n";
	my @MyColumns;
	my @columns;
	my @errors2;
	my $observation;
	foreach my $line(`cat $file`){
	    if ( !defined $title && ($line =~ /process/) ){
		$title = $line;
		$title =~ s/\n//;
		my @titleColumns = split(/\s+/,$title); 
		@titleColumns = FilterColumns(@titleColumns);
		push @titleColumns, "\$\\sum\$Bkg";
		push @titleColumns, "Data";
		#my $separators = "l | c c | ".("c "x(@titleColumns - 5))." | c c";
		#$title = join(" & ",@titleColumns);

                #my title columns
                my @MyTitleColumns;
#                 #No Top
#                 push @MyTitleColumns, "Higgs Mass";
#                 push @MyTitleColumns, "Signal";
#                 push @MyTitleColumns, "WW";
#                 push @MyTitleColumns, "Drell-Yan";
#                 push @MyTitleColumns, "W+Jets";
#                 push @MyTitleColumns, "Other";
#                 push @MyTitleColumns, "\$\\sum\$Bkg";
# 		push @MyTitleColumns, "Data";
#                 #No DY
#                 push @MyTitleColumns, "Higgs Mass";
#                 push @MyTitleColumns, "Signal";
#                 push @MyTitleColumns, "WW";
#                 push @MyTitleColumns, "Top";
#                 push @MyTitleColumns, "W+Jets";
#                 push @MyTitleColumns, "Other";
#                 push @MyTitleColumns, "\$\\sum\$Bkg";
# 		push @MyTitleColumns, "Data";
                #No WJets
                push @MyTitleColumns, "Higgs Mass";
                push @MyTitleColumns, "Signal";
                push @MyTitleColumns, "WW";
                push @MyTitleColumns, "Drell-Yan";
                push @MyTitleColumns, "Top";
                push @MyTitleColumns, "Other";
                push @MyTitleColumns, "\$\\sum\$Bkg";
		push @MyTitleColumns, "Data";


                my $separators = "l | c | c c c c | c c";
		$title = join(" & ",@MyTitleColumns);

		$title .= " \\\\";
		my $tab = $tableBegin;
		$tab =~ s/SEPARATORS/$separators/m;
		$tab =~ s/TITLE/$title/m;
		print OUT $tab;
            }
	    $line =~ s/\n//;
	    if ($line =~ /rate/){
		$line =~ s/rate/$mass/;
		@columns = FilterColumns(split(/\s+/,$line));
# 		foreach my $i(1..$#columns){
# 		    $columns[$i] = sprintf("%.1f", $columns[$i]);
# 		}
	    }
	    if ($line =~ /lnN/){
		my @errors = FilterColumns(split(/\s+/,$line));
		die "Errors and yeild don't match\n" unless (scalar @errors == scalar @columns);
		if (! @errors2){
		    @errors2 = @errors;
		    foreach my $i(1..$#errors2){
			my $err = $errors2[$i];
			$err = 1 if ($err eq "-");
			$errors2[$i] = ($err-1)*($err-1);
		    }
		} else {
		    foreach my $i(1..$#errors2){
			my $err = $errors[$i];
			$err = 1 if ($err eq "-");
			$errors2[$i] += ($err-1)*($err-1);
		    }
		}
	    }   
	    if ($line =~ /Observation\s+(\d+)/){
		$observation = $1;
	    }
	}
        my $signalYield = 0;
        my $signalYieldErr2 = 0;
        my $wwYield = 0;
        my $wwYieldErr2 = 0;
        my $dyYield = 0;
        my $dyYieldErr2 = 0;
        my $wjetsYield = 0;
        my $wjetsYieldErr2 = 0;
        my $topYield = 0;
        my $topYieldErr2 = 0;
        my $otherWithTopYield = 0;
        my $otherWithTopYieldErr2 = 0;
        my $otherWithDYYield = 0;
        my $otherWithDYYieldErr2 = 0;
        my $otherWithWJetsYield = 0;
        my $otherWithWJetsYieldErr2 = 0;
	my $sum = 0;
	my $sumErr2 = 0;

	foreach my $i(1..$#columns){
            #add signals
            if ($i==1 || $i ==2) {
                $signalYield += $columns[$i];
		$signalYieldErr2 += $errors2[$i]*$columns[$i]*$columns[$i];
            }
            #add ww
            if ($i==3 || $i ==4) {
                $wwYield += $columns[$i];
		$wwYieldErr2 += $errors2[$i]*$columns[$i]*$columns[$i];
            }
            #add dy
            if ($i==7 ) {
                $dyYield += $columns[$i];
		$dyYieldErr2 += $errors2[$i]*$columns[$i]*$columns[$i];
            }
            #add wjets
            if ($i==8 ) {
                $wjetsYield += $columns[$i];
		$wjetsYieldErr2 += $errors2[$i]*$columns[$i]*$columns[$i];
            }
            #add top
            if ($i==6 ) {
                $topYield += $columns[$i];
		$topYieldErr2 += $errors2[$i]*$columns[$i]*$columns[$i];
            }
            #otherWithTop
            if ($i==5 || $i==6 ||$i==9 ||$i==10 ) {
                $otherWithTopYield += $columns[$i];
		$otherWithTopYieldErr2 += $errors2[$i]*$columns[$i]*$columns[$i];
            }
            #otherWithDY
            if ($i==5 || $i==7 ||$i==9 ||$i==10 ) {
                $otherWithDYYield += $columns[$i];
		$otherWithDYYieldErr2 += $errors2[$i]*$columns[$i]*$columns[$i];
            }
            #otherWithWJets
            if ($i==5 || $i==8 ||$i==9 ||$i==10 ) {
                $otherWithWJetsYield += $columns[$i];
		$otherWithWJetsYieldErr2 += $errors2[$i]*$columns[$i]*$columns[$i];
            }
	    if ($i>2){ # skip first two columns, which are qqH and ggH
		$sum += $columns[$i];
		$sumErr2 += $errors2[$i]*$columns[$i]*$columns[$i];
	    }
	    $columns[$i] = sprintf("\$%0.1f\\pm%0.1f\$", $columns[$i], $columns[$i]*sqrt($errors2[$i]));
	}
	push @columns, sprintf("\$%0.1f\\pm%0.1f\$", $sum, sqrt($sumErr2));
	push @columns, "$observation";
	#print OUT join(" & ",@columns)." \\\\\n";

#         #make my custom columns (No Top)
#         push @MyColumns, sprintf("%s", $columns[0]);
#         push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $signalYield, sqrt($signalYieldErr2));
#         push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $wwYield, sqrt($wwYieldErr2));
#         push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $dyYield, sqrt($dyYieldErr2));
#         push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $wjetsYield, sqrt($wjetsYieldErr2));
#         push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $otherWithTopYield, sqrt($otherWithTopYieldErr2));
#         push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $sum, sqrt($sumErr2));
#         push @MyColumns, "$observation";
#         #make my custom columns (No DY)
#         push @MyColumns, sprintf("%s", $columns[0]);
#         push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $signalYield, sqrt($signalYieldErr2));
#         push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $wwYield, sqrt($wwYieldErr2));
#         push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $topYield, sqrt($topYieldErr2));
#         push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $wjetsYield, sqrt($wjetsYieldErr2));
#         push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $otherWithDYYield, sqrt($otherWithDYYieldErr2));
#         push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $sum, sqrt($sumErr2));
#         push @MyColumns, "$observation";
        #make my custom columns (No WJets)
        push @MyColumns, sprintf("%s", $columns[0]);
        push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $signalYield, sqrt($signalYieldErr2));
        push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $wwYield, sqrt($wwYieldErr2));
        push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $dyYield, sqrt($dyYieldErr2));
        push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $topYield, sqrt($topYieldErr2));
        push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $otherWithWJetsYield, sqrt($otherWithWJetsYieldErr2));
        push @MyColumns, sprintf("\$%0.1f\\pm%0.1f\$", $sum, sqrt($sumErr2));
        push @MyColumns, "$observation";



        print OUT join(" & ",@MyColumns)." \\\\\n";

    }
    my $tab = $tableEnd;
    my $label = $names[$card];
    $label =~ s/_/\\_/g;
    $tab =~ s/TEXT/Summary of card $label/m;
    print OUT $tab;
}
print OUT '\end{document}'."\n";
close OUT;
exit
