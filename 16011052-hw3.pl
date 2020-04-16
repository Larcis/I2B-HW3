use strict;
use warnings; 
use Scalar::Util qw(looks_like_number);
use Win32::Console;
my $CONSOLE = Win32::Console->new(STD_OUTPUT_HANDLE);
my $attr = $CONSOLE->Attr(); # Get current console colors

open(INPUT, "<input.txt") or die "Could not open input.txt!";
my %total_background_nuc_counts = ('A', 0, 'C', 0, 'G', 0, 'T', 0);
my $total_lines = 0;
my $length = 0;
my $is_line_valid = 0;
my $window_size = 0;
#my @motif_start_indexes = (11, 3, 14, 6, 11, 1);
my @motif_start_indexes = ();
my @motifs = ();
my @dnas = ();
while((my $line  = <INPUT>)){
    chomp $line;
    if(length($line) < 2){
        print "\nThis line does not contain enough nucleotids, discarding it!\n";
        $is_line_valid = 0;
    } elsif($total_lines == 0){
        $length = length($line);
        print "Please enter a valid window size between [1, $length] : ";
        do{
            $window_size = <STDIN>;
        }while($window_size < 1 || $window_size > $length);
        $is_line_valid = 1;
    } elsif(length($line) != $length){
        print "\nThis line does not have same length with the first one, discarding it!\n";
        $is_line_valid = 0;
    } else {
        $is_line_valid = 1;
    }
    if($is_line_valid){
        my $motif_start_idx = int(rand($length - $window_size + 1)); #$motif_start_indexes[$total_lines]; #
        #print "\nmotif starts $motif_start_idx \n";
        push @dnas, $line;
        push @motif_start_indexes, $motif_start_idx;
        push @motifs, ""; 
        my @line_words = split //, $line;
        for(my $i = 0; $i < length($line); $i++){
            #print $line_words[$i] , " $i \n";
            if($i < $motif_start_idx || $i >= ($motif_start_idx + $window_size)){
                $total_background_nuc_counts{$line_words[$i]}++;
            }else{
                $motifs[$total_lines].= $line_words[$i];
                #print $line_words[$i] 
            }
        }
        $total_lines++; 
    }
}

for(my $i = 0; $i < $total_lines; $i++){
    my @line_words = split //, $dnas[$i];
    my $motif_start_idx = $motif_start_indexes[$i];
    for(my $j = 0; $j < length($dnas[$i]); $j++){
        if($j < $motif_start_idx || $j >= ($motif_start_idx + $window_size)){
            print $line_words[$j];
        }else{
            $| = 1;
            $CONSOLE->Attr($FG_RED | $BG_BLACK);
            print $line_words[$j];
            $CONSOLE->Attr($attr); 
        }
    }
    print "\n";
}

my %freqs_in_motifs = (
    'A' => [],
    'C' => [],
    'G' => [],
    'T' => []
);
for(my $i = 0; $i < $window_size; $i++){
    $freqs_in_motifs{'A'}[$i]= 0;
    $freqs_in_motifs{'C'}[$i]= 0;
    $freqs_in_motifs{'G'}[$i]= 0;
    $freqs_in_motifs{'T'}[$i]= 0;
}
for(my $i = 0; $i < $total_lines; $i++){
    my @motif_words = split //, $motifs[$i];
    for(my $j = 0; $j < length($motifs[$i]); $j++){
        $freqs_in_motifs{$motif_words[$j]}[$j]++;
    }
}

my @table_1;
my @table_2;
my @table_3;
my $val = 0;
my @nucs = keys %total_background_nuc_counts;
print "\n";
for(my $i = 0; $i < 5 ; $i++){ 
    for(my $j = 0; $j < $window_size + 2; $j++){
        if($i == 0 && $j == 0){
            $val = "N"
        }elsif($j == 0){
            $val = $nucs[$i-1];
        }elsif($i == 0){
            $val = $j - 1;
        } elsif($j == 1){
            $val = $total_background_nuc_counts{$nucs[$i-1]};
        } else {
            $val = $freqs_in_motifs{$nucs[$i-1]}[$j-2];
        }
        $table_1[$i][$j] = $val;
        $table_2[$i][$j] = $val;
        $table_3[$i][$j] = $val;
    } 
}

sub print_table {
    my ($name, $window_size, @table ) = @_;
    print "\nTABLE $name: \n";
    for(my $k= 1; $k < (2+$window_size) * 10+1; $k++){
        print "_";
    }
    print "\n";
    for(my $i = 0; $i < 5 ; $i++){ 
        for(my $j = 0; $j < $window_size + 2; $j++){
            if(looks_like_number($table[$i][$j]) && ($table[$i][$j] - int($table[$i][$j]))){
                print sprintf("%.2f%6s", $table[$i][$j], "|");
            } else {
                print sprintf("%4s%6s", $table[$i][$j], "|");
            }
        } 
        print "\n";
        for(my $k= 1; $k < (2+$window_size) * 10+1; $k++){
            if($k % 10 == 0){
                print "|";
            } else {
                print "_";
            }
        }
        print "\n";
        
    }
}

print_table("1",$window_size, @table_1);

my @sums = [];
for(my $i = 1; $i < $window_size + 2 ; $i++){ 
    my $sum = 0;
    for(my $j = 1; $j < 5; $j++){
        $sum += $table_1[$j][$i];
    } 
    push @sums, $sum;
}

for(my $i = 1; $i < 5 ; $i++){ 
    for(my $j = 1; $j < $window_size + 2; $j++){
        my $value = $table_1[$i][$j];
        $table_2[$i][$j] = $value / $sums[$j];
    } 
}
print_table("2",$window_size, @table_2);

for(my $i = 1; $i < 5 ; $i++){ 
    for(my $j = 1; $j < $window_size + 2; $j++){
        my $value = $table_1[$i][$j];
        $table_3[$i][$j] = ($value + 0.25) / ($sums[$j] + 1);
    } 
}
print_table("3",$window_size, @table_3);