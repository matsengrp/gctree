#!/usr/bin/env perl
#Kenneth B. Hoehn (kenneth.hoehn@oriel.ox.ac.uk)

# Script written by Kenneth B. Hoehn and used in article:
# A Phylogenetic Codon Substitution Model for Antibody Lineages
# http://www.genetics.org/content/early/2017/03/15/genetics.116.196303
=pod
@article{hoehn2017phylogenetic,
  title={A phylogenetic codon substitution model for antibody lineages},
  author={Hoehn, Kenneth B and Lunter, Gerton and Pybus, Oliver G},
  journal={Genetics},
  volume={206},
  number={1},
  pages={417--427},
  year={2017},
  publisher={Genetics Soc America}
}
=cut


use strict;
use warnings;
use PDL;
use PDL::LinearAlgebra::Trans;
my $SMALL = 1e-200;

sub Fill_upp{
    my $node = $_[0];
    my $seqs = $_[1];
    my $Qs = $_[2]; #array reference of q matrixes
    my $partition = $_[3]; #array of partition indexes
    my $nparts = $_[4]; #number of unique partitions
  
    my @keys = keys %$seqs;
    my $length = length($seqs->{$keys[0]})/3;
    my @bigempty;
    for(my $i=0; $i < $length; $i++){
        my @empty = (-200)x61;
        push(@bigempty,\@empty);
    }
    $node->{"uppmat"}=\@bigempty;

   if($node->{"level"}==0){#if at root node
    #don't do anything, because this isn't a real node :)
   }elsif($node->{"level"}==1){ #one below root
    my $other; #other node to either the left or right
    if($node->{"up"}->{"left"} eq $node){$other="right"}
    elsif($node->{"up"}->{"right"} eq $node){$other="left"}
    else{die("something weird happened")}
    my @Pxz;
      for(my $i=0;$i<$nparts;$i++){
        push(@Pxz,mexp($Qs->[$i]*$node->{"up"}->{$other}->{"dist"}));
      }
      for(my $i=0; $i < $length;$i++){
          for(my $j=0;$j<61;$j++){
              my $sumxz;
              if ($Pxz[$partition->[$i]]->at(0,$j) > 0) {
                $sumxz = log($Pxz[$partition->[$i]]->at(0,$j))+$node->{"up"}->{$other}->{"mat"}->[$i][0];
              } else {
                $sumxz = log($SMALL)+$node->{"up"}->{$other}->{"mat"}->[$i][0];
              }
              for(my $k=1;$k<61;$k++){
                  if($Pxz[$partition->[$i]]->at($k,$j)==0){print $node->{"dist"}." $k $j\n"}
                  my $pxz;
                  if ($Pxz[$partition->[$i]]->at($k,$j) > 0) {
                    $pxz = log($Pxz[$partition->[$i]]->at($k,$j)) + $node->{"up"}->{$other}->{"mat"}->[$i][$k];
                  } else {
                    $pxz = log($SMALL) + $node->{"up"}->{$other}->{"mat"}->[$i][$k];
                  }
                  $sumxz = $sumxz + log(1+exp($pxz-$sumxz));
              }
              $node->{"uppmat"}->[$i][$j] = $sumxz;
          }
      }
   }else{#if not at the root node
    my $other; #other node to either the left or right
    if($node->{"up"}->{"left"} eq $node){$other="right"}
    elsif($node->{"up"}->{"right"} eq $node){$other="left"}
    else{die("something weird happened")}
    my @Pxy; my @Pyv;
    for(my $i=0;$i<$nparts;$i++){
      push(@Pxy,mexp($Qs->[$i]*$node->{"up"}->{"dist"}));
      push(@Pyv,mexp($Qs->[$i]*$node->{"up"}->{$other}->{"dist"}));
    }
    #pxy
    for(my $i=0; $i < $length;$i++){
        for(my $j=0;$j<61;$j++){
            my $sumxy;
            if ($Pxy[$partition->[$i]]->at($j,0) > 0) {
              $sumxy = log($Pxy[$partition->[$i]]->at($j,0))+$node->{"up"}->{"uppmat"}->[$i][0];
            } else {
              $sumxy = log($SMALL)+$node->{"up"}->{"uppmat"}->[$i][0];
            }
            for(my $k=1;$k<61;$k++){
                if($Pxy[$partition->[$i]]->at($k,$j)==0){print $node->{"up"}->{"dist"}." $k $j\n"}
                my $pxy;
                if ($Pxy[$partition->[$i]]->at($j,$k) > 0) {
                  $pxy = log($Pxy[$partition->[$i]]->at($j,$k)) + $node->{"up"}->{"uppmat"}->[$i][$k];
                } else {
                  $pxy = log($SMALL) + $node->{"up"}->{"uppmat"}->[$i][$k];
                }
                $sumxy = $sumxy + log(1+exp($pxy-$sumxy));
            }
            $node->{"uppmat"}->[$i][$j] = $sumxy;
        }
    }
    #pyv
    for(my $i=0; $i < $length;$i++){
        for(my $j=0;$j<61;$j++){
            my $sumyv;
            if ($Pyv[$partition->[$i]]->at(0,$j) > 0) {
              $sumyv = log($Pyv[$partition->[$i]]->at(0,$j))+$node->{"up"}->{$other}->{"mat"}->[$i][0];
            } else {
              $sumyv = log($SMALL)+$node->{"up"}->{$other}->{"mat"}->[$i][0];
            }
            for(my $k=1;$k<61;$k++){
                if($Pyv[$partition->[$i]]->at($k,$j)==0){print $node->{"up"}->{$other}->{"dist"}." $k $j ".$Pyv[$partition->[$i]]->at($k,$j)."\n"}
                my $pyv;
                if ($Pyv[$partition->[$i]]->at($k,$j) > 0) {
                  $pyv = log($Pyv[$partition->[$i]]->at($k,$j)) + $node->{"up"}->{$other}->{"mat"}->[$i][$k];
                } else {
                  $pyv = log($SMALL) + $node->{"up"}->{$other}->{"mat"}->[$i][$k];
                }
                $sumyv = $sumyv + log(1+exp($pyv-$sumyv));
            }
            $node->{"uppmat"}->[$i][$j] += $sumyv;
        }
    }
   }
   #tally up stuff
   my $upphood = 0;
   for(my $i=0; $i < $length;$i++){
    my $sumyv;
        for(my $j=0;$j<61;$j++){
            if($j==0){$sumyv=$node->{"uppmat"}->[$i][$j];}
            else{$sumyv = $sumyv + log(1+exp($node->{"uppmat"}->[$i][$j]-$sumyv));}
        }
        $upphood += $sumyv;
    }
    print "Upphood\t".$node->{"subtaxa"}."\t$upphood\n";

   if(exists($node->{"left"})){ #recurse!
    Fill_upp($node->{"right"},$seqs,$Qs,$partition,$nparts);
    Fill_upp($node->{"left"},$seqs,$Qs,$partition,$nparts);
   }  
}

sub Marginal_ASR{
  my $node = $_[0];
  my $seqs = $_[1];
  my $Qs = $_[2]; #array reference of q matrixes
  my $partition = $_[3]; #array of partition indexes
  my $nparts = $_[4]; #number of unique partitions

  my @keys = keys %$seqs;
  my $length = length($seqs->{$keys[0]})/3;

  if($node->{"level"} != 0){
    my @charmat;
    my @Pyv;
    for(my $i=0;$i<$nparts;$i++){
      push(@Pyv,mexp($Qs->[$i]*$node->{"dist"}));
  }
  my $lhood = 0;
  for(my $i=0; $i < $length;$i++){
    my @sitemat;#relative lhoods of each v at the site
    my $sitelhood;
    my $maxchar;
    my $maxlhood;
    for(my $v=0;$v<61;$v++){
        my $lhoodv;
        for(my $y=0;$y<61;$y++){
          my $val;
          if ($Pyv[$partition->[$i]]->at($v,$y) > 0) {
            $val = $node->{"uppmat"}->[$i][$y]+log($Pyv[$partition->[$i]]->at($v,$y))+$node->{"mat"}->[$i][$v];
          } else {
            $val = $node->{"uppmat"}->[$i][$y]+log($SMALL)+$node->{"mat"}->[$i][$v];
          }
          if($y==0){$lhoodv=$val;}
          else{$lhoodv = $lhoodv + log(1+exp($val-$lhoodv));}
       }
       push(@sitemat,$lhoodv);
       if($v == 0){$sitelhood = $lhoodv;}
       else{$sitelhood=$sitelhood + log(1+exp($lhoodv-$sitelhood));}
    }
    push(@charmat,\@sitemat);
    $lhood += $sitelhood;
  }
  $node->{"Codon_lhoods"}=\@charmat;
  print "Anc recon:\t".$node->{"subtaxa"}."\tLikelihood:\t$lhood\n";
  }
  if(exists($node->{"left"})){
    Marginal_ASR($node->{"right"},$seqs,$Qs,$partition,$nparts);
    Marginal_ASR($node->{"left"},$seqs,$Qs,$partition,$nparts);
  }
}

sub Pruning_Lhood{
    my $node = $_[0];
    my $seqs = $_[1];
    my $codoni = $_[2];
    my $Qs = $_[3]; #array reference of q matrixes
    my $partition = $_[4]; #array of partition indexes
    my $nparts = $_[5]; #number of unique partitions
    my $ambig_char = $_[6];

    my @keys = keys %$seqs;
    my $length = length($seqs->{$keys[0]})/3;
    my @bigempty;
    for(my $i=0; $i < $length; $i++){
        my @empty = (-200)x61;
        push(@bigempty,\@empty);
    }
    $node->{"mat"}=\@bigempty;

    if(!exists($seqs->{$node->{"id"}})){ #internal node
        my $r = Pruning_Lhood($node->{"right"},$seqs,$codoni,$Qs,$partition,$nparts,$ambig_char);
        my $l = Pruning_Lhood($node->{"left"},$seqs,$codoni,$Qs,$partition,$nparts,$ambig_char);
        my @Prs; my @Pls;
        for(my $i=0;$i<$nparts;$i++){
          push(@Prs,mexp($Qs->[$i]*$node->{"right"}->{"dist"}));
          push(@Pls,mexp($Qs->[$i]*$node->{"left"}->{"dist"}));
        }
        for(my $i=0; $i < $length;$i++){
            for(my $j=0;$j<61;$j++){
               my $sumr;
               my $suml;
               if ($Prs[$partition->[$i]]->at(0,$j) > 0) {
                 $sumr = log($Prs[$partition->[$i]]->at(0,$j))+$node->{"right"}->{"mat"}->[$i][0];
               } else {
                 $sumr = log($SMALL)+$node->{"right"}->{"mat"}->[$i][0];
               }
               if ($Pls[$partition->[$i]]->at(0,$j) > 0) {
                 $suml = log($Pls[$partition->[$i]]->at(0,$j))+$node->{"left"}->{"mat"}->[$i][0];
               } else {
                 $suml = log($SMALL)+$node->{"left"}->{"mat"}->[$i][0];
               }
                for(my $k=1;$k<61;$k++){
                    if($Prs[$partition->[$i]]->at($k,$j)==0){print $node->{"dist"}." $k $j\n"}
                      my $pr;
                      my $pl;
                      if ($Prs[$partition->[$i]]->at($k,$j) > 0) {
                        $pr = log($Prs[$partition->[$i]]->at($k,$j)) + $node->{"right"}->{"mat"}->[$i][$k];
                      } else {
                        $pr = log($SMALL) + $node->{"right"}->{"mat"}->[$i][$k];
                      }
                      if ($Pls[$partition->[$i]]->at($k,$j) > 0) {
                        $pl = log($Pls[$partition->[$i]]->at($k,$j)) + $node->{"left"}->{"mat"}->[$i][$k];
                      } else {
                        $pl = log($SMALL) + $node->{"left"}->{"mat"}->[$i][$k];
                      }
                    $sumr = $sumr + log(1+exp($pr-$sumr));
                    $suml = $suml + log(1+exp($pl-$suml));
                }
                $node->{"mat"}->[$i][$j] = $sumr+$suml;
            }
        }
    }else{ #external node
        my @s = split("",$seqs->{$node->{"id"}});
        my @t = @{transarrayCodon(\@s,$codoni)};
        for(my $i=0; $i < scalar(@{$node->{"mat"}});$i++){
            my $val=log(1);
            if($t[$i] ne "NA"){ #adjust for ambiguous sites
                $node->{"mat"}->[$i][$t[$i]]=$val;
            }else{
                #fill in equilibrium frequencies for ambiguous state
                for(my $j = 0; $j < 61; $j++){
                  if(!exists($ambig_char->{$node->{"id"}}->{$i}->[$j])){
                    print $node->{"id"}."\t$i\t$j\n";
                    print $ambig_char->{$node->{"id"}}."\n";
                    print $ambig_char->{$node->{"id"}}->{$i}."\n";
                    print $ambig_char->{$node->{"id"}}->{$i}->[$j]."\n";
                    die();
                  }
                  my $val = $ambig_char->{$node->{"id"}}->{$i}->[$j];
                  if($val == 0){$val=-200}
                  else{$val = log($val)}
                  $node->{"mat"}->[$i][$j] =  $val;
                }
            }
        }
    }
   my $lhood=0;
   for(my $i=0; $i < $length;$i++){
       my $slhood = $node->{"mat"}->[$i]->[0];
       for(my $j=1;$j<61;$j++){
           $slhood = $slhood + log(1+exp($node->{"mat"}->[$i]->[$j]-$slhood));
       }
       $lhood += $slhood;
   }
   if($node->{"level"} != 0){
       print $node->{"id"}."\t".$node->{"up"}->{"id"}."\t$lhood\n";
     }else{
       print $node->{"id"}."\t"."NONE"."\t$lhood\n";
    }
    return($lhood);
}

#Make Q matrix for HLP16
sub getQmat_HLP16{
    my $bij = $_[0];
    my $kappa =$_[1];
    my $omega = $_[2];
    my %freqs = %{$_[3]};
    my @codons = @{$_[4]};
    my $print = $_[5];
    my %tr = %{codonTable()};
    my %q;
    for(my $i=0; $i < scalar(@codons); $i++){
        my $from = $codons[$i];
        for(my $j=0; $j < scalar(@codons); $j++){
            my $to = $codons[$j];
            if($from eq $to){$q{$from.$to}=0;next;}
            my @diff = @{diffPos($from,$to)};
            if(scalar(@diff) > 1){
                $q{$from.$to}=$SMALL;
            }else{
                my @dc = sort {$a cmp $b} (substr($from,$diff[0],1),substr($to,$diff[0],1));
                if("@dc" eq "a g" || "@dc" eq "c t"){
                    if($tr{$from} eq $tr{$to}){
                        $q{$from.$to} = $kappa*$freqs{$to}*(1+$bij->[$i*61+$j]);
                    }else{
                        $q{$from.$to} = $omega*$kappa*$freqs{$to}*(1+$bij->[$i*61+$j]);
                    }
                }
                else{
                    if($tr{$from} eq $tr{$to}){
                        $q{$from.$to} = $freqs{$to}*(1+$bij->[$i*61+$j]);
                    }else{
                        $q{$from.$to} = $omega*$freqs{$to}*(1+$bij->[$i*61+$j]);
                    }
                }
            }
        }
    }
    my $trate = 0;
    foreach my $from (@codons){
        my $fsum = 0;
        foreach my $to (@codons){
            $fsum += $q{$from.$to};
        }
        $q{$from.$from}=-$fsum;
        $trate += $freqs{$from}*$fsum;
    }
    if($print){ print "$trate\n";}
    #check that rows sum to 0
    foreach my $from (@codons){
        my $fsum = 0;
        foreach my $to (@codons){
            $fsum += $q{$from.$to};
        }
        if($fsum > 0.0001){die("fsum > 0.0001f irst $from $fsum\n")}
    }
    #scale to a relative rate of 1
    my $rate = 0;
    foreach my $from (@codons){
        foreach my $to (@codons){
            $q{$from.$to} = $q{$from.$to}/$trate;
        }
        $rate -= $freqs{$from}*$q{$from.$from};
    }
    if($print){print "Mean rate: ".$rate."\n";} 
    #check that rows sum to 0
    foreach my $from (@codons){
        my $fsum = 0;
        foreach my $to (@codons){
            $fsum += $q{$from.$to};
        }
        if($fsum > 0.0001){die("$from $fsum\n")}
    }
 
    my @pdls;
    foreach my $from (@codons){
        my @row;
        foreach my $to (@codons){
            push(@row,$q{$from.$to})
        }
        push(@pdls,pdl([@row]));
    }
    my $Q = $pdls[0];
    for(my $i = 1; $i < 61; $i++){
        $Q = $Q->glue(1,$pdls[$i]);
    }
    return($Q);
}


#Read in a fasta file
sub getfasta{
  my($filename) = @_;
  chomp($filename);
  my $trying = "Unable to open $filename.\n";
  open(INPUTFILE, $filename) or die($trying);  # Create a new file
  my @filedata = <INPUTFILE>;
  close INPUTFILE;
  my %seqs;
  my $seqname;
  my $seq = '';
  if(scalar(@filedata)==0){return \%seqs}
  foreach my $line (@filedata) {
    if($line =~/^\s*$/) { next; }     # Ignore blank lines
    elsif($line =~/^\s*#/) { next; }  # Ignore comments
    elsif($line =~ /^\>\s*(.*)\s*$/) {
      my $temp = $1;
      $temp =~ s/\s//sg;
      if(length($seq) > 0) {
        $seq =~ s/\s//sg;
        $seqs{$seqname} = $seq;
        $seq = '';
      }
      $seqname = $temp;
      next;
    } else { $seq .= $line; }
  }
  $seq =~ s/\s//sg;
  $seqs{$seqname} = $seq;
  return \%seqs;
}

#Read in a rooted Newick tree
sub rootedNewick {
  my $in = $_[0];
  my $node = $_[1];
  my $level = $_[2];
  my $first; my $id;
  if($in =~ /(,|\(|\)|;)/){
    $first = $&;
    $node->{"id"} = $`;
    $node->{"level"} = $level;
    $in = $';
    if($first eq ","){$in = $first.$in}
  }else{
    die($in);
  }
  if($first eq "("){#left
    my %n;
    $node->{"left"} = \%n;
    $in = rootedNewick($in,\%n,$level+1);
  }
  elsif($first  eq ","){#up
    return($in);
  }
  elsif($first  eq ")"){#up
    return($in);
  }
  elsif($first  eq ";"){#up
    return($in);
  }
  my $second;
  my $pre;
  if($in =~ /(,|\(|\)|;)/){
    $second = $&;
    $in = $';
    $pre = $`;
  }else{
    die($in);
  }
  if($second eq ","){#right
    my %n;
    $node->{"right"} = \%n;
    $in = rootedNewick($in,\%n,$level+1);
  }elsif($second  eq ")"){#up
  }
  if($in =~ /(,|\(|\)|;)/){
    $node->{"id"} = $`;
    $in = $';
    if($& eq ","){$in = $&.$'}
  }else{
    die($in);
  }
  return($in);
}

#read in a rooted newick tree
sub readInRootedNewick{
  my $in = $_[0];
  my $printtree = $_[1];

   my %root; #set up root node
   my $oin = $in;
   $root{"dist"}=0;
   rootedNewick($in,\%root,0);
   getdists(\%root); #parse distance
   getdivergence(\%root,0); #get divergences

   my $t = printtreestring(\%root,"").";";
   if($t ne $oin){print "Tree read in incorrectly!Probably not an issue if you had zero branch lengths\n$oin\n$t\n";}
   else{print "Tree read in correctly\n";}
  
   if($printtree){
     $t = printtreestring(\%root,"").";";
     print("$oin\n$t\n");
   }
   return(\%root);
}

sub relevel{
  my $node = $_[0];
  my $increase = $_[1];
  if(exists($node->{"right"})){
    relevel($node->{"right"},$increase);
    relevel($node->{"left"},$increase);
  }
  $node->{"level"}=$node->{"level"}+$increase;
}


#Once tree is read in, need to get the branch lengths
#Once tree is read in, need to get the branch lengths
sub getdists{
  my $node = $_[0];
  if(exists($node->{"left"})){
    getdists($node->{"left"});
    getdists($node->{"right"});
  }
  if(!exists($node->{"id"})){
    die("Node ID doens't exist!");
  }else{
    if($node->{"id"}=~/\:/){
      $node->{"dist"} = $';
      $node->{"id"} = $`;
      if($node->{"dist"} == 0){
        $node->{"dist"} = 0.0000000000000000000001;
        print "zero length branch length caught\n";
      }
    }else{
      if($node->{"level"} != 0){
        die($node->{"id"}." level ".$node->{"level"}." is formatted incorrectly!");
      }
    }
  }
}

#Get the divergences for each node in the tree
sub getdivergence{
  my $node = $_[0];
  my $div = $_[1];
  $node->{"divergence"} = $div + $node->{"dist"};
  if(exists($node->{"left"})){
    getdivergence($node->{"left"},$node->{"divergence"});
    getdivergence($node->{"right"},$node->{"divergence"});
  }
}

#Print out the tree in Newick format to a string
sub printtreestring{
  my $node = $_[0];
  my $string = $_[1];
  if(exists($node->{"left"})){
    $string = $string."(";
    $string=printtreestring($node->{"left"},$string);
    $string = $string.",";
    $string=printtreestring($node->{"right"},$string);
    $string=$string.")";
  }
  if($node->{"level"} != 0){
    $string = $string.$node->{"id"}.":".$node->{"dist"};
  }
  return($string);
}


#translate codons to indexes in Q matrix
sub transarrayCodon{
  my @in = @{$_[0]};
  my %codoni = %{$_[1]};
  my @trans;
  for(my $i = 0; $i < scalar(@in); $i+=3){
    if(!exists($codoni{lc $in[$i].$in[$i+1].$in[$i+2]})){
      push(@trans,"NA");
    }else{
      push(@trans,$codoni{lc $in[$i].$in[$i+1].$in[$i+2]});
    }
  }
  return(\@trans);
}

#re-translate indexes to codons
sub untransarrayCodon{
  my @in = @{$_[0]};
  my @codons = @{$_[1]};
  my @trans;
  for(my $i = 0; $i < scalar(@in); $i++){
    if($in[$i] ne "NA"){
      push(@trans,$codons[$in[$i]]);
    }else{
      print "NA found $i\n";
      push(@trans,"NNN");
    }
  }
  return(\@trans);
}

#Print ML codon sequence at all nodes
sub printML_codon{
  my $node = $_[0];
  my $string = $_[1];
  my $codons = $_[2];
  if(exists($node->{"left"})){
    $string = printML_codon($node->{"left"},$string,$codons);
    $string = printML_codon($node->{"right"},$string,$codons);
  }

  if($node->{"level"}!=0){
  my @sequence;
  for(my $i=0; $i < scalar(@{$node->{"Codon_lhoods"}});$i++){
    my $max = -inf;
    my $maxchar = -1;
    for(my $j=0; $j < 61; $j++){
      my $n = $node->{"Codon_lhoods"}->[$i][$j];
      if($n > $max){
        $max = $n; 
        $maxchar=$j;
      }
    }
    push(@sequence,$maxchar);
  }
  $node->{"sequence"} = \@sequence;
  
  my @seq = @{untransarrayCodon($node->{"sequence"},$codons)};
  my $sequence = "";
  foreach my $s (@seq){$sequence = $sequence.$s;}
  $string = $string.">".$node->{"level"}.";".$node->{"subtaxa"}.";".$node->{"divergence"}."\n$sequence\n";
  }

  return($string);
}

#Print ML amino acid sequence at all nodes
sub printML_aa{
  my $node = $_[0];
  my $string = $_[1];
  my $codons = $_[2];
  my $aatable1 = $_[3];
  my $aatable3 = $_[4];
  if(exists($node->{"left"})){
    $string = printML_aa($node->{"left"}, $string,$codons,$aatable1,$aatable3);
    $string = printML_aa($node->{"right"},$string,$codons,$aatable1,$aatable3);
  }

  if($node->{"level"}!=0){
  my @sequence;
    for(my $i=0; $i < scalar(@{$node->{"Codon_lhoods"}});$i++){
      my %aas;
      for(my $j=0; $j < 61; $j++){ #tally up relative likelihoods of amino acids
        my $n = $node->{"Codon_lhoods"}->[$i][$j];
        my $aa = $aatable1->{uc $aatable3->{$codons->[$j]}}; #single letter amino acid
        if(exists($aas{$aa})){
          $aas{$aa} = $aas{$aa} + log(1+exp($n - $aas{$aa}))
        }else{
          $aas{$aa} = $n
        }
      }

      my $max = -inf;
      my $maxchar = -1;
      foreach my $key (keys %aas){
        my $n = $aas{$key};
        if($n > $max){
          $max = $n; 
          $maxchar=$key;
        }
      }
      push(@sequence,$maxchar);
    }
    #$node->{"sequence"} = \@sequence;
    my $sequence = "";
    foreach my $s (@sequence){$sequence = $sequence.$s;}
    $string = $string.">".$node->{"level"}.";".$node->{"subtaxa"}.";".$node->{"divergence"}."\n$sequence\n";
    }
  return($string);
}

#make subtaxa labels
sub getSubTaxa{
  my $node = $_[0];
  #if a tip
  if(!exists($node->{"right"})){
    $node->{"subtaxa"}=$node->{"id"};
  }else{
    #if internal node
    getSubTaxa($node->{"right"});
    getSubTaxa($node->{"left"});
    my @l = split(",",$node->{"left"}->{"subtaxa"});
    my @r = split(",",$node->{"right"}->{"subtaxa"});
    my @total = (@l,@r);
    @total = sort @total;
    $node->{"subtaxa"}=join(",",@total);
  }
}

#assign parents to each node
sub assignparents{
  my $node = $_[0];
  my $parent = $_[1];
  if($node->{"level"} != 0){
    $node->{"up"}=$parent;
  }
  if(exists($node->{"left"})){
    assignparents($node->{"left"},$node);
    assignparents($node->{"right"},$node);
  }
}


#Position of differences between two strings
sub diffPos{
  my $s1 = $_[0];
  my $s2 = $_[1];
  if(length($s1) != length($s2)){die("$s1 $s2 no same length")}
  my @diffpos;
  for(my $i = 0; $i < length($s1); $i++){
    if(substr($s1,$i,1) ne substr($s2,$i,1)){
      push(@diffpos,$i);
    }
  }
  return(\@diffpos);
}

#return codon translation table
sub codonTable{
  my %codons = (
"ttt" =>  "Phe",
"ttc" =>  "Phe",
"tta" =>  "Leu",
"ttg" =>  "Leu",
"ctt" =>  "Leu",
"ctc" =>  "Leu",
"cta" =>  "Leu",
"ctg" =>  "Leu",
"att" =>  "Ile",
"atc" =>  "Ile",
"ata" =>  "Ile",
"atg" =>  "Met",
"gtt" =>  "Val",
"gtc" =>  "Val",
"gta" =>  "Val",
"gtg" =>  "Val",
"tct" =>  "Ser",
"tcc" =>  "Ser",
"tca" =>  "Ser",
"tcg" =>  "Ser",
"cct" =>  "Pro",
"ccc" =>  "Pro",
"cca" =>  "Pro",
"ccg" =>  "Pro",
"act" =>  "Thr",
"acc" =>  "Thr",
"aca" =>  "Thr",
"acg" =>  "Thr",
"gct" =>  "Ala",
"gcc" =>  "Ala",
"gca" =>  "Ala",
"gcg" =>  "Ala",
"tat" =>  "Tyr",
"tac" =>  "Tyr",
"taa" =>  "STOP",
"tag" =>  "STOP",
"cat" =>  "His",
"cac" =>  "His",
"caa" =>  "Gln",
"cag" =>  "Gln",
"aat" =>  "Asn",
"aac" =>  "Asn",
"aaa" =>  "Lys",
"aag" =>  "Lys",
"gat" =>  "Asp",
"gac" =>  "Asp",
"gaa" =>  "Glu",
"gag" =>  "Glu",
"tgt" =>  "Cys",
"tgc" =>  "Cys",
"tga" =>  "STOP",
"tgg" =>  "Trp",
"cgt" =>  "Arg",
"cgc" =>  "Arg",
"cga" =>  "Arg",
"cgg" =>  "Arg",
"agt" =>  "Ser",
"agc" =>  "Ser",
"aga" =>  "Arg",
"agg" =>  "Arg",
"ggt" =>  "Gly",
"ggc" =>  "Gly",
"gga" =>  "Gly",
"ggg" =>  "Gly"
);
  return \%codons;
}

#return codon translation table
sub codonTableSingle{
  my %codons = (
"ALA" => "A",
"CYS" => "C",
"ASP" => "D",
"GLU" => "E",
"PHE" => "F",
"GLY" => "G",
"HIS" => "H",
"ILE" => "I",
"LYS" => "K",
"LEU" => "L",
"MET" => "M",
"ASN" => "N",
"PRO" => "P",
"GLN" => "Q",
"ARG" => "R",
"SER" => "S",
"THR" => "T",
"VAL" => "V",
"TRP" => "W",
"TYR" => "Y");
  return \%codons;
}


my $nsim;
my $kappa;
my @omegas;
my @motifs;
my $partfile;
my $treefile;
my @hs;
my $freqs;
my $length;
my $outdir;
my $igphyml;
my $seqfile;
my $context=0;
my $rooted=0;
my $rootid;
my $ancstate=0;
my $statsfile;
my $stem;
my $statslhood;
my $ambigfile;

my $bstats;


#Read in parameters from congif file
# my $config = $ARGV[0];
# open(C,$config) or die("Couldn't open config file ($config)");
if (open(my $C, $ARGV[0])) {
  print "Found readable config file at first command line argument. Now reading it...\n";
  while(defined(my $line = <$C>)){
    chomp($line);
    if($line =~ /nsim\s+(\S+)/){
      $nsim = $1;
    }
    if($line =~ /omegas\s+(\S+)/){
      @omegas = split(",",$1);
    }
    if($line =~ /kappa\s+(\S+)/){
      $kappa = $1;
    }
    if($line =~ /motifs\s+(\S+)/){
      @motifs = split(",",$1);
    }
    if($line =~ /hs\s+(\S+)/){
      @hs = split(",",$1);
    }
    if($line =~ /freqs\s+(\S+)/){
      $freqs = $1;
    }
    if($line =~ /tree\s+(\S+)/){
      $treefile = $1;
    }
    if($line =~ /fullcontext\s+(\S+)/){
      $context=$1;
    }
    if($line =~ /outdir\s+(\S+)/){
      $outdir=$1;
    }
    if($line =~ /rooted\s+(\S+)/){
      $rooted=$1;
    }
    if($line =~ /length\s+(\S+)/){
      $length=$1;
    }
    if($line =~ /rootid\s+(\S+)/){
      $rootid=$1;
    }
    if($line =~ /part\s+(\S+)/){
      $partfile=$1;
    }
    if($line =~ /igphyml\s+(\S+)/){
      $igphyml=$1;
    }
    if($line =~ /seqfile\s+(\S+)/){
      $seqfile=$1;
    }
    if($line =~ /stats\s+(\S+)/){
      $statsfile=$1;
    }
    if($line =~ /stem\s+(\S+)/){
      $stem=$1;
    }
    if($line =~ /ambigfile\s+(\S+)/){
      $ambigfile=$1;
    }
  }
close $C;
} else {
  print "Did not find a readable config file at first command line argument.\n";
}



#check to see if stats file was specified in command line
for(my $i = 1; $i < scalar(@ARGV); $i++){
  my $line = $ARGV[$i];
  if($line =~ /-stats/){
    $statsfile=$ARGV[$i+1];
  }
}

#Read in igphyml stats file, if specified
if(defined $statsfile && $statsfile ne "N"){
  open(STATS,$statsfile)or die("Couldn't open $statsfile\n");
  my @stats = <STATS>;
  @motifs = (0)x0;
  @omegas = (0)x0;
  @hs = (0)x0;
  $freqs = "stats";
  foreach my $l (@stats){
    chomp($l);
    if($l =~ /Motif:\s+(\S+)\s+\d\s+\d\s+(\S+)/){
      push(@motifs,$1);
      push(@hs,$2);
      print "Read motif h $1 = $2 from $statsfile\n";
    }
    if($l =~ /\. Omega\s+(\d+)\s+\S+:\s+(\S+)/){
      $omegas[$1]=$2;
      print "Read omega $1 = $2 from $statsfile\n";
    }
    if($l =~ /. Nonsynonymous\/synonymous ratio:\s+(\S+)/){
      $omegas[0]=$1;
      print "Read old school omega 0 = $1 from $statsfile\n";
    }
    if($l =~ /\. Transition\/transversion ratio:\s+(\S+)/){
      $kappa=$1;
      print "Read kappa $kappa from $statsfile\n";
    }
    if($l =~ /\. Log-likelihood:\s+(\S+)/){
      $statslhood=$1;
      print "Read stats lhood $statslhood from $statsfile\n";
    }
    $bstats .= " $l";
  }
  print "Reading in frequency parameters from IgPhyML stats file\n";
}

#check command line args to see if any should be over-ridden
for(my $i = 1; $i < scalar(@ARGV); $i++){
  my $line = $ARGV[$i];
  if($line =~ /-nsim/){
    $nsim = $ARGV[$i+1];
  }
  if($line =~ /-omegas/){
    @omegas = split(",",$ARGV[$i+1]);
  }
  if($line =~ /-kappa/){
    $kappa = $ARGV[$i+1];
  }
  if($line =~ /-part/){
    $partfile = $ARGV[$i+1];
  }
  if($line =~ /-motifs/){
    @motifs = split(",",$ARGV[$i+1]);
  }
  if($line =~ /-hs/){
    @hs = split(",",$ARGV[$i+1]);
  }
  if($line =~ /-freqs/){
    $freqs = $ARGV[$i+1];
  }
  if($line =~ /-tree/){
    $treefile = $ARGV[$i+1];
  }
  if($line =~ /-fullcontext/){
    $context=$ARGV[$i+1];
  }
  if($line =~ /-outdir/){
    $outdir=$ARGV[$i+1];
  }
  if($line =~ /-rooted/){
    $rooted=$ARGV[$i+1];
  }
  if($line =~ /-length/){
    $length=$ARGV[$i+1];
  }
  if($line =~ /-rootid/){
    $rootid=$ARGV[$i+1];
  }
  if($line =~ /-igphyml/){
    $igphyml=$ARGV[$i+1];
  }
  if($line =~ /-seqfile/){
    $seqfile=$ARGV[$i+1];
  }
  if($line =~ /-stats/){
    $statsfile=$ARGV[$i+1];
  }
  if($line =~ /-stem/){
    $stem=$ARGV[$i+1];
  }
  if($line =~ /ambigfile/){
    $ambigfile=$1;
  }
}

#check that all necessary parameters are specified
if(!defined $kappa){die("kappa needs to be specified")}
if(scalar(@omegas)==0){die("omegas needs to be specified")}
if(scalar(@motifs)==0){die("motifs needs to be specified")}
if(scalar(@hs)==0){die("hs needs to be specified")}
if(!defined $freqs){die("freqs needs to be specified")}
if(!defined $outdir){die("outdir needs to be specified")}
if(!defined $rootid){die("rootid needs to be specified")}
if(!defined $igphyml){die("igphyml needs to be specified")}
if(!defined $stem){die("stem needs to be specified")}
if(!defined $seqfile){die("seqfile needs to be specified");}
my $seqs = getfasta($seqfile);
if(!defined $length){die("Length needs to be specified or set explicitly to default by \"-length D\”")}
if ($length eq 'D') {
  my @keys = keys %$seqs;
  $length = length($seqs->{$keys[0]})/3;
}
if(!defined $statsfile){$statsfile="N";}
if(!defined $partfile){$partfile="N";}
if(!defined $ambigfile){$ambigfile="N";}

print "\nReconstruction Settings\n";
print "kappa\t$kappa\n";
print "omegas\t@omegas\n";
print "motifs\t@motifs\n";
print "hs\t@hs\n";
print "freqs\t$freqs\n";
print "outdir\t$outdir\n";
print "length\t$length\n";
print "rootid\t$rootid\n";
print "seqfile\t$seqfile\n";
print "stats\t$statsfile\n";
print "ambigfile\t$ambigfile\n";
print "outfile format: $outdir/$stem\_\n";
print "\n";

if(scalar(@hs) ne scalar(@motifs)){die(scalar(@hs)." h values but ".scalar(@motifs)." motifs!\n")}

#Read in tree
open(TREE,$treefile)or die("Couldn't open $treefile");
my $tree = <TREE>;
chomp($tree);
my %root;
if($rooted==1){
  %root = %{readInRootedNewick($tree,0)};
}else{
  %root = %{readInUnrootedNewick($tree,$rootid,0)};
}

#read in ambiguous character file
my %ambig;
  if($ambigfile ne "N"){
  open(AM,$ambigfile) or die();
  while(<AM>){
    my $line = $_;
    chomp($line);
    my @in = split(" ",$line);
    if(!exists($ambig{$in[0]})){
      my %new;
      $ambig{$in[0]} = \%new;
    }
    if(!exists($ambig{$in[0]}->{$in[1]})){
      my @new = (0)x61;
      $ambig{$in[0]}->{$in[1]} = \@new;
    }
    $ambig{$in[0]}->{$in[1]}->[$in[2]]=$in[3];
  }
}

#Set up partition model from file
my @part = ((-1)x$length);
my $nparts=0;
if($partfile ne "N"){
  open(P,$partfile) or die("Couldn't open $partfile");
  my $h = <P>;
  while(<P>){
    my $line = $_;
    chomp($line);
    my @in1 = split(":",$line);
    my @in2 = split(",",$in1[1]);
    for(my $i=0;$i<scalar(@in2);$i++){
      my @in3=split("\\.\\.",$in2[$i]);
      for(my $j=$in3[0];$j<=$in3[1];$j++){
        if($j >= $length){die("Partition $nparts extends beyond specified sequnece length!")}
        $part[$j]=$nparts;
      }
    }
    $nparts++;
  }
}else{ #default of a single partition
  @part = ((0)x$length);
  $nparts=1;
}
if($nparts != scalar(@omegas)){die("$nparts partitions, but ".(scalar(@omegas)." omegas!"))}
print "Partition index: ";
for(my$i=0;$i<$length;$i++){
  if($part[$i] == -1){
    die("Position $i unspecified in partition file.");
  }else{
    print "$part[$i] ";
  }
}
print "\n";

#Make codon indexes and get frequencies
my @transfreq;
my @codons;
my %codoni;
my $index = 0;
my @chars = ("t","c","a","g");
my %freqs;
my $fsum=0;
foreach my $a (@chars){
  foreach my $b (@chars){
    foreach my $c (@chars){
      if($a.$b.$c ne "tga" && $a.$b.$c ne "taa" && $a.$b.$c ne "tag"){
        push(@codons,$a.$b.$c);
        $codoni{$a.$b.$c}=$index;
        my $match = uc $a.$b.$c;
        if(defined $statsfile && $statsfile ne "N"){
          if($bstats =~ /f\($match\)=(\d+\.*\d*)/){
            $freqs{lc $match} = $1;
            if($freqs{lc $match} == 0){
              print "Zero frequency caught\n";
              $freqs{lc $match}=1e-10;
            }
          }else{die($match)}
        }elsif($freqs eq "uniform"){
          $freqs{lc $match} = 1/61;
        }else{
          die("freqs not set properly\n");
        }
        $fsum += $freqs{$a.$b.$c};
        $index++;
      }
    }
  }
}
foreach my $k (keys %freqs){
  $freqs{$k} = $freqs{$k}/$fsum;
  $transfreq[$codoni{$k}]=$freqs{$k};
}
print "Codon frequencies: @transfreq\n";

#Make B matrix
print "Setting up B matrix\n";
my @Bmat = (0)x(61*61);
my $fi;my $ti;my $li;my $ri;

for(my $mi = 0; $mi < scalar(@motifs); $mi++){
  my $motif = $motifs[$mi];
  print "Reading in $motif table\n";
  my @htable;
  open(HTABLE,"$igphyml/src/motifs/HTABLE_$motif") or die("Couldn't open $igphyml/src/motifs/HTABLE_$motif\n");
  while(<HTABLE>){
    my $l = $_;
    chomp($l);
    push(@htable,$l);
  }
  close(HTABLE);

  for($fi=0;$fi<61;$fi++){
    for($ti=0;$ti<61;$ti++){
      my @htotals = (0)x(1);
      for($li=0;$li<61;$li++){
        for($ri=0;$ri<61;$ri++){
          $htotals[0] += $transfreq[$li]*$transfreq[$ri]*$htable[$fi*61*61*61+$ti*61*61+$li*61+$ri];
        }
      }
      my $hsum = $hs[$mi]*$htotals[0];
      $Bmat[61*$fi+$ti]+=$hsum;
    }
  }
  if(scalar(@Bmat) != (61*61)){die("@Bmat")}
}

#Make Q matrices
my @Qs;
for(my $i=0;$i<scalar(@omegas);$i++){
  print "Making Q matrix $i, omega: $omegas[$i]\n";
  push(@Qs,getQmat_HLP16(\@Bmat,$kappa,$omegas[$i],\%freqs,\@codons,0))
}

#Check to see that the root sequence is in good order
my @rs;
if(defined $seqfile || $seqfile ne "N"){
  if(!exists($seqs->{$rootid})){
    die("$rootid not found in sequence file: $seqfile\n");
  }
  if(length($seqs->{$rootid}) ne $length*3){
    die("Specified root sequence is not specified length: $length ".(length($seqs->{$rootid})/3)."\n");
  }
}else{
  die("Need sequences and root!");
}

#make site/tip matrix for ambiguous characters
my @tipr;
my @keyst = keys %$seqs;
for(my $i = 0; $i < length($seqs->{$keyst[0]})/3; $i++){
  my @temp = (1e-10)x61;
  my $count=61e-10;
  foreach my $k (keys %$seqs){
    my @s = split("",$seqs->{$k});
    my @t = @{transarrayCodon(\@s,\%codoni)};
    if($t[$i] ne "NA"){
      $temp[$t[$i]]++;
      $count++;
    }
  }
  for(my $j=0;$j<61;$j++){
    $tipr[61*$i+$j] = $temp[$j]/$count
  }
}

#collect subtaxa of each node
getSubTaxa(\%root);

#assign node parents
assignparents(\%root,"NA");

print "About to do reconstructions\n";

#Do felsenstein's pruning algorithm
print "\nDoing lower partial likelihoods\n";
my $lhood = Pruning_Lhood(\%root,$seqs,\%codoni,\@Qs,\@part,$nparts,\%ambig);

print "\nDoing upper partial likelihoods\n";
Fill_upp(\%root,$seqs,\@Qs,\@part,$nparts);

print "\nCalculating marginal ASRs at each node\n";
#Do marginal ancestral state reconstructions
Marginal_ASR(\%root,$seqs,\@Qs,\@part,$nparts);

my $mcodons = printML_codon(\%root,"",\@codons);
#print "$mcodons\n";

my $codonTableSingle = codonTableSingle();
my $codonTableTriple = codonTable();

my $maa = printML_aa(\%root,"",\@codons,$codonTableSingle,$codonTableTriple);
#print "$maa\n";

print "Likelihood comparison - StatsFile: $statslhood, Reconstruction: $lhood\n";

# Die if the likelihoods don't match:
sub isnan { ! defined( $_[0] <=> 9**9**9 ) }
if (isnan($statslhood) or isnan($lhood)) {
  print "*** Failed to reconstruct the likelihood. See the line above. ***\n";
  #system("rm $outdir/gctree.simulation.fasta")
  die("*** Failed to reconstruct the likelihood. ***\n");
}
my @sort_array = ($lhood, $statslhood);
my $min_lhood = (sort {$a <=> $b} @sort_array)[0];
if (abs($lhood - $statslhood) / $min_lhood > 0.01) {
  print "*** Failed to reconstruct the likelihood with sufficient precision (>1%). See the line above. ***\n";
  #system("rm $outdir/gctree.simulation.fasta")
  die("*** Failed to reconstruct the likelihood with sufficient precision (>1%). ***\n");
}

open(OUT,">$outdir/$stem.MLcodons.fa") or die();
print OUT uc "$mcodons\n";
close(OUT);

open(OUT,">$outdir/$stem.MLaas.fa") or die();
print OUT uc "$maa\n";
close(OUT);


print "\n\n".'###############################################################################
#                            It worked - congrats!                            #
###############################################################################'."\n\n";
