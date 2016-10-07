#chain_filter.py

This is just a simple python script to filter out chained alignments. It consists of two simple python scripts
The first one 'filter_self.py' removes exact contig-contig alignments identified when mapping a set of contigs back to itself.
This step is unnecessary if you are identifying matches in distinct reference/query sequences.

After you've identified your matches, filtered and listed the coordinates, you can run synteny_filter.py on the *.coords file
to identify chains longer than a specified number of links.

The basic flow is:

### Identify matches between sequences:
$ nucmer -p out ref.fasta query.fasta

### If you are mapping a set of contigs back to itself, you must filter out self hits prior to using delta-filter
$ filter_self.py out.delta

### Filter repetitive elements from delta file
$ delta-filter -r -q out_noself.delta >out_noself_filtered.delta

### List coordinates of identified sequence matches
$ show-coords -H -T -L 1000 out_noself_filtered.delta >out_noself_filtered.coords

### Parse coordinates file for chained alignments longer than a specified number of links (default=5)
$ synteny_filter.py out_noself_filtered.coords --chain_length 5


Outputs:
shared_homologies.coords - parsed version of input file containing only homologies of certain chain length and over
shared_homologies.out - summary file listing the number of bases / chained homologies for each pair of alignments
