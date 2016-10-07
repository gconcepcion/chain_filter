#chain_filter.py

This is just a simple python script to filter out chained alignments from a nucmer *.coords file.

It consists of two simple python scripts:

1) filter_self.py - removes exact contig-contig alignments identified when mapping a set of contigs back to itself.
2) chain_filter.py - filters a chained alignments from a *.coords file.

The filter_self.py script is necessary if you've aligned one contig set against itself and you need to filter out the
exact 1:1 contig hits that prevent the delta-filter from working properly. You can omit this step if you used two
distinct input fasta files for the initial nucmer search.

##The basic flow is:

1) Identify matches using `nucmer`
2) Filter self hits using `filter_self.py` (optional if input fastas were distinct)
3) Filter repetitive hits using `delta-filter`
4) List coordinates of alignments using `show-coords`
5) Filter chains greater than specified length using `chain_filter.py`


```bash
  $ nucmer -p out ref.fasta query.fasta
  $ filter_self.py out.delta
  $ delta-filter -r -q out_noself.delta >out_noself_filtered.delta
  $ show-coords -H -T -L 1000 out_noself_filtered.delta >out_noself_filtered.coords*
  $ chain_filter.py out_noself_filtered.coords --chain_length 5
```

###Outputs:
  * shared_homologies.coords - parsed version of input file containing only homologies of certain chain length and over
  * shared_homologies.out - summary file listing the number of bases / chained homologies for each pair of alignments
