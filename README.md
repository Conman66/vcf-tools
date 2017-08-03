VCF Tools
-----------------------

To find the shared SNPs between two different VCF files and exclude the SNPs found in a third, use VCF_Match.py:

```
./VCF_Match.py VCF_1 VCF_2 -r remove_VCF_3 > shared_SNPs.txt
```

To convert positions from one namespace to another with a delta file, use map_vcf.py:

```
./map_vcf.py shared_SNPs.txt delta_file
```

To weed out the SNPs that aren't within core regions, use in-core.py:

```
./in-core.py -i shared_SNPs.txt.map
```

Finally, to get the genes that these SNPs are in, use in-gene.pl:

```
perl in-gene.pl ref_file shared_SNPs.txt.map.core > genes.txt
```

