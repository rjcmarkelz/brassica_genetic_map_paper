
Combine Chromosome and Scaffold bin genotypes:

```sh
cd data/bins/
cp bin-genotypes-2015-06-11 bin-genotypes.scaffolds-chromosomes.2015-07-13
sed '1d' bin-genotypes_ref1.5_v0.1.1.tsv >> bin-genotypes.scaffolds-chromosomes.2015-07-13
```
