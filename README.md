# UnseenEst

UnseenEst is an algorithm that estimates the frequency distribution of genetic variations in a population based on empirically observed site-frequency spectrum (also called fingerprint) in a sub-sample of individuals.

### Getting started

UnseenEst is written in Python2/3 
**Requirements**:
-[cxvopt](http://cvxopt.org/)
-[networkx](https://networkx.github.io/)
-[pandas](http://pandas.pydata.org/)
-[statsmodels](http://statsmodels.sourceforge.net/)
-[scipy](https://www.scipy.org/)
-[matplotlib](http://matplotlib.org/)

**Inputs**:
- The empirical site-frequency spectrum (SFS). This should be a text file where the 1st row contains the number of variants observed once in the dataset (i.e. found in one individual) and the ith row is the number of variants observed i times in the dataset (found in i individuals). 
- The total number of alleles in the dataset. This should be the number of sequenced individuals times 2. 

**Outputs**:
- UnseenEst outputs a text file with two columns. The first column indicates the frequency bin and the second column indicates the number of variants with population frequency in that bin. 

The command to run UnseenEst is
```
python unseen_est input_SFS_name alleles_numbers output_filename
```

### Example
The example SFS is generated from a dataset of 118396 alleles. 
```
python unseen_est.py ../example/example_input.txt 118396 ../example/example_output.txt
``` 

### Reference

Please cite:

Zou, J. et al. ***Quantifying the unobserved protein-coding variants in human populations provides a roadmap for large-scale sequencing projects.*** Nature Communications (2016).

### License

This project is licensed under the [MIT license](https://opensource.org/licenses/MIT).

