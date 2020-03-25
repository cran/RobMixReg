# Robust Mixture Regression
r package

### Robust Mixture Regression Plot (with outliers)
<p align="center">
  <img width="400"  src="https://github.com/changwn/RobMixReg/blob/master/pic1.png">
</p>

### Add Regresseion Line
<p align="center">
  <img width="350"  src="https://github.com/changwn/RobMixReg/blob/master/pic2.png">
</p>

# Install from CRAN
```
install.packages("RobMixReg)
library("RobMixReg")
```

# Install from github for most updated package. 
#### Please report the bug as the description in the Question&Problem.
```
library("devtools")
devtools::install_github("changwn/RobMixReg")
```

# Example
```
library(RobMixReg)
#library(robust)
library(flexmix)
library(robustbase)
library(MASS)
library(gtools)

# gaussData
x=(gaussData$x);y=as.numeric(gaussData$y);
formula01=as.formula("y~x")
example_data01=data.frame(x,y)

res_rmr = rmr(lr.method='flexmix', formula=formula01, data=example_data01)
res_rmr = rmr(lr.method='TLE', formula=formula01, data=example_data01)
res_rmr = rmr(lr.method='CTLE', formula=formula01, data=example_data01)
res_rmr = rmr(lr.method='mixbi', formula=formula01, data=example_data01)
res_rmr = rmr(lr.method='mixLp', formula=formula01, data=example_data01)

# simuData
example_data02 <- simuData[,1:3]
formula02=as.formula("y~X1+X2")

res_rmr = rmr(lr.method='flexmix', formula=formula01, data=example_data01, nc=3)
res_rmr = rmr(lr.method='TLE', formula=formula01, data=example_data01, nc=3,tRatio=0.05)
res_rmr = rmr(lr.method='CTLE', formula=formula01, data=example_data01, nc=3)
res_rmr = rmr(lr.method='mixbi', formula=formula01, data=example_data01, nc=3)
res_rmr = rmr(lr.method='mixLp', formula=formula01, data=example_data01, nc=3)

```
# Questions & Problems

If you have any questions or problems, please feel free to open a new issue [here](https://github.com/changwn/RMR/issues). We will fix the new issue ASAP.  You can also email the maintainers and authors below.

- [Wennan Chang](https://zcslab.github.io/people/wennan/)
(wnchang@iu.edu)

PhD candidate at BDR group, Indiana University School of Medicine

- [Sha Cao](https://medicine.iu.edu/faculty/38873/cao-sha/)
(shacao@iu.edu)

Assistant Professor

Department of Biostatistics, Indiana University School of Medicine
