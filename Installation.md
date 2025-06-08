# The package depends on R(>= 3.5.0)
--------------
## 1.Introduction
This package includes a function which reads a genotype file,a truth file, processes the data, and detects change points.
## 2. Installation
- First of all, you need install dependencies required by GenoChangePoint under R, if you do not have them yet. Open a **R terminal**,
```R
install.package("changepoint")
```
- Secondly, you should intsall GenoCOPoint-main.zip from [GenoCOPoint](https://github.com/StellarMech/GenoCOPoint)
  - way: unzip in your folder and directly install it
    ```R
     install.package("devtools")
     library(devtools)
     install()
    ```