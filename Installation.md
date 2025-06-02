# The package depends on R(>= 3.5.0)
--------------
## 1.Introduction
This package includes a function which reads a genotype file, processes the data, and detects change points.
## 2. Installation
- First of all, you need install dependencies required by GenoChangePoint under R, if you do not have them yet. Open a **R terminal**,
```R
install.package("changepoint")
```
- Secondly, you have two ways to intsall GenoCOPoint
  - way 1 : install with a local copy
    - 1.1 download a zipped GenoCOPoint and create a tar.gz
      ```linux
      cd path/to/target/folder
      wget https://github.com/StellarMech/GenoCOPoint/archive/master.zip
      unzip master.zip
      tar -czvf GenoCOPoint-master.tar.gz GenoCOPoint-master/
        ```
    - 1.2 install under R
      ```R
       install.packages("GenoCOPoint-master.tar.gz")
       q("no")
      ```
  - way 2 : directly install from github
    ```R
     install.packages("devtools")
     devtools::install_github("StellarMech/GenoCOPoint")
    ```