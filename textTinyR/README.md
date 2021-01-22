

## textTinyR
<br>

The *textTinyR* package consists of text processing functions for small or big data files. More details on the functionality of textTinyR can be found in [blog-post1](http://mlampros.github.io/2017/01/05/textTinyR_package/) and [blog-post2](http://mlampros.github.io/2018/04/04/extending_textTinyR_package/). The R package can be installed, in the following Operating Systems: Linux, Mac and Windows. However, there is one limitation : *chinese*, *japanese*, *korean*, *thai* or *languages with ambiguous word boundaries* are not supported.


<br>


**UPDATE 01-04-2018** : *boost-locale* is no longer a system requirement for the textTinyR package.


<br>


### **Installation of the textTinyR package (CRAN, Github)**

<br>

To install the package from CRAN use, 

```R

install.packages('textTinyR')


```
<br>

and to download the latest version from Github use the *install_github* function of the devtools package,
<br><br>

```R

devtools::install_github(repo = 'mlampros/textTinyR')


```
<br>
Use the following link to report bugs/issues,
<br><br>

[https://github.com/mlampros/textTinyR/issues](https://github.com/mlampros/textTinyR/issues)
