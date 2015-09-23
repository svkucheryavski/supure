Sprectra Unmixing with Pure Variables
===========================================

supure is an R package implementing methods for spectral unmixing (curve resolution) based on
the pure variables approach. It includes SIMPLSIMA and a simple implementation of MCR-ALS with 
angle conctraint. The package supplements a chapter in a book RESOLVING SPECTRAL MIXTURES, DATA 
HANDLING IN SCIENCE AND TECHNOLOGY, edited by Cyril Ruckebusch which will be published by Elsevier
in 2016. See `?supure` for all details.

How to install
--------------
To get the latest release plase use GitHub sources. You can [download](https://github.com/svkucheryavski/supure/releases) a zip-file with source package and 
install it using the `install.packages` command, e.g. if the downloaded file is 
`supure_0.0.1.tar.gz` and it is located in a current working directory, just run the following:

```
install.packages('mdatools_0.0.1.tar.gz')
```

If you have `devtools` package installed, the following command will install the latest release 
from the master branch of GitHub repository (do not forget to load the `devtools` package first):

```
install_github('supure', username = 'svkucheryavski')
```
