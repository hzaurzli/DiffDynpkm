# DiffDynpkm
fpkm/rpkm/cpm/tpm Calculate differential expression and predic dynamic expression in different conditions

# Install package
```
devtools::install_github('hzaurzli/DiffDynpkm')
library(DiffDynpkm)
```

# How to use DiffDynpkm
## Calculate differential expression
### 1.without NB model to calculate logFoldChange
#### [1].glm model p value calculation
```
data("gene_fpkm")
# treat is the treatment group
treat = gene_fpkm[,1:3]
# control is the control group
control = gene_fpkm[,4:6]
result = cal_diff(treat,control,method = "glm",shrink = F)

result

          Treat_1    Treat_2   Treat_3 Control_1  Control_2 Control_3    logFC        p_val        p,adj
Gene_1   9.168998 40.6542223  6.494375  7.444491 6.66855797 6.7189040 1.434788 5.395155e-02 1.000000e+00
Gene_2  19.655950  3.4254479 52.739200  7.618138 2.92421019 4.9873463 2.287560 1.325098e-02 2.650196e-01
Gene_3  72.357121 26.7231457 68.143043  3.603424 5.53719624 8.6450711 3.232987 2.247679e-09 4.495359e-08
Gene_4   4.426875 91.5598305 56.573256  9.084924 1.37014205 6.7992116 3.144350 2.788986e-03 5.577971e-02
Gene_5   3.485587 55.6122392 81.369367  3.899194 8.22648960 0.5721399 3.467580 2.741380e-03 5.482761e-02
Gene_6  81.364895 58.3585172 29.646583  5.893894 6.75695019 4.7783573 3.280600 1.601472e-10 3.202945e-09
Gene_7  15.648395 91.3456218 25.503929  7.391459 3.91020792 7.6312400 2.807002 1.936954e-04 3.873908e-03
Gene_8  51.295896 96.9578617 31.766921  2.902030 0.48393840 8.8162287 3.882950 2.127887e-07 4.255773e-06
Gene_9  73.491901 95.1658370 52.630370  9.224134 4.93290758 8.4202829 3.292979 1.637155e-17 3.274309e-16
Gene_10 35.228322 17.0702782 16.839530  7.279689 1.88318049 9.8737986 1.860701 4.552902e-04 9.105804e-03
Gene_11 25.636701 43.8891836 32.841739  2.142500 4.23089056 8.6638516 2.767148 2.317405e-10 4.634809e-09
Gene_12 75.291490 62.2583791 82.054535  3.784741 0.10510219 2.5552558 5.090561 1.055756e-18 2.111512e-17
Gene_13  8.608052 32.1660879 73.823465  1.734581 5.53250048 5.3283339 3.185606 1.989149e-04 3.978299e-03
Gene_14 72.793265 98.0866878 46.147087  9.046380 0.01719581 1.4083363 4.373278 1.843331e-08 3.686663e-07
Gene_15 94.116855 21.9074829 29.551971  2.216388 0.48820646 5.8018264 4.097079 2.924627e-06 5.849253e-05
Gene_16 27.323464 11.8033823 51.398272  6.677323 7.31623466 7.9694762 2.043241 4.405745e-04 8.811489e-03
Gene_17 66.353787 52.8402758 62.794364  4.493450 1.10471320 7.5514195 3.790758 3.536471e-20 7.072942e-19
Gene_18 77.868354  0.8073661 31.499290  5.529022 2.05673142 1.2623550 3.638284 5.188957e-03 1.037791e-01
Gene_19 41.133018 41.6323200 62.676640  9.832506 5.18450921 1.1238172 3.171657 4.052762e-13 8.105523e-12
Gene_20 50.613347  9.4147402 79.222878  3.229073 9.21405144 9.2609417 2.681650 7.120254e-04 1.424051e-02
```
#### [2].wlicox p value calculation
```
data("gene_fpkm")
# treat is the treatment group
treat = gene_fpkm[,1:3]
# control is the control group
control = gene_fpkm[,4:6]
result = cal_diff(treat,control,method = "wilcox",shrink = F)

result

          Treat_1    Treat_2   Treat_3 Control_1  Control_2 Control_3    logFC p_val p,adj
Gene_1   9.168998 40.6542223  6.494375  7.444491 6.66855797 6.7189040 1.434788  0.50     1
Gene_2  19.655950  3.4254479 52.739200  7.618138 2.92421019 4.9873463 2.287560  0.25     1
Gene_3  72.357121 26.7231457 68.143043  3.603424 5.53719624 8.6450711 3.232987  0.25     1
Gene_4   4.426875 91.5598305 56.573256  9.084924 1.37014205 6.7992116 3.144350  0.50     1
Gene_5   3.485587 55.6122392 81.369367  3.899194 8.22648960 0.5721399 3.467580  0.50     1
Gene_6  81.364895 58.3585172 29.646583  5.893894 6.75695019 4.7783573 3.280600  0.25     1
Gene_7  15.648395 91.3456218 25.503929  7.391459 3.91020792 7.6312400 2.807002  0.25     1
Gene_8  51.295896 96.9578617 31.766921  2.902030 0.48393840 8.8162287 3.882950  0.25     1
Gene_9  73.491901 95.1658370 52.630370  9.224134 4.93290758 8.4202829 3.292979  0.25     1
Gene_10 35.228322 17.0702782 16.839530  7.279689 1.88318049 9.8737986 1.860701  0.25     1
Gene_11 25.636701 43.8891836 32.841739  2.142500 4.23089056 8.6638516 2.767148  0.25     1
Gene_12 75.291490 62.2583791 82.054535  3.784741 0.10510219 2.5552558 5.090561  0.25     1
Gene_13  8.608052 32.1660879 73.823465  1.734581 5.53250048 5.3283339 3.185606  0.25     1
Gene_14 72.793265 98.0866878 46.147087  9.046380 0.01719581 1.4083363 4.373278  0.25     1
Gene_15 94.116855 21.9074829 29.551971  2.216388 0.48820646 5.8018264 4.097079  0.25     1
Gene_16 27.323464 11.8033823 51.398272  6.677323 7.31623466 7.9694762 2.043241  0.25     1
Gene_17 66.353787 52.8402758 62.794364  4.493450 1.10471320 7.5514195 3.790758  0.25     1
Gene_18 77.868354  0.8073661 31.499290  5.529022 2.05673142 1.2623550 3.638284  0.50     1
Gene_19 41.133018 41.6323200 62.676640  9.832506 5.18450921 1.1238172 3.171657  0.25     1
Gene_20 50.613347  9.4147402 79.222878  3.229073 9.21405144 9.2609417 2.681650  0.25     1
```

## 2.with NB model to calculate logFoldChange
#### [1].glm model p value calculation
```
data("gene_fpkm")
# treat is the treatment group
treat = gene_fpkm[,1:3]
# control is the control group
control = gene_fpkm[,4:6]
result = cal_diff(treat,control,method = "glm",shrink = T)

result

          Treat_1    Treat_2   Treat_3 Control_1  Control_2 Control_3    logFC        p_val        p,adj
Gene_1   9.168998 40.6542223  6.494375  7.444491 6.66855797 6.7189040 1.836622 5.395155e-02 1.000000e+00
Gene_2  19.655950  3.4254479 52.739200  7.618138 2.92421019 4.9873463 2.793592 1.325098e-02 2.650196e-01
Gene_3  72.357121 26.7231457 68.143043  3.603424 5.53719624 8.6450711 3.235614 2.247679e-09 4.495359e-08
Gene_4   4.426875 91.5598305 56.573256  9.084924 1.37014205 6.7992116 3.407558 2.788986e-03 5.577971e-02
Gene_5   3.485587 55.6122392 81.369367  3.899194 8.22648960 0.5721399 3.486260 2.741380e-03 5.482761e-02
Gene_6  81.364895 58.3585172 29.646583  5.893894 6.75695019 4.7783573 3.436710 1.601472e-10 3.202945e-09
Gene_7  15.648395 91.3456218 25.503929  7.391459 3.91020792 7.6312400 3.234502 1.936954e-04 3.873908e-03
Gene_8  51.295896 96.9578617 31.766921  2.902030 0.48393840 8.8162287 3.463808 2.127887e-07 4.255773e-06
Gene_9  73.491901 95.1658370 52.630370  9.224134 4.93290758 8.4202829 3.332764 1.637155e-17 3.274309e-16
Gene_10 35.228322 17.0702782 16.839530  7.279689 1.88318049 9.8737986 1.720784 4.552902e-04 9.105804e-03
Gene_11 25.636701 43.8891836 32.841739  2.142500 4.23089056 8.6638516 2.552892 2.317405e-10 4.634809e-09
Gene_12 75.291490 62.2583791 82.054535  3.784741 0.10510219 2.5552558 4.801232 1.055756e-18 2.111512e-17
Gene_13  8.608052 32.1660879 73.823465  1.734581 5.53250048 5.3283339 3.319636 1.989149e-04 3.978299e-03
Gene_14 72.793265 98.0866878 46.147087  9.046380 0.01719581 1.4083363 3.484497 1.843331e-08 3.686663e-07
Gene_15 94.116855 21.9074829 29.551971  2.216388 0.48820646 5.8018264 3.913389 2.924627e-06 5.849253e-05
Gene_16 27.323464 11.8033823 51.398272  6.677323 7.31623466 7.9694762 2.426139 4.405745e-04 8.811489e-03
Gene_17 66.353787 52.8402758 62.794364  4.493450 1.10471320 7.5514195 3.455286 3.536471e-20 7.072942e-19
Gene_18 77.868354  0.8073661 31.499290  5.529022 2.05673142 1.2623550 4.328954 5.188957e-03 1.037791e-01
Gene_19 41.133018 41.6323200 62.676640  9.832506 5.18450921 1.1238172 3.004767 4.052762e-13 8.105523e-12
Gene_20 50.613347  9.4147402 79.222878  3.229073 9.21405144 9.2609417 2.888848 7.120254e-04 1.424051e-02
```
#### [2].wilcox p value calculation
```
data("gene_fpkm")
# treat is the treatment group
treat = gene_fpkm[,1:3]
# control is the control group
control = gene_fpkm[,4:6]
result = cal_diff(treat,control,method = "wilcox",shrink = T)

result

          Treat_1    Treat_2   Treat_3 Control_1  Control_2 Control_3    logFC p_val p,adj
Gene_1   9.168998 40.6542223  6.494375  7.444491 6.66855797 6.7189040 1.836622  0.50     1
Gene_2  19.655950  3.4254479 52.739200  7.618138 2.92421019 4.9873463 2.793592  0.25     1
Gene_3  72.357121 26.7231457 68.143043  3.603424 5.53719624 8.6450711 3.235614  0.25     1
Gene_4   4.426875 91.5598305 56.573256  9.084924 1.37014205 6.7992116 3.407558  0.50     1
Gene_5   3.485587 55.6122392 81.369367  3.899194 8.22648960 0.5721399 3.486260  0.50     1
Gene_6  81.364895 58.3585172 29.646583  5.893894 6.75695019 4.7783573 3.436710  0.25     1
Gene_7  15.648395 91.3456218 25.503929  7.391459 3.91020792 7.6312400 3.234502  0.25     1
Gene_8  51.295896 96.9578617 31.766921  2.902030 0.48393840 8.8162287 3.463808  0.25     1
Gene_9  73.491901 95.1658370 52.630370  9.224134 4.93290758 8.4202829 3.332764  0.25     1
Gene_10 35.228322 17.0702782 16.839530  7.279689 1.88318049 9.8737986 1.720784  0.25     1
Gene_11 25.636701 43.8891836 32.841739  2.142500 4.23089056 8.6638516 2.552892  0.25     1
Gene_12 75.291490 62.2583791 82.054535  3.784741 0.10510219 2.5552558 4.801232  0.25     1
Gene_13  8.608052 32.1660879 73.823465  1.734581 5.53250048 5.3283339 3.319636  0.25     1
Gene_14 72.793265 98.0866878 46.147087  9.046380 0.01719581 1.4083363 3.484497  0.25     1
Gene_15 94.116855 21.9074829 29.551971  2.216388 0.48820646 5.8018264 3.913389  0.25     1
Gene_16 27.323464 11.8033823 51.398272  6.677323 7.31623466 7.9694762 2.426139  0.25     1
Gene_17 66.353787 52.8402758 62.794364  4.493450 1.10471320 7.5514195 3.455286  0.25     1
Gene_18 77.868354  0.8073661 31.499290  5.529022 2.05673142 1.2623550 4.328954  0.50     1
Gene_19 41.133018 41.6323200 62.676640  9.832506 5.18450921 1.1238172 3.004767  0.25     1
Gene_20 50.613347  9.4147402 79.222878  3.229073 9.21405144 9.2609417 2.888848  0.25     1
```


##Predic dynamic expression
```
data('exp')

       heart-0_1 heart-5_1 heart-10_1 heart-0_2 heart-5_2   heart-10_2 brain-0_1 brain-5_1 brain-10_1 brain-0_2 brain-5_2 brain-10_2
gene_1 0.7943423 0.4398317  0.7544752 0.6292211 0.7101824 0.0006247733 0.4753166 0.2201189  0.3798165 0.6127710 0.3517979  0.1111354
gene_2 0.3181810 0.2316258  0.1428000 0.4145463 0.4137243 0.3688454509 0.1524447 0.1388061  0.2330341 0.4659625 0.2659726  0.8578277
gene_3 0.4667790 0.5115055  0.5999890 0.3328235 0.4886130 0.9544738275 0.4829024 0.8903502  0.9144382 0.6087350 0.4106898  0.1470947
gene_4 0.1419069 0.6900071  0.6192565 0.8913941 0.6729991 0.7370777379 0.5211357 0.6598384  0.8218055 0.7862816 0.9798219  0.4394315
gene_5 0.4447680 0.2179907  0.5022996 0.3539046 0.6499852 0.3747139566 0.3554454 0.5336879  0.7403344 0.2211029 0.4127461  0.2656867
          gut-0_1   gut-5_1   gut-10_1   gut-0_2   gut-5_2  gut-10_2  lung-0_1   lung-5_1 lung-10_1  lung-0_2  lung-5_2  lung-10_2
gene_1 0.24361947 0.6680556 0.41764678 0.7881958 0.1028646 0.4348927 0.9849570 0.89305111 0.8864691 0.1750527 0.1306957 0.65310193
gene_2 0.04583117 0.4422001 0.79892485 0.1218993 0.5609480 0.2065314 0.1275317 0.75330786 0.8950454 0.3744628 0.6651152 0.09484066
gene_3 0.93529980 0.3012289 0.06072057 0.9477269 0.7205963 0.1422943 0.5492847 0.95409124 0.5854834 0.4045103 0.6478935 0.31982062
gene_4 0.31170220 0.4094750 0.01046711 0.1838495 0.8427293 0.2311618 0.2391000 0.07669117 0.2457237 0.7321352 0.8474532 0.49752727
gene_5 0.62997305 0.1838285 0.86364411 0.7465680 0.6682846 0.6180179 0.3722381 0.52983569 0.8746823 0.5817501 0.8397678 0.31244816
       liver-0_1 liver-5_1 liver-10_1 liver-0_2 liver-5_2 liver-10_2
gene_1 0.3435165 0.6567581  0.3203732 0.1876911 0.7822943 0.09359499
gene_2 0.3839696 0.2743836  0.8146400 0.4485163 0.8100644 0.81238951
gene_3 0.3077200 0.2197676  0.3694889 0.9842192 0.1542023 0.09104400
gene_4 0.3879090 0.2464490  0.1110965 0.3899944 0.5719353 0.21689276
gene_5 0.7082903 0.2650178  0.5943432 0.4812898 0.2650327 0.56459043

data('design)

   group_name rep tissue time
1   heart-0_1   1      3    0
2   heart-5_1   2      3    5
3  heart-10_1   3      3   10
4   heart-0_2   1      3    0
5   heart-5_2   2      3    5
6  heart-10_2   3      3   10
7   brain-0_1   4      1    0
8   brain-5_1   5      1    5
9  brain-10_1   6      1   10
10  brain-0_2   4      1    0
11  brain-5_2   5      1    5
12 brain-10_2   6      1   10
13    gut-0_1   7      2    0
14    gut-5_1   8      2    5
15   gut-10_1   9      2   10
16    gut-0_2   7      2    0
17    gut-5_2   8      2    5
18   gut-10_2   9      2   10
19   lung-0_1  10      5    0
20   lung-5_1  11      5    5
21  lung-10_1  12      5   10
22   lung-0_2  10      5    0
23   lung-5_2  11      5    5
24  lung-10_2  12      5   10
25  liver-0_1  13      4    0
26  liver-5_1  14      4    5
27 liver-10_1  15      4   10
28  liver-0_2  13      4    0
29  liver-5_2  14      4    5
30 liver-10_2  15      4   10

cal_dyn(exp = exp,design = design,k=c(3,3),gene = 'gene_1',bs = 'cr')
```
