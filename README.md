# Non-Local Low-Rank Normal Filtering for Mesh Denoising
by Xianzhi Li, [Lei Zhu](http://appsrv.cse.cuhk.edu.hk/~lzhu/), [Chi-Wing Fu](https://www.cse.cuhk.edu.hk/~cwfu/) and [Pheng-Ann Heng](http://www.cse.cuhk.edu.hk/~pheng/).

### Introduction
This repository is for our Pacific Graphics 2018 paper 'Non-Local Low-Rank Normal Filtering for Mesh Denoising'. In this paper, we present a non-local low-rank normal filtering method for mesh denoising. By exploring the geometric similarity between local surface patches on 3D meshes in the form of normal fields, we devise a low-rank recovery model that filters normal vectors by means of patch groups.

In this repository, we release demo (exe files), code (C++), and data. 

### Citation
If you find our work useful in your research, please consider citing:
> @article{xianzhi2018nllr,
title={Non-local low-rank normal filtering for mesh denoising},
author={Li, Xianzhi and Lei, Zhu and Fu, Chi-Wing and Heng, Pheng-Ann},
journal={Computer Graphics Forum (Pacific Graphics)},
note={to appear}
}

### Usage
To try our method for mesh denoising, you can directly run the 'NLLR.exe' inside `demo.rar`.
```
upzip the demo.rar
copy noisy mesh into the demo folder (e.g., child_n3.off)
run: NLLR.exe child_n3.off 0.39 10 10
```
'child_n3.off' is the input noisy mesh, 0.39 is the sigma_M in our paper, 10 is the N_k and 10 is the number of vertex updating. You will see the denoised mesh inside the same folder. For the 'NLLR.exe', you can refer to `NLLR` folder for the source code.

To evaluate the denoising performance, you can directly run the 'evaluation.exe' inside `demo.rar`.
```
upzip the demo.rar
copy ground truth mesh into the demo folder (e.g., child.off)
run: evaluation.exe child.off denoised_result.off
```
'child.off' is the ground truth mesh, 'denoised_result.off' is the denoised mesh. You will see the mean square angle error (MSAE).

### Questions
Please constact 'lixianzhi123@gmail.com'

