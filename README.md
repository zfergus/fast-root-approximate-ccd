# Fast Root Approximate in Cubic CCD

[![Build](https://github.com/zfergus/fast-approximate-root-ccd/actions/workflows/continuous.yml/badge.svg)](https://github.com/zfergus/fast-approximate-root-ccd/actions/workflows/continuous.yml)

Implementation of "Fast Root Approximate in Cubic CCD" from "Penetration-free Projective Dynamics on the GPU" [Lan et al. 2022]

## False Negatives (i.e., Missed Collision)

After testing this method can produce false negatives (i.e., miss collisions).

Note that while this is a faithful implementation of the paper's algorithmic description, technical implementation details may be missing from the paper that are required for the algorithm to work correctly. However, without the original authors' implementation, it is difficult to know.

### Example: Point-Triangle

This point-triangle query produces a false negative:
```c++
p_t0  <<  -5066329 / 1048576.0,  3747703 /  8388608.0,  -9110393 / 4194304.0;
t0_t0 <<  -5097349 / 1048576.0, 15565867 / 33554432.0,  -8419809 / 4194304.0;
t1_t0 << -10108719 / 2097152.0,  4441449 /  8388608.0,  -2398085 / 1048576.0;
t2_t0 << -10030265 / 2097152.0,  5321347 /  8388608.0,  -4791705 / 2097152.0;
p_t1  <<  -9984255 / 2097152.0,  8516467 / 16777216.0,  -7675421 / 4194304.0;
t0_t1 <<  -5031053 / 1048576.0,  7841285 / 16777216.0, -14946117 / 8388608.0;
t1_t1 <<  -4854525 / 1048576.0,  4763077 /  8388608.0,   -269775 /  131072.0;
t2_t1 <<  -9736347 / 2097152.0, 10062471 / 16777216.0,  -4408415 / 2097152.0;
```


### Example: Edge-Edge

This edge-edge query produces a false negative:
```c++
ea0_t0 << 1, 0.50803125, 2.10835646075301e-18;
ea1_t0 << -2.38233935445388e-18, 0.50803125, 1;
eb0_t0 << -4.99999999958867e-07, 0.5, 0;
eb1_t0 << -4.99999999958867e-07, 0.5, 1;
ea0_t1 << 1, 0.47124375, 4.11078309465837e-18;
ea1_t1 << -2.8526707189104e-18, 0.47124375, 1;
eb0_t1 << -4.99999999958867e-07, 0.5, 0;
eb1_t1 << -4.99999999958867e-07, 0.5, 1;
```
