# Device-Agnostic Programming with Chapel: Exercises

## Getting Started

We will use NCI's Intel-based Gadi system for all exercises. We will use both the Cascade-Lake-based `normal` nodes and the GPU-accelerated `gpuvolta` nodes.

To set up your environment for Chapel development, run the following command:

```
source /scratch/vp91/chapel-2.1/setup.bash
```

If you use Visual Studio Code as your editor, you may wish to install the [Chapel Language Extension for VS Code](https://marketplace.visualstudio.com/items?itemName=chpl-hpe.chapel-vscode).


## Numerical Computation of PI

The file [pi.chpl](pi.chpl) contains a sequential code that numerically computes the integral of $`4/(1+x*x)`$ over the interval $`[0..1)`$, which should equal $\pi$. Review the code so you understand what it is doing. Note the line that controls the number of integration steps:

```chapel
config const num_steps = 100000000;
```

Build a CPU-only executable using `make pi_cpu`. Run the executable on the login nodes with small numbers of steps to see how increasing the number of steps improves the accuracy of integration, e.g.

```
./pi_cpu -snum_steps 4
```

### Parallel Pi

As provided, the program computes $\pi$ in a sequential `for` loop. Modify the code so that it uses Chapel's features for [data-parallelism](https://chapel-lang.org/docs/language/spec/data-parallelism.html) to compute the integral.

On the login nodes, you can test your changes using small numbers of threads by changing the number of [worker threads](https://chapel-lang.org/docs/usingchapel/tasks.html#controlling-the-number-of-threads) that the Chapel runtime creates. For example:

```
CHPL_RT_NUM_THREADS_PER_LOCALE=4 ./pi_cpu
```

Now run the CPU-only version on a Gadi Cascade Lake compute node using the provided jobscript:

```
qsub job_pi_cpu.sh
```

### GPU-Accelerated PI

The [`CHPL_LOCALE_MODEL` environment variable](https://chapel-lang.org/docs/usingchapel/chplenv.html#readme-chplenv-chpl-locale-model) determines whether to compile for GPU, or CPU only. You can check the value of this environment variable in Chapel code using the [`ChplConfig` module](https://chapel-lang.org/docs/modules/standard/ChplConfig.html). For example, the following code sets the value of the `targetLoc` locale to be the first GPU sub-locale if compiling for GPUs; otherwise, it sets the value of `targetLoc` to be the current (CPU) locale.

```chapel
use ChplConfig;
const targetLoc = if CHPL_LOCALE_MODEL == "gpu" then here.gpus[0] else here;
```

Modify `pi.chpl` so that it works on either CPU or GPU, depending on how it is compiled.

Build the GPU version using `make pi_gpu`. What happens if you run it on the (CPU-only) login node?

Run the GPU version on a Gadi GPU Volta compute node using the provided jobscript:

```
qsub job_pi_gpu.sh`
```

### Diagnostics and Profiling

You may wonder: how does the Chapel code translate into kernel launches and data movement? Chapel provides a variety of [diagnostic utilities](https://chapel-lang.org/docs/technotes/gpu.html#diagnostics-and-utilities) to help count and trace kernel launches, data movement, and memory allocations - try adding these diagnostics to `pi.chpl`.

How does performance compare with the CPU version? What factors might be contributing to the relative performance of each version? You may wish to conduct [GPU profiling using `nvprof` or Visual Profiler](https://docs.nvidia.com/cuda/profiler-users-guide/index.html) to better understand the performance of the GPU code.

## Heat Equation Solver

The file [heat.chpl](heat.chpl) contains a sequential code that numerically solves the 2D heat equation using an explicit finite difference discretization.

### Parallel Heat

Modify `heat.chpl` to parallelize the solver as much as possible, making sure that correctness (as measured by `Error (L2norm)`) is maintained.

Once you are happy with your parallel solver, consider also parallelizing the initialization and solution check code.

Run your parallel solver using the provided jobscript:

```
qsub job_heat_cpu.sh
```

### GPU-Accelerated Heat

Modify `heat.chpl` so that it works on either CPU or GPU, depending on how it is compiled.

Run your GPU solver using the provided jobscript:

```
qsub job_heat_gpu.sh
```

How does the performance compare to the CPU version? Can you use Chapel GPU diagnostics or profiling (e.g. `nvprof`) to understand and improve the performance of your code?

### Inspecting the Generated Code

If you are comfortable with reading [PTX code](https://docs.nvidia.com/cuda/parallel-thread-execution/index.html), you can inspect the PTX that the Chapel compiler has created from your data-parallel loops. Add the compile option `--savec tmp` to the `CHPL_FLAGS` variable in the Makefile, to instruct the Chapel compiler to save all intermediate generated binary code to the directory `tmp`. You should find the generated PTX in the file `tmp/chpl__gpu.s`.

In the PTX, each generated kernel is named after the file and line number of the Chapel code that it is generated from. For example, if your heat file contains a `forall` data-parallel loop on line 147, then the PTX should contain a generated kernel starting with the following line:

```
	// .globl	chpl_gpu_kernel_heat_line_147_ // -- Begin function chpl_gpu_kernel_heat_line_147_
```
