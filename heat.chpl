
/*
** PROGRAM: heat equation solve
**
** PURPOSE: This program will explore use of an explicit
**          finite difference method to solve the heat
**          equation under a method of manufactured solution (MMS)
**          scheme. The solution has been set to be a simple 
**          function based on exponentials and trig functions.
**
**          A finite difference scheme is used on a 1000x1000 cube.
**          A total of 0.5 units of time are simulated.
**
**          The MMS solution has been adapted from
**          G.W. Recktenwald (2011). Finite difference approximations
**          to the Heat Equation. Portland State University.
**
**
** USAGE:   Run with two arguments:
**          First is the number of cells.
**          Second is the number of timesteps.
**
**          For example, with 100x100 cells and 10 steps:
**
**          ./heat 100 10
**
**
** This Chapel program is translated from a C version originally
** written by Tom Deakin, Oct 2018
**
*/

use Time;
use Math;

param LINE = "--------------------\n"; // A line for fancy output

// Problem size, forms an nxn grid
config const n = 1000;
if n < 0 then halt("Error: n must be positive");

// Number of timesteps
config const nsteps = 10;
if nsteps < 0 then halt("Error: nsteps must be positive");

// Start the total program runtime timer
const start = timeSinceEpoch().totalSeconds();

//
// Set problem definition
//
param alpha = 0.1;          // heat equation coefficient
param length = 1000.0;      // physical size of domain: length x length square
const dx = length / (n+1);  // physical size of each cell (+1 as don't simulate boundaries as they are given)
const dt = 0.5 / nsteps;    // time interval (total time of 0.5s)

// Stability requires that dt/(dx^2) <= 0.5,
const r = alpha * dt / (dx*dx);

// Print message detailing runtime configuration
writef("\n");
writef(" MMS heat equation\n\n");
writef(LINE);
writef("Problem input\n\n");
writef(" Grid size: %i x %i\n", n, n);
writef(" Cell width: %er\n", dx);
writef(" Grid length: %dr x %dr\n", length, length);
writef("\n");
writef(" Alpha: %er\n", alpha);
writef("\n");
writef(" Steps: %i\n", nsteps);
writef(" Total time: %er\n", dt*nsteps);
writef(" Time step: %er\n", dt);
writef(LINE);

// Stability check
writef("Stability\n\n");
writef(" r value: %dr\n", r);
if r > 0.5 then
  writef(" Warning: unstable\n");
writef(LINE);

// Allocate two nxn grids
const outerDom: domain(2) = {0..#n, 0..#n};
var arr1: [outerDom] real;
var arr2: [outerDom] real;
ref u = arr1;
ref u_tmp = arr2;

// Set the initial value of the grid under the MMS scheme
initialValue(n, dx, length, u);

//
// Run through timesteps under the explicit scheme
//

// Start the solve timer
const tic = timeSinceEpoch().totalSeconds();

for 0..#nsteps {
  // Call the solve kernel
  // Computes u_tmp at the next timestep
  // given the value of u at the current timestep
  solve(n, alpha, dx, dt, u, u_tmp);

  // reference swap
  u <=> u_tmp;
}

// Stop solve timer
const toc = timeSinceEpoch().totalSeconds();

//
// Check the L2-norm of the computed solution
// against the *known* solution from the MMS scheme
//
const norm = l2Norm(n, u, nsteps, dt, alpha, dx, length);

// Stop total timer
const stop = timeSinceEpoch().totalSeconds();

// Print results
writef("Results\n\n");
writef("Error (L2norm): %er\n", norm);
writef("Solve time (s): %dr\n", toc-tic);
writef("Total time (s): %dr\n", stop-start);
writef(LINE);

// Sets the mesh to an initial value, determined by the MMS scheme
proc initialValue(const n: int, const dx: real, const length: real, ref u: [?outerDom] real) {
  var y = dx; // Physical y position
  for j in 0..#n {
    var x = dx; // Physical x position
    for i in 0..#n {
      u[i,j] = sin(Math.pi * x / length) * sin(Math.pi * y / length);
      x += dx;
    }
    y += dx;
  }
}

// Compute the next timestep, given the current timestep
proc solve(const n: int, const alpha: real, const dx: real, const dt: real, const ref u: [?outerDom] real, ref u_tmp: [outerDom] real) {
  // Finite difference constant multiplier
  const r = alpha * dt / (dx*dx);
  const r2 = 1.0 - 4.0*r;

  // Loop over the nxn grid
  for (i,j) in outerDom {
    // Update the 5-point stencil, using boundary conditions on the edges of the domain.
    // Boundaries are zero because the MMS solution is zero there.
    u_tmp[i,j] =  r2 * u[i,j] +
    r * (if i < n-1 then u[i+1,j] else 0.0) +
    r * (if i > 0   then u[i-1,j] else 0.0) +
    r * (if j < n-1 then u[i,j+1] else 0.0) +
    r * (if j > 0   then u[i,j-1] else 0.0);
  }
}

// True answer given by the manufactured solution
proc solution(const t: real, const x: real, const y: real, const alpha: real, const length: real) {
  return exp(-2.0*alpha*Math.pi*Math.pi*t/(length*length)) * sin(Math.pi*x/length) * sin(Math.pi*y/length);
}

// Computes the L2-norm of the computed grid and the MMS known solution
// The known solution is the same as the boundary function.
proc l2Norm(const n: int, const ref u: [?outerDom] real, const nsteps: int, 
  const dt: real, const alpha: real, const dx: real, const length: real) {

  // Final (real) time simulated
  const time = dt * nsteps;

  // L2-norm error
  var l2norm = 0.0;

  // Loop over the grid and compute difference of computed and known solutions as an L2-norm
  var y = dx;
  for j in 0..#n {
    var x = dx;
    for i in 0..#n {
      const answer = solution(time, x, y, alpha, length);
      const err = u[i,j] - answer;
      l2norm += err*err;

      x += dx;
    }
    y += dx;
  }

  return sqrt(l2norm);
}
