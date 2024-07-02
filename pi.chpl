/*
This program will numerically compute the integral of

                  4/(1+x*x)

from 0 to 1.  The value of this integral is pi -- which
is great since it gives us an easy way to check the answer.

This Chapel program is translated from a C version originally
written by Tim Mattson, 11/99.
*/

use Time;

config const num_steps = 100000000;

const start_time = timeSinceEpoch().totalSeconds();

const step: real = 1.0 / num_steps;
var sum: real = 0.0;
for i in 1..#num_steps {
    const x = (i - 0.5) * step;
    sum = sum + 4.0 / (1.0 + x * x);
}

const pi = step * sum;
const run_time = timeSinceEpoch().totalSeconds() - start_time;
writef("pi with %i steps is %dr in %dr seconds\n ", num_steps, pi, run_time);