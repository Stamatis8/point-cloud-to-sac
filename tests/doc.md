# Tests

## test-1.cpp

Verify validity of `SectionalAreaXwiseYsymmetrical()` by comparing it with the analytic sectional area curve of the Wigley hull

## test-2.cpp

Test first derivative of sectional are curve generated through `SectionalAreaXwiseYsymmetrical()` at boundary points, by comparing it with the analytic version of the Wigley hull

## test-3.cpp

Evaluate cross sectional area curve for KCSsim hull (see `KCSsimModeler.hpp`) using `SectionalAreaXwiseYsymmetrical()`. Test for validity by integating and comparing with model volume