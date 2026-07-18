# Floor Function Paper - An Alternative Analytic Representation

## Overview

This repository contains the implementation code and test suite for the paper **"An Alternative Analytic Representation of the Floor Function and Related Discrete Operations"** by Mohammad Rakibul Islam.

The core contribution is an alternative closed-form expression for the floor function using elementary trigonometric functions, derivatives, and arithmetic operations:

$$F(x) = x - g\!\left(x + \frac{g'^{+}(x) - 1}{2}\right)$$

where $g(x) = \pi^{-1}\arccos(\cos(\pi x))$ and $g'^{+}(x)$ is computed via numerical sampling.

## Features

### Core Functions
- **Floor Function**: Alternative analytic representation avoiding piecewise definitions
- **Derived Operators**: Ceiling, rounding, precision rounding, modulo
- **Divisibility Testing**: Indicator functions for divisibility and whole number checks
- **Block Extraction**: Both left-to-right and right-to-left variants (anchored at decimal point)
- **Base Conversion**: Decimal-to-target and target-to-decimal for bases 2–9
- **Array Encoding**: Single integer encoding for sequences (Gödel-style)
- **Logical Predicates**: Equality and order comparisons (≥, ≤, >, <)
- **Number Theory**: Divisor counting and primality testing
- **Pythagorean Parameterization**: Universal parameterization valid for complex arguments

### Test Suite
- **250/251 test cases passed (99.6%)**
- Comprehensive coverage of all formulas
- Edge case handling documented
- Full verification log available

## Installation

### Requirements
- Python 3.8 or higher
- No external dependencies (uses only standard library `math` and `cmath`)

### Setup
bash

Clone the repository
git clone https://github.com/rakibulofficial1017/floor-function-paper.git

Navigate to project directory
cd floor-function-paper

Usage
Running Tests
python test.py
This runs the comprehensive test suite covering all implemented formulas. Expected output:

Total tests passed: 250
Total tests failed: 1
Success rate: 99.6%
Viewing Verification Log
# View the complete verification results
cat verification.log

# Or run tests with verbose output
python test.py --verbose 2>&1 | tee verification.log
The verification.log file contains detailed output from all 251 test cases, including:

Individual test pass/fail status
Actual vs. expected values
Error messages for failed cases
Timing information (if enabled)
Individual Function Usage
from test import F, E_right, encode_array, is_prime, O

# Floor function
print(F(3.7))           # Output: 3

# Block extraction (right-to-left)
print(E_right(123456, 1, 2))  # Output: 56

# Array encoding
N, d = encode_array([56, 34, 12], d=2)
print(N)                 # Output: 123456

# Primality test
print(is_prime(17))      # Output: True

# Base conversion
print(O(13, 2))          # Output: 1101 (binary digits as decimal)
File Structure
floor-function-paper/
├── test.py                        # Main implementation and test suite
├── verification.log               # Complete test execution log (251 cases)
├── README.md                      # This file
├── LICENSE                        # MIT License
└── floor-function-paper.tex       # LaTeX manuscript
Known Limitations
IEEE 754 Precision Constraints
The base-10 array encoding scheme encounters precision issues for encoded integers exceeding ~10¹⁵ (15–16 significant digits). For example:

'Hello' (15 digits) ✓ Works correctly
'Test123' (20 digits) ✗ Precision loss
This is a numerical implementation constraint, not a formulaic flaw. For larger datasets:

Use arbitrary-precision libraries (Python's int, GMP, mpmath)
Split large arrays into chunks below the precision threshold
Implement symbolic computation where possible
See Section 5 of the paper for detailed discussion.

Complex Domain Behavior
There is system-dependent behavior in complex evaluation:

Python (cmath): Shows residue term −2bi in lower half-plane
Wolfram Alpha: Returns ⌊a⌋ for all complex inputs
This discrepancy requires further rigorous analysis. See Section 2.4 of the paper.

Adaptive Perturbation
For large numbers (> 10¹¹), perturbation scaling must be adaptive:

eps = max(1e-7, abs(x) * 1e-10)  # Ensures above machine epsilon
Academic Reference
Paper Citation
@misc{islam2026floormathrepresentation,
  author = {Islam, Mohammad Rakibul},
  title = {An Alternative Analytic Representation of the Floor Function and Related Discrete Operations},
  year = {2026},
  note = {High school research project},
  url = {https://github.com/rakibulofficial1017/floor-function-paper}
}
Related Work
Gödel, K. (1931). Über formal unentscheidbare Sätze der Principia Mathematica. Monatshefte für Mathematik und Physik, 38, 173–198.
Stack Exchange. Formula for the floor function. https://math.stackexchange.com/questions/2807610/formula-for-the-floor-function
Verification Log
The verification.log file provides reproducible evidence of all test results. It documents:

250 passed test cases
1 expected failure (precision limit on large strings)
Exact values for debugging and replication
License
This project is licensed under the MIT License. See the LICENSE file for details.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software.

Contributing & Feedback
This work was completed independently by the author. Feedback and discussion are welcome through GitHub Issues or the author's contact information below.

Note: All mathematical content, formulas, proofs, and theoretical developments are the original work of the author. AI tools were used for language editing and formatting assistance during manuscript preparation.

Contact
Email: rakibulofficial1017@gmail.com
GitHub: @rakibulofficial1017
Author: Mohammad Rakibul Islam
Date: July 2026
Classification: High School Research Project
License: MIT

"The floor function is traditionally defined through piecewise cases... This paper presents an alternative closed-form expression"
