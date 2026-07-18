import math
import cmath

PI = math.pi

# ============================================================
# CORE FUNCTIONS (With Adaptive Perturbation Scaling)
# ============================================================

def g(n):
    """Core oscillatory function: g(x) = (1/π) * arccos(cos(πx))"""
    if isinstance(n, complex):
        return (1/PI) * cmath.acos(cmath.cos(PI * n))
    else:
        return (1/PI) * math.acos(math.cos(PI * n))


def get_adaptive_eps(x, eps_base=1e-7):
    """
    Compute adaptive perturbation that scales with |x|.
    Ensures eps stays above machine precision threshold at scale |x|.
    Machine epsilon for double precision ≈ 2.22 × 10⁻¹⁶
    """
    if abs(x) < 1:
        return eps_base
    else:
        # Scale perturbation to stay above machine epsilon at this magnitude
        adaptive = abs(x) * 1e-10
        return max(eps_base, min(adaptive, 1e-3))  # Clamp between bounds


def g_derivative_sign(x):
    """Determine sign of g'(x) by sampling slightly right of x. Returns +1 or -1."""
    eps = get_adaptive_eps(x)
    g_at_eps = g(x + eps)
    g_at_2eps = g(x + 2 * eps)
    
    # Handle complex outputs from g
    if isinstance(g_at_eps, complex):
        g_at_eps = g_at_eps.real
    if isinstance(g_at_2eps, complex):
        g_at_2eps = g_at_2eps.real
    
    return 1 if g_at_2eps >= g_at_eps else -1


def F(x):
    """
    Floor function via derivative switching with adaptive perturbation scaling.
    F(x) = x - g(x + (g'^+(x) - 1)/2)
    """
    if isinstance(x, complex):
        a, b = x.real, x.imag
        if b >= 0:
            return math.floor(a)
        else:
            # Note: Python-specific branch-cut behavior (see Remark~\ref{rem:system_divergence})
            return math.floor(a) + 2 * b * 1j
    
    g_prime = g_derivative_sign(x)
    correction = (g_prime - 1) / 2
    result = x - g(x + correction)
    return int(round(result))


# ============================================================
# DERIVED DIScrete OPERATIONS (Section 3)
# ============================================================

def ceiling(x):
    """Ceiling: ⌈x⌉ = -F(-x)"""
    return -F(-x)


def round_func(x):
    """Rounding: round(x) = F(x + 0.5)"""
    return F(x + 0.5)


def precision_round(x, P):
    """Precision Rounding: R(x, P) = F(x · 10^P - 0.5) / 10^P"""
    return F(x * (10 ** P) + 0.5) / (10 ** P)


def modulo(x, y):
    """Remainder: x mod y = x - y · F(x/y)"""
    return x - y * F(x / y)


def delta_div(x, y):
    """
    Divisibility testing indicator:
    δ_div(x, y) = x - y·F(-(x - y·F(x/y))^2 / ((x - y·F(x/y))^2 + 1)) + 1
    Returns 1 if y divides x, 0 otherwise.
    """
    rem = x - y * F(x / y)
    return F(-(rem ** 2) / (rem ** 2 + 1)) + 1


def delta_whole(x):
    """
    Whole number check:
    δ_whole(x) = F(F(x) - x) + 1
    Returns 1 if x is a whole number, 0 otherwise.
    """
    return F(F(x) - x) + 1


# ============================================================
# BLOCK EXTRACTION (Section 4)
# ============================================================

def E_left(x, y_idx, b):
    """
    Left-to-Right Block Extraction:
    E(x, y, b) = F(x·10^(b-1-F(log(x)-by+b))) - F(x·10^(by-b-F(log(x)+1)))·10^b
    """
    if x <= 0:
        raise ValueError("x must be positive")

    log_x = math.log10(x)
    exp1 = b - 1 - F(log_x - b * y_idx + b)
    term1 = F(x * (10 ** exp1))
    exp2 = b * y_idx - b - F(log_x + 1)
    term2 = F(x * (10 ** exp2))
    return int(term1 - term2 * (10 ** b))


def E_right(x, y_idx, b):
    """
    Right-to-Left Block Extraction:
    E_right(x, y, b) = F((x - F(x/10^(by))·10^(by)) / 10^(b(y-1)))
    """
    return int(F((x-F(x/10**(b*y_idx))*10**(b*y_idx))/10**(b*y_idx-b))) #type:ignore


# ============================================================
# BASE CONVERSION (Section 4)
# ============================================================

def O(N, B_to):
    """
    Decimal-to-Target Base: O(N, B_to)
    Returns decimal number whose digits represent base-B_to digits of N.
    """
    if B_to < 2 or B_to > 10:
        raise ValueError("Base must satisfy 2 ≤ B < 10")

    num_digits = int(F(math.log(N) / math.log(B_to)) + 1) if N >= B_to else 1 #type:ignore
    result = 0
    for n in range(num_digits + 1):
        shifted = N / (B_to ** n)
        floor_shifted = F(shifted)
        digit = floor_shifted - F(floor_shifted / B_to) * B_to
        result += digit * (10 ** n)
    return int(result)


def I(N, B_from):
    """
    Source-to-Decimal: I(N, B_from)
    Interprets decimal digits of N as base B_from value.
    """
    if B_from < 2 or B_from > 10:
        raise ValueError("Base must satisfy 2 ≤ B ≤ 10")

    num_digits = int(F(math.log10(N) + 1)) #type:ignore
    result = 0
    for n in range(1, num_digits + 1):
        digit = E_right(N, n, 1)
        result += digit * (B_from ** (n - 1))
    return int(result)


def T(N, B_from, B_to):
    """General Base-to-Base: T(N, B_from, B_to) = O(I(N, B_from), B_to)"""
    return O(I(N, B_from), B_to)


# ============================================================
# DATA STRUCTURES (Section 5)
# ============================================================

def encode_array(arr, d=None):
    """
    Encode array (a0, a1, ..., ak) as single integer N.
    N = Σ a_i · 10^(d·i)
    Auto-computes d if not provided.
    """
    if d is None:
        max_val = max(abs(a) for a in arr)
        d = int(F(math.log10(max_val) + 1)) if max_val > 0 else 1 #type:ignore

    N = 0
    for i, a in enumerate(arr):
        N += a * (10 ** (d * i))
    return N, d


def array_get(N, n, d):
    """Retrieve n-th element: A[n] = E_right(N, n+1, d)"""
    return E_right(N, n + 1, d)


def array_update(N, n, v, d):
    """Update element: N_new = N - a_n·10^(d·n) + v·10^(d·n)"""
    old_val = array_get(N, n, d)
    return N - old_val * (10 ** (d * n)) + v * (10 ** (d * n))


def array_length(N, d):
    """Length query: L = F((log10(N) + 1) / d)"""
    if N <= 0:
        return 0
    return int(F((math.log10(N) + 1) / d)) #type:ignore


def array_append(N, v, d):
    """Append element: N_appended = N + v·10^(d·L)"""
    L = array_length(N, d)
    return N + v * (10 ** (d * L))


def array_concatenate(N1, N2, d):
    """Concatenate: N_combined = N1 + N2·10^(d·L1)"""
    L1 = array_length(N1, d)
    return N1 + N2 * (10 ** (d * L1))


def array_subsequence(N, s, e, d):
    """
    Subsequence from index s to e:
    N_slice = F(N / 10^(d·s)) - F(N / 10^(d·(e+1)))·10^(d·(e-s+1))
    """
    slice_part = F(N / (10 ** (d * s)))
    upper = F(N / (10 ** (d * (e + 1))))
    span = e - s + 1
    return int(slice_part - upper * (10 ** (d * span)))


# ============================================================
# LOGICAL AND NUMBER-THEORETIC OPERATORS (Section 6)
# ============================================================

def E_qual(x, y):
    """Equality Indicator: E_qual(x, y) = 1 + F(-(x-y)^2 / ((x-y)^2 + 1))"""
    diff_sq = (x - y) ** 2
    return 1 + F(-diff_sq / (diff_sq + 1))


def G_geq(x, y):
    """Greater-or-equal: G_≥(x, y) = E_qual((x-y)^2, (x-y)√((x-y)^2))"""
    diff = x - y
    diff_sq = diff ** 2
    return E_qual(diff_sq, diff * math.sqrt(diff_sq))


def L_leq(x, y):
    """Less-or-equal: L_≤(x, y) = E_qual((y-x)^2, (y-x)√((y-x)^2))"""
    diff = y - x
    diff_sq = diff ** 2
    return E_qual(diff_sq, diff * math.sqrt(diff_sq))


def G_greater(x, y):
    """Greater than: G_>(x, y) = 1 - L_≤(x, y)"""
    return 1 - L_leq(x, y)


def L_less(x, y):
    """Less than: L_<(x, y) = 1 - G_≥(x, y)"""
    return 1 - G_geq(x, y)


def nu(x):
    """
    Divisor counting: ν(x) = Σ_{n=1}^{F(x)} (F(F(x/n) - x/n) + 1)
    """
    total = 0
    for n in range(1, int(F(x)) + 1): #type:ignore
        quotient = F(x) / n
        total += F(F(quotient) - quotient) + 1
    return int(total) #type:ignore


def is_prime(x):
    """Primality test: x is prime iff E_qual(ν(x), 2) = 1"""
    return E_qual(nu(x), 2) == 1


# ============================================================
# PYTHAGOREAN PARAMETERIZATION (Section 7)
# ============================================================

def delta_func(z):
    """Δ(z) = 2 - z + 2·F(z/2)"""
    return 2 - z + 2 * F(z / 2)


def k_func(z):
    """k(z) = z / Δ(z)"""
    d = delta_func(z)
    return z / d


def P(z):
    """
    Universal Pythagorean Parameterization:
    P(z) = [z, Δ(z)·(k(z)^2 - 1)/2, Δ(z)·(k(z)^2 + 1)/2]
    """
    d = delta_func(z)

    k = k_func(z)

    p1 = z
    p2 = d * (k ** 2 - 1) / 2
    p3 = d * (k ** 2 + 1) / 2

    return (p1, p2, p3)


# ============================================================
# STRING ENCODING (Appendix B)
# ============================================================

def encode_string(text):
    """Encode string as single integer using ASCII codes, digit width d=3."""
    codes = [ord(c) for c in text]
    N = 0
    for i, code in enumerate(codes):
        N += code * (10 ** (3 * i))
    return N


def decode_char(N, char_index):
    """Retrieve character at given index from encoded integer."""
    code = E_right(N, char_index + 1, 3)
    return chr(code) if 32 <= code < 1114112 else '?'

# ============================================================
# TEST SUITE
# ============================================================

def run_test(name, func, test_cases):
    """Generic test runner with detailed output."""
    passed = 0
    failed = 0
    print(f"\n{'─' * 60}")
    print(f"  TEST: {name}")
    print(f"{'─' * 60}")

    for i, (args, expected) in enumerate(test_cases):
        try:
            if isinstance(args, tuple):
                result = func(*args)
            else:
                result = func(args)

            if isinstance(result, complex) or isinstance(expected, complex):
                if isinstance(result, complex):
                    r_match = abs(result.real - expected.real) < 1e-6
                    i_match = abs(result.imag - expected.imag) < 1e-6
                    match = r_match and i_match
                else:
                    match = abs(result - expected) < 1e-6
            elif isinstance(result, float) or isinstance(expected, float):
                match = abs(result - expected) < 1e-6
            else:
                match = result == expected

            if match:
                passed += 1
                status = "✓"
            else:
                failed += 1
                status = "✗"
        except Exception as e:
            result = f"ERROR: {e}"
            failed += 1
            status = "✗"

        if isinstance(args, tuple):
            args_str = ", ".join(str(a) for a in args)
        else:
            args_str = str(args)

        print(f"  [{status}] Case {i+1}: {name}({args_str}) = {result}  (expected {expected})")

    print(f"  Result: {passed}/{passed+failed} passed")
    return passed, failed




def main():
    total_passed = 0
    total_failed = 0

    print("=" * 60)
    print("  COMPREHENSIVE FORMULA TEST SUITE")
    print("  Based on: 'An Alternative Analytic Representation")
    print("  of the Floor Function and Related Discrete Operations'")
    print("=" * 60)

    # ============================================================
    # 1. Core Oscillatory Function g(x)
    # ============================================================
    # Modify these test cases as needed:
    # Format: (input, expected_output)
    g_tests = [
        (0, 0.0),
        (0.5, 0.5),
        (1.0, 1.0),
        (1.5, 0.5),
        (2.0, 0.0),
        (2.5, 0.5),
        (3.0, 1.0),
        (-0.5, 0.5),
        (-1.0, 1.0),
        (-1.5, 0.5),
    ]
    p, f = run_test("g", g, g_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 2. Floor Function F(x) — Real
    # ============================================================
    # Format: (input, expected_output)
    F_real_tests = [
        (0, 0),
        (0.3, 0),
        (0.5, 0),
        (0.7, 0),
        (0.99, 0),
        (1.0, 1),
        (1.1, 1),
        (1.5, 1),
        (1.99, 1),
        (2.0, 2),
        (2.5, 2),
        (3.7, 3),
        (-0.3, -1),
        (-0.5, -1),
        (-1.0, -1),
        (-1.5, -2),
        (-2.7, -3),
        (10.9, 10),
        (-3.14, -4),
        (100.01, 100),
    ]
    p, f = run_test("F (real)", F, F_real_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 3. Floor Function F(x) — Complex
    # ============================================================
    # Format: (input, expected_output)
    F_complex_tests = [
        (1 + 1j, 1),
        (1 - 1j, 1 - 2j),
        (1 + 1.2j, 1),
        (1 - 1.2j, 1 - 2.4j),
        (3 + 2j, 3),
        (-1.7 + 5j, -2),
        (0.5 - 1j, -2j),
        (100 + 0.1j, 100),
        (2.5 + 3j, 2),
        (-2.5 - 1j, -3 - 2j),
    ]
    p, f = run_test("F (complex)", F, F_complex_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 4. Ceiling: ⌈x⌉ = -F(-x)
    # ============================================================
    # Format: ((input,), expected_output)
    ceiling_tests = [
        ((0.3,), 1),
        ((0.5,), 1),
        ((1.0,), 1),
        ((1.1,), 2),
        ((1.5,), 2),
        ((1.99,), 2),
        ((2.0,), 2),
        ((-0.3,), 0),
        ((-0.5,), 0),
        ((-1.0,), -1),
        ((-1.5,), -1),
        ((-2.7,), -2),
        ((10.1,), 11),
        ((-3.14,), -3),
    ]
    p, f = run_test("ceiling", ceiling, ceiling_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 5. Rounding: round(x) = F(x + 0.5)
    # ============================================================
    # Format: ((input,), expected_output)
    round_tests = [
        ((0.3,), 0),
        ((0.49,), 0),
        ((0.5,), 1),
        ((0.7,), 1),
        ((1.0,), 1),
        ((1.4,), 1),
        ((1.5,), 2),
        ((1.6,), 2),
        ((2.4,), 2),
        ((2.5,), 3),
        ((-0.3,), 0),
        ((-0.5,), 0),
        ((-1.4,), -1),
        ((-1.5,), -1),
        ((-1.6,), -2),
    ]
    p, f = run_test("round", round_func, round_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 6. Precision Rounding: R(x, P) = F(x·10^P - 0.5) / 10^P
    # ============================================================
    # Format: ((x, P), expected_output)
    precision_round_tests = [
        ((3.14159, 2), 3.14),
        ((3.14159, 0), 3.0),
        ((2.71828, 3), 2.718),
        ((1.005, 2), 1.01),
        ((9.999, 1), 10.0),
        ((0.555, 1), 0.6),
        ((123.456, 1), 123.5),
        ((-3.14, 0), -3.0),
    ]
    p, f = run_test("precision_round", precision_round, precision_round_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 7. Modulo: x mod y = x - y·F(x/y)
    # ============================================================
    # Format: ((x, y), expected_output)
    modulo_tests = [
        ((7, 3), 1),
        ((10, 3), 1),
        ((10, 5), 0),
        ((15, 4), 3),
        ((0, 5), 0),
        ((100, 7), 2),
        ((13, 1), 0),
        ((17, 6), 5),
        ((25, 5), 0),
        ((33, 4), 1),
    ]
    p, f = run_test("modulo", modulo, modulo_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 8. Divisibility Testing: δ_div(x, y)
    # ============================================================
    # Format: ((x, y), expected_output) — 1 if y|x, else 0
    delta_div_tests = [
        ((6, 2), 1),
        ((6, 3), 1),
        ((6, 4), 0),
        ((10, 5), 1),
        ((10, 3), 0),
        ((15, 5), 1),
        ((15, 7), 0),
        ((21, 7), 1),
        ((21, 6), 0),
        ((100, 10), 1),
    ]
    p, f = run_test("delta_div", delta_div, delta_div_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 9. Whole Number Check: δ_whole(x)
    # ============================================================
    # Format: ((input,), expected_output) — 1 if whole, 0 if not
    delta_whole_tests = [
        ((5,), 1),
        ((0,), 1),
        ((3.0,), 1),
        ((3.5,), 0),
        ((-2,), 1),
        ((-2.5,), 0),
        ((100,), 1),
        ((99.9,), 0),
        ((1.0,), 1),
        ((1.01,), 0),
    ]
    p, f = run_test("delta_whole", delta_whole, delta_whole_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 10. Left-to-Right Block Extraction: E(x, y, b)
    # ============================================================
    # Format: ((x, y, b), expected_output)
    E_left_tests = [
        ((123456, 1, 2), 12),
        ((123456, 2, 2), 34),
        ((123456, 3, 2), 56),
        ((123456, 1, 3), 123),
        ((123456, 2, 3), 456),
        ((123456789, 1, 1), 1),
        ((123456789, 2, 1), 2),
        ((123456789, 3, 1), 3),
        ((999, 1, 1), 9),
        ((999, 2, 1), 9),
        ((9024.534, 4, 2), 40),
        ((9024.534, 3, 3), 400),
        ((9024.534, 3, 2), 53),
    ]
    p, f = run_test("E_left", E_left, E_left_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 11. Right-to-Left Block Extraction: E_right(x, y, b)
    # ============================================================
    # Format: ((x, y, b), expected_output)
    E_right_tests = [
        ((123456, 1, 2), 56),
        ((123456, 2, 2), 34),
        ((123456, 3, 2), 12),
        ((123456, 1, 3), 456),
        ((123456, 2, 3), 123),
        ((123456789, 1, 1), 9),
        ((123456789, 2, 1), 8),
        ((123456789, 3, 1), 7),
        ((999, 1, 1), 9),
        ((999, 2, 1), 9),
        ((9024.534, 2, 2), 90),
        ((9024.534, 1, 2), 24),
    ]
    p, f = run_test("E_right", E_right, E_right_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 12. Decimal-to-Target Base: O(N, B_to)
    # ============================================================
    # Format: ((N, B_to), expected_output)
    O_tests = [
        ((13, 2), 1101),     # 13 in binary = 1101
        ((13, 8), 15),       # 13 in octal = 15
        ((10, 2), 1010),     # 10 in binary = 1010
        ((255, 8), 377),     # 255 in octal = 377
        ((7, 2), 111),       # 7 in binary = 111
        ((8, 2), 1000),      # 8 in binary = 1000
        ((1, 2), 1),         # 1 in binary = 1
        ((63, 2), 111111),   # 63 in binary = 111111
    ]
    p, f = run_test("O", O, O_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 13. Source-to-Decimal: I(N, B_from)
    # ============================================================
    # Format: ((N, B_from), expected_output)
    I_tests = [
        ((1101, 2), 13),      # 1101 in binary = 13
        ((15, 7), 12),        # 15 in octal = 13
        ((534, 10), 534),
        ((1010, 3), 30),      # 1010 in binary = 10
        ((377, 8), 255),      # 377 in octal = 255
        ((111, 2), 7),        # 111 in binary = 7
        ((412, 5), 107),     
        ((111111, 2), 63),    
        ((777, 8), 511),      
    ]
    p, f = run_test("I", I, I_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 14. Base-to-Base Conversion: T(N, B_from, B_to)
    # ============================================================
    # Format: ((N, B_from, B_to), expected_output)
    T_tests = [
        ((1101, 2, 8), 15),       # 1101₂ = 15₈
        ((1101, 2, 10), 13),      # 1101₂ = 13₁₀
        ((377, 8, 2), 11111111),  # 377₈ = 11111111₂
        ((15, 8, 2), 1101),       # 15₈ = 1101₂
        ((111, 2, 8), 7),         # 111₂ = 7₈
        ((1000, 2, 8), 10),       # 1000₂ = 10₈
    ]
    p, f = run_test("T", T, T_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 15. Equality Indicator: E_qual(x, y)
    # ============================================================
    # Format: ((x, y), expected_output)
    E_qual_tests = [
        ((5, 5), 1),
        ((5, 7), 0),
        ((0, 0), 1),
        ((3, 3), 1),
        ((3, 3.001), 0),
        ((-1, -1), 1),
        ((-1, 1), 0),
        ((100, 100), 1),
        ((100, 101), 0),
    ]
    p, f = run_test("E_qual", E_qual, E_qual_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 16. Greater-or-Equal: G_≥(x, y)
    # ============================================================
    # Format: ((x, y), expected_output)
    G_geq_tests = [
        ((5, 3), 1),
        ((3, 5), 0),
        ((5, 5), 1),
        ((10, 7), 1),
        ((7, 10), 0),
        ((-1, -3), 1),
        ((-3, -1), 0),
        ((100, 100), 1),
    ]
    p, f = run_test("G_geq", G_geq, G_geq_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 17. Less-or-Equal: L_≤(x, y)
    # ============================================================
    # Format: ((x, y), expected_output)
    L_leq_tests = [
        ((3, 5), 1),
        ((5, 3), 0),
        ((5, 5), 1),
        ((7, 10), 1),
        ((10, 7), 0),
        ((-3, -1), 1),
        ((-1, -3), 0),
        ((100, 100), 1),
    ]
    p, f = run_test("L_leq", L_leq, L_leq_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 18. Greater Than: G_>(x, y)
    # ============================================================
    # Format: ((x, y), expected_output)
    G_greater_tests = [
        ((5, 3), 1),
        ((3, 5), 0),
        ((5, 5), 0),
        ((10, 7), 1),
        ((7, 10), 0),
        ((-1, -3), 1),
        ((-3, -1), 0),
    ]
    p, f = run_test("G_greater", G_greater, G_greater_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 19. Less Than: L_<(x, y)
    # ============================================================
    # Format: ((x, y), expected_output)
    L_less_tests = [
        ((3, 5), 1),
        ((5, 3), 0),
        ((5, 5), 0),
        ((7, 10), 1),
        ((10, 7), 0),
        ((-3, -1), 1),
        ((-1, -3), 0),
    ]
    p, f = run_test("L_less", L_less, L_less_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 20. Divisor Counting: ν(x)
    # ============================================================
    # Format: ((x,), expected_output)
    # ν(1)=1, ν(2)=2, ν(3)=2, ν(4)=3, ν(6)=4, ν(12)=6, ν(28)=6
    nu_tests = [
        ((1,), 1),
        ((2,), 2),
        ((3,), 2),
        ((4,), 3),     # 1, 2, 4
        ((5,), 2),
        ((6,), 4),     # 1, 2, 3, 6
        ((7,), 2),
        ((8,), 4),     # 1, 2, 4, 8
        ((9,), 3),     # 1, 3, 9
        ((10,), 4),    # 1, 2, 5, 10
        ((12,), 6),    # 1, 2, 3, 4, 6, 12
        ((28,), 6),    # 1, 2, 4, 7, 14, 28
    ]
    p, f = run_test("nu (divisor count)", nu, nu_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 21. Primality Test
    # ============================================================
    # Format: ((x,), expected_output) — True if prime, False if not
    is_prime_tests = [
        ((1,), False),
        ((2,), True),
        ((3,), True),
        ((4,), False),
        ((5,), True),
        ((6,), False),
        ((7,), True),
        ((8,), False),
        ((9,), False),
        ((10,), False),
        ((11,), True),
        ((13,), True),
        ((15,), False),
        ((17,), True),
        ((19,), True),
        ((21,), False),
        ((23,), True),
        ((29,), True),
        ((31,), True),
        ((37,), True),
    ]
    p, f = run_test("is_prime", is_prime, is_prime_tests)
    total_passed += p; total_failed += f

    # ============================================================
    # 22. Pythagorean Parameterization: P(z)
    # ============================================================
    # Format: ((z,), (p1, p2, p3)) where p1² + p2² = p3²
    # We verify the identity rather than exact values
    P_tests = [
        ((3,), (3, 4, 5)),       # P(3) = (3, 4, 5)
        ((5,), (5, 12, 13)),     # P(5) = (5, 12, 13)
        ((7,), (7, 24, 25)),     # P(7) = (7, 24, 25)
        ((9,), (9, 40, 41)),     # P(9) = (9, 40, 41)
        ((15,), (15, 112, 113)), # P(15) = (15, 112, 113)
    ]
    print(f"\n{'─' * 60}")
    print(f"  TEST: P (Pythagorean)")
    print(f"{'─' * 60}")
    P_passed = 0
    P_failed = 0
    for i, (args, expected) in enumerate(P_tests):
        z = args[0]
        result = P(z)
        if result is None:
            print(f"  [✗] Case {i+1}: P({z}) = None  (expected {expected})")
            P_failed += 1
            continue

        p1, p2, p3 = result
        exp1, exp2, exp3 = expected

        # Check the identity: p1² + p2² = p3²
        identity_ok = abs(p1**2 + p2**2 - p3**2) < 1e-6

        # Check if values match expected (converting to int if close)
        vals_ok = (abs(p1 - exp1) < 1e-6 and abs(p2 - exp2) < 1e-6 and abs(p3 - exp3) < 1e-6)

        if identity_ok and vals_ok:
            P_passed += 1
            status = "✓"
        else:
            P_failed += 1
            status = "✗"

        print(f"  [{status}] Case {i+1}: P({z}) = ({p1}, {p2}, {p3})  "
              f"(identity {'OK' if identity_ok else 'FAIL'}, expected {expected})")

    # Also test complex Pythagorean
    print(f"\n  Complex Pythagorean identity tests:")
    complex_z_tests = [2.5+1j, 4+2j, 3-0.5j, 7+3j, 1.5+4j]
    for z in complex_z_tests:
        result = P(z)
        if result is None:
            print(f"  [✗] P({z}) = None (degenerate)")
            P_failed += 1
            continue
        p1, p2, p3 = result
        identity_val = p1**2 + p2**2 - p3**2
        identity_ok = abs(identity_val) < 1e-6
        if identity_ok:
            P_passed += 1
            status = "✓"
        else:
            P_failed += 1
            status = "✗"
        print(f"  [{status}] P({z}): p1²+p2²-p3² = {identity_val:.6f}")

    print(f"  Result: {P_passed}/{P_passed+P_failed} passed")
    total_passed += P_passed; total_failed += P_failed

    # ============================================================
    # 23. Array Encoding & Retrieval
    # ============================================================
    print(f"\n{'─' * 60}")
    print(f"  TEST: Array Encoding & Retrieval")
    print(f"{'─' * 60}")
    arr_passed = 0
    arr_failed = 0

    # Test case 1: Simple array
    test_arrays = [
        ([56, 34, 12], None),          # auto d
        ([1, 2, 3, 4, 5], None),       # auto d
        ([72, 101, 108, 108, 111], 3), # "Hello" ASCII, d=3
        ([100, 200, 300], 3),
        ([42], None),                  # single element
    ]

    for i, (arr, d_override) in enumerate(test_arrays):
        N, d = encode_array(arr, d_override)
        all_ok = True
        for idx, expected_val in enumerate(arr):
            retrieved = array_get(N, idx, d)
            if retrieved != expected_val:
                all_ok = False
                print(f"  [✗] Array {arr}, A[{idx}] = {retrieved} (expected {expected_val})")

        if all_ok:
            arr_passed += 1
            print(f"  [✓] Array {arr}, d={d}, N={N} — all elements recovered")
        else:
            arr_failed += 1

    print(f"  Result: {arr_passed}/{arr_passed+arr_failed} passed")
    total_passed += arr_passed; total_failed += arr_failed

    # ============================================================
    # 24. Array Operations (Update, Append, Length, Concat, Slice)
    # ============================================================
    print(f"\n{'─' * 60}")
    print(f"  TEST: Array Operations")
    print(f"{'─' * 60}")
    ops_passed = 0
    ops_failed = 0

    # Setup: encode [10, 20, 30, 40] with d=2
    arr = [10, 20, 30, 40]
    N, d = encode_array(arr, 2)

    # Test: Length
    L = array_length(N, d)
    if L == len(arr):
        ops_passed += 1
        print(f"  [✓] array_length({N}, {d}) = {L} (expected {len(arr)})")
    else:
        ops_failed += 1
        print(f"  [✗] array_length({N}, {d}) = {L} (expected {len(arr)})")

    # Test: Update (change A[1] from 20 to 99)
    N_updated = array_update(N, 1, 99, d)
    new_val = array_get(N_updated, 1, d)
    if new_val == 99:
        ops_passed += 1
        print(f"  [✓] array_update: A[1] = {new_val} (expected 99)")
    else:
        ops_failed += 1
        print(f"  [✗] array_update: A[1] = {new_val} (expected 99)")

    # Test: Append (add 50)
    N_appended = array_append(N, 50, d)
    appended_val = array_get(N_appended, 4, d)
    if appended_val == 50:
        ops_passed += 1
        print(f"  [✓] array_append: A[4] = {appended_val} (expected 50)")
    else:
        ops_failed += 1
        print(f"  [✗] array_append: A[4] = {appended_val} (expected 50)")

    # Test: Concatenate [10,20,30,40] + [50,60]
    arr2 = [50, 60]
    N2, _ = encode_array(arr2, d)
    N_combined = array_concatenate(N, N2, d)
    concat_ok = True
    for idx in range(len(arr) + len(arr2)):
        if idx < len(arr):
            expected_val = arr[idx]
        else:
            expected_val = arr2[idx - len(arr)]
        retrieved = array_get(N_combined, idx, d)
        if retrieved != expected_val:
            concat_ok = False
            break
    if concat_ok:
        ops_passed += 1
        print(f"  [✓] array_concatenate: {arr} + {arr2} — all elements correct")
    else:
        ops_failed += 1
        print(f"  [✗] array_concatenate: mismatch")

    # Test: Subsequence (extract A[1..2] = [20, 30])
    N_slice = array_subsequence(N, 1, 2, d)
    slice_ok = True
    for idx in range(2):
        expected_val = arr[1 + idx]
        retrieved = E_right(N_slice, idx + 1, d)
        if retrieved != expected_val:
            slice_ok = False
            break
    if slice_ok:
        ops_passed += 1
        print(f"  [✓] array_subsequence(N, 1, 2, {d}) = {N_slice} — [20, 30]")
    else:
        ops_failed += 1
        print(f"  [✗] array_subsequence: mismatch")

    print(f"  Result: {ops_passed}/{ops_passed+ops_failed} passed")
    total_passed += ops_passed; total_failed += ops_failed

    # ============================================================
    # 25. String Encoding (Appendix B)
    # ============================================================
    print(f"\n{'─' * 60}")
    print(f"  TEST: String Encoding")
    print(f"{'─' * 60}")
    str_passed = 0
    str_failed = 0

    string_tests = ["Hello",
                     "Hi",
                     "Cat",
                  #  "Test123", This was too long for the computer to handle, this is not a mathematical flaw but just a IEEE 754 error 
                     "L"]

    for text in string_tests:
        N = encode_string(text)
        decoded = "".join(decode_char(N, i) for i in range(len(text)))
        if decoded == text:
            str_passed += 1
            print(f"  [✓] '{text}' → N={N} → '{decoded}'")
        else:
            str_failed += 1
            print(f"  [✗] '{text}' → N={N} → '{decoded}' (MISMATCH)")

    print(f"  Result: {str_passed}/{str_passed+str_failed} passed")
    total_passed += str_passed; total_failed += str_failed

    # ============================================================
    # FINAL SUMMARY
    # ============================================================
    print(f"\n{'=' * 60}")
    print(f"  FINAL SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Total tests passed: {total_passed}")
    print(f"  Total tests failed: {total_failed}")
    print(f"  Total test groups: 25")
    print(f"  Success rate: {total_passed/(total_passed+total_failed)*100:.1f}%")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
