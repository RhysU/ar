# Doxygen Warnings Summary

**Total Warnings: 110**

Generated on: 2026-01-12
Doxygen Version: 1.9.8

---

## Warning Categories

### 1. Broken Cross-References (9 warnings)

#### Unresolved `\ref` Commands
- **collomb2009.cpp:81** (3 occurrences): Unable to resolve reference to 'BurgAlgorithm'
- **faber1986.cpp:104** (3 occurrences): Unable to resolve reference to 'faber1986'

#### Unresolved `\link` Commands (ar.hpp)
- **ar.hpp:100**: Explicit link request to 'define' could not be resolved
- **ar.hpp:111**: Explicit link request to 'define' could not be resolved
- **ar.hpp:121**: Explicit link request to 'define' could not be resolved
- **ar.hpp:132**: Explicit link request to 'define' could not be resolved
- **ar.hpp:142**: Explicit link request to 'define' could not be resolved

---

### 2. Incorrect Parameter Documentation (4 warnings)

**ar.hpp:743** - `ar::burg_method()` function has mismatched parameter names:
- Parameter 'f' documented but not in function signature
- Parameter 'b' documented but not in function signature
- Parameter 'Ak' documented but not in function signature
- Parameter 'ac' documented but not in function signature

---

### 3. Undocumented Members (97 warnings)

#### ar-python.cpp (17 warnings)
**Macros:**
- Line 19: NPY_ALIGNED
- Line 24: NPY_ELEMENTSTRIDES
- Line 29: NPY_NOTSWAPPED
- Line 40: ARSEL_VERSION
- Line 44: DEFAULT_SUBMEAN
- Line 45: DEFAULT_ABSRHO
- Line 46: DEFAULT_CRITERION
- Line 47: DEFAULT_MINORDER
- Line 48: DEFAULT_MAXORDER
- Line 49: STRINGIFY(x)
- Line 50: STRINGIFY_HELPER(x)
- Line 413: GETSTATE(m)
- Line 437: INITERROR

**Functions:**
- Line 125: ar_arsel()
- Line 439: PyInit_ar()

#### ar6.cpp (3 warnings)
- Line 31: OptionIndex (enumeration)
- Line 34: usage[] (variable)
- Line 68: main()

#### arsel.cpp (7 warnings)
**Macros:**
- Line 26: STRINGIFY_HELPER(x)
- Line 27: STRINGIFY(x)
- Line 30: ARSEL_VERSION
- Line 34: ARSEL_CXXFLAGS

**Other:**
- Line 45: OptionIndex (enumeration)
- Line 48: usage[] (variable)
- Line 75: main()

#### example.cpp (1 warning)
- Line 23: main()

#### lorenz.cpp (6 warnings)
**Macros:**
- Line 21: STRINGIFY_HELPER(x)
- Line 22: STRINGIFY(x)

**Enumerations:**
- Line 34: OptionIndex
- Line 38: SchemeType

**Other:**
- Line 41: usage[] (variable)
- Line 228: main()

#### real.hpp (2 warnings)
- Line 15: REAL (macro)
- Line 17: real (typedef)

#### test.cpp (6 warnings)
**Macros:**
- Line 22: STRINGIFY_HELPER(x)
- Line 23: STRINGIFY(x)

**Other:**
- Line 26: OptionIndex (enumeration)
- Line 29: usage[] (variable)
- Line 45: pdiff()
- Line 52: main()

#### zohar.cpp (1 warning)
- Line 21: main()

#### ar.hpp - empirical_variance_function (4 warnings)
- Line 1630: first_argument_type (typedef)
- Line 1631: second_argument_type (typedef)
- Line 1632: result_type (typedef)
- Line 1634: operator()()

#### ar.hpp - empirical_variance_generator (3 warnings)
- Line 1660: result_type (typedef)
- Line 1662: empirical_variance_generator() (constructor)
- Line 1664: operator()()

#### ar.hpp - empirical_variance_iterator (24 warnings)
**Typedefs:**
- Line 1692: iterator_category
- Line 1693: value_type
- Line 1694: difference_type
- Line 1695: pointer
- Line 1696: reference

**Operators:**
- Line 1710: operator++()
- Line 1713: operator++(int)
- Line 1716: operator+()
- Line 1719: operator+=()
- Line 1724: operator--()
- Line 1727: operator--(int)
- Line 1730: operator-(difference_type)
- Line 1733: operator-=()
- Line 1738: operator-(iterator)
- Line 1752: operator==()
- Line 1763: operator!=()
- Line 1768: operator<()
- Line 1771: operator<=()
- Line 1774: operator>()
- Line 1777: operator>=()
- Line 1782: operator*()
- Line 1790: operator[]()

#### ar.hpp - predictor (5 warnings)
- Line 823: iterator_category (typedef)
- Line 824: value_type (typedef)
- Line 825: difference_type (typedef)
- Line 826: pointer (typedef)
- Line 827: reference (typedef)

#### ar.hpp - strided_adaptor (23 warnings)
**Typedefs:**
- Line 2655: iterator_category
- Line 2656: value_type
- Line 2657: difference_type
- Line 2658: pointer
- Line 2659: reference

**Constructors/Operators:**
- Line 2663: strided_adaptor(const strided_adaptor &o)
- Line 2665: strided_adaptor(Iterator x, difference_type n)
- Line 2667: operator++()
- Line 2674: operator++(int)
- Line 2682: operator+=()
- Line 2689: operator--()
- Line 2696: operator--(int)
- Line 2704: operator-=()
- Line 2711: operator[]()
- Line 2716: operator*()
- Line 2721: operator==()
- Line 2726: operator!=()
- Line 2731: operator<()
- Line 2736: operator-(const strided_adaptor &o)
- Line 2741: operator+(difference_type x)

---

## Summary by File

| File | Warning Count |
|------|--------------|
| ar.hpp | 64 |
| ar-python.cpp | 17 |
| arsel.cpp | 7 |
| lorenz.cpp | 6 |
| test.cpp | 6 |
| collomb2009.cpp | 3 |
| faber1986.cpp | 3 |
| ar6.cpp | 3 |
| real.hpp | 2 |
| example.cpp | 1 |
| zohar.cpp | 1 |

---

## Recommendations

### Priority 1 - Broken References
Fix the broken cross-references as these represent documentation errors:
- Resolve or remove `\ref BurgAlgorithm` references in collomb2009.cpp
- Resolve or remove `\ref faber1986` references in faber1986.cpp
- Fix the `\link define` references in ar.hpp

### Priority 2 - Parameter Mismatches
Fix the parameter documentation mismatch in ar.hpp:743 for the `burg_method()` function

### Priority 3 - Undocumented Public API
Document the public API members in ar.hpp, especially:
- Iterator classes and their standard typedefs/operators
- Public functions and classes

### Priority 4 - Example/Test Code
Consider whether example and test file members need documentation or should be excluded from documentation generation
