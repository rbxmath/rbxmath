# LuaMath

A general purpose Lua math library (with a heavy emphasis on polynomial algebra!).

# Documentation

## SymbolicAlgebra

symbolicAlgebra is a library that adds symbolic polynomials to Lua. What does it mean for something to be symbolic? I am glad you asked. In normal lua, if you create a variable `local x = a` you are creating a location in memory "called" `x`, and then storing the *current* value of `a` within it. Thus, unless `a` is a table, if you change the value of `a`, you will not change the value of `x`. This is not how variables work in mathematics. In mathematics, if I have a variable `x = a`, then `x` ***is*** `a`. If `a` changes, `x` changes and visa versa. Symbolic variables are somewhat different from both of these, but they are much closer to the math version. 

Suppose I write the foilow piece of code

```lua
local poly = require(symbolicAlgebra).polynomial

local lin = poly.new({1,1}) --lin is now the polynomial 1 + x
local quad = lin * lin --quad is now the polynomial 1 + 2x + x^2
```

As mentioned in the comments, `lin` is the polynomial `1 + x` and `quad` is the polynomial `1 + 2x + x^2`. Now, if I were to change `lin`, `quad` would remain unchanged, so the symbolic variable is not quite a mathematical variable, but, we can do interesting things with this. Suppose I were to define a new variable `local liny = poly.replace(lin, "x", "y")`. What would `lin * liny` be then? Well, in algebra, we would say that `(1 + x) * (1 + y) = 1 + x + y + xy`, and using symbolic polynomials, we get this result! Namely

```lua
local poly = require(symbolicAlgebra).polynomial

local lin = poly.new({1,1}) --lin is now the polynomial 1 + x
local quad = lin * lin --quad is now the polynomial 1 + 2x + x^2
local quadish = lin * liny --quadish is the polynomial 1 + x + y + xy
```

Thus, we can actually do real, mathematical algebra using these symbolic polynomials. Currently, only polynomial addition, subtraction, and multplication are implemented (division is tricky because you can't just divide one polynomial by another polynomial in general) so you can't do things like factor, but this will be implemented shortly (perhaps by April 10th, 2022).

## symbolicAlgebra.polynomial

Polynomials are the fundamental object within the symbolic algebra library. "Under the hood," polynomials are arrays of monomials, which will be introduced latter. They are always sorted in ascending order (the monomial ordering is graded reverse lexicographic though this might change to graded lexicographic because `table.sort` sorts ascending.) Polynomials can be multivariable or single variable.

The current functionality is as follows:

`+, -, *`

Polynomials can be multiplied, added and subtracted with each other.

`symbolicAlgebra.polynomial.new(array)`:

Creates a new polynomial based on `array`. If `array` is an array of numbers, this function creates a polynomial with symbol `x` and coefficients given by `array`. Otherwise, if `array` is an array of monomials, then the function returns the polynomial obtained by "adding" all of the monomials in the array.

`symbolicAlgebra.polynomial.copy(polynomial)`

Creates a shallow copy of `polynomial`.

`symbolicAlgebra.polynomial.eval(polynomial, rules)`

Evaluates the polynomial at some list of values. `rules` should be of the form `{x = a, y = b, ...}`. For instance, the polynomial `1 + x` evaluated at `1` is `2`.

`symbolicAlgebra.polynomial.replace(polynomial, symbolToBeReaplaced, symbol)`

Returns a polynomial where `symbolToBeReplaced` has been replaced by `symbol`. Note that replacing `y` with `x` in `xy` would yield `x^2`, so this process is not always reversable.

`symbolicAlgebra.polynomial.scale(c, polynomial)`

Returns the product of the constant `c` and `polynomial`.

`symbolicAlgebra.polynomial.simplify(polynomial)`

Returns a cleaned-up version of polynomial.

### Coming soon...

`symbolicAlgebra.polynomial.derivative(polynomial)`

Computes the formal derivative of `polynomial`.

`symbolicAlgebra.polynomial.integral(polynomial, c)`

Compute the formal antiderivative of `polynomial` plus the constant `c`.

## symbolicAlgebra.monomial

Monomials are the building blocks of polynomials. A monomial is a constant multiple times some product of variables. For instance, `2xy` is a monomial, as is `xy^2 z`, but `1 + x` is not. "Under the hood," a monomial is an array with two entries and one field `monomial.coeff`. The first entry is a set of symbols and the second is the array of exponents that each symbol should be raised to. `monomial.coeff` is the coefficient of the monomial.

The current functionality is as follows:

`+, -, *`

Monomials can be multiplied, added and subtracted with each other.

`symbolicAlgebra.monomial.new(symbols, exponents, coefficient)`:

Creates a new monomial with symbols given by `symbols`, exponents given by `exponents`, and coefficient given by `coefficient`. For example, `symbolicAlgebra.monomial.new({"x","y"}, {2,3}, 3)` creates the monomial `3x^2 y^3`

`symbolicAlgebra.monomial.copy(monomial)`

Creates a shallow copy of `monomial`.

`symbolicAlgebra.monomial.eval(monomial, rules)`

Evaluates the monomial at some list of values. `rules` should be of the form `{x = a, y = b, ...}`.

`symbolicAlgebra.monomial.replace(monomial, symbolToBeReaplaced, symbol)`

Returns a monomial where `symbolToBeReplaced` has been replaced by `symbol`. Note that replacing `y` with `x` in `xy` would yield `x^2`, so this process is not always reversable.

`symbolicAlgebra.monomial.scale(c, monomial)`

Returns a monomial which only differs from `monomial` by the coefficient. The coefficient of the returned monomial will be `c * monomial.coeff`.

`symbolicAlgebra.monomial.simplify(monomial)`

Simplifies the monomial by removing zero exponents and monomials equal to `0`.

## MatrixAlgebra

TODO
