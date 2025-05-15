"""
    SymbolicSigns

This package provides types to work with symbolic signs that often come up
in differential homological algebra and other graded contexts. Specifically,
given any commutative ring `R` and any type `T`, the package allows to define
a new commutative ring `SymbolicSignRing{T,R}` that extends `R` by expressions
of the form `(-1)^|x|` and `(-1)^(|x|*|y|)` where `|x|` and `|y|` represent
the symbolic degrees of `x::T` and `y::T`, respectively.

The symbolic signs `(-1)^(|x|*|y|)` and `(-1)^(|y|*|x|)` are treated as equal,
and `(-1)^(|x|*|x|)` is simplified to `(-1)^|x|`.

See also [`SymbolicSignRing`](@ref), [`DegTerm`](@ref).
"""
module SymbolicSigns

export DegTerm, SignExp, Sign, signexp, SymbolicSignRing

using StructEqualHash, LinearCombinations, Modulo2

import Base: show, ==, hash, one, isone, iszero, iseven, isodd, +, -, *, ^,
    convert, length, iterate, eltype, copy, promote_rule

using LinearCombinations: linear_convert
using LinearCombinations: Sign as LCSign
import LinearCombinations: Zero, termcoeff, signed, sign_type, show_summand

#
# DegTerm
#

"""
    DegTerm{T}

    DegTerm{T}() where T
    DegTerm(x::T) where T
    DegTerm(x::T, y::T) where T

This type represents a *degree term* representing the degree of an element from `T`
or the product of the degrees of two elements from `T`. There is also a version
with no elements from `T` to represent the constant `1` as a degree term.

The degree terms `DegTerm(x, y)` and `DegTerm(y, x)` are treated as equal,
and `DegTerm(x, x)` is simplified to `DegTerm(x)`.

See also [`SignExp`](@ref).

# Examples
```jldoctest
julia> dx = DegTerm(:x)
|x|

julia> typeof(dx)
DegTerm{Symbol}

julia> dy = DegTerm(:y)
|y|

julia> dx*dy
|x||y|

julia> dx*dy == DegTerm(:x, :y)
true

julia> dx*dy == DegTerm(:y, :x)   # order doesn't matter
true

julia> dx*dx   # note the simplification
|x|

julia> dz = DegTerm(:z); dx*dy*dz
ERROR: cannot multiply two DegTerm's with more than 2 factors in total
[...]

julia> d1 = DegTerm{Symbol}()
1

julia> d1 == one(dx)
true
```
"""
struct DegTerm{T}
    n::Int8
    x::T
    y::T
    DegTerm{T}(x, y) where T = x == y ? new{T}(1, x) : new{T}(2, x, y)
    # do we need this test for equality in applications?
    DegTerm{T}(x) where T = new{T}(1, x)
    DegTerm{T}() where T = new{T}(0)
end

DegTerm(x::T...) where T = DegTerm{T}(x...)

function show(io::IO, t::T) where T <: DegTerm
    if t.n == 2
        xr = string(t.x)
        yr = string(t.y)
        if xr <= yr
            print(io, '|', xr, "||", yr, '|')
        else
            print(io, '|', yr, "||", xr, '|')
        end
    elseif t.n == 1
        print(io, '|', t.x, '|')
    else
        print(io, '1')
    end
end

length(t::DegTerm)::Int = t.n

eltype(t::DegTerm{T}) where T = T

function iterate(t::DegTerm, i = 1)
    if i > length(t)
        nothing
    elseif i == 1
        t.x, 2
    elseif i == 2
        t.y, 3
    end
end

function ==(t::DegTerm{T}, u::DegTerm{T}) where T
    if t.n != u.n
        false
    elseif t.n == 2
        (t.x == u.x && t.y == u.y) || (t.x == u.y && t.y == u.x)
    elseif t.n == 1
        t.x == u.x
    else
        true
    end
end

function hash(t::T, h::UInt) where T <: DegTerm
    h = hash(T, h)
    if t.n == 2
        hash(t.x, h) * hash(t.y, h)
    elseif t.n == 1
        hash(t.x, h)
    else
        h
    end
end

isone(t::DegTerm) = iszero(length(t))

one(::Type{T}) where T <: DegTerm = T()
one(::T) where T <: DegTerm = one(T)

function *(t::DegTerm{T}, u::DegTerm{T}) where T
    if t.n == 1 && u.n == 1
        DegTerm(t.x, u.x)
    elseif t.n == 0
        u
    elseif u.n == 0
        t
    else
        error("cannot multiply two DegTerm's with more than 2 factors in total")
    end
end

@linear_broadcastable DegTerm

#
# SignExp
#

"""
    SignExp{T} == Linear{DegTerm{T}, ZZ2}

A *sign exponent* is a linear combination of degree terms with coefficients in `ZZ2`.

See also [`LinearCombinations.Linear`](https://matthias314.github.io/LinearCombinations.jl/stable/linear/#LinearCombinations.Linear),
`Modulo2.ZZ2`, [`DegTerm`](@ref), [`Sign`](@ref).

# Examples
```jldoctest
julia> se = dx + dx*dy + 1
|x||y|+|x|+1

julia> typeof(se) == SignExp{Symbol}
true

julia> se == SignExp(dx => 1, dx*dy => -1, DegTerm{Symbol}() => -1)
true

julia> SignExp(dx) == SignExp(dx => 1)   # note the short form
true
```
"""
const SignExp{T} = Linear{DegTerm{T}, ZZ2}

SignExp(t::Pair{DegTerm{T}}...) where T = SignExp{T}(t)
SignExp(t::DegTerm{T}...) where T = SignExp{T}(map(u -> u => one(ZZ2), t))
SignExp{T}(n::Number) where T = Linear(one(DegTerm{T}) => ZZ2(n))

+(t::DegTerm{T}, u::DegTerm{T}) where T =
    t == u ? zero(SignExp{T}) : SignExp{T}(t => one(ZZ2), u => one(ZZ2))

function +(t::DegTerm, n::Number)
    u = one(t)
    t == u ? SignExp(t => ZZ2(n+1)) : SignExp(t => one(ZZ2), u => ZZ2(n))
end

+(e::SignExp{T}, n::Number) where T = addmul(e, one(DegTerm{T}), ZZ2(n))

+(n::Number, t::DegTerm) = t+n
+(n::Number, e::SignExp) = e+n

-(t::DegTerm{T}, u::DegTerm{T}) where T = t+u
-(t::DegTerm, n::Number) = t+n
-(e::SignExp, n::Number) = e+n
-(n::Number, t::DegTerm) = t+n
-(n::Number, e::SignExp) = e+n

*(::DegTerm, ::Zero) = Zero()
*(::Zero, ::DegTerm) = Zero()
*(::SignExp, ::Zero) = Zero()
*(::Zero, ::SignExp) = Zero()
*(c::Number, t::DegTerm) = SignExp(t => c)
*(t::DegTerm, c::Number) = c*t

signed(e::Union{SignExp,DegTerm}, c) = Sign(e)*c
signed(::Union{SignExp, DegTerm}, c::ZZ2) = c

convert(::Type{SignExp{T}}, n::Number) where T = SignExp{T}(n)

#
# Sign
#

"""
    Sign{T}

This is a wrapper around `SignExp{T}` that allows to use a multiplicative notation.

See also [`SignExp`](@ref), [`SymbolicSignRing`](@ref).

# Examples
```jldoctest
julia> s1 = Sign(dx + dx*dy)
(-1)^(|x||y|+|x|)

julia> s2 = Sign(dy + 1)
(-1)^(|y|+1)

julia> s1*s2
(-1)^(|y|+|x||y|+|x|+1)

julia> isone(s1^2)
true

julia> convert(ZZ2, s1)
1
```
"""
struct Sign{T}
    e::SignExp{T}
end

Sign(args...) = Sign(SignExp(args...))

function show(io::IO, s::Sign)
    if isone(s)
        print(io, '1')
    else
        e = signexp(s)
        print(io, "(-1)^")
        length(e) > 1 && print(io, '(')
        show(IOContext(io, :compact => true), MIME"text/plain"(), e)
        length(e) > 1 && print(io, ')')
    end
end

@struct_equal_hash Sign{T} where T

signexp(s::Sign) = s.e

copy(s::Sign) = Sign(copy(s.e))

length(s::Sign) = length(signexp(s))

eltype(s::Sign) = eltype(signexp(s))

iterate(s::Sign, state...) = iterate(signexp(s), state...)

iszero(::Sign) = false

+(s::Sign) = s
# do we need this?
# -(s::Sign) = SymbolicSignRing(s => -LinearCombinations.ONE)

isone(s::Sign) = iszero(signexp(s))

one(::Type{Sign{T}}) where T = Sign(zero(SignExp{T}))
one(::T) where T <: Sign = one(T)

*(s1::Sign{T}, s2::Sign{T}) where T = Sign(signexp(s1) + signexp(s2))

function mul!(s1::Sign{T}, s2::Sign{T}) where T
    add!(signexp(s1), signexp(s2))
    s1
end

^(s::Sign, n::Integer) = isodd(n) ? s : one(s)

convert(::Type{ZZ2}, ::Sign) = one(ZZ2)

iseven(::Sign) = true
isodd(::Sign) = false

@linear_broadcastable Sign

function show_summand(io::IO, s::Sign, cs)
    show(io, MIME"text/plain"(), cs)
    isone(s) && return
    print(io, '*')
    show(io, MIME"text/plain"(), s)
end

#
# SymbolicSignRing
#

"""
    LinearCombinations.termcoeff((s, c)::Pair{Sign{T}}) where T -> Pair{Sign{T}}

Convert the term-coefficient pair `(s, c)` to one where the exponent of the sign
has no constant term. If the exponent of `s` has a constant term, then it is
removed and instead the sign of `c` changed.

This is important to make linear combinations of signs unique.

See also [`LinearCombinations.termcoeff`](https://matthias314.github.io/LinearCombinations.jl/stable/extensions/#LinearCombinations.termcoeff).

# Examples
```jldoctest
julia> using LinearCombinations: termcoeff

julia> dx = DegTerm(:x)
|x|

julia> termcoeff(Sign(dx) => 2)
(-1)^(|x|) => 2

julia> termcoeff(Sign(dx+1) => 2)
(-1)^(|x|) => -2
```
"""
function termcoeff(sc::Pair{Sign{T}}) where T
    s, c = sc
    x = DegTerm{T}()
    if isone(s.e[x])
        Sign(addmul(s.e, x, one(ZZ2))) => -c
    else
        sc
    end
end

"""
    SymbolicSignRing{T,R} = Linear{Sign{T},R}

This type extends the commutative ring `R` by elements of type `Sign{T}`.
Elements are formal linear combinations with terms of type `Sign{T}` and
coefficients of type `R`.

See also [`LinearCombinations.Linear`](https://matthias314.github.io/LinearCombinations.jl/stable/linear/#LinearCombinations.Linear),
[`Sign`](@ref).

# Examples
```jldoctest
julia> dx, dy = DegTerm(:x), DegTerm(:y)
(|x|, |y|)

julia> s1, s2 = Sign(dx + dx*dy), Sign(dy + 1)
((-1)^(|x||y|+|x|), (-1)^(|y|+1))

julia> a = Linear(s1 => 2, s2 => -1)
(-1)^(|y|)+2*(-1)^(|x||y|+|x|)

julia> convert(ZZ2, a)
1

julia> 3.0*a
6.0*(-1)^(|x|+|x||y|)+3.0*(-1)^(|y|)

julia> ZZ2(1)*a
1
```
Note that the extended ring has zero divisors:
```jldoctest
julia> a = Linear(Sign(dx) => 1, Sign(dy) => 1)
(-1)^(|y|)+(-1)^(|x|)

julia> b = Linear(Sign(dx) => 1, Sign(dy) => -1)
-(-1)^(|y|)+(-1)^(|x|)

julia> a*b
0
```
"""
const SymbolicSignRing{T,R} = Linear{Sign{T},R}

# needed because "SymbolicSignRing" means "SymbolicSignRing{T,R} where {T,R}"
SymbolicSignRing(x::Pair{S}...) where S <: Sign = Linear{S}(x...)

*(s::Sign{T}, c::R) where {T,R<:Number} = SymbolicSignRing{T,R}(s => c)
*(c::R, s::Sign{T}) where {T,R<:Number} = s*c

*(a::SymbolicSignRing{T}, s::Sign{T}) where T = LinearCombinations.mul(a, s)

for op in (:(+), :(-), :(*))
    @eval $op(a::SymbolicSignRing, b::SymbolicSignRing) = invoke($op, Tuple{AbstractLinear, AbstractLinear}, a, b)
    if op in (:(+), :(-))
        @eval $op(a::SymbolicSignRing{T,R}, s::Sign{T}) where {T,R} = addmul(a, s, $op(one(R)))
        @eval $op(s::Sign{T}, a::SymbolicSignRing{T}) where T = add!($op(a), s)
    end
    @eval function $op(a::SymbolicSignRing{T}, c) where T
        S = typeof(c)
        if has_char2(S)
            $op(convert(ZZ2, a), c)
        elseif $op in (+, -) && c isa Number
            addmul(a, one(Sign{T}), $op(c))
        else
            R = coefftype(a)
            invoke($op, Tuple{AbstractLinear{Sign{T},R}, S}, a, c)
        end
    end
end

*(a::SymbolicSignRing, s::LCSign) = isone(s) ? a : -a

function -(c, a::SymbolicSignRing{T,R}) where {T,R}
    S = typeof(c)
    if has_char2(S)
        convert(ZZ2, a) + c
    elseif c isa Number
        addmul!((-one(R)*one(c))*a, one(Sign{T}), c)
    else
        invoke(-, Tuple{S, AbstractLinear{Sign{T},R}}, c, a)
    end
end

function convert(::Type{ZZ2}, a::SymbolicSignRing)
    # error("converting: $a::$(typeof(a))")
    sum(ZZ2, coeffs(a); init = zero(ZZ2))
end

iseven(a::SymbolicSignRing) = iseven(ZZ2(a))
isodd(a::SymbolicSignRing) = isodd(ZZ2(a))

# many conversion methods are needed to avoid ambiguities
convert(::Type{SymbolicSignRing{T,R}}, a::SymbolicSignRing{T}) where {T,R} =
    linear_convert(SymbolicSignRing{T,R}, a)
# convert(::Type{S}, a::SymbolicSignRing) where S <: SymbolicSignRing = linear_convert(S, a)
convert(::Type{SymbolicSignRing{T,R}}, s::Sign{T}) where {T,R} = SymbolicSignRing{T,R}(s => one(R))
convert(::Type{SymbolicSignRing{T,R}}, c) where {T,R} = SymbolicSignRing{T,R}(one(Sign{T}) => c)

# promote_rule(::Type{SymbolicSignRing{T,R}}, ::Type{R}) where {T,R} = SymbolicSignRing{T,R}
# promote_rule(::Type{SymbolicSignRing{T,R}}, ::Type{Sign{T}}) where {T,R} = SymbolicSignRing{T,R}

promote_rule(::Type{<:Sign}, ::Type{ZZ2}) = ZZ2
promote_rule(::Type{SymbolicSignRing{T,R}}, ::Type{ZZ2}) where {T,R} = promote_rule(R, ZZ2)

sign_type(::Type{SignExp{T}}) where T = SymbolicSignRing{T,LCSign}
sign_type(::Type{DegTerm{T}}) where T = SymbolicSignRing{T,LCSign}

end # module
