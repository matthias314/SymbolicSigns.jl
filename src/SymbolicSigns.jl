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

struct Sign{T}
    e::SignExp{T}
end

Sign(args...) = Sign(SignExp(args...))

show(io::IO, s::Sign) = print(io, isone(s) ? "1" : "(-1)^($(repr(signexp(s))))")

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

show_summand(io::IO, s::Sign, cs) = isone(s) ? print(io, cs) : print(io, cs, '*', s)

#
# SymbolicSignRing
#

function termcoeff(sc::Pair{Sign{T}}) where T
    s, c = sc
    x = DegTerm{T}()
    if isone(s.e[x])
        Sign(addmul(s.e, x, one(ZZ2))) => -c
    else
        sc
    end
end

const SymbolicSignRing{T,R} = Linear{Sign{T},R}

# SymbolicSignRing{T,R}(itr...) where {T,R} = Linear(Dict{Hashed{Sign{T}},R}(hashedtermcoeff(xc) for xc in itr...))

# the next two definitions are necessary because "SymbolicSignRing" means "SymbolicSignRing{T,R} where {T,R}"
SymbolicSignRing{T,R}(xc::Pair...) where {T,R} = SymbolicSignRing{T,R}(xc)
SymbolicSignRing(x...) = Linear(x...)

+(a::SymbolicSignRing{T,R}, c::R) where {T,R} = addmul(a, one(Sign{T}), c)
+(c::R, a::SymbolicSignRing{T,R}) where {T,R} = a+c

-(a::SymbolicSignRing{T,R}, c::R) where {T,R} = a + (-c)
-(c::R, a::SymbolicSignRing{T,R}) where {T,R} = addmul!((-a), one(Sign{T}), c)

*(s::Sign{T}, c::R) where {T,R<:Number} = SymbolicSignRing{T,R}(s => c)
*(c::R, s::Sign{T}) where {T,R<:Number} = s*c

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