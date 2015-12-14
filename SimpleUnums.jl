import Base:
    sign, repr, show,
    prevfloat, nextfloat

const global esizesize = 2
const global fsizesize = 2

function set_environment(ess, fss)
    global esizesize = ess
    global fsizesize = fss
end

"""Object representing a unum"""
immutable Unum
    s::Int      # sign
    e::Int      # exponent
    f::Int      # fraction
    u::Int      # ubit
    esm1::Int   # exponent size - 1
    fsm1::Int   # fraction size - 1

    ########    # the following are not really required
    es::Int     # exponent size
    fs::Int     # fraction size
    bias::Int   # exponent bias

    Unum(s, e, f, u, esm1, fsm1, bias) = new(s, e, f, u, esm1, fsm1, esm1+1, fsm1+1, (1<<bias)-1)

end

Unum(s, e, f, u, esm1, fsm1) = Unum(s, e, f, u, esm1, fsm1, esm1)  # standard bias is esm1

sign(u::Unum) = u.s == 0 ? +1 : -1

esizesize,fsizesize = 2, 2  # utag = [esizesize, fsizesize]

"Convert Unum to Float64"
# Float64(u)
function call(::Type{Float64}, u::Unum)
    denominator = 1 << u.fs

    if u.e == 0
        sign(u) * 2.0^(1 - u.bias) * u.f / denominator  # assumes u.s is 0 or 1

    elseif u.e == 1<<u.es - 1 && u.f == 1<<u.fs - 1 && u.es == (1<<esizesize) && u.fs == (1<<fsizesize)
        # not efficient!
        copysign(Inf, u)

    else
        copysign(2.0^(u.e - u.bias) * (1.0 + u.f/denominator), sign(u))
    end

end

Base.convert(::Type{Float64}, u::Unum) = Float64(u)


#all_unums(esizesize, fsizesize) = all_unums(esizesize, fsizesize, esizesize-1)

"List all unums with given parameters"
# function all_unums(esizesize, fsizesize, bias)
#
#     es = 1 << esizesize
#     fs = 1 << fsizesize
#
#     unums = vec([Unum(0, e, f, 0, es-1, fs-1, bias) for e=0:2^es-1, f=0:2^fs-1])
#
#     nums = map(Float64, unums)
#     permutation = sortperm(nums)
#
#     unums = unums[permutation]
#     nums = nums[permutation]
#
#     unums, nums
# end

function all_unums(esizesize, fsizesize)  # standard bias es-1

    es = 1 << esizesize
    fs = 1 << fsizesize

    unums = vec([Unum(0, e, f, 0, es-1, fs-1) for e=0:2^es-1, f=0:2^fs-1])

    nums = map(Float64, unums)
    permutation = sortperm(nums)

    unums = unums[permutation]
    nums = nums[permutation]

    unums, nums
end

# last n bits
Base.bits(x, n) = bits(x)[end-n+1:end]

doc"""Extract the parts of the Float64 x.
Returns the sign, exponent and fractional part
(mantissa or significand).
"""
function extract_parts(x::Float64)

    s = Int(signbit(x))
    x = copysign(x, 1)  # make positive

    e = floor(log2(x))

    y = x / 2^e
    f = Int64( (y - 1) * 2^52 )

    s, e, f #y, m, s
end


"Minimal number of bits required to represent the integer `x`"
bits_to_represent(x::Integer) = ceil(Int, log2(x))
# (much) more efficient: 64 - leading_zeros(x)

doc"""Minimal number of bits required to represent the significand / mantissa m
represented as an integer (coming from extract_parts)
"""
function significand_bits(m::Int64)
    # assuming that x is exact
    # e.g. fails for float(pi)

    52 - trailing_zeros(m)
    # more efficient: 52 - trailing_zeros(reinterpret(Int64, x))
end

"""Find the representation of a given `Float64` as a `Unum` with given environment (`esizesize` and `fsizesize`).

Note that *any* `Float64` may be represented via a (possibly inexact) Unum in any environment.

Exponent=0 is special.
"""
function unum_representation(x::Float64, ess=esizesize, fss=fsizesize)

    # special cases:
    if isinf(x)
        # plus and minus infinity

    elseif isnan(x)
        # NaN
    end

    s, e, f = extract_parts(x)
    #@show s, e, f


    e_bits = ceil(Int, log2(abs(e)+1)) + 1
    f_bits = 52 - trailing_zeros(f) #significand_bits(x)

    bias = 2^(e_bits - 1) - 1
    e += bias

    max_exponent_size = 2^ess
    bias = 2^(max_exponent_size-1)  # largest bias

    min_exponent = 0
    max_exponent = 2^max_exponent_size

    min_exponent -= bias
    max_exponent -= bias

    if e == min_exponent
        # subnormal

    elseif e < min_exponent
        # underflow; return (0, minreal)

    elseif e == max_exponent
        # may need to be overflow if one of the special cases

    elseif e > max_exponent
        # overflow; return (maxreal, ∞)

    end

    # if everything OK:

    ubit = 0
    if f_bits > 2^fss
        ubit = 1
        f_bits = 2^fss
    end

    #@show ubit, e_bits, f_bits

    #m = m / 2^(52-f_bits)
    f = f >> (52 - f_bits)

    Unum(s, e, f, ubit, e_bits-1, f_bits-1)
end

unum_representation(x::Integer) = unum_representation(Float64(x))

function special_values(esizesize, fsizesize)

    # max exponent is special
    # min exponent is for subnormals

    ## utagsize =
    ## maxubits
    # posinfu
    ## neginfu
    # qNaNu
    # sNaNu
    # maxrealu
    # negbigu
    # maxreal
    # smallsubnormalu
    # smallsubnormal
end

function nextfloat(u::Unum)
    # NEEDS IMPROVING!

    if u.u == 0  # u is exact; return the next open interval
        return Unum(u.s, u.e, u.f, 1, u.esm1, u.fsm1)
    end

    # u.u is 1, so move to next representable float
    new_e = u.e
    new_es = u.es
    new_fs = u.fs

    #new_f = u.f / (2^u.fs) * 2^fsizesize
    new_fs = 2^fsizesize
    new_f = u.f * 2^(new_fs - u.fs)
    new_f += 1

    #@show new_f


    while (new_f % 2) == 0  # if divisible by 2 #  new_f & 1
        new_f >>= 1
        new_fs -= 1
    end

    if new_fs > 2^fsizesize
        new_f = 0
        new_e += 1
        new_fs = 1
    end

    # else
    #     new_fs += 1
    #     new_f += 1  # changed number of bits
    # end


    if new_e == 2^u.es
        new_es += 1
        # CHECK FOR OVERFLOW
    end


    Unum(u.s, new_e, new_f, 0, new_es-1, new_fs-1)
end

function repr(u::Unum)
    comma = ","
    string("Unum(", u.s, comma, u.e, comma, u.f, comma, u.u, comma, u.esm1, comma, u.fsm1, ")")
end

function show(io::IO, u::Unum)
    space = " "

    print(io, repr(u), "\n")

    print_with_color(:red, io, "s: ", bits(u.s, 1), "; ")
    print_with_color(:blue, io, "e: ", bits(u.e, u.esm1+1), "; ")
    print_with_color(:black, io, "f: ", bits(u.f, u.fsm1+1), "; ")
    print_with_color(:magenta, io, "u: ", bits(u.u, 1), "; ")
    print_with_color(:green, io, "esm1: ", bits(u.esm1, esizesize), "; ")
    print_with_color(:grey, io, "fsm1: ", bits(u.fsm1, fsizesize))

    println(io)

    print(io, "  ", Float64(u))

    if u.u == 1
        println(io, "...")
    else
        println(io, "↓")
    end

    if u.u == 1
        u_next = nextfloat(u)
        f1 = Float64(u)
        f2 = Float64(u_next)

        println(io, "  (", f1, ", ", f2, ")")

        mid = (f1 + f2) / 2.
        rad = f2 - mid

        print(io, "  ", mid, " ± ", rad)
    end


end

#=
TODO:
- Output with + sign to denote that more digits available (use big)
- Interval output for inexact
=#
