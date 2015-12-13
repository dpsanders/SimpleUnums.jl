import Base:
    nextfloat, repr, show


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

Base.sign(u::Unum) = u.s == 0 ? +1 : -1

esizesize,fsizesize = 2, 2  # utag = [esizesize, fsizesize]

"Convert Unum to Float64"
# Float64(u)
function call(::Type{Float64}, u::Unum)
    denominator = 1 << u.fs

    if u.e == 0
        sign(u) * 2.0^(1 - u.bias) * u.f / denominator  # assumes u.s is 0 or 1

    elseif u.e == 1<<u.es - 1 && u.f == 1<<u.fs - 1 && u.es == (1<<esizesize) && u.fs == (1<<fsizesize)
        # not efficient!
        sign(u) * Inf

    else
        sign(u) * 2.0^(u.e - u.bias) * (1.0 + u.f/denominator)
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

function extract_parts(x::Float64)
    e = floor(log2(x))
    y = x / 2^e
    m = Int64( (y - 1) * 2^52 )
    s = Int(signbit(x))

    e, y, m, s
end

# function inrange(x::Integer, sizesize)
#     # 0 exponent has a special meaning
#     bias = 2^(sizesize-1) - 1
#     max_exponent = 2^sizesize
#     min_exponent = 0
#
#     min_exponent = max_exponent - bias
#     max_exponent -= bias
#
#     max_fractional = 2^sizesize - 1
#
# end

function representation(x::Float64, esizesize, fsizesize)
    e, y, m, s = extract_parts(x)
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
function unum_representation(x::Float64, esizesize, fsizesize)

    # special cases:
    if isinf(x)
        # plus and minus infinity

    elseif isnan(x)
        # NaN
    end

    e, y, m, s = extract_parts(x)
    #@show e, y, m, s


    e_bits = ceil(Int, log2(abs(e))) + 1
    f_bits = 52 - trailing_zeros(m) #significand_bits(x)


    max_exponent_size = 2^esizesize
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
    if f_bits > fsizesize
        ubit = 1
        f_bits = fsizesize
    end

    #@show ubit, e_bits, f_bits

    #m = m / 2^(52-f_bits)
    m = m >> (52-f_bits)

    Unum(s, e, m, ubit, e_bits-1, f_bits-1)


end

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
    u = Unum(u.s, u.e, u.f+1, u.u, u.esm1, u.fsm1)
end

function repr(u::Unum)
    comma = ","
    string("Unum(", u.s, comma, u.e, comma, u.f, comma, u.u, comma, u.esm1, comma, u.fsm1, ")")
end

function show(io::IO, u::Unum)
    space = " "

    print(io, repr(u), "\n")

    print_with_color(:red, io, "s: ", bits(u.s, 1))
    print_with_color(:blue, io, "; e: ", bits(u.e, u.esm1+1))
    print_with_color(:black, io, "; f: ", bits(u.f, u.fsm1+1))
    print_with_color(:magenta, io, "; u: ", bits(u.u, 1))
    print(io, "\n")

    print(io, "value: ", Float64(u))
    if u.u == 1
        print(io, "...")
    end

end
