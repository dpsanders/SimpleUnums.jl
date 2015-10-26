immutable Unum
    s::Int    # sign
    e::Int    # exponent
    f::Int    # fraction
    u::Int    # ubit
    esm1::Int   # exponent size - 1
    fsm1::Int   # fraction size - 1
    ########
    es::Int
    fs::Int
    bias::Int

    Unum(s, e, f, u, esm1, fsm1, bias) = new(s, e, f, u, esm1, fsm1, esm1+1, fsm1+1, (1<<bias)-1)

end

Unum(s, e, f, u, esm1, fsm1) = Unum(s, e, f, u, esm1, fsm1, esm1)  # bias is esm1

Base.sign(u::Unum) = u.s == 0 ? +1 : -1

esizesize,fsizesize = 2, 2  # utag = [esizesize, fsizesize]

function call(::Type{Float64}, u::Unum)
    denominator = 1 << u.fs

    if u.e == 0
        sign(u) * 2.0^(1-u.bias) * u.f / denominator  # assumes u.s is 0 or 1

    elseif u.e == 1<<u.es-1 && u.f==1<<u.fs-1 && u.es==(1<<esizesize) && u.fs==(1<<fsizesize)
        # not efficient!
        sign(u) * Inf

    else
        sign(u) * 2.0^(u.e-u.bias) * (1.0 + u.f/denominator)
    end

end

Base.convert(::Type{Float64}, u::Unum) = Float64(u)


#all_unums(esizesize, fsizesize) = all_unums(esizesize, fsizesize, esizesize-1)

function all_unums(esizesize, fsizesize, bias)

    es = 1 << esizesize
    fs = 1 << fsizesize

    unums = vec([Unum(0, e, f, 0, es-1, fs-1, bias) for e=0:2^es-1, f=0:2^fs-1])

    nums = map(Float64, unums)
    permutation = sortperm(nums)

    unums = unums[permutation]
    nums = nums[permutation]

    unums, nums
end

function all_unums(esizesize, fsizesize)

    es = 1 << esizesize
    fs = 1 << fsizesize

    unums = vec([Unum(0, e, f, 0, es-1, fs-1) for e=0:2^es-1, f=0:2^fs-1])

    nums = map(Float64, unums)
    permutation = sortperm(nums)

    unums = unums[permutation]
    nums = nums[permutation]

    unums, nums
end

function extract_parts(x::Float64)
    e = floor(log2(x))
    y = x/ 2^e
    m = Int((y-1) * 2^52)
    s = copysign(1, x)

    e, y, m, s
end

function inrange(x::Integer, sizesize)
    # 0 exponent has a special meaning
    bias = 2^(sizesize-1) - 1
    max_value = 2^sizesize-1
    min_value = max_value - bias
    max_value -= bias


end

function representation(x::Float64, esizesize, fsizesize)
    e, y, m, s = extract_parts(x)
end
