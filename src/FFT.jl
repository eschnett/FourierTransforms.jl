module FFT

################################################################################

function direct_ft!(X::AbstractVector{T}, x::AbstractVector{T}) where {T<:Complex}
    N = length(x)
    @assert length(X) == N
    RT = typeof(real(zero(T)))
    for k in 1:N
        s = zero(T)
        for n in 1:N
            s += cispi(-2 * (n - 1) * (k - 1) / RT(N)) * x[n]
        end
        X[k] = s
    end
    return X
end

direct_ft(x::AbstractVector) = direct_ft!(similar(x), x)

################################################################################

function direct_ft_1!(X::AbstractVector{T}, x::AbstractVector{T}) where {T<:Complex}
    @assert length(X) == length(x) == 1
    X[1] = x[1]
    return X
end

direct_ft_1(x::AbstractVector) = direct_ft_1!(similar(x), x)

################################################################################

function direct_ft_2!(X::AbstractVector{T}, x::AbstractVector{T}) where {T<:Complex}
    @assert length(X) == length(x) == 2
    X[1] = x[1] + x[2]
    X[2] = x[1] - x[2]
    return X
end

direct_ft_2(x::AbstractVector) = direct_ft_2!(similar(x), x)

################################################################################

function direct_ft_4!(X::AbstractVector{T}, x::AbstractVector{T}) where {T<:Complex}
    @assert length(X) == length(x) == 4
    X[1] = (x[1] + x[3]) + (x[2] + x[4])
    X[2] = (x[1] - x[3]) - im * (x[2] - x[4])
    X[3] = (x[1] + x[3]) - (x[2] + x[4])
    X[4] = (x[1] - x[3]) + im * (x[2] - x[4])
    return X
end

direct_ft_4(x::AbstractVector) = direct_ft_4!(similar(x), x)

################################################################################

# <https://en.wikipedia.org/w/index.php?title=Cooley%E2%80%93Tukey_FFT_algorithm&oldid=1056450535>
function ditfft2!(X::AbstractVector{T}, x::AbstractVector{T}) where {T<:Complex}
    N = length(x)
    @assert length(X) == N
    RT = typeof(real(zero(T)))
    if N == 0
        # do nothing
    elseif N == 1
        X[1] = x[1]
    else
        # Ensure N is even
        @assert N % 2 == 0
        N2 = N ÷ 2
        ditfft2!((@view X[(0 * N2 + 1):(1 * N2)]), (@view x[1:2:end]))
        ditfft2!((@view X[(1 * N2 + 1):(2 * N2)]), (@view x[2:2:end]))
        for k in 1:N2
            ϕ = cispi(-2 * (k - 1) / RT(N))
            p = X[k]
            q = ϕ * X[k + N ÷ 2]
            X[k] = p + q
            X[k + N ÷ 2] = p - q
        end
    end
    return X
end

ditfft2(x::AbstractVector) = ditfft2!(similar(x), x)

################################################################################

function radix_fft!(X::AbstractVector{T}, x::AbstractVector{T}; radix::Int=2) where {T<:Complex}
    N = length(x)
    @assert length(X) == N
    RT = typeof(real(zero(T)))
    if N == 0
        # do nothing
    elseif N == 1
        X[1] = x[1]
    else
        N₁ = radix              # decimation in time
        @assert N % N₁ == 0
        N₂ = N ÷ N₁
        for n₁ in 1:N₁
            fft!((@view X[(1 + (n₁ - 1) * N₂):(n₁ * N₂)]), (@view x[n₁:N₁:end]); radix=radix)
            for n₂ in 1:N₂
                X[(n₁ - 1) * N₂ + n₂] *= cispi(-2 * (n₁ - 1) * (n₂ - 1) / RT(N))
            end
        end
        for n₂ in 1:N₂
            error("CONTINUE HERE")
            fft!((@view X[(1 + (n₁ - 1) * N₂):(n₁ * N₂)]), (@view x[n₁:N₁:end]); radix=radix)
        end
    end
    return X
end

radix_fft(x::AbstractVector; radix=radix) = radix_fft!(similar(x), x; radix=radix)

################################################################################

export fft!, fft, inv_fft!, inv_fft

fft!(X::AbstractVector, x::AbstractVector) = ditfft2!(X, x)

function inv_fft!(X::AbstractVector{T}, x::AbstractVector{T}) where {T<:Complex}
    # TODO: Don't modify `x`
    x .= conj(x)
    fft!(X, x)
    x .= conj(x)
    X .= conj(X) / length(X)
    return X
end

fft(x::AbstractVector) = fft!(similar(x), x)
inv_fft(x::AbstractVector) = inv_fft!(similar(x), x)

end
