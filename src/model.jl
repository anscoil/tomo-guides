using LinearAlgebra

struct TiltedPropagator1D{T,N}
    G
    pfft
    pifft
    
    function TiltedPropagator1D(A::AType, λ::Real, dx::Real, z::Real,
                                θx::Real) where {T <: Real,
                                                 AType <: AbstractVector{Complex{T}}}
        n = length(A)
        fx = fftfreq(n, 1/dx)
        f0 = sin(θx)/λ
        f0c = cos(θx)/λ
        G = AType([exp(im*T(2*π*abs(z))*(sqrt(Complex{T}(1/λ^2-(fx+f0)^2))-T(f0c)))
                   for fx in fx])
        if z < 0; G .= conj.(G) end
        A_tmp = similar(G)
        pfft = plan_fft!(A_tmp)
        pifft = plan_ifft!(A_tmp)
        new{T,1}(G, pfft, pifft)
    end

    function TiltedPropagator1D(AType::Type{A}, ::Type{T},
                                nx::Integer, λ::Real, dx::Real, z::Real,
                                θv::AbstractVector{U}) where {T, U <: Real,
                                                              A <: AbstractArray}
        nθ = length(θv)
        fx = fftfreq(nx, 1/dx)
        f0 = sin.(θv)/λ
        f0c = cos.(θv)/λ
        G = AType([exp(im*T(2*π*abs(z))*(sqrt(Complex{T}(1/λ^2-(fx+f0[j])^2))-T(f0c[j])))
                   for fx in fx, j in 1:nθ])
        if z < 0; G .= conj.(G) end
        A_tmp = similar(G)
        pfft = plan_fft!(A_tmp, 1)
        pifft = plan_ifft!(A_tmp, 1)
        new{T,2}(G, pfft, pifft)
    end

    function TiltedPropagator1D(A::AbstractArray{Complex{T},2},
                                λ::Real, dx::Real, z::Real,
                                θv::AbstractVector{U}) where {T, U <: Real}
        nx, nθ = size(A)
        @assert length(θv) == nθ
        TiltedPropagator1D(typeof(A).name.wrapper, T, nx, λ, dx, z, θv)
    end
end

function propagate!(A::AbstractArray{Complex{T},N}, p;
                    direction::Symbol=:forward) where {T, N}
    p.pfft * A
    if direction == :forward
        A .*= p.G
    elseif direction == :backward
        A .*= conj.(p.G)
    else
        error("Wrong propagation direction")
    end
    p.pifft * A
end

function make_poly(fxv::AbstractVector{T},
                   f0v::AbstractVector{T},
                   λ::Real, NA::Real, p::Integer) where {T <: Real}
    [if abs(fx+f0) <= abs(NA/λ); T(abs(fx+f0)*λ/NA)^p; else T(0) end
     for fx in fxv, f0 in f0v]
end

struct Wg_bpm
    p_bpm
    p_bpm_half
    p_medium
    cosθv
    k0dz
    λ; n0; nx; nz; dx; dz
    A_u
    uf
    Ith

    function Wg_bpm(λ::Real, n0::Real,
                    nx::Integer, nz::Integer,
                    dx::Real, dz::Real, z_medium::Real,
                    Ith::AbstractArray{T,2},
                    θv::AbstractVector{U}) where {T <: Real,
                                                  U <: Real}
        nθ = length(θv)
        @assert size(Ith) == (nx,nθ)
        AType = typeof(Ith).name.wrapper
        p_bpm = TiltedPropagator1D(AType, T, nx, λ/n0, dx, dz, θv)
        p_bpm_half = TiltedPropagator1D(AType, T, nx, λ/n0, dx, dz/2, θv)
        p_medium = TiltedPropagator1D(AType, T, nx, λ/n0, dx, z_medium, θv)
        cosθv = similar(p_bpm.G, T, nθ)
        copyto!(cosθv, cos.(T.(θv)))
        k0dz = T(2π/λ*dz)
        A_u = similar(p_bpm.G, (nx,nθ,nz))
        uf = similar(p_bpm.G)
        new(p_bpm, p_bpm_half, p_medium, cosθv, k0dz,
            T(λ), T(n0), nx, nz, T(dx), T(dz), A_u, uf, Ith)
    end
end

struct Wg_ϕ{Wg_T}
    wg :: Wg_T
    uf
    Ith
    fft_uf
    kxv
    ϕ
    ∇ϕ
    pfft
    pifft
    Zkv
    NA
    
    function Wg_ϕ{Wg_T}(λ::Real, n0::Real,
                        nx::Integer, nz::Integer,
                        dx::Real, dz::Real, z_medium::Real, NA::Real,
                        Ith::AbstractArray{T,2},
                        θv::AbstractVector{U}; nzk::Integer=1) where {Wg_T <: Wg_bpm,
                                                                  T <: Real,
                                                                  U <: Real}
        nθ = length(θv)
        wg = Wg_T(λ, n0, nx, nz, dx, dz, z_medium, Ith, θv)
        AType = typeof(Ith).name.wrapper
        fft_uf = similar(wg.uf)
        ϕ = similar(wg.uf, T)
        ∇ϕ = similar(ϕ)
        A_tmp = similar(wg.uf)
        pfft = plan_fft!(A_tmp, 1)
        pifft = plan_ifft!(A_tmp, 1)
        fxv = fftfreq(nx, T(1/dx))
        kxv = fftfreq(nx, T(2π/dx))
        f0v = T.(sin.(θv).*(n0/λ))
        pupil = AType([if abs(fx+f0) <= abs(NA/λ); T(1) else T(0) end
                       for fx in fxv, f0 in f0v])
        Zkv = similar(ϕ, (nx,nθ,1+nzk))
        poly_defocus = [if abs(fx+f0) <= abs(NA/λ);
                        T(sqrt(1/(λ/n0)^2-(fx+f0)^2) - n0/λ);
                        else T(0) end for fx in fxv, f0 in f0v]
        copyto!(@view(Zkv[:,:,1]), poly_defocus)
        f_zk = p -> make_poly(fxv, f0v, λ, NA, p)
        for p in 1:nzk
            copyto!(@view(Zkv[:,:,1+p]), f_zk(2*p))
        end
        new{Wg_T}(wg, wg.uf, wg.Ith, fft_uf, kxv, ϕ, ∇ϕ, pfft, pifft, Zkv, T(NA))
    end
end

Wg_bpm_ϕ = Wg_ϕ{Wg_bpm}

function compute_fwd_wg!(wg::Wg_bpm, RI)
    nz = size(RI, 2)
    wg.uf .= 1
    propagate!(wg.uf, wg.p_bpm_half)
    @views wg.uf .*= exp.(im.*wg.k0dz.*(RI[:,1] .- n0)./wg.cosθv')
    # The fields are saved just after crossing the phase masks forward
    wg.A_u[:,:,1] .= wg.uf
    @inbounds for l in 2:nz
        propagate!(wg.uf, wg.p_bpm)
        @views wg.uf .*= exp.(im.*wg.k0dz.*(RI[:,l] .- n0)./wg.cosθv')
        wg.A_u[:,:,l] .= wg.uf
    end
    propagate!(wg.uf, wg.p_bpm_half)
    propagate!(wg.uf, wg.p_medium)
end

function get_compute_fwd_wg!(nx::Integer, dx::Real, λ::Real,
                             θv::AbstractVector{T}, α::Real) where {T <: Real}
    fxv = fftfreq(nx, 1/dx)
    function f!(wg::Wg_bpm, RI)
        uf = compute_fwd_wg!(wg, RI)
        fft!(uf,1)
        uf .*= exp.(im*2π*α*(fxv .+ sin.(θv)'/λ).^2)
        ifft!(uf,1)
    end
    f!
end

function compute_err_wg!(wg::Wg_bpm, uf)
    nθ = size(uf,2)
    sum((abs.(uf).^2 .- wg.Ith).^2)/nθ
end

function compute_err_vector_wg!(wg::Wg_bpm, uf)
    sum((abs.(uf).^2 .- wg.Ith).^2, dims=1)[1,:]
end

function compute_grad_wg_err!(wg::Wg_bpm, RI, ∇RI, uerr)
    nz = size(RI,2)
    propagate!(uerr, wg.p_medium, direction=:backward)
    propagate!(uerr, wg.p_bpm_half, direction=:backward)
    @views ∇RI[:,end] .= wg.k0dz*sum(imag.(conj(wg.A_u[:,:,end])
                                           .*uerr./wg.cosθv'), dims=2)[:,1]
    @views uerr .*= exp.(-im.*wg.k0dz.*(RI[:,end] .- n0))
    @inbounds for l in nz-1:-1:1
        propagate!(uerr, wg.p_bpm, direction=:backward)
        @views ∇RI[:,l] .= wg.k0dz*sum(imag.(conj(wg.A_u[:,:,l])
                                             .*uerr./wg.cosθv'), dims=2)[:,1]
        @views uerr .*= exp.(-im.*wg.k0dz.*(RI[:,l] .- n0))
    end
    propagate!(uerr, wg.p_bpm_half, direction=:backward)
    ∇RI
end

function compute_grad_wg!(wg::Wg_bpm, RI, ∇RI, uf)
    uf .*= 4 .*(abs.(uf).^2 .- wg.Ith)
    compute_grad_wg_err!(wg, RI, ∇RI, uf)
end

function get_compute_grad_wg!(nx::Integer, dx::Real, λ::Real,
                              θv::AbstractVector{T}, α::Real) where {T <: Real}
    fxv = fftfreq(nx, 1/dx)
    function f!(wg::Wg_bpm, RI, ∇I, uf)
        uf .*= 4 .*(abs.(uf).^2 .- wg.Ith)
        fft!(uf,1)
        uf .*= exp.(-im*2π*α*(fxv .+ sin.(θv)'/λ).^2)
        ifft!(uf,1)
        compute_grad_wg_err!(wg, RI, ∇RI, uf)
    end
    f!
end

function compute_fwd_wg!(wg::Wg_ϕ,
                         RIϕ::Tuple{AType,AType}) where {T <: Real,
                                                         AType <: AbstractArray{T,2}}
    RI, ϕ = RIϕ
    uf = compute_fwd_wg!(wg.wg, RI)
    wg.pfft * uf
    uf .*= exp.(im.*ϕ)
    wg.fft_uf .= uf
    wg.pifft * uf
    uf
end

function compute_err_wg!(wg::Wg_ϕ, uf)
    compute_err_wg!(wg.wg, uf)
end

function compute_grad_wg!(wg::Wg_ϕ,
                          RIϕ::Tuple{AType,AType},
                          ∇RIϕ::Tuple{AType,AType},
                          uf) where {T <: Real,
                                     AType <: AbstractArray{T,2}}
    RI, ϕ = RIϕ
    ∇RI, ∇ϕ = ∇RIϕ
    uf .*= 4 .*(abs.(uf).^2 .- wg.Ith)
    wg.pfft * uf
    @views ∇ϕ .= imag.(conj.(wg.fft_uf).*uf)
    uf .*= exp.(-im.*ϕ)
    wg.pifft * uf
    compute_grad_wg_err!(wg.wg, RI, ∇RI, uf)
    ∇RI
end

function compute_fwd_wg!(wg::Wg_ϕ,
                         RIX::Tuple{AType,VType}) where {T <: Real,
                                                         AType <: AbstractArray{T,2},
                                                         VType <: AbstractVector{T}}
    RI, Xv = RIX
    wg.ϕ .= wg.kxv.*Xv'
    compute_fwd_wg!(wg, (RI, wg.ϕ))
end

function compute_grad_wg!(wg::Wg_ϕ,
                          RIX::Tuple{AType,VType},
                          ∇RIX::Tuple{AType,VType},
                          uf) where {T <: Real,
                                     AType <: AbstractArray{T,2},
                                     VType <: AbstractVector{T}}
    RI, Xv = RIX
    ∇RI, ∇Xv = ∇RIX
    wg.ϕ .= wg.kxv.*Xv'
    compute_grad_wg!(wg, (RI, wg.ϕ), (∇RI, wg.∇ϕ), uf)
    @views ∇Xv .= sum(wg.kxv .* wg.∇ϕ, dims=1)[1,:]
    ∇RIX
end

function compute_fwd_wg_Zk!(wg::Wg_ϕ,
                            RIZ::Tuple{AType,VType}) where {T <: Real,
                                                            AType <: AbstractArray{T,2},
                                                            VType <: AbstractVector{T}}
    RI, zkv = RIZ
    nzk = length(zkv)
    @assert nzk == size(wg.Zkv,3)
    zkvr = reshape(zkv, (1,1,nzk))
    @views wg.ϕ .= sum(wg.Zkv .* zkvr, dims=3)[:,:,1]
    compute_fwd_wg!(wg, (RI, wg.ϕ))
end

function compute_grad_wg_Zk!(wg::Wg_ϕ,
                             RIZ::Tuple{AType,VType},
                             ∇RIZ::Tuple{AType,VType},
                             uf) where {T <: Real,
                                        AType <: AbstractArray{T,2},
                                        VType <: AbstractVector{T}}
    RI, zkv = RIZ
    ∇RI, ∇zkv = ∇RIZ
    nzk = length(zkv)
    @assert nzk == size(wg.Zkv,3)
    zkvr = reshape(zkv, (1,1,nzk))
    @views wg.ϕ .= sum(wg.Zkv .* zkvr, dims=3)[:,:,1]
    compute_grad_wg!(wg, (RI, wg.ϕ), (∇RI, wg.∇ϕ), uf)
    @views ∇zkv .= sum(wg.Zkv .* wg.∇ϕ, dims=(1,2))[1,1,:]
    ∇RIZ
end

function get_aberration(wg::Wg_bpm_ϕ, zkv)
    λ = wg.wg.λ
    nx = wg.wg.nx
    dx = wg.wg.dx
    n0 = wg.wg.n0
    NA = wg.NA
    
    pmax = length(zkv)
    fxv = fftfreq(nx, 1/dx)
    defocus = [if abs(fx) <= abs(NA/λ);
               sqrt(1/(λ/n0)^2-fx^2) - n0/λ;
               else 0.0 end for fx in fxv]
    l = [if abs(fx) <= abs(NA/λ); (abs(fx)*λ/NA)^(2*p); else 0.0 end
         for fx in fxv, p in 1:pmax-1]
    A = sum(l .* zkv[2:end]', dims=2)[:,1]
    fftshift(fxv)*λ, fftshift(A .+ zkv[1]*defocus)
end
