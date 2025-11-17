using BufferedStreams, Match, Nemo

function parse_int(io)
  @assert !eof(io)
  coeff = 0
  while true
    d = peek(io, Char)
    isdigit(d) || break
    coeff = coeff * 10 + Int(read(io, Char) - '0')
    eof(io) && break
  end
  return coeff
end

function parse_monomial!(io, exp::Vector{Int})
  while true
    @assert !eof(io)
    var = Int(read(io, Char) - 'a') + 1

    if !eof(io) && peek(io, Char) == '^'
      read(io, Char)
      exp[var] += parse_int(io)
    else
      exp[var] += 1
    end

    !eof(io) && peek(io, Char) == '*' || break
    read(io, Char)
  end
end

function parse_term!(io, exp::Vector{Int})
  @assert !eof(io)
  @match peek(io, Char) begin
    'a':'z' => begin
      parse_monomial!(io, exp)
      return 1
    end
    '0':'9' => begin
      coeff = parse_int(io)
      if !eof(io) && peek(io, Char) == '*'
        read(io, Char)
        parse_monomial!(io, exp)
        return coeff
      else
        return coeff
      end
    end
  end
end

function parse_polynomial(io, ctx::MPolyBuildCtx, exp::Vector{Int}, S::Ring)
  empty = true
  while !eof(io)
    @match peek(io, Char) begin
      '+' => begin
        read(io, Char)
        coeff = parse_term!(io, exp)
        push_term!(ctx, S(coeff), exp)
        exp .= 0
        empty = false
      end
      '-' => begin
        read(io, Char)
        coeff = parse_term!(io, exp)
        push_term!(ctx, S(-coeff), exp)
        exp .= 0
        empty = false
      end
      _ => begin
        empty || return finish(ctx)

        coeff = parse_term!(io, exp)
        push_term!(ctx, S(coeff), exp)
        exp .= 0
        empty = false
      end
    end
  end

  error("expected polynomial")
end

function parse_ideal(io, ctx::MPolyBuildCtx, exp::Vector{Int}, R::MPolyRing, s::Int)
  I = Vector{elem_type(R)}(undef, s)
  S = base_ring(R)

  for j in 1:s
    I[j] = parse_polynomial(io, ctx, exp, S)
    read(io, Char)  # skip delimiter | or line break at end
  end

  return I
end

function load_distribution(filename::String, n::Int, d::Int, s::Int)
  K = GF(32003)
  R, _ = polynomial_ring(K, n)
  D = Vector{Vector{elem_type(R)}}()

  let io = BufferedInputStream(open(filename))
    title = readline(io)
    @assert title == "Ideal"

    ctx = MPolyBuildCtx(R)
    exp = zeros(Int, n)
    while !eof(io)
      I = parse_ideal(io, ctx, exp, R, s)
      push!(D, I)
    end

    close(io)
  end

  return D
end

function load_stats(filename::String)
end
