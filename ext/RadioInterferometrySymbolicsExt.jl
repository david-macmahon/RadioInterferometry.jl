# Extra functionality for when Symbolics is being used.
module RadioInterferometrySymbolicsExt

if isdefined(Base, :get_extension)
    import RadioInterferometry: signclean
    import Symbolics
else
    import ..RadioInterferometry: signclean
    import ..Symbolics
end

export signclean

"""
A Symbolics rewriter for use with `simplify` that uses the trig identities:
    sin(-x) == -sin(x)
    cos(-x) ==  cos(x)
to rewrite sin/cos expressions without unary minus operations inside the
arguments to `sin` snd `cos`.
"""
signclean_rewriter = Symbolics.Rewriters.Prewalk(Symbolics.Rewriters.Chain([
    Symbolics.@rule       sin(-1*(~x))       =>        -sin(~x)
    Symbolics.@rule       cos(-1*(~x))       =>         cos(~x)
]))

"""
Cleanup any "negative" arguments in Symbolic `sin`/`cos` expressions in `r3`
(typically rotation matrix) using trig identities.  Returns a new Array.
"""
function signclean(r3::AbstractArray{<:Real})
    Symbolics.simplify.(r3; rewriter=signclean_rewriter)
end

end # module RadioInterferometrySymbolicsExt