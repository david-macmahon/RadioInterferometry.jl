# Extra functionality for when Symbolics is being used.

using Symbolics

export signclean

"""
A Symbolics rewriter for use with `simplify` that uses the trig identities:
   sin(-x) == -sin(x)
   cos(-x) ==  cos(x)
to rewrite sin/cos expressions without unary minus operations inside the
arguments to `sin` snd `cos`.
"""
signclean_rewriter = Rewriters.Prewalk(Rewriters.Chain([
    @rule         sin(-1*(~x))         =>          -sin(~x)
    @rule         cos(-1*(~x))         =>           cos(~x)
    @rule (~w...)*sin(-1*(~x))*(~y...) => (~w...)*(-sin(~x))*(~y...)
    @rule (~w...)*cos(-1*(~x))*(~y...) => (~w...)*( cos(~x))*(~y...)
]))

"""
Cleanup any "negative" arguments in Symbolic `sin`/`cos` expressions in `r3`
(typically rotation matrix) using trig identities.  Returns a new Array.
"""
function signclean(r3=Array{Real})
    simplify.(r3; rewriter=signclean_rewriter)
end
