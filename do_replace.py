#!/usr/bin/env python3

# Read the file
with open('BMSZ2-revised.tex', 'r') as f:
    content = f.read()

lines = content.split('\n')

# We want lines 7705-7715 (1-indexed), which is 7704-7714 (0-indexed)
old_section_lines = lines[7704:7715]  # 7715 is exclusive in 0-indexed (so 7714 is last)
print(f"Old section has {len(old_section_lines)} lines")
print(f"First line: {repr(old_section_lines[0])}")
print(f"Last line: {repr(old_section_lines[-1])}")

old_section = '\n'.join(old_section_lines) + '\n'
print(f"\nOld section length: {len(old_section)}")

# Check if this is in the content
if old_section in content:
    print("Found old section in content!")
else:
    print("Old section NOT found in content")

# Now create the new section
new_section = """\\\smallskip\\noindent\\textbf{Normalization to $X_{S,T}$.}
Using the $\\GL_D(E)$-part of the Levi action (9), we may arrange
\\[
\\ker\\bar C=L'_0,\\qquad
\\bar C|_{L'_1}:L'_1\\xrightarrow{\\sim}V_0/\\operatorname{Im}X_0.
\\]
After identifying $V_0/\\operatorname{Im}X_0$ with the chosen complement $V_+$, this gives an isomorphism $S:L'_1\\xrightarrow{\\sim}V_+$.

For the $C$-block: since $C-C_S$ takes values in $\\operatorname{Im}X_0$, choose $R:E'\\to V_0$ such that $X_0R=C-C_S$.
By the unipotent formula (10), conjugating by $u_R$ replaces $C$ by $C_S$.

For the $B$-block: keep $C=C_S$ fixed.  To preserve this $C$-block, use only unipotents $u_R$ with $R(E')\\subset\\ker X_0$.
Then $X_0R=0$, and (10) becomes
\\[
u_RX_{C_S,B}u_R^{-1}=X_{C_S,B+R^\\vee C_S-C_S^\\vee R}.
\\tag{13}
\\]
Write $R=R_1\\oplus R_0$ with $R_1:L'_1\\to\\ker X_0$ and $R_0:L'_0\\to\\ker X_0$.
The correction term in (14) shows that the term $-S^\\vee R_0$ kills the $L'_0\\to L_1$ cross-block of $B$, and the term $R_1^\\vee S-S^\\vee R_1$ kills the $L'_1\\to L_1$ block.
Thus the only surviving block comes from some $T\\in\\mathscr T$, and the normal form is $X_{S,T}$.

Every point of the open locus is $H_{X_0}$-conjugate to one with $C$-block equal to $C_S$ and $B$-block corresponding to $T$.

\\\smallskip\\noindent\\textbf{The residual invariant.}
The stabilizer of the normalized $C$-block acts on $T$ by change of coordinates:
\\[
T\\longmapsto (\\gamma^{-1})^\\vee T\\gamma,
\\tag{15}
\\]
which is exactly the isometry class of the form $B_T(u,v)=\\langle Tu,v\\rangle$.
Consequently, the remaining $H_{X_0}$-orbits in the normalized slice are classified by $[B_T]\\in\\operatorname{Iso}_{r,-\\epsilon}.

Inside the normalized slice, the open part requires the residual form to have maximal possible rank, namely $T\\in\\operatorname{Form}_r(L'_0)$.
By (8), the assignment $[B_T]\\mapsto G\\cdot X_{S,T}$ gives the bijection between $\\operatorname{Iso}_{r,-\\epsilon}$ and the open $G$-orbits in $\\Ind_P^G(\\cO_0)$.

\\\smallskip\\noindent\\textbf{The multiplicity.}
Let $x=X_{S,T}$ with $T\\in\\operatorname{Form}_r(L'_0)$.
From (16)--(17) we have
\\[
\\ker x=E\\oplus\\ker T,\\qquad
\\Im x=(L_1\\oplus\\operatorname{Im}T)\\oplus V_0.
\\]

If $B_T$ is non-degenerate, then $\\ker T=0$, so $\\ker x=E$.
Every element of $Z_G(x)$ preserves $E$, hence lies in $P$, so $Z_G(x)\\subseteq P$ and $m(G\\cdot x,P)=1$.

If $B_T$ has one-dimensional radical (which occurs only when $D=\\mathbb R$, $\\epsilon=+1$, $l$ odd), then $F:=\\ker T$ is one-dimensional.
The quotient $\\ker x/(\\ker x\\cap\\operatorname{Im}x)\\simeq\\ell\\oplus F$ is a hyperbolic plane, whose two isotropic lines give $[Z_G(x):Z_G(x)\\cap P]=2$, hence $m(G\\cdot x,P)=2$.

\\end proof}
\\subsection{The generalized good descent case}\\label{sec:app-gd}
"""

# Replace
if old_section in content:
    new_content = content.replace(old_section, new_section, 1)
    with open('BMSZ2-revised.tex', 'w') as f:
        f.write(new_content)
    print("\nReplacement successful!")
else:
    print("\nCould not find old section for replacement")