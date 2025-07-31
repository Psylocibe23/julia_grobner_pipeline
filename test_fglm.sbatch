#!/bin/bash
#SBATCH --job-name=fglm_test
#SBATCH --output=logs/fglm_test_%j.out
#SBATCH --partition=turing-long
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00

module load sagemath/10.6-singularity

sage -c "
F = GF(2)
R = PolynomialRing(F, ['x0','x1','x2','x3'], order='deglex')
G_drl = [
    R('x0^3 + x0^2*x2 + x0*x3^2 + x0 + x1^3 + x1^2*x2 + x1^2*x3 + x1 + x2^3 + x2*x3^2 + x2 + x3^4 + x3^3 + x3^2 + x3'),
    R('x0^2*x3 + x0^2 + x0*x1^2 + x0 + x1^2*x3 + x1^2 + x1*x2^2 + x1*x3^2 + x1 + x2^4 + x2^2*x3 + x2*x3^2 + x3^2 + x3 + 1'),
    R('x0^3 + x0*x1^2 + x0*x2^2 + x0*x3^2 + x0 + x1^4 + x1^2*x2 + x1^2*x3 + x1^2 + x1*x3^2 + x2^3 + x3^2 + x3 + 1'),
    R('x0^4 + x0^2*x1 + x0^2*x3 + x0*x1^2 + x0*x3^2 + x0 + x1^3 + x1^2*x2 + x1^2*x3 + x1*x2^2 + x1 + x2^2 + x2 + x3^3 + x3^2 + 1')
]
I = R.ideal(G_drl)
Rlex = PolynomialRing(F, ['x0','x1','x2','x3'], order='lex')
G_lex = [Rlex(str(p)) for p in G_drl]
Ilex = Rlex.ideal(G_lex)
import time
start = time.time()
G_lex_fglm = Ilex.groebner_basis(algorithm='singular:stdfglm')
end = time.time()
print(f'FGLM finished in {end-start:.3f} seconds')
for g in G_lex_fglm: print(g)
"
