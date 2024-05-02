using MathObjects;

namespace Solver;

public class LU_LOS : ISolver
{
    private const int _maxIter = 100_000;

    private const double _eps = 1E-15;

    private void PartitialLU(GlobalMatrix A)
    {
        for (int i = 0; i < A.Size; i++)
        {
            for (int j = A._ig[i]; j < A._ig[i + 1]; j++)
            {
                int jCol = A._jg[j];
                int jk = A._ig[jCol];
                int k = A._ig[i];

                int sdvig = A._jg[A._ig[i]] - A._jg[A._ig[jCol]];

                if (sdvig > 0)
                    jk += sdvig;
                else
                    k -= sdvig;

                double sumL = 0.0;
                double sumU = 0.0;

                for (; k < j && jk < A._ig[jCol + 1]; k++, jk++)
                {
                    sumL += A._al[k] * A._au[jk];
                    sumU += A._au[k] * A._al[jk];
                }

                A._al[j] -= sumL;
                A._au[j] -= sumU;
                A._au[j] /= A._diag[jCol];
            }

            double sumD = 0.0;
            for (int j = A._ig[i]; j < A._ig[i + 1]; j++)
                sumD += A._al[j] * A._au[j];

            A._diag[i] -= sumD;
        }
    }

    private static GlobalVector Forward(GlobalMatrix Matrix, GlobalVector b)
    {
        var result = new GlobalVector(b);

        for (int i = 0; i < Matrix.Size; i++)
        {
            for (int j = Matrix._ig[i]; j < Matrix._ig[i + 1]; j++)
            {
                var aaa = Matrix._al[j];
                var aaaa = result[Matrix._jg[j]];
                var aaaaa = result[i];
                if (double.IsNaN(aaa) || double.IsNaN(aaaa) || double.IsNaN(aaaaa) || double.IsNaN(result[i]) ||
                    double.IsInfinity(aaaa) || double.IsInfinity(aaaaa) || double.IsInfinity(aaa) || double.IsInfinity(result[i]))
                    Console.WriteLine();
            
                if (double.IsInfinity(Matrix._al[j] * result[Matrix._jg[j]]))
                    Console.WriteLine();           
                result[i] -= Matrix._al[j] * result[Matrix._jg[j]];
                if (double.IsNaN(aaa) || double.IsNaN(aaaa) || double.IsNaN(aaaaa) || double.IsNaN(result[i]) ||
                    double.IsInfinity(aaaa) || double.IsInfinity(aaaaa) || double.IsInfinity(aaa))
                    Console.WriteLine();
            }
            var a = Matrix._diag[i];
            var aa = result[i];
            result[i] /= Matrix._diag[i];
            if (double.IsNaN(result[i]))
                Console.WriteLine();
        }
        return result;
    }

    private static GlobalVector Backward(GlobalMatrix A, GlobalVector b)
    {
        var result = new GlobalVector(b);
        for (int i = A.Size - 1; i >= 0; i--)
            for (int j = A._ig[i + 1] - 1; j >= A._ig[i]; j--)
            {
                var a = A._au[j];
                var aa = result[i];
                var aaa = A._au[j] * result[i];
                result[A._jg[j]] -= A._au[j] * result[i];
            }
        return result;
    }

    public (GlobalVector, GlobalVector) Solve(GlobalMatrix A, GlobalVector b)
    {
        for (int i = 0; i < b.Size; i++)
            if (double.IsNaN(b[i]) || double.IsInfinity(b[i]))
                Console.WriteLine(i);

        GlobalVector x = new(b);
        GlobalVector x_ = new(b.Size);

        GlobalMatrix LU = new(A);

        PartitialLU(LU);

        for (int i = 0; i < LU._al.Length; i++)
            if (double.IsNaN(LU._al[i]) || double.IsInfinity(LU._al[i]))
                Console.WriteLine(i);

        for (int i = 0; i < LU._au.Length; i++)
            if (double.IsNaN(LU._au[i]) || double.IsInfinity(LU._au[i]))
                Console.WriteLine(i);

        GlobalVector r = Forward(LU, b - A * x);
        for (int i = 0; i < r.Size; i++)
            if (double.IsNaN(r[i]) || double.IsInfinity(r[i]))
                Console.WriteLine(i);

        GlobalVector z = Backward(LU, r);
        for (int i = 0; i < z.Size; i++)
            if (double.IsNaN(z[i]) || double.IsInfinity(z[i]))
                Console.WriteLine(i);

        GlobalVector p = Forward(LU, A * z);
        for (int i = 0; i < p.Size; i++)
            if (double.IsNaN(p[i]) || double.IsInfinity(p[i]))
                Console.WriteLine(i);

        GlobalVector tmp = new(b.Size);
        for (int i = 0; i < tmp.Size; i++)
            if (double.IsNaN(tmp[i]) || double.IsInfinity(tmp[i]))
                Console.WriteLine(i);

        GlobalVector r_ = new(b.Size);
        GlobalVector z_ = new(b.Size);
        GlobalVector p_ = new(b.Size);

        double alph = 0.0D;
        double beta = 0.0D;

        int iter = 0;
        do
        {
            x_ = x;
            z_ = z;
            r_ = r;
            p_ = p;

            var ds = p_ * p_;

            alph = (p_ * r_) / (p_ * p_);

            x = x_ + alph * z_;
            r = r_ - alph * p_;

            tmp = Forward(LU, A * Backward(LU, r));
            beta = -1.0D * (p_ * tmp) / (p_ * p_);

            z = Backward(LU, r) + beta * z_;
            p = tmp + beta * p_;

            iter++;
            Console.WriteLine($"{r.Norma() / b.Norma():E15}");
        } while (iter < _maxIter && r.Norma() / b.Norma() >= _eps);

        Console.WriteLine(
        $@"Computing finished!
Total iterations: {iter}
Relative residuality: {r.Norma() / b.Norma():E15}");
        return (x, x_);
    }
}