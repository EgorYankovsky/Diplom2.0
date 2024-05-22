using MathObjects;

namespace Solver;

public class LOS : ISolver
{
    private const int _maxIter = 100000;

    private const double _eps = 1E-15;

    public (GlobalVector, GlobalVector) Solve(GlobalMatrix A, GlobalVector b)
    {
        GlobalVector x = new(b.Size);
        GlobalVector x_ = new(b.Size);

        GlobalVector r = new(b.Size);
        GlobalVector r_ = new(b.Size);

        GlobalVector z = new(b.Size);
        GlobalVector z_ = new(b.Size);

        GlobalVector p = new(b.Size);
        GlobalVector p_ = new(b.Size);

        double alph = 0.0D;
        double beta = 0.0D;

        r = b - A * x;
        z = r;
        p = A * r;

        int iter = 0;
        do
        {
            x_ = x;
            z_ = z;
            r_ = r;
            p_ = p;
         
            alph = (p_ * r_) / (p_ * p_);

            x = x_ + alph * z_;
            r = r_ - alph * p_;

            beta = -1.0D * (p_ * (A * r)) / (p_ * p_);
            z = r + beta * z_;
            p = A * r + beta * p_;

            iter++;
            if (iter % 10 == 0)
                Console.WriteLine($"{r.Norma() / b.Norma():E15}");
        } while (iter < _maxIter && r.Norma() / b.Norma() >= _eps);

        Console.WriteLine(
        $@"Computing finished!
Total iterations: {iter}
Relative residuality: {r.Norma() / b.Norma():E15}");
        return (x, x_);
    }
}