using MathObjects;

namespace Solver;

public class LOS : ISolver
{
    private const int _maxIter = 100000;

    private const double _eps = 1E-15;

    public GlobalVector Solve(GlobalMatrix A, GlobalVector b)
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

        r_ = b - A * x_;
        z_ = r_;
        p_ = A * r_;

        int iter = 0;
        do
        {
            alph = (p_ * r_) / (p_ * p_);

            x = x_ + alph * z_;
            r = r_ - alph * p_;

            beta = -1.0D * (p_ * (A * r)) / (p_ * p_);
            z = r + beta * z_;
            p = A * r + beta * p_;

            iter++;

            x_ = x;
            z_ = z;
            r_ = r;
            p_ = p;
            Console.WriteLine($"{r.Norma() / b.Norma():E15}");
        } while (iter < _maxIter && r.Norma() / b.Norma() >= _eps);

        Console.WriteLine(
        $@"Computing finished!
Total iterations: {iter}
Relative residuality: {r.Norma() / b.Norma():E15}");
        return x;
    }
}