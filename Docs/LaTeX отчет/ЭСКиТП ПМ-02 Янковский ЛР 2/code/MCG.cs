using MathObjects;

namespace Solver;

public class MCG : ISolver
{
    private const int _maxIter = 10000;

    private const double _eps = 1E-15;

    public GlobalVector Solve(GlobalMatrix A, GlobalVector b)
    {
        GlobalVector x = new(b.Size);
        GlobalVector x_ = new(b.Size);

        GlobalVector r = new(b.Size);
        GlobalVector r_ = new(b.Size);

        GlobalVector z = new(b.Size);
        GlobalVector z_ = new(b.Size);

        double alph = 0.0D;
        double beta = 0.0D;

        r_ = b - A * x_;
        z_ = r_;

        int iter = 0;
        do
        {
            alph = (r_ * r_) / ((A * z_) * z_);

            x = x_ + alph * z_;
            r = r_ - alph * (A * z_);

            beta = (r * r) / (r_ * r_);
            z = r + beta * z_;

            iter++;

            x_ = x;
            z_ = z;
            r_ = r;
            Console.WriteLine($"{r.Norma() / b.Norma():E15}");
        } while (iter < _maxIter && r.Norma() / b.Norma() >= _eps);

        Console.WriteLine(
        $@"Computing finished!
Total iterations: {iter}
Relative residuality: {r.Norma() / b.Norma():E15}");
        return x;
    }
}