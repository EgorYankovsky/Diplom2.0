using MathObjects;

namespace Solver;

public class OptimalMCG : ISolver
{
    private const int _maxIter = 10000;

    private const double _eps = 1E-15;

    public GlobalVector Solve(GlobalMatrix A, GlobalVector b)
    {
        GlobalVector x = new(b.Size);
        
        return x;
    }
}