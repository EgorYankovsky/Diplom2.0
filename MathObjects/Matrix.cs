namespace MathObjects;

public abstract class Matrix : IMathObject
{
    public abstract double this[int i, int j]
    {
        get;
        set;
    }
}