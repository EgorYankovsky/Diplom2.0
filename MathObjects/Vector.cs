namespace MathObjects;

public abstract class Vector : IMathObject
{
    public double[]? _values;

    public virtual double this[int i]
    {
        get => i < _values.Length ? _values[i] : throw new IndexOutOfRangeException();
        set => _values[i] = value;
    }

    public int Size => _values.Length;

}