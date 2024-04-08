using System.Collections;

namespace DataStructs;

public enum TypeOfRib
{
    Inner, 
    Outer,
    BoundaryI,
    BoundaryII,
    BoundaryIII,
    NotStated
}

public record Rib
{
    public Point a;
    public Point b;

    public Rib()
    {
        a = new Point(0.0D, 0.0D, 0.0D);
        b = new Point(1.0D, 1.0D, 1.0D);
    }

    public TypeOfRib typeOfRib;

    public Rib(Point a, Point b)
    {
        this.a = a;
        this.b = b;
        typeOfRib = (a.Type, b.Type) switch
        {
            (Location.Inside, not (Location.OutSide or Location.NotStated)) or
            (not (Location.OutSide or Location.NotStated), Location.Inside) => TypeOfRib.Inner,
            (Location.BoundaryII, Location.BoundaryII) => TypeOfRib.BoundaryII,
            (Location.BoundaryI, _) or (_, Location.BoundaryI) => TypeOfRib.BoundaryI,
            _ => TypeOfRib.NotStated,
        };

    }

    public double Length => Math.Sqrt((b.X - a.X) * (b.X - a.X) + (b.Y - a.Y) * (b.Y - a.Y) + (b.Z - a.Z) * (b.Z - a.Z));
}

public class ArrayOfRibs : IEnumerable
{
    private List<Rib> _ribs;

    public ArrayOfRibs(int size)
    {
        _ribs = new(size);
    }

    public int Count => _ribs.Count;

    public void Add(Rib rib) => _ribs.Add(rib);

    public void Remove(int i)
    {
        if (i >= _ribs.Count) throw new ArgumentOutOfRangeException("Out of range ArrayOfRibs");
        _ribs.RemoveAt(i);
    }

    public Rib this[int i] => _ribs[i];

    public IEnumerator GetEnumerator() => _ribs.GetEnumerator();
}