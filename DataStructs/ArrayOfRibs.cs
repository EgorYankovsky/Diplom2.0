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
    public Point3D a;
    public Point3D b;

    public Rib()
    {
        a = new Point3D(0.0D, 0.0D, 0.0D);
        b = new Point3D(1.0D, 1.0D, 1.0D);
        typeOfRib = TypeOfRib.NotStated;
    }

    public TypeOfRib typeOfRib;

    public Rib(Point3D a, Point3D b)
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

public class ArrayOfRibs(int size) : IEnumerable
{
    private List<Rib> _ribs = new(size);

    public int Count => _ribs.Count;

    public void Add(Rib rib) => _ribs.Add(rib);

    public void Remove(int i)
    {
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(i, _ribs.Count);
        _ribs.RemoveAt(i);
    }

    public Rib this[int i] => _ribs[i];

    public IEnumerator GetEnumerator() => _ribs.GetEnumerator();
}