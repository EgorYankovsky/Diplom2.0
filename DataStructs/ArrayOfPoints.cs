using System.Collections;

namespace DataStructs;

public abstract class ArrayOfPoints : IEnumerable
{
    public abstract int GetLength();

    public abstract IPoint? this[int i] { get; }

    public abstract IEnumerator GetEnumerator();

    public abstract void Append(IPoint p);
}

public class ArrayOfPoints2D : ArrayOfPoints
{
    private readonly List<Point2D> _list;

    public ArrayOfPoints2D(string path)
    {
        _list = [];
        var data = File.ReadAllText(path).Split("\n");
        int pointsAmount = int.Parse(data[0]);
        for (int i = 0; i < pointsAmount; i++)
            _list.Add(new Point2D([.. data[i + 1].Split()]));
    }

    public ArrayOfPoints2D(int pointsAmount) => _list = new(pointsAmount);

    public override Point2D this[int i] => _list[i];

    public override void Append(IPoint p) => _list.Add((Point2D)p);
    
    public override IEnumerator GetEnumerator() => _list.GetEnumerator();
    
    public override int GetLength() => _list.Count;
}

public class ArrayOfPoints3D : ArrayOfPoints
{
    private readonly List<Point3D> _list;

    public ArrayOfPoints3D(string path)
    {
        _list = [];
        var data = File.ReadAllText(path).Split("\n");
        int pointsAmount = int.Parse(data[0]);
        for (int i = 0; i < pointsAmount; i++)
            _list.Add(new Point3D([.. data[i + 1].Split()]));
    }
    
    public ArrayOfPoints3D(int pointsAmount) => _list = new(pointsAmount);

    public override Point3D this[int i] => _list[i];

    public override void Append(IPoint p) => _list.Add((Point3D)p);

    public override IEnumerator GetEnumerator() => _list.GetEnumerator();

    public override int GetLength() => _list.Count;
}