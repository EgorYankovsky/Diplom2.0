namespace DataStructs;

public class ArrayOfPoints
{
    // Путь из которого берем данные.
    private string _path = Path.GetFullPath("../../../../Data/Subtotals/Points.dat");

    // Массив точек.
    private List<Point> _list;

    // Длина массива точек.
    public int Length => _list.Count;

    /// <summary>
    /// Метод, возвращающий i-ую точку.
    /// </summary>
    /// <param name="i">Итератор.</param>
    /// <returns>Точка.</returns>
    public Point? this[int i] => _list[i];

    public void Append(Point p) => _list.Add(p);

    /// <summary>
    /// Конструктор класса.
    /// </summary>
    public ArrayOfPoints()
    {
        _list = new();
        using var sr = new StreamReader(_path);
        int length = int.Parse(sr.ReadLine() ?? "0");
        for (int i = 0; i < length; i++)
            _list.Add(new Point(sr.ReadLine().Split().ToList()));
    }

    public ArrayOfPoints(int pointsAmount)
    {
        _list = new(pointsAmount);
    }
}