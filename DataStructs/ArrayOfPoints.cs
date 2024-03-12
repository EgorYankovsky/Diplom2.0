namespace DataStructs;

public class ArrayOfPoints
{
    // Путь из которого берем данные.
    private string _path = Path.GetFullPath("../../../../Data/Subtotals/Points.dat");

    // Массив точек.
    private readonly List<Point> _list;

    // Длина массива точек.
    public int Length { get; set; }

    /// <summary>
    /// Метод, возвращающий i-ую точку.
    /// </summary>
    /// <param name="i">Итератор.</param>
    /// <returns>Точка.</returns>
    public Point? this[int i] => _list[i];

    /// <summary>
    /// Конструктор класса.
    /// </summary>
    public ArrayOfPoints()
    {
        _list = new();
        using var sr = new StreamReader(_path);
        Length = int.Parse(sr.ReadLine() ?? "0");
        for (int i = 0; i < Length; i++)
            _list.Add(new Point(sr.ReadLine().Split().ToList()));
    }
}