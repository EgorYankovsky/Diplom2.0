namespace DataStructs;

public class Elem
{
    // Массив всех ссылок на элементы.
    public int[] Arr = new int[8];

    /// <summary>
    /// Метод, возвращающий элемент в строчном формате.
    /// </summary>
    /// <returns>Строка в которой перечень номеров узлов в одном элементе.</returns>
    public override string ToString() => $"{Arr[0]} {Arr[1]} {Arr[2]} {Arr[3]} {Arr[4]} {Arr[5]} {Arr[6]} {Arr[7]}";

    /// <summary>
    /// Конструктор класса Elem.
    /// </summary>
    /// <param name="a"></param>
    /// <param name="b"></param>
    /// <param name="c"></param>
    /// <param name="d"></param>
    /// <param name="e"></param>
    /// <param name="f"></param>
    /// <param name="g"></param>
    /// <param name="h"></param>
    public Elem(int a, int b, int c, int d,
                int e, int f, int g, int h)
    {
        Arr[0] = a; Arr[1] = b; Arr[2] = c; Arr[3] = d;
        Arr[4] = e; Arr[5] = f; Arr[6] = g; Arr[7] = h;
    }
}