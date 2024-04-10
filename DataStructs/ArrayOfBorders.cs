namespace DataStructs;

public class ArrayOfBorders
{
    // Относительный путь до файла с данными.
    private string _path = Path.GetFullPath("../../../../Data/Subtotals/Borders.dat");

    // Массив с содержимым ссылок на границы.
    public List<List<int>> Arr;

    // Длина массива.
    public int Length { get; set; }

    /// <summary>
    /// Метод, возвращающий i-ую границу.
    /// </summary>
    /// <param name="i">Индекс i.</param>
    /// <returns>Граница.</returns>
    public List<int> this[int i] => Arr[i];

    public MyEnumerator GetEnumerator() => new(this);

    public class MyEnumerator
    {  
        int nIndex;  
        ArrayOfBorders collection;  
        public MyEnumerator(ArrayOfBorders coll)
        {  
            collection = coll;  
            nIndex = -1;  
        }  
    
        public bool MoveNext()
        {  
            nIndex++;  
            return nIndex < collection.Arr.Count;  
        }  
    
        public List<int> Current => collection.Arr[nIndex];
    }

    /// <summary>
    /// Метод, сортирующий границы в порядке убывания краевого условия.
    /// </summary>
    private void SortList()
    {
        for (int i = 1; i < Length; i++)
            for (int j = 0; j < Length - i; j++)
                if (Arr[j][0] < Arr[j + 1][0])
                    (Arr[j], Arr[j + 1]) = (Arr[j + 1], Arr[j]);
    }

    /// <summary>
    /// Конструктор массива границ.
    /// </summary>
    public ArrayOfBorders()
    {
        Arr = new();
        using var sr = new StreamReader(_path);
        Length = int.Parse(sr.ReadLine() ?? "0");
        for (int i = 0; i < Length; i++)
            Arr.Add(sr.ReadLine().Split().Select(int.Parse).ToList());
        SortList();
    }

    public ArrayOfBorders(object o)
    {
        Arr = new();
    }
}
