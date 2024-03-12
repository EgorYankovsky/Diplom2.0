namespace DataStructs;

public class ArrayOfElems
{
    private string _path = Path.GetFullPath("../../../../Data/Subtotals/Elems.dat");

    public List<List<int>> Arr;
    
    public List<double> mu0i;

    public int Length { get; set; }

    public MyEnumerator GetEnumerator() => new(this);

    public List<int> this[int i] => Arr[i];

    public class MyEnumerator
    {  
        int nIndex;  
        ArrayOfElems collection;  
        public MyEnumerator(ArrayOfElems coll)
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

    public ArrayOfElems()
    {
        Arr = new();
        mu0i = new();
        using var sr = new StreamReader(_path);
        Length = int.Parse(sr.ReadLine() ?? "0");
        for (int i = 0; i < Length; i++)
        {
            var info = sr.ReadLine().Split().ToList();
            Arr.Add(new List<int> {int.Parse(info[0]), int.Parse(info[1]), int.Parse(info[2]), int.Parse(info[3])});
            mu0i.Add(double.Parse(info[4]));
        }
    }
}